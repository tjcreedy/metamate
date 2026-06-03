#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the pure-Python logic in ``metamate.core`` that backs
``filter-adaptive`` (and the validation-enforcement step shared by ``dump`` and
``filter-adaptive``).

These tests are deliberately independent of the metaMATE production scripts:
they build tiny, hand-computable fixtures in a temporary directory and call the
existing public functions in :mod:`metamate.core` without modifying them. They
need only ``numpy`` and ``biopython`` (already required by metaMATE) and run
with the standard library test runner::

    python -m unittest discover -s tests -v

What is covered:
* ``enforce_validation_fasta`` -- authentic ASVs re-added, non-authentic ASVs
  removed, and (regression) the ``refpass`` n ``nontarget`` overlap is *not*
  re-added. See ``metamate_filter_adaptive_bug_report.md``.
* ``adaptive_filter`` -- per-sample threshold application, target "rescue" under
  ``enforcevalidation`` on/off, the no-non-authentics fallback threshold, and
  the summary-table statistics.
* ``write_adaptive_csv`` -- the abundance table reflects the caller-supplied
  retained-ASV set (the column set), not just whatever is in the counts dict.
* ``load_validation_from_control`` / ``load_validation_from_otu_summary`` -- the
  control file produces overlapping target/nontarget sets (the upstream cause of
  the enforce bug) while the OTU summary produces disjoint sets.
"""

import os
import sys
import csv
import tempfile
import unittest

# Make the repo root importable regardless of the working directory.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import numpy

from metamate import core


class _Rec:
    """Minimal stand-in for a Bio.SeqRecord: ``adaptive_filter`` only reads
    ``record.seq`` and wraps it with ``str()``."""

    def __init__(self, seq):
        self.seq = seq


def _write_fasta(path, seqs):
    """seqs: dict {header: sequence}."""
    with open(path, "w") as f:
        for head, seq in seqs.items():
            f.write(f">{head}\n{seq}\n")


def _read_fasta(path):
    """Return {header: sequence} from a simple (one-line-per-seq) FASTA."""
    out = {}
    head = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                head = line[1:]
                out[head] = ""
            elif head is not None:
                out[head] += line
    return out


def _read_csv(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


class EnforceValidationFastaTests(unittest.TestCase):
    """``core.enforce_validation_fasta`` post-processes an output FASTA: it must
    add any missing authentics and remove all non-authentics."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        # Source contains the full ASV set with distinct sequences.
        self.source = os.path.join(self.tmp, "source.fasta")
        _write_fasta(self.source, {
            "a1": "AAAA",
            "a2": "CCCC",
            "a3": "GGGG",
            "a4": "TTTT",
            "a5": "ACGT",
        })
        self.out = os.path.join(self.tmp, "out.fasta")

    def test_disjoint_sets_add_authentics_remove_nonauthentics(self):
        # OTU-like case: target and nontarget are disjoint.
        # Output is missing authentic a2 and wrongly contains non-authentic a3.
        _write_fasta(self.out, {"a1": "AAAA", "a3": "GGGG"})
        target = {"a1", "a2"}
        nontarget = {"a3", "a4"}

        core.enforce_validation_fasta(self.out, self.source, target, nontarget)

        result = _read_fasta(self.out)
        self.assertEqual(set(result), {"a1", "a2"})       # a2 re-added, a3 removed
        self.assertEqual(result["a2"], "CCCC")            # sequence pulled from source
        self.assertNotIn("a3", result)
        self.assertNotIn("a4", result)

    def test_overlap_is_not_readded_regression(self):
        # Regression for metamate_filter_adaptive_bug_report.md:
        # a5 is in BOTH target (refpass) and nontarget (length/stop fail). After
        # enforcement it must be removed and must NOT be re-added.
        _write_fasta(self.out, {"a1": "AAAA", "a3": "GGGG", "a5": "ACGT"})
        target = {"a1", "a2", "a5"}     # a5 = refpass n nontarget overlap
        nontarget = {"a3", "a4", "a5"}

        core.enforce_validation_fasta(self.out, self.source, target, nontarget)

        result = _read_fasta(self.out)
        self.assertEqual(set(result), {"a1", "a2"})
        self.assertNotIn("a5", result,
                         "overlap (refpass n nontarget) must not be re-added")
        self.assertEqual(result["a2"], "CCCC")

    def test_noop_when_already_correct(self):
        _write_fasta(self.out, {"a1": "AAAA", "a2": "CCCC"})
        core.enforce_validation_fasta(self.out, self.source,
                                      target={"a1", "a2"}, nontarget=set())
        self.assertEqual(_read_fasta(self.out), {"a1": "AAAA", "a2": "CCCC"})


class AdaptiveFilterTests(unittest.TestCase):
    """``core.adaptive_filter`` with the default ``verified_removed`` criteria.

    Fixture (percentile 0.5 so thresholds are the medians of per-library
    non-authentic counts)::

        target    = {t_low, t_high}
        nontarget = {n1, n2}

        S1: t_low=1 t_high=9 n1=2 n2=4 x_keep=10 x_drop=3   -> median(2,4)=3.0
        S2: t_low=1 x3=5                                    -> no n* -> fallback
        S3: n1=7 n2=8 x4=1                                  -> median(7,8)=7.5

        fallback = mean(valid thresholds) = mean(3.0, 7.5) = 5.25
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.fasta = os.path.join(self.tmp, "adaptive.fasta")
        self.summary = os.path.join(self.tmp, "adaptive_summary.csv")

        self.target = {"t_low", "t_high"}
        self.nontarget = {"n1", "n2"}
        self.librarycounts = {
            "S1": {"t_low": 1, "t_high": 9, "n1": 2, "n2": 4, "x_keep": 10, "x_drop": 3},
            "S2": {"t_low": 1, "x3": 5},
            "S3": {"n1": 7, "n2": 8, "x4": 1},
        }
        # adaptive_filter reads asvs[name].seq for the retained sequences.
        all_names = {a for c in self.librarycounts.values() for a in c}
        self.asvs = {n: _Rec(n.upper()) for n in all_names}

        self.expected_t1 = numpy.percentile([2, 4], 50)     # 3.0
        self.expected_t3 = numpy.percentile([7, 8], 50)     # 7.5
        self.expected_fallback = numpy.mean([self.expected_t1, self.expected_t3])  # 5.25

    def _run(self, enforcevalidation):
        return core.adaptive_filter(
            self.asvs, self.librarycounts, self.target, self.nontarget,
            percentile=0.5, criteria="verified_removed",
            output_fasta=self.fasta, output_summary=self.summary,
            enforcevalidation=enforcevalidation,
        )

    def test_per_sample_thresholds_and_fallback(self):
        _, sample_thresholds, fallback, _ = self._run(enforcevalidation=True)
        self.assertAlmostEqual(sample_thresholds["S1"], self.expected_t1)
        self.assertAlmostEqual(sample_thresholds["S3"], self.expected_t3)
        # S2 has no non-authentic ASVs, so it inherits the fallback threshold.
        self.assertAlmostEqual(fallback, self.expected_fallback)
        self.assertAlmostEqual(sample_thresholds["S2"], self.expected_fallback)

    def test_threshold_applied_strictly_greater_than(self):
        filtered, _, _, _ = self._run(enforcevalidation=True)
        # S1 threshold 3.0: keep counts strictly > 3 (or rescued targets).
        # n2=4 survives the threshold; x_drop=3 does not (3 is not > 3).
        self.assertIn("n2", filtered["S1"])
        self.assertNotIn("x_drop", filtered["S1"])
        self.assertNotIn("n1", filtered["S1"])      # n1=2 below threshold

    def test_target_rescued_when_enforcevalidation_on(self):
        filtered, _, _, retained = self._run(enforcevalidation=True)
        # t_low=1 is below every threshold but is authentic -> rescued.
        self.assertIn("t_low", filtered["S1"])
        self.assertIn("t_low", filtered["S2"])
        self.assertIn("t_low", retained)
        self.assertEqual(
            filtered["S1"], {"t_low": 1, "t_high": 9, "n2": 4, "x_keep": 10})
        self.assertEqual(retained, {"t_low", "t_high", "n2", "x_keep"})

    def test_target_not_rescued_when_enforcevalidation_off(self):
        filtered, _, _, retained = self._run(enforcevalidation=False)
        # Without rescue, the below-threshold authentic t_low is dropped.
        self.assertNotIn("t_low", retained)
        self.assertNotIn("t_low", filtered["S1"])
        self.assertEqual(filtered["S1"], {"t_high": 9, "n2": 4, "x_keep": 10})
        self.assertEqual(retained, {"t_high", "n2", "x_keep"})

    def test_fasta_contains_exactly_retained_set(self):
        _, _, _, retained = self._run(enforcevalidation=True)
        written = _read_fasta(self.fasta)
        self.assertEqual(set(written), retained)
        # Sequences come from the supplied ASV records.
        self.assertEqual(written["t_high"], "T_HIGH")

    def test_summary_statistics(self):
        self._run(enforcevalidation=True)
        with open(self.summary, newline="") as fh:
            rows = {r["Sample"]: r for r in csv.DictReader(fh)}
        s1 = rows["S1"]
        self.assertEqual(int(s1["Total_ASVs_Pre"]), 6)
        self.assertEqual(int(s1["Verified_Authentic_Pre"]), 2)
        self.assertEqual(int(s1["Verified_NonAuthentic_Pre"]), 2)
        self.assertEqual(int(s1["Total_ASVs_Post"]), 4)   # t_low,t_high,n2,x_keep
        # Post-authentic counts *threshold survivors only* (t_low was rescued,
        # so it does not count here) -> only t_high.
        self.assertEqual(int(s1["Verified_Authentic_Post"]), 1)
        # n2=4 survived the threshold and is non-authentic.
        self.assertEqual(int(s1["Verified_NonAuthentic_Post"]), 1)

    def test_estimated_removed_picks_threshold(self):
        # 'estimated_removed' estimates each library's true authentic/non-authentic
        # composition (from the verified retention rates) and searches the observed
        # counts for the threshold whose estimated non-authentic removal proportion
        # is closest to --percentile. It needs >=10 ASVs and >=1 verified authentic
        # and >=1 verified non-authentic in the library.
        target = {"t1", "t2", "t3"}
        nontarget = {"n1", "n2", "n3", "n4"}
        # 12 ASVs: 3 authentic (high), 4 non-authentic (1..4), 5 unclassified (5..9).
        counts = {"t1": 10, "t2": 20, "t3": 30,
                  "n1": 1, "n2": 2, "n3": 3, "n4": 4,
                  "u1": 5, "u2": 6, "u3": 7, "u4": 8, "u5": 9}
        lib = {"BIG": counts}
        asvs = {n: _Rec(n) for n in counts}

        _, sample_thresholds, fallback, _ = core.adaptive_filter(
            asvs, lib, target, nontarget, percentile=0.5,
            criteria="estimated_removed",
            output_fasta=os.path.join(self.tmp, "est.fasta"),
            output_summary=os.path.join(self.tmp, "est.csv"),
            enforcevalidation=True)

        # A real threshold is computed (not the all-skipped fallback). At t=2,
        # half the verified non-authentics (n3,n4) survive -> 50% removed,
        # matching percentile 0.5.
        self.assertEqual(sample_thresholds["BIG"], 2)
        survivors = [c for a, c in counts.items() if a in nontarget and c > 2]
        self.assertAlmostEqual(1 - len(survivors) / 4, 0.5)

    def test_estimated_removed_falls_back_when_library_too_small(self):
        # Fewer than 10 ASVs -> estimation is skipped -> fallback threshold.
        target, nontarget = {"t1"}, {"n1"}
        lib = {"S": {"t1": 9, "n1": 3, "u1": 4}}
        asvs = {n: _Rec(n) for n in ("t1", "n1", "u1")}
        _, sample_thresholds, fallback, _ = core.adaptive_filter(
            asvs, lib, target, nontarget, percentile=0.5,
            criteria="estimated_removed",
            output_fasta=os.path.join(self.tmp, "est2.fasta"),
            output_summary=os.path.join(self.tmp, "est2.csv"),
            enforcevalidation=True)
        self.assertEqual(sample_thresholds["S"], fallback)

    def test_all_samples_skipped_falls_back_to_one(self):
        # No non-authentic ASVs anywhere -> no valid thresholds -> fallback 1.
        target = {"t1"}
        nontarget = {"nX"}     # present in no library
        lib = {"S1": {"t1": 1, "a": 2, "b": 5}}
        asvs = {n: _Rec(n) for n in ("t1", "a", "b")}
        filtered, sample_thresholds, fallback, retained = core.adaptive_filter(
            asvs, lib, target, nontarget, percentile=0.5,
            criteria="verified_removed",
            output_fasta=os.path.join(self.tmp, "f2.fasta"),
            output_summary=os.path.join(self.tmp, "s2.csv"),
            enforcevalidation=True)
        self.assertEqual(fallback, 1)
        self.assertEqual(sample_thresholds["S1"], 1)
        # keep counts > 1 or rescued target: t1(rescued), b=5; a=2>1 kept.
        self.assertEqual(filtered["S1"], {"t1": 1, "a": 2, "b": 5})
        self.assertEqual(retained, {"t1", "a", "b"})


class WriteAdaptiveCsvTests(unittest.TestCase):
    """``core.write_adaptive_csv`` -- the table columns are driven by the
    caller-supplied retained-ASV set, which lets ``metamate.py`` exclude
    enforced non-authentics (``retained - ev_nontarget``) so the table matches
    the enforced FASTA."""

    def test_columns_follow_retained_set_not_counts(self):
        tmp = tempfile.mkdtemp()
        out = os.path.join(tmp, "table.csv")
        filtered = {"S1": {"a": 5, "c": 2}, "S2": {"b": 3}}
        retained = {"a", "b"}      # "c" deliberately excluded (e.g. an enforced nontarget)
        thresholds = {"S1": 3.0, "S2": 5.25}

        core.write_adaptive_csv(filtered, retained, thresholds,
                                fallback_threshold=5.25, output_csv=out)

        rows = _read_csv(out)
        self.assertEqual(rows[0], ["Sample", "Threshold", "a", "b"])  # sorted, no "c"
        body = {r[0]: r for r in rows[1:]}
        self.assertEqual(body["S1"], ["S1", "3.0", "5", "0"])  # b absent -> 0
        self.assertEqual(body["S2"], ["S2", "5.25", "0", "3"])


class LoadValidationTests(unittest.TestCase):

    def test_control_file_sets_overlap(self):
        # The control file's refpass rows hold the full reference-match set, so
        # target and nontarget can overlap (the upstream cause of the enforce
        # bug). load_validation_from_control must reproduce that faithfully.
        tmp = tempfile.mkdtemp()
        ctrl = os.path.join(tmp, "x_control.txt")
        with open(ctrl, "w") as f:
            f.write("lengthfail\tL1\n")
            f.write("lengthfail\tBOTH\n")
            f.write("stopfail\tS1\n")
            f.write("refpass\tR1\n")
            f.write("refpass\tR2\n")
            f.write("refpass\tBOTH\n")     # overlap: also lengthfail

        target, nontarget = core.load_validation_from_control(ctrl)
        self.assertEqual(target, {"R1", "R2", "BOTH"})
        self.assertEqual(nontarget, {"L1", "S1", "BOTH"})
        self.assertEqual(target & nontarget, {"BOTH"},
                         "control-derived sets overlap; enforce step must handle it")

    def test_otu_summary_sets_disjoint(self):
        tmp = tempfile.mkdtemp()
        summ = os.path.join(tmp, "otu_summary.csv")
        with open(summ, "w") as f:
            f.write("ASV,OTU,ASV_Status,OTU_Status\n")
            f.write("asv1,otu1,RefMatch,Authentic\n")
            f.write("asv2,otu1,NonAuthentic,Authentic\n")
            f.write("asv3,otu2,NonAuthentic,Non-Authentic\n")
            f.write("asv4,otu3,Unclassified,Unclassified\n")

        target, nontarget = core.load_validation_from_otu_summary(summ)
        self.assertEqual(target, {"otu1"})
        self.assertEqual(nontarget, {"otu2"})
        self.assertEqual(target & nontarget, set(), "OTU sets must be disjoint")
        self.assertNotIn("otu3", target | nontarget)   # Unclassified -> neither


if __name__ == "__main__":
    unittest.main()
