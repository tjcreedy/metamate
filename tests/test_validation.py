#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the two control-group determinants used by ``find`` and
``filter-adaptive``: length filtering (``metamate.filterlength``) and
translation/stop-codon filtering (``metamate.filtertranslate``).

Independent of the production scripts; run with::

    python -m unittest discover -s tests -p 'test_*.py' -v
"""

import os
import sys
import types
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from metamate import filterlength, filtertranslate


class _Parser:
    """Stub argparse parser: ``error`` raises SystemExit like the real one."""

    def error(self, message):
        raise SystemExit(f"parser.error: {message}")


def _len_args(**kw):
    base = dict(expectedlength=None, minimumlength=0, maximumlength=float("inf"),
                basesvariation=None, percentvariation=None, codonsvariation=None,
                onlyvarybycodon=False)
    base.update(kw)
    return types.SimpleNamespace(**base)


class CheckLengthTests(unittest.TestCase):
    """``check_length(seq, [min, max], expected, varbycodon)`` passes a sequence
    if it is within [min, max] and (optionally) differs from ``expected`` by a
    whole number of codons."""

    def test_within_range(self):
        self.assertTrue(filterlength.check_length("A" * 60, [54, 66], 60))

    def test_too_short(self):
        self.assertFalse(filterlength.check_length("A" * 30, [54, 66], 60))

    def test_too_long(self):
        self.assertFalse(filterlength.check_length("A" * 70, [54, 66], 60))

    def test_only_vary_by_codon_rejects_frameshift_length(self):
        # 61 bp is in range but differs from 60 by 1 (not a multiple of 3).
        self.assertFalse(filterlength.check_length("A" * 61, [54, 66], 60, varbycodon=True))

    def test_only_vary_by_codon_accepts_codon_multiple(self):
        # 63 bp differs from 60 by exactly 3.
        self.assertTrue(filterlength.check_length("A" * 63, [54, 66], 60, varbycodon=True))

    def test_multi_pass_and_fail_lists(self):
        seqs = {"ok60": "A" * 60, "short": "A" * 30, "ok66": "A" * 66}
        args = _len_args(expectedlength=60, onlyvarybycodon=False)
        self.assertEqual(set(filterlength.check_length_multi(seqs, [54, 66], args)),
                         {"ok60", "ok66"})
        self.assertEqual(filterlength.check_length_multi(seqs, [54, 66], args, fail=True),
                         ["short"])


class ResolveLengthSpecTests(unittest.TestCase):
    """``resolve_length_spec`` turns the various length arguments into a single
    [minimum, maximum] range (mirroring how metaMATE's CLI supplies them, with
    ``None`` defaults)."""

    def test_percent_variation(self):
        args, lset = filterlength.resolve_length_spec(
            _len_args(expectedlength=60, percentvariation=0), _Parser())
        self.assertEqual((args.minimumlength, args.maximumlength), (60, 60))
        self.assertTrue(lset)

    def test_bases_variation(self):
        args, _ = filterlength.resolve_length_spec(
            _len_args(expectedlength=60, basesvariation=6), _Parser())
        self.assertEqual((args.minimumlength, args.maximumlength), (54, 66))

    def test_codons_variation_matches_bases(self):
        # Regression: --codonsvariation 2 must equal --basesvariation 6 (README),
        # i.e. 3 * codons. Previously this raised TypeError (used basesvariation,
        # which defaults to None).
        args, _ = filterlength.resolve_length_spec(
            _len_args(expectedlength=60, codonsvariation=2), _Parser())
        self.assertEqual((args.minimumlength, args.maximumlength), (54, 66))

    def test_minmax_only(self):
        args, lset = filterlength.resolve_length_spec(
            _len_args(minimumlength=400, maximumlength=420), _Parser())
        self.assertEqual((args.minimumlength, args.maximumlength), (400, 420))
        self.assertTrue(lset)

    def test_expected_with_minmax_errors(self):
        with self.assertRaises(SystemExit):
            filterlength.resolve_length_spec(
                _len_args(expectedlength=60, minimumlength=50), _Parser())

    def test_two_variations_errors(self):
        with self.assertRaises(SystemExit):
            filterlength.resolve_length_spec(
                _len_args(expectedlength=60, basesvariation=6, percentvariation=10), _Parser())

    def test_only_vary_by_codon_without_variation_errors(self):
        with self.assertRaises(SystemExit):
            filterlength.resolve_length_spec(
                _len_args(expectedlength=60, onlyvarybycodon=True), _Parser())


def _rec(seq):
    return SeqRecord(Seq(seq), id="x")


class TranslationTests(unittest.TestCase):
    """``stopcount`` counts stop codons per reading frame; ``min_stopcount``
    returns the minimum across all six frames; ``check_stops_multi`` flags a
    sequence as non-authentic only if *every* frame contains a stop."""

    def test_stopcount_per_frame(self):
        # Table 5 (invertebrate mito): TAA is a stop, AAA (Lys) is not.
        # "TAAAAAAAA" -> frame1 TAA AAA AAA = 1 stop; frames 2 and 3 = 0.
        counts = filtertranslate.stopcount(_rec("TAAAAAAAA"), 5, frame=(1, 2, 3))
        self.assertEqual(counts, [1, 0, 0])

    def test_stopcount_single_frame_returns_scalar(self):
        self.assertEqual(filtertranslate.stopcount(_rec("TAAAAAAAA"), 5, frame=1), 1)

    def test_min_stopcount_zero_when_a_clean_frame_exists(self):
        # A clean forward ORF has a stop-free frame, so the 6-frame minimum is 0.
        self.assertEqual(filtertranslate.min_stopcount(_rec("ATG" + "AAA" * 9), 5), 0)

    def test_check_stops_multi_clean_sequence_not_flagged(self):
        seqs = {"clean": _rec("ATG" + "AAA" * 9)}
        args = types.SimpleNamespace(table=5)
        # fail=False returns sequences WITHOUT stops in their best frame.
        self.assertEqual(filtertranslate.check_stops_multi(seqs, args), ["clean"])
        # fail=True (the non-authentic set) returns nothing for a clean sequence.
        self.assertEqual(filtertranslate.check_stops_multi(seqs, args, fail=True), [])


if __name__ == "__main__":
    unittest.main()
