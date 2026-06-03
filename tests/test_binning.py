#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the OTU/library parsing helpers in ``metamate.binning`` that
``filter-adaptive`` and OTU mode rely on. Independent of the production scripts:
fixtures are built in a temp directory and the existing functions are called
as-is.

Run with::

    python -m unittest discover -s tests -v
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from metamate import binning


class ParseUcTests(unittest.TestCase):
    """``binning.parse_uc`` maps every ASV to its OTU centroid: ``S`` (seed)
    lines map a centroid to itself, ``H`` (hit) lines map a member to its
    centroid, and ``C`` (cluster summary) lines are ignored."""

    def test_seed_hit_and_cluster_lines(self):
        tmp = tempfile.mkdtemp()
        uc = os.path.join(tmp, "clusters.uc")
        with open(uc, "w") as f:
            # type cluster size %id strand - - cigar query target
            f.write("S\t0\t250\t*\t+\t*\t*\t*\tasv1\t*\n")
            f.write("H\t0\t250\t99.2\t+\t0\t0\t250M\tasv2\tasv1\n")
            f.write("H\t0\t250\t98.0\t+\t0\t0\t250M\tasv3\tasv1\n")
            f.write("S\t1\t250\t*\t+\t*\t*\t*\tasv4\t*\n")
            f.write("C\t0\t3\t*\t*\t*\t*\t*\tasv1\t*\n")   # ignored

        mapping = binning.parse_uc(uc)
        self.assertEqual(mapping, {
            "asv1": "asv1",   # centroid -> itself
            "asv2": "asv1",   # member  -> centroid
            "asv3": "asv1",
            "asv4": "asv4",
        })


class DetectAlignedTests(unittest.TestCase):
    """``binning.detect_aligned`` reports a fasta as aligned only when *every*
    sequence shares the same length. It must scan all records, not just the
    first few -- an earlier version only looked at the first ``n`` sequences
    and so missed divergence that appeared later in the file."""

    def _write_fasta(self, records):
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, "asvs.fasta")
        with open(path, "w") as f:
            for name, seq in records:
                f.write(f">{name}\n{seq}\n")
        return path

    def test_equal_length_sequences_are_aligned(self):
        path = self._write_fasta([
            ("a", "ACGT--ACGT"),
            ("b", "ACGTACGTAC"),
            ("c", "AC--GTACGT"),
        ])
        self.assertTrue(binning.detect_aligned(path))

    def test_divergence_only_at_end_is_not_aligned(self):
        # 50 identical-length records, then one longer record at the very end.
        # This is exactly the case the old first-n-only check would miss.
        records = [(f"s{i}", "ACGT" * 5) for i in range(50)]
        records.append(("last", "ACGT" * 6))
        path = self._write_fasta(records)
        self.assertFalse(binning.detect_aligned(path))

    def test_empty_file_is_not_aligned(self):
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, "empty.fasta")
        open(path, "w").close()
        self.assertFalse(binning.detect_aligned(path))


class ParseReadmapTests(unittest.TestCase):
    """``binning.parse_readmap`` reads a sample-by-ASV count table, auto-detects
    the separator and orientation, drops zero counts, and returns
    ``(counts_by_library, {'total': per_asv_totals})``."""

    def test_asvs_as_columns_tab_separated(self):
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, "readmap.tsv")
        with open(path, "w") as f:
            f.write("sample\tasv1\tasv2\tasv3\n")   # leading label cell is tolerated
            f.write("S1\t5\t0\t2\n")
            f.write("S2\t0\t3\t4\n")

        master = {"asv1": None, "asv2": None, "asv3": None}
        counts_by_library, totals = binning.parse_readmap(master, path)

        self.assertEqual(counts_by_library["S1"], {"asv1": 5, "asv3": 2})  # 0s dropped
        self.assertEqual(counts_by_library["S2"], {"asv2": 3, "asv3": 4})
        self.assertEqual(totals["total"], {"asv1": 5, "asv2": 3, "asv3": 6})


if __name__ == "__main__":
    unittest.main()
