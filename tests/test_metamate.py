#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for two pieces of ``metamate.metamate``:

* ``subset_table`` -- filtering an abundance table to a FASTA's IDs (the
  dump / filter-adaptive table<->FASTA consistency step), including separator
  detection and per-sample cell zeroing.
* ``getcliargs`` -- the required-input validation rules for each mode, i.e.
  exactly what every mode demands as input.

Independent of the production scripts; run with::

    python -m unittest discover -s tests -p 'test_*.py' -v
"""

import os
import sys
import csv
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from metamate import metamate


def _write(path, text):
    with open(path, "w", newline="") as f:
        f.write(text)


def _read_rows(path, sep):
    with open(path, newline="") as f:
        return list(csv.reader(f, delimiter=sep))


class SubsetTableTests(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.fasta = os.path.join(self.tmp, "kept.fasta")
        _write(self.fasta, ">keep1\nACGT\n>keep2\nACGT\n")
        self.out = os.path.join(self.tmp, "out.txt")

    def test_keeps_only_fasta_ids(self):
        table = os.path.join(self.tmp, "t.tsv")
        _write(table, "id\tS1\tS2\nkeep1\t5\t0\ndrop1\t9\t9\nkeep2\t3\t7\n")
        metamate.subset_table(self.fasta, table, self.out)
        rows = _read_rows(self.out, "\t")
        self.assertEqual(rows[0], ["id", "S1", "S2"])
        self.assertEqual({r[0] for r in rows[1:]}, {"keep1", "keep2"})  # drop1 removed

    def test_detects_comma_separator(self):
        table = os.path.join(self.tmp, "t.csv")
        _write(table, "id,S1,S2\nkeep1,5,0\ndrop1,9,9\nkeep2,3,7\n")
        metamate.subset_table(self.fasta, table, self.out)
        rows = _read_rows(self.out, ",")
        self.assertEqual(len(rows), 3)               # header + keep1 + keep2
        self.assertEqual(rows[0], ["id", "S1", "S2"])

    def test_zeroes_cells_not_in_filtered_counts(self):
        table = os.path.join(self.tmp, "t.tsv")
        _write(table, "id\tS1\tS2\nkeep1\t5\t0\nkeep2\t3\t7\n")
        # keep2 survived only in S2; its S1 cell must be zeroed.
        filtered = {"S1": {"keep1": 5}, "S2": {"keep2": 7}}
        metamate.subset_table(self.fasta, table, self.out,
                              filtered_library_counts=filtered)
        body = {r[0]: r for r in _read_rows(self.out, "\t")[1:]}
        self.assertEqual(body["keep1"], ["keep1", "5", "0"])
        self.assertEqual(body["keep2"], ["keep2", "0", "7"])   # S1 zeroed


class GetCliArgsValidationTests(unittest.TestCase):
    """Each mode's required-input rules. Error paths raise SystemExit (argparse
    ``parser.error``); the success paths confirm a valid minimal invocation."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _out(self, name="out"):
        return os.path.join(self.tmp, name)

    def test_filter_adaptive_requires_references(self):
        with self.assertRaises(SystemExit):
            metamate.getcliargs(["filter-adaptive", "-A", "a.fasta", "-M", "m.tsv",
                                 "--expectedlength", "60", "--percentvar", "0",
                                 "-o", self._out()])

    def test_libraries_and_readmap_mutually_exclusive(self):
        with self.assertRaises(SystemExit):
            metamate.getcliargs(["filter-adaptive", "-A", "a.fasta",
                                 "-L", "lib.fastq", "-M", "m.tsv",
                                 "-R", "r.fasta", "--expectedlength", "60",
                                 "--percentvar", "0", "-o", self._out()])

    def test_neither_libraries_nor_readmap_errors(self):
        with self.assertRaises(SystemExit):
            metamate.getcliargs(["filter-adaptive", "-A", "a.fasta",
                                 "-R", "r.fasta", "--expectedlength", "60",
                                 "--percentvar", "0", "-o", self._out()])

    def test_partial_otu_arguments_error(self):
        # Supplying one OTU argument requires all three.
        with self.assertRaises(SystemExit):
            metamate.getcliargs(["filter-adaptive", "-A", "a.fasta",
                                 "-R", "r.fasta", "--expectedlength", "60",
                                 "--percentvar", "0", "--uc", "x.uc",
                                 "-o", self._out()])

    def test_dump_resultcache_requires_resultindex(self):
        with self.assertRaises(SystemExit):
            metamate.getcliargs(["dump", "-A", "a.fasta", "-o", self._out(),
                                 "-C", "cache"])

    def test_filter_adaptive_minimal_success(self):
        args = metamate.getcliargs(["filter-adaptive", "-A", "a.fasta", "-M", "m.tsv",
                                    "-R", "r.fasta", "--expectedlength", "60",
                                    "--percentvar", "0", "-o", self._out("ok")])
        self.assertEqual(args.mode, "filter-adaptive")
        self.assertFalse(args.otu_mode)
        self.assertEqual((args.minimumlength, args.maximumlength), (60, 60))
        self.assertEqual(args.refmatchlength, int(0.8 * 60))   # default = 80% of min length
        self.assertTrue(os.path.isdir(self._out("ok")))        # output dir created

    def test_otu_mode_success_without_library_source(self):
        # In OTU mode, per-library counts come from --otu_table, so -L/-M aren't
        # required.
        args = metamate.getcliargs(["filter-adaptive", "-A", "a.fasta",
                                    "-R", "r.fasta", "--expectedlength", "60",
                                    "--percentvar", "0",
                                    "--uc", "x.uc", "--otu_fasta", "o.fasta",
                                    "--otu_table", "o.csv", "-o", self._out("otu")])
        self.assertTrue(args.otu_mode)

    def test_optional_second_reference_defaults_to_none(self):
        args = metamate.getcliargs(["filter-adaptive", "-A", "a.fasta", "-M", "m.tsv",
                                    "-R", "r.fasta", "--expectedlength", "60",
                                    "--percentvar", "0", "-o", self._out("noref2")])
        self.assertIsNone(args.references2)

    def test_optional_second_reference_is_parsed(self):
        args = metamate.getcliargs(["filter-adaptive", "-A", "a.fasta", "-M", "m.tsv",
                                    "-R", "r.fasta", "--references2", "local.fasta",
                                    "--expectedlength", "60", "--percentvar", "0",
                                    "-o", self._out("ref2")])
        self.assertEqual(args.references2, "local.fasta")


class CombineReferencesTests(unittest.TestCase):
    """``filterreference.combine_references`` merges the primary reference with
    an optional second (e.g. more local) reference. A single usable path is
    returned unchanged; two or more are concatenated into one fasta in the
    working directory so downstream matching treats them as one set."""

    def _write_fasta(self, name, records):
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, name)
        with open(path, "w") as f:
            for rid, seq in records:
                f.write(f">{rid}\n{seq}\n")
        return path

    def test_single_reference_returned_unchanged(self):
        from metamate import filterreference
        r1 = self._write_fasta("r1.fasta", [("big1", "ACGTACGT")])
        wd = tempfile.mkdtemp()
        # A None second reference (the unset optional) is ignored.
        self.assertEqual(filterreference.combine_references([r1, None], wd), r1)

    def test_two_references_are_concatenated(self):
        from metamate import filterreference
        from Bio import SeqIO
        r1 = self._write_fasta("r1.fasta", [("big1", "ACGTACGT"), ("big2", "TTTTGGGG")])
        r2 = self._write_fasta("r2.fasta", [("local1", "CCCCAAAA")])
        wd = tempfile.mkdtemp()
        combined = filterreference.combine_references([r1, r2], wd)
        self.assertEqual(os.path.dirname(combined), wd)        # written into wd
        ids = [rec.id for rec in SeqIO.parse(combined, "fasta")]
        self.assertEqual(ids, ["big1", "big2", "local1"])      # all records, in order

    def test_missing_reference_exits(self):
        from metamate import filterreference
        r1 = self._write_fasta("r1.fasta", [("big1", "ACGTACGT")])
        wd = tempfile.mkdtemp()
        with self.assertRaises(SystemExit):
            filterreference.combine_references([r1, os.path.join(wd, "nope.fasta")], wd)


if __name__ == "__main__":
    unittest.main()
