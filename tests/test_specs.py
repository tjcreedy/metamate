#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the filtering-specification parser in ``metamate.core``
(``find`` and ``dump`` modes). Expectations are taken from the worked examples
in the top-level README ("Specifications" section).

Independent of the production scripts; run with::

    python -m unittest discover -s tests -p 'test_*.py' -v
"""

import os
import sys
import types
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from metamate import core


class ResolveRangesTests(unittest.TestCase):
    """``resolve_ranges`` expands a threshold string into a flat list of floats.
    ``start-stop/nsteps`` is a linspace; comma-separated parts are concatenated."""

    def test_single_value(self):
        self.assertEqual(core.resolve_ranges("5"), [5.0])

    def test_range_is_linspace(self):
        # README: "1-2/5" -> 1, 1.25, 1.5, 1.75, 2
        self.assertEqual(core.resolve_ranges("1-2/5"), [1.0, 1.25, 1.5, 1.75, 2.0])

    def test_mixed_values_and_ranges(self):
        # README: "1,2,3,4-10/4" -> 1, 2, 3, 4, 6, 8, 10
        self.assertEqual(core.resolve_ranges("1,2,3,4-10/4"),
                         [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0])

    def test_integer_range(self):
        # README: "2-10/5" -> 2, 4, 6, 8, 10
        self.assertEqual(core.resolve_ranges("2-10/5"), [2.0, 4.0, 6.0, 8.0, 10.0])


class ResolveSpecTests(unittest.TestCase):
    """``resolve_spec`` parses one ``[category(/ies); metric; thresholds]`` term.
    (Called after spaces have been stripped, as ``parse_specs`` does.)"""

    def test_simple_term(self):
        name, cats, metric, thresholds = core.resolve_spec("[total;n;2-10/5]")
        self.assertEqual(name, "total_n")
        self.assertEqual(cats, ["total"])
        self.assertEqual(metric, "n")
        self.assertEqual(thresholds, [2.0, 4.0, 6.0, 8.0, 10.0])

    def test_partitioned_term(self):
        name, cats, metric, thresholds = core.resolve_spec("[library|taxon;p;0.4,0.6]")
        self.assertEqual(name, "library|taxon_p")
        self.assertEqual(cats, ["library", "taxon"])
        self.assertEqual(metric, "p")
        self.assertEqual(thresholds, [0.4, 0.6])

    def test_bad_metric_errors(self):
        with self.assertRaises(SystemExit):
            core.resolve_spec("[total;x;2]")

    def test_bad_primary_category_errors(self):
        # First category must be 'library' or 'total'.
        with self.assertRaises(SystemExit):
            core.resolve_spec("[clade;n;2]")

    def test_bad_secondary_category_errors(self):
        # Secondary categories must be 'clade' or 'taxon'.
        with self.assertRaises(SystemExit):
            core.resolve_spec("[total|library;n;2]")

    def test_wrong_number_of_parts_errors(self):
        with self.assertRaises(SystemExit):
            core.resolve_spec("[total;n]")


def _args(mode, specification, taxgroups=None):
    return types.SimpleNamespace(mode=mode, specification=specification, taxgroups=taxgroups)


class ParseSpecsTests(unittest.TestCase):

    def _write_spec(self, text):
        path = os.path.join(tempfile.mkdtemp(), "spec.txt")
        with open(path, "w") as f:
            f.write(text)
        return path

    def test_dump_single_term(self):
        specs, termset, nterm, nthresh, thresholds = core.parse_specs(
            _args("dump", ["[total;n;5]"]))
        self.assertEqual(specs["name"], ["total_n"])
        self.assertEqual(nterm, 1)
        self.assertEqual(nthresh, 1)
        self.assertEqual(list(thresholds), [(0, ["total_n"], 5.0)])

    def test_dump_rejects_multiple_thresholds(self):
        # dump allows only a single threshold set.
        with self.assertRaises(SystemExit):
            core.parse_specs(_args("dump", ["[total;n;1,2]"]))

    def test_find_additive_terms_count(self):
        # 3 + 2 = 5 threshold sets across two additive terms.
        path = self._write_spec("[total;n;1,2,3] + [library;p;0.1,0.2]")
        specs, termset, nterm, nthresh, thresholds = core.parse_specs(_args("find", path))
        self.assertEqual(nthresh, 5)
        self.assertEqual(nterm, 2)               # 'total' and 'library'
        self.assertEqual(specs["name"], ["total_n", "library_p"])

    def test_find_multiplicative_terms_count(self):
        # 1 * 2 = 2 combined threshold sets.
        path = self._write_spec("[library;p;0.1] * [total;n;2,5]")
        _, _, _, nthresh, _ = core.parse_specs(_args("find", path))
        self.assertEqual(nthresh, 2)

    def test_find_ignores_comments_and_blank_lines(self):
        path = self._write_spec("# a comment\n\n[total;n;2-10/5]\n")
        _, _, _, nthresh, _ = core.parse_specs(_args("find", path))
        self.assertEqual(nthresh, 5)

    def test_filter_adaptive_returns_empty_specs(self):
        specs, termset, nterm, nthresh, thresholds = core.parse_specs(
            _args("filter-adaptive", None))
        self.assertEqual(specs, {"name": []})
        self.assertEqual(nthresh, 0)
        self.assertEqual(list(thresholds), [])


if __name__ == "__main__":
    unittest.main()
