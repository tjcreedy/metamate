# Tests

## Unit tests

Fast, dependency-light unit tests for the pure-Python logic behind
`filter-adaptive`, OTU mode, and the validation-enforcement step shared by
`dump` and `filter-adaptive`. They only need `numpy` and `biopython` (already
required by metaMATE) — no BBMap/MAFFT/R, no network — and they build their own
fixtures in temporary directories.

Run from the repository root:

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
```

* [`test_core.py`](test_core.py) — `enforce_validation_fasta` (including the
  `refpass ∩ nontarget` overlap **regression test**), `adaptive_filter`
  (per-sample thresholds, target rescue with `--enforcevalidation` on/off, the
  no-non-authentics fallback, `estimated_removed`, summary statistics),
  `write_adaptive_csv`, and the `load_validation_from_control` /
  `load_validation_from_otu_summary` loaders.
* [`test_binning.py`](test_binning.py) — `parse_uc` and `parse_readmap`.
* [`test_specs.py`](test_specs.py) — the `find`/`dump` specification parser:
  `resolve_ranges`, `resolve_spec`, and `parse_specs` (additive/multiplicative
  term counting, malformed-spec errors), checked against the README examples.
* [`test_validation.py`](test_validation.py) — the two control-group
  determinants: `filterlength.check_length` / `resolve_length_spec` (including a
  regression test for `--codonsvariation`) and `filtertranslate.stopcount` /
  `min_stopcount` / `check_stops_multi`.
* [`test_metamate.py`](test_metamate.py) — `subset_table` (table↔FASTA
  consistency, separator detection, per-sample cell zeroing) and `getcliargs`
  required-input validation for each mode.

## End-to-end fixtures

* [`data/filter_adaptive/`](data/filter_adaptive/) — a tiny deterministic
  dataset with runnable `filter-adaptive` commands (ASV mode and OTU mode) and
  documented expected outputs, including a control file that lets the run
  *resume* validation without BBMap. See its
  [README](data/filter_adaptive/README.md).
* [`data/`](data/) — the larger metabarcoding example dataset used by the
  `find`/`dump` examples in the top-level README.
