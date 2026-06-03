# `filter-adaptive` test fixtures

A tiny, fully deterministic dataset for exercising metaMATE's `filter-adaptive`
mode end-to-end, in both **ASV mode** and **OTU mode** and with both
`--criteria` settings, *without* needing BBMap / MAFFT / R for the validation
step.

All files are produced by [`make_fixtures.py`](make_fixtures.py)
(`python make_fixtures.py`). Edit that script and re-run it to change the
fixtures; do not hand-edit the data files.

## The twelve ASVs

The dataset is built so every classification outcome is represented, including
the `refpass ∩ nontarget` overlap that caused the FASTA/table mismatch bug
(see [`../../../metamate_filter_adaptive_bug_report.md`](../../../metamate_filter_adaptive_bug_report.md)):

| ASV | length | classification | control rows |
|---|---|---|---|
| `ref1`, `ref2` | 60 bp | authentic (reference match) | `refpass` |
| `both1` | 63 bp | **authentic *and* a length outlier (overlap)** | `refpass` + `lengthfail` |
| `short1` | 30 bp | non-authentic (too short) | `lengthfail` |
| `stop1` | 60 bp | non-authentic (in-frame stop) | `stopfail` |
| `noise1`…`noise7` | 60 bp | unclassified (threshold-only) | — |

The seven `noise` ASVs give samples SA/SB at least 10 ASVs each, which is the
minimum the `estimated_removed` criteria needs to compute a threshold.

## Files

| File | Purpose |
|---|---|
| `asvs.fasta` | the 12 input ASVs (`-A`) |
| `readmap.tsv` | per-sample read counts, samples as rows, ASVs as columns (`-M`) |
| `references.fasta` | reference sequences matching `ref1`, `ref2`, `both1` (`-R`) |
| `asvs_control.txt` | **pre-computed validation** (`lengthfail`/`stopfail`/`refpass`); copy into the output dir to let metaMATE *resume* validation instead of running BBMap/translation |
| `asv2otu.uc` | ASV→OTU centroid map for OTU mode (`--uc`) |
| `otus.fasta` | OTU centroid sequences (`--otu_fasta`) |
| `otu_table.csv` | per-sample OTU counts, OTUs as rows (`--otu_table`) |

## Running without BBMap/MAFFT/R (resume path)

`filter-adaptive` calls the same validation step as `find`. If a
`<asv-basename>_control.txt` file already exists in the output directory (and
`--overwrite` is not given), metaMATE *resumes* from it and skips reference
matching and length/translation checks. The only external tools still needed
are the ones `metamate` always probes for at startup (`mafft`, `Rscript`,
`bbmap.sh` must be on `PATH`, but they are not actually invoked on this path).

### ASV mode

```bash
OUT=./out_asv
mkdir -p "$OUT"
cp asvs_control.txt "$OUT/asvs_control.txt"      # enables resume (no BBMap)
metamate filter-adaptive \
  -A asvs.fasta -M readmap.tsv -R references.fasta \
  --expectedlength 60 --percentvar 0 -s 5 \
  --percentile 0.5 \
  -o "$OUT"
```

Swap `--criteria estimated_removed` in to use the estimation-based threshold
(see below). For OTU mode add `--uc asv2otu.uc --otu_fasta otus.fasta
--otu_table otu_table.csv` and drop `-M readmap.tsv` (counts come from the OTU
table).

To run the *full* pipeline instead (recomputing validation with BBMap etc.),
omit the `cp` step (or pass `--overwrite`) and ensure BBMap/MAFFT/R are
installed.

## Expected results (`--percentile 0.5`, validation enforced — the default)

Per-sample thresholds come from the **non-authentic** read counts present in
that sample; a sample with no non-authentic ASVs inherits the mean of the other
samples' thresholds (the *fallback*). SC has no non-authentic ASVs in any of
these runs, so it always uses the fallback.

### ASV mode, `--criteria verified_removed` (default)

Threshold = the requested percentile of the verified non-authentic counts.
Thresholds: `SA=5.0`, `SB=6.0`, `SC=5.5` (fallback).

Retained FASTA / table columns: **`ref1 ref2 noise1 noise2 noise3 noise4 noise5 noise6 noise7`**

This demonstrates, in one run:
* authentic ASVs are **kept even below threshold** (`ref2`=2 in SA, `ref1`=2 in SB);
* non-authentic ASVs are **removed even above threshold** (`both1`=40 in SA, `short1`=8 in SB);
* the overlap `both1` (refpass ∩ nontarget) is **removed and not re-added** — the regression-tested bug;
* unclassified ASVs are kept only if they pass the per-sample threshold (e.g. `noise2`=1 dropped in SA);
* the FASTA and the abundance table contain exactly the same IDs.

`filter-adaptive_summary.csv` records pre/post counts per sample.
`Verified_Authentic_Post` counts only authentic ASVs that passed the threshold
*on their own* (rescued ones are not counted), so for SA it is 2 of 3.

### ASV mode, `--criteria estimated_removed`

Threshold = the observed count whose *estimated* fraction of true non-authentics
removed is closest to the percentile (see the top-level README). Thresholds:
`SA=5`, `SB=2`, `SC=3.5` (fallback = mean of 5 and 2). SB's lower threshold
keeps more unclassified ASVs than `verified_removed` did. The enforced
authentic/non-authentic handling is identical, so `both1`/`short1`/`stop1` are
still removed.

### OTU mode, `--criteria verified_removed`

OTUs are classified from their member ASVs (any authentic member ⇒ Authentic;
all members non-authentic ⇒ Non-Authentic; else Unclassified):

| OTU (centroid) | members | status |
|---|---|---|
| `ref1` | ref1, ref2 | Authentic |
| `both1` | both1, short1 | Authentic (contains authentic `both1`) |
| `stop1` | stop1 | Non-Authentic |
| `noise1` | noise1…noise7 | Unclassified |

Thresholds: `SA=3.0`, `SB=6.0`, `SC=4.5` (fallback).

Retained FASTA / table columns: **`ref1 both1 noise1`** — the Non-Authentic OTU
`stop1` is removed, the Unclassified OTU `noise1` is kept by threshold, and the
Authentic OTUs are kept. `otu_summary.csv` lists the per-ASV/per-OTU statuses.
