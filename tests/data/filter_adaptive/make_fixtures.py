#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the tiny, deterministic fixture set used to exercise metaMATE's
``filter-adaptive`` mode end-to-end (both ASV mode and OTU mode), without
needing BBMap/MAFFT/R.

Run from this directory::

    python make_fixtures.py

It (re)writes the data files described in this directory's README.md. The ASV
classification is *fixed* by the supplied control file (``asvs_control.txt``),
which lets ``filter-adaptive`` resume validation instead of recomputing it with
external tools. The eight ASVs are chosen to cover every behaviour worth
testing:

    ref1, ref2  -> authentic (reference match)           [refpass]
    both1       -> authentic AND a length outlier         [refpass + lengthfail]
                   (the refpass n nontarget *overlap*)
    short1      -> non-authentic, too short               [lengthfail]
    stop1       -> non-authentic, in-frame stop codon      [stopfail]
    noise1-3    -> unclassified (subject only to thresholds)

Expected length is 60 bp (percentvar 0): ref/noise/stop sequences are 60 bp,
``short1`` is 30 bp, and ``both1`` is 63 bp so a full run would also flag it as
a length outlier, matching the control file.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _w(name, text):
    with open(os.path.join(HERE, name), "w", newline="\n") as f:
        f.write(text)


# --- Sequences (60 bp unless noted) ---------------------------------------
# Seven "noise" (unclassified) ASVs give SA/SB >= 10 present ASVs, which is the
# minimum the 'estimated_removed' criteria needs to engage (it otherwise falls
# back). SC keeps the non-authentics absent so it always uses the fallback.
def _seq60(i):
    # 60 bp, distinct per i and distinct from the period-4 ref sequences so a
    # full (non-resume) run would not accidentally match noise to a reference.
    prefix = ["AAAA", "CCCC", "GGGG", "TTTT", "ACAC", "GTGT", "TATA"][i - 1]
    return (prefix + "TTGCA" * 12)[:60]

SEQS = {
    "ref1":   "ACGTACGTACGTACGTACGT" * 3,                 # 60
    "ref2":   "ACGTACGTACGTACGTACGA" * 3,                 # 60
    "stop1":  "TAAACGTACGTACGTACGTA" * 3,                 # 60, contains TAA
    "short1": "ACGTACGTACGTACGTACGTACGTACGTAC",           # 30 (length fail)
    "both1":  "ACGTACGTACGTACGTACGT" * 3 + "TGA",         # 63 (length fail) + ref match
}
for _i in range(1, 8):
    SEQS[f"noise{_i}"] = _seq60(_i)

# References that ASVs ref1, ref2 and both1 match (identical sequences).
REFS = {"refA": SEQS["ref1"], "refB": SEQS["ref2"], "refC": SEQS["both1"]}

NOISE = [f"noise{i}" for i in range(1, 8)]
ASV_ORDER = ["ref1", "ref2", "both1", "short1", "stop1"] + NOISE

# --- Per-sample, per-ASV read counts (samples as rows, ASVs as columns) ---
# Chosen so per-sample thresholds and the final retained set are hand-checkable
# (see README.md). SC has no non-authentic ASVs, exercising the fallback.
READMAP = {
    #     ref1 ref2 both1 short1 stop1  n1   n2  n3  n4  n5  n6  n7
    "SA": [50, 2, 40, 5, 3, 100, 1, 8, 12, 2, 6, 4],
    "SB": [2, 30, 1, 8, 6, 4, 60, 2, 7, 11, 3, 9],
    "SC": [20, 15, 0, 0, 0, 2, 9, 30, 5, 1, 7, 2],
}

# --- OTU mapping: ASV -> OTU centroid -------------------------------------
# otuA (ref1):  ref1, ref2          -> Authentic
# otuB (both1): both1, short1       -> Authentic (both1 is authentic)
# otuC (stop1): stop1               -> Non-Authentic
# otuD (noise1):noise1, noise2, noise3 -> Unclassified
OTU_MEMBERS = {
    "ref1":  ["ref1", "ref2"],
    "both1": ["both1", "short1"],
    "stop1": ["stop1"],
    "noise1": NOISE,
}
OTU_TABLE = {
    #          SA   SB   SC
    "ref1":  [52, 32, 35],
    "both1": [45, 2, 0],
    "stop1": [3, 6, 0],
    "noise1": [109, 66, 34],
}


def main():
    # asvs.fasta
    _w("asvs.fasta", "".join(f">{n}\n{SEQS[n]}\n" for n in ASV_ORDER))

    # references.fasta
    _w("references.fasta", "".join(f">{n}\n{s}\n" for n, s in REFS.items()))

    # readmap.tsv (samples as rows, ASVs as columns)
    lines = ["Sample\t" + "\t".join(ASV_ORDER)]
    for s, row in READMAP.items():
        lines.append(s + "\t" + "\t".join(str(v) for v in row))
    _w("readmap.tsv", "\n".join(lines) + "\n")

    # asvs_control.txt (fixed validation; copy into the output dir to resume)
    ctrl = []
    for a in ("short1", "both1"):
        ctrl.append(f"lengthfail\t{a}")
    ctrl.append("stopfail\tstop1")
    for a in ("ref1", "ref2", "both1"):
        ctrl.append(f"refpass\t{a}")
    _w("asvs_control.txt", "\n".join(ctrl) + "\n")

    # asv2otu.uc (vsearch-style: S = seed/centroid, H = hit/member)
    uc = []
    for centroid, members in OTU_MEMBERS.items():
        uc.append(f"S\t0\t{len(SEQS[centroid])}\t*\t+\t*\t*\t*\t{centroid}\t*")
        for m in members:
            if m == centroid:
                continue
            uc.append(f"H\t0\t{len(SEQS[m])}\t99.0\t+\t0\t0\t*\t{m}\t{centroid}")
    _w("asv2otu.uc", "\n".join(uc) + "\n")

    # otus.fasta (centroid sequences)
    _w("otus.fasta", "".join(f">{c}\n{SEQS[c]}\n" for c in OTU_MEMBERS))

    # otu_table.csv (OTUs as rows, samples as columns)
    rows = ["OTU,SA,SB,SC"]
    for otu, row in OTU_TABLE.items():
        rows.append(otu + "," + ",".join(str(v) for v in row))
    _w("otu_table.csv", "\n".join(rows) + "\n")

    print("Wrote fixtures to", HERE)


if __name__ == "__main__":
    main()
