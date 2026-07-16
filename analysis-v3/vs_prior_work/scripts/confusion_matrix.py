"""Confusion matrix: vcfdist-v3's benchmark decision vs another caller's.

For every evaluated variant, keyed by standard-VCF (CHROM, POS, REF, ALT), we take the
benchmark decision made by the reference tool (vcfdist-v3) and by one other caller, then
cross-tabulate over labels {TP, FP, FN, N}. N = the variant is absent from that tool's
evaluated set (a representation difference — e.g. a repeat indel the two tools left-align
differently). The diagonal is agreement; off-diagonal cells are variants the two tools
classified differently. One matrix is produced per other caller.

Labels are read in a common coordinate system so the two callers join:
- vcfdist: from summary.vcf (standard VCF). Per record, label = TP if the QUERY sample's
  BD is TP, else FP if QUERY BD is FP, else FN if the TRUTH sample's BD is FN.
- vcfeval: from its classified VCFs — tp.vcf (TP) / fp.vcf (FP) on the query side, fn.vcf
  (FN) on the truth side.
"""
import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from parse_common import open_text, size_class

LABELS = ["TP", "FP", "FN", "N"]
REF_TOOL = "vcfdist-v3"


def load_vcfdist(prefix):
    """key (chrom, pos, ref, alt) -> label, from {prefix}summary.vcf (standard VCF)."""
    m = {}
    with open_text(prefix + "summary.vcf") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, fmt = c[0], c[1], c[3], c[4], c[8].split(":")
            if "BD" not in fmt:
                continue
            bi = fmt.index("BD")
            tvals, qvals = c[9].split(":"), c[10].split(":")
            tbd = tvals[bi] if bi < len(tvals) else "."
            qbd = qvals[bi] if bi < len(qvals) else "."
            if qbd == "TP":
                lab = "TP"
            elif qbd == "FP":
                lab = "FP"
            elif tbd == "FN":
                lab = "FN"
            else:
                continue
            for a in alts.split(","):
                if a in (".", "*"):
                    continue
                key = (chrom, pos, ref, a)
                # collapse the per-haplotype rows; a recovered (TP) call wins.
                if key not in m or lab == "TP":
                    m[key] = lab
    return m


def load_vcfeval(dirpath):
    """key (chrom, pos, ref, alt) -> label, from vcfeval's classified VCFs."""
    m = {}
    spec = (("tp.vcf", "TP"), ("fp.vcf", "FP"), ("fn.vcf", "FN"))
    for base, lab in spec:
        path = None
        for cand in (os.path.join(dirpath, base + ".gz"), os.path.join(dirpath, base)):
            if os.path.exists(cand):
                path = cand
                break
        if path is None:
            continue
        with open_text(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                chrom, pos, ref, alts = c[0], c[1], c[3], c[4]
                for a in alts.split(","):
                    if a in (".", "*"):
                        continue
                    key = (chrom, pos, ref, a)
                    if key not in m or lab == "TP":
                        m[key] = lab
    return m


def load_tool(tool, outdir, ds):
    """Dispatch to the right loader for a tool/dataset."""
    if tool.startswith("vcfdist"):
        return load_vcfdist(os.path.join(outdir, tool, f"{ds}."))
    if tool == "vcfeval":
        return load_vcfeval(os.path.join(outdir, "vcfeval", ds))
    raise SystemExit(f"unknown tool: {tool}")


def build_cm(ref_map, other_map, sv_threshold=50):
    """counts[size_class][ref_label][other_label] over the union of variant keys."""
    counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for key in set(ref_map) | set(other_map):
        _chrom, _pos, ref, alt = key
        sc = size_class(ref, alt, sv_threshold) or "OTHER"
        rl = ref_map.get(key, "N")
        ol = other_map.get(key, "N")
        rl = rl if rl in LABELS else "N"
        ol = ol if ol in LABELS else "N"
        counts[sc][rl][ol] += 1
    return counts


def _matrix(counts, size_classes=None):
    """Aggregate counts (optionally restricted to size_classes) into a 4x4 array."""
    mat = np.zeros((len(LABELS), len(LABELS)), dtype=int)
    for sc, byref in counts.items():
        if size_classes is not None and sc not in size_classes:
            continue
        for rl, byother in byref.items():
            for ol, n in byother.items():
                mat[LABELS.index(rl)][LABELS.index(ol)] += n
    return mat


def write_tsv(path, per_dataset, other_tool):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["dataset", "other_tool", "size_class",
                    f"{REF_TOOL}_label", f"{other_tool}_label", "count"])
        for ds, counts in per_dataset.items():
            for sc in sorted(counts):
                for rl in LABELS:
                    for ol in LABELS:
                        n = counts[sc][rl][ol]
                        if n:
                            w.writerow([ds, other_tool, sc, rl, ol, n])


def plot_cm(per_dataset, other_tool, out_pdf):
    datasets = sorted(per_dataset)
    fig, ax = plt.subplots(1, len(datasets), figsize=(3.4 * len(datasets), 3.4), squeeze=False)
    for di, ds in enumerate(datasets):
        a = ax[0][di]
        mat = _matrix(per_dataset[ds])
        # log-scaled color so the off-diagonal (rare) cells stay visible.
        a.imshow(np.log10(mat + 1), cmap="Blues", aspect="equal")
        for i in range(len(LABELS)):
            for j in range(len(LABELS)):
                v = mat[i][j]
                frac = v / mat.sum() if mat.sum() else 0
                a.text(j, i, f"{v}", ha="center", va="center", fontsize=8,
                       color="white" if frac > 0.25 else "black")
        a.set_xticks(range(len(LABELS))); a.set_xticklabels(LABELS, fontsize=8)
        a.set_yticks(range(len(LABELS))); a.set_yticklabels(LABELS, fontsize=8)
        a.set_title(ds, fontsize=9)
        if di == 0:
            a.set_ylabel(f"{REF_TOOL}", fontsize=8)
        a.set_xlabel(other_tool, fontsize=8)
    fig.suptitle(f"Per-variant decision: {REF_TOOL} (rows) vs {other_tool} (cols)", fontsize=10)
    fig.text(0.5, 0.01, "Per variant (standard-VCF key); N = absent from that tool "
                        "(e.g. representation difference). Diagonal = agreement.",
             ha="center", fontsize=6, style="italic")
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.savefig(out_pdf)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--other-tool", required=True)
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument("--sv-threshold", type=int, default=50)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-pdf", required=True)
    a = ap.parse_args(argv)

    per_dataset = {}
    for ds in a.datasets:
        ref_map = load_tool(REF_TOOL, a.outdir, ds)
        other_map = load_tool(a.other_tool, a.outdir, ds)
        per_dataset[ds] = build_cm(ref_map, other_map, a.sv_threshold)

    write_tsv(a.out_tsv, per_dataset, a.other_tool)
    plot_cm(per_dataset, a.other_tool, a.out_pdf)


if __name__ == "__main__":
    main()
