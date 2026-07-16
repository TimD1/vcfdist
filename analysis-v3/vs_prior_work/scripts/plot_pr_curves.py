"""Overlay precision-recall curves for all four tools, per size class x dataset."""
import argparse
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

TOOLS = ["vcfdist-v3", "vcfdist-v2", "vcfeval"]  # "happy" disabled for now
TOOL_LABELS = {"vcfdist-v3": "vcfdist v3", "vcfdist-v2": "vcfdist v2.6.4",
               "vcfeval": "vcfeval", "happy": "hap.py (xcmp)"}
COLORS = {"vcfdist-v3": "#2e7d32", "vcfdist-v2": "#73c375",
          "vcfeval": "#fa6949", "happy": "#6aadd5"}
SIZES = ["SNP", "INDEL", "SV"]
FOOTNOTE = ("All curves are QUAL-swept (per MIN_QUAL); counts are hard TP/FP/FN — vcfdist "
            "thresholds credit at --credit-threshold, vcfeval/hap.py use their own matching. "
            "Axes zoomed per panel to the data range; dots mark plotted points.")


def load_curve(path):
    tmp = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tmp[(r["tool"], r["dataset"], r["size_class"])].append(
                (float(r["recall"]), float(r["precision"])))
    return {k: sorted(v) for k, v in tmp.items()}


def main(curve_tsv, out_pdf):
    data = load_curve(curve_tsv)
    datasets = sorted({k[1] for k in data})
    fig, ax = plt.subplots(len(SIZES), len(datasets),
                           figsize=(3.2 * len(datasets), 3 * len(SIZES)),
                           squeeze=False)
    for si, size in enumerate(SIZES):
        for di, ds in enumerate(datasets):
            a = ax[si][di]
            allpts = []
            for tool in TOOLS:
                pts = data.get((tool, ds, size))
                if not pts:
                    continue
                allpts.extend(pts)
                rec = [p[0] for p in pts]; prec = [p[1] for p in pts]
                a.plot(rec, prec, marker="o", markersize=3, linestyle="-",
                       color=COLORS[tool], label=TOOL_LABELS[tool], linewidth=1.0)
            if si == 0:
                a.set_title(ds, fontsize=9)
            if di == 0:
                a.set_ylabel(f"{size}\nPrecision", fontsize=8)
            if si == len(SIZES) - 1:
                a.set_xlabel("Recall", fontsize=8)
            # Dynamic upper-right zoom: bound the plotted points, pad below, cap near 1.
            if allpts:
                rmin = min(p[0] for p in allpts)
                pmin = min(p[1] for p in allpts)
                xr = max(1.0 - rmin, 1e-3)
                yr = max(1.0 - pmin, 1e-3)
                a.set_xlim(rmin - 0.05 * xr, 1.0 + 0.02 * xr)
                a.set_ylim(pmin - 0.05 * yr, 1.0 + 0.02 * yr)
            a.grid(True, alpha=0.3)
    handles = [Line2D([0], [0], color=COLORS[t], lw=1.2, label=TOOL_LABELS[t]) for t in TOOLS]
    ax[0][0].legend(handles=handles, fontsize=6, loc="lower left")
    fig.text(0.5, 0.005, FOOTNOTE, ha="center", fontsize=6, style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_pdf)   # format inferred from extension (.pdf / .png)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--curve", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    main(a.curve, a.out)
