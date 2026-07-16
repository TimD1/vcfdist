"""Plot FNR and FDR bar charts (SNP/INDEL/SV panels) for all four tools."""
import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TOOLS = ["vcfdist-v3", "vcfdist-v2", "vcfeval"]  # "happy" disabled for now
TOOL_LABELS = {"vcfdist-v3": "vcfdist v3", "vcfdist-v2": "vcfdist v2.6.4",
               "vcfeval": "vcfeval", "happy": "hap.py (xcmp)"}
COLORS = {"vcfdist-v3": "#2e7d32", "vcfdist-v2": "#73c375",
          "vcfeval": "#fa6949", "happy": "#6aadd5"}
SIZES = ["SNP", "INDEL", "SV"]


def load_counts(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["tool"], r["dataset"], r["size_class"])
            out[key] = {k: int(r[k]) for k in ("tp_query", "tp_truth", "fp", "fn")}
    return out


def fnr(row):
    denom = row["fn"] + row["tp_truth"]
    return row["fn"] / denom if denom else 0.0


def fdr(row):
    denom = row["fp"] + row["tp_query"]
    return row["fp"] / denom if denom else 0.0


def _datasets(counts):
    return sorted({k[1] for k in counts})


def _plot(counts, metric_fn, ylabel, out_pdf):
    datasets = _datasets(counts)
    fig, ax = plt.subplots(1, len(SIZES), figsize=(9, 3), squeeze=False)
    x = np.arange(len(datasets))
    w = 0.8 / len(TOOLS)
    for si, size in enumerate(SIZES):
        a = ax[0][si]
        for ti, tool in enumerate(TOOLS):
            vals = []
            for ds in datasets:
                row = counts.get((tool, ds, size))
                vals.append(metric_fn(row) if row else 0.0)
            a.bar(x + (ti - (len(TOOLS) - 1) / 2) * w, vals, w,
                  color=COLORS[tool], label=TOOL_LABELS[tool])
        a.set_title(f"{size}", fontsize=9)
        a.set_xticks(x); a.set_xticklabels(datasets, fontsize=7, rotation=20)
        a.set_ylabel(ylabel if si == 0 else "", fontsize=8)
    ax[0][0].legend(fontsize=6, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)


def main(counts_tsv, out_fnr_pdf, out_fdr_pdf):
    counts = load_counts(counts_tsv)
    _plot(counts, fnr, "False Negative Rate (FN / (FN+TP))", out_fnr_pdf)
    _plot(counts, fdr, "False Discovery Rate (FP / (FP+TP))", out_fdr_pdf)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--fnr-out", required=True)
    ap.add_argument("--fdr-out", required=True)
    a = ap.parse_args()
    main(a.counts, a.fnr_out, a.fdr_out)
