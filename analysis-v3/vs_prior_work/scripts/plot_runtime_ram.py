"""Plot wall-clock runtime and peak RAM per tool x dataset."""
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


def load_runtime(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["tool"], r["dataset"])] = {
                "wall_seconds": float(r["wall_seconds"]),
                "max_rss_kb": int(r["max_rss_kb"]),
            }
    return out


def _panel(ax, data, datasets, value_fn, ylabel):
    x = np.arange(len(datasets))
    w = 0.8 / len(TOOLS)
    for ti, tool in enumerate(TOOLS):
        vals = [value_fn(data.get((tool, ds))) if data.get((tool, ds)) else 0.0
                for ds in datasets]
        ax.bar(x + (ti - (len(TOOLS) - 1) / 2) * w, vals, w,
               color=COLORS[tool], label=TOOL_LABELS[tool])
    ax.set_xticks(x); ax.set_xticklabels(datasets, fontsize=7, rotation=20)
    ax.set_ylabel(ylabel, fontsize=8)


def main(runtime_tsv, out_pdf):
    data = load_runtime(runtime_tsv)
    datasets = sorted({k[1] for k in data})
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
    _panel(ax[0], data, datasets, lambda r: r["wall_seconds"], "Wall-clock time (s)")
    _panel(ax[1], data, datasets, lambda r: r["max_rss_kb"] / 1e6, "Peak RAM (GB)")
    ax[0].legend(fontsize=6, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runtime", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    main(a.runtime, a.out)
