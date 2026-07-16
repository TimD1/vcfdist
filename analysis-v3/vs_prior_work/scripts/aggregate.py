"""Per-(tool,dataset) aggregation: dispatch to parsers, write unified TSV shards."""
import argparse
import csv

import parse_vcfdist
import parse_vcfeval
import parse_happy


def parse_time_log(path):
    """Return (wall_seconds, max_rss_kb) from a GNU `/usr/bin/time -v` log."""
    wall, rss = None, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("Elapsed (wall clock) time"):
                t = line.rsplit(": ", 1)[1]
                parts = t.split(":")
                parts = [float(x) for x in parts]
                if len(parts) == 3:
                    wall = parts[0] * 3600 + parts[1] * 60 + parts[2]
                else:
                    wall = parts[0] * 60 + parts[1]
            elif line.startswith("Maximum resident set size"):
                rss = int(line.rsplit(": ", 1)[1])
    return wall, rss


def write_tsv(path, rows, columns):
    """Write rows (list of dicts) to a TSV with the given column order."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _dispatch(tool, inputs):
    if tool in ("vcfdist-v3", "vcfdist-v2"):
        counts = parse_vcfdist.parse_counts(inputs["summary"])
        curve = parse_vcfdist.parse_curve(inputs["curve"])
        return counts, curve
    if tool == "vcfeval":
        return parse_vcfeval.parse_dir(inputs["dir"])
    if tool == "happy":
        return parse_happy.parse_vcf(inputs["vcf"])
    raise SystemExit(f"unknown tool: {tool}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--summary")
    ap.add_argument("--curve")
    ap.add_argument("--dir")
    ap.add_argument("--vcf")
    ap.add_argument("--time-log", required=True)
    ap.add_argument("--counts-out", required=True)
    ap.add_argument("--curve-out", required=True)
    ap.add_argument("--runtime-out", required=True)
    a = ap.parse_args(argv)

    counts, curve = _dispatch(a.tool, {"summary": a.summary, "curve": a.curve,
                                       "dir": a.dir, "vcf": a.vcf})
    for r in counts:
        r["tool"] = a.tool; r["dataset"] = a.dataset
    for r in curve:
        r["tool"] = a.tool; r["dataset"] = a.dataset

    write_tsv(a.counts_out, counts,
              ["tool", "dataset", "size_class", "tp_query", "tp_truth", "fp", "fn"])
    write_tsv(a.curve_out, curve,
              ["tool", "dataset", "size_class", "min_qual", "precision", "recall"])

    wall, rss = parse_time_log(a.time_log)
    write_tsv(a.runtime_out,
              [{"tool": a.tool, "dataset": a.dataset, "wall_seconds": wall, "max_rss_kb": rss}],
              ["tool", "dataset", "wall_seconds", "max_rss_kb"])


if __name__ == "__main__":
    main()
