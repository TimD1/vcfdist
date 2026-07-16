"""Parse vcfdist (v2 and v3) precision-recall outputs into unified rows."""
import csv
from parse_common import open_text

SIZE_CLASSES = {"SNP", "INDEL", "SV"}


def parse_counts(summary_tsv_path):
    """All-records (THRESHOLD==NONE) TP/FP/FN counts per size class."""
    out = []
    with open_text(summary_tsv_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["THRESHOLD"] != "NONE":
                continue
            if row["VAR_TYPE"] not in SIZE_CLASSES:
                continue
            out.append({
                "size_class": row["VAR_TYPE"],
                "tp_query": int(row["QUERY_TP"]),
                "tp_truth": int(row["TRUTH_TP"]),
                "fp": int(row["QUERY_FP"]),
                "fn": int(row["TRUTH_FN"]),
            })
    return out


def parse_curve(curve_tsv_path):
    """Full precision/recall sweep per size class from precision-recall.tsv."""
    out = []
    with open_text(curve_tsv_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["VAR_TYPE"] not in SIZE_CLASSES:
                continue
            out.append({
                "size_class": row["VAR_TYPE"],
                "min_qual": float(row["MIN_QUAL"]),
                "precision": float(row["PREC"]),
                "recall": float(row["RECALL"]),
            })
    return out
