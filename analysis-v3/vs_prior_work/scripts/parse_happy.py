"""Parse hap.py annotated VCF into unified counts and a QUAL-swept PR curve."""
from parse_common import open_text, size_class, compute_pr_curve


def _fmt_field(fmt_keys, sample_field, key):
    values = sample_field.split(":")
    for k, v in zip(fmt_keys, values):
        if k == key:
            return v
    return "."


def parse_vcf(happy_vcf_path, sv_threshold=50):
    """Return (counts_rows, curve_rows) from a hap.py two-sample annotated VCF."""
    counts = {}
    query_records = []
    truth_totals = {}

    def bump(sc, key, n=1):
        counts.setdefault(sc, {"size_class": sc, "tp_query": 0, "tp_truth": 0, "fp": 0, "fn": 0})
        counts[sc][key] += n

    with open_text(happy_vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt, qual_s, fmt = f[3], f[4], f[5], f[8].split(":")
            truth_bd = _fmt_field(fmt, f[9], "BD")
            query_bd = _fmt_field(fmt, f[10], "BD") if len(f) > 10 else "."
            query_qq = _fmt_field(fmt, f[10], "QQ") if len(f) > 10 else "."
            # hap.py carries ROC quality in QQ (QUAL is often '.'); prefer QQ, else QUAL.
            qual_src = query_qq if query_qq not in (".", "") else qual_s
            qual = 0.0 if qual_src in (".", "") else float(qual_src)
            for a in alt.split(","):
                sc = size_class(ref, a, sv_threshold)
                if sc is None:
                    continue
                if query_bd == "TP":
                    bump(sc, "tp_query"); query_records.append((sc, qual, True))
                elif query_bd == "FP":
                    bump(sc, "fp"); query_records.append((sc, qual, False))
                if truth_bd == "TP":
                    bump(sc, "tp_truth"); truth_totals[sc] = truth_totals.get(sc, 0) + 1
                elif truth_bd == "FN":
                    bump(sc, "fn"); truth_totals[sc] = truth_totals.get(sc, 0) + 1

    curve = compute_pr_curve(query_records, truth_totals)
    return list(counts.values()), curve
