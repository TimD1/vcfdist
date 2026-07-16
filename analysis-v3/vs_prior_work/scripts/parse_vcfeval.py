"""Parse vcfeval classified VCFs into unified counts and a QUAL-swept PR curve."""
import os
from parse_common import open_text, size_class, allele_count, compute_pr_curve


def _iter_records(path, sv_threshold=50):
    """Yield (size_class, qual, allele_ct) for each usable record; qual None if '.'."""
    if not os.path.exists(path):
        return
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt, qual_s = f[3], f[4], f[5]
            for a in alt.split(","):
                sc = size_class(ref, a, sv_threshold)
                if sc is None:
                    continue
                gt = f[9].split(":")[0] if len(f) > 9 else "0/1"
                qual = None if qual_s in (".", "") else float(qual_s)
                yield sc, qual, allele_count(gt)


def _find(eval_dir, base):
    for ext in (".vcf", ".vcf.gz"):
        p = os.path.join(eval_dir, base + ext)
        if os.path.exists(p):
            return p
    return os.path.join(eval_dir, base + ".vcf")


def parse_dir(eval_dir, sv_threshold=50):
    """Return (counts_rows, curve_rows) for a vcfeval output directory."""
    tp = list(_iter_records(_find(eval_dir, "tp"), sv_threshold))
    fp = list(_iter_records(_find(eval_dir, "fp"), sv_threshold))
    tp_base = list(_iter_records(_find(eval_dir, "tp-baseline"), sv_threshold))
    fn = list(_iter_records(_find(eval_dir, "fn"), sv_threshold))

    counts = {}
    def bump(key, sc, weight):
        counts.setdefault(sc, {"size_class": sc, "tp_query": 0, "tp_truth": 0, "fp": 0, "fn": 0})
        counts[sc][key] += weight
    for sc, _q, w in tp:      bump("tp_query", sc, w)
    for sc, _q, w in fp:      bump("fp", sc, w)
    for sc, _q, w in tp_base: bump("tp_truth", sc, w)
    for sc, _q, w in fn:      bump("fn", sc, w)

    # PR curve (one record per call; missing QUAL treated as 0)
    query_records = ([(sc, q or 0.0, True) for sc, q, _w in tp]
                     + [(sc, q or 0.0, False) for sc, q, _w in fp])
    truth_totals = {}
    for sc, _q, _w in tp_base + fn:
        truth_totals[sc] = truth_totals.get(sc, 0) + 1
    curve = compute_pr_curve(query_records, truth_totals)
    return list(counts.values()), curve
