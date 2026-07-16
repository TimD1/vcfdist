"""Shared helpers for parsing tool outputs into unified benchmarking tables."""
import gzip

SNP = "SNP"
INDEL = "INDEL"
SV = "SV"


def open_text(path):
    """Open a text or gzip file for reading, transparently by extension."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def size_class(ref, alt, sv_threshold=50):
    """Classify a REF/ALT pair as SNP/INDEL/SV, or None if not evaluated.

    None for symbolic/missing alts and for same-length non-SNP records.
    """
    if alt in (".", "*") or alt.startswith("<") or "[" in alt or "]" in alt:
        return None
    if len(ref) == 1 and len(alt) == 1:
        return SNP
    diff = abs(len(ref) - len(alt))
    if diff == 0:
        return None
    return INDEL if diff < sv_threshold else SV


def allele_count(gt):
    """Count non-reference, non-missing alleles in a genotype string."""
    total = 0
    for a in gt.replace("/", "|").split("|"):
        if a not in (".", "0", ""):
            total += 1
    return total


def compute_pr_curve(query_records, truth_totals):
    """Sweep QUAL thresholds to produce PR points per size class.

    query_records: iterable of (size_class, qual, is_tp) for query-side calls.
    truth_totals: {size_class: total truth count in class}.
    Returns list of {"size_class","min_qual","precision","recall"}.
    """
    by_class = {}
    for sc, qual, is_tp in query_records:
        by_class.setdefault(sc, []).append((float(qual), bool(is_tp)))

    rows = []
    for sc, recs in by_class.items():
        truth_total = truth_totals.get(sc, 0)
        thresholds = sorted({0.0} | {q for q, _ in recs})
        for t in thresholds:
            tpq = sum(1 for q, tp in recs if tp and q >= t)
            fp = sum(1 for q, tp in recs if not tp and q >= t)
            precision = tpq / (tpq + fp) if (tpq + fp) else 1.0
            recall = tpq / truth_total if truth_total else 0.0
            rows.append({"size_class": sc, "min_qual": t,
                         "precision": precision, "recall": recall})
    return rows
