import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import parse_common as pc

def test_size_class_snp():
    assert pc.size_class("A", "T") == pc.SNP

def test_size_class_indel_below_threshold():
    assert pc.size_class("A", "ACGTG") == pc.INDEL  # size-diff 4

def test_size_class_sv_at_threshold():
    assert pc.size_class("A", "A" + "C" * 50) == pc.SV  # size-diff 50

def test_size_class_same_length_nonsnp_skipped():
    assert pc.size_class("AC", "GT") is None  # size-diff 0

def test_size_class_star_allele_skipped():
    assert pc.size_class("A", "*") is None

def test_allele_count():
    assert pc.allele_count("0|1") == 1
    assert pc.allele_count("1/1") == 2
    assert pc.allele_count("0/0") == 0
    assert pc.allele_count(".|1") == 1

def test_compute_pr_curve_monotone_thresholds():
    # 2 TP (qual 10, 20), 1 FP (qual 5); 4 truth SNPs total
    recs = [(pc.SNP, 10.0, True), (pc.SNP, 20.0, True), (pc.SNP, 5.0, False)]
    rows = pc.compute_pr_curve(recs, {pc.SNP: 4})
    by_q = {r["min_qual"]: r for r in rows}
    # at threshold 0: TPq=2, FP=1 -> prec 2/3, recall 2/4
    assert abs(by_q[0.0]["precision"] - 2 / 3) < 1e-9
    assert abs(by_q[0.0]["recall"] - 0.5) < 1e-9
    # at threshold 20: TPq=1, FP=0 -> prec 1.0, recall 1/4
    assert abs(by_q[20.0]["precision"] - 1.0) < 1e-9
    assert abs(by_q[20.0]["recall"] - 0.25) < 1e-9

def test_open_text_plain(tmp_path):
    p = tmp_path / "a.txt"; p.write_text("hello\nworld\n")
    with pc.open_text(str(p)) as fh:
        assert fh.read() == "hello\nworld\n"

def test_open_text_gzip(tmp_path):
    import gzip
    p = tmp_path / "a.txt.gz"
    with gzip.open(str(p), "wt") as fh:
        fh.write("gz-content\n")
    with pc.open_text(str(p)) as fh:
        assert fh.read() == "gz-content\n"

def test_compute_pr_curve_zero_truth_total_recall_zero():
    # SNP present in query records but absent from truth_totals -> recall 0.0
    recs = [(pc.SNP, 10.0, True), (pc.SNP, 5.0, False)]
    rows = pc.compute_pr_curve(recs, {})
    assert rows and all(r["recall"] == 0.0 for r in rows)
    by_q = {r["min_qual"]: r for r in rows}
    assert abs(by_q[0.0]["precision"] - 0.5) < 1e-9
