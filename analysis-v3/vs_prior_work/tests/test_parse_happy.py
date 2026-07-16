import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import parse_happy as ph

HDR = ("##fileformat=VCFv4.2\n"
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
       '##FORMAT=<ID=BD,Number=1,Type=String,Description="decision">\n'
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n")

def _happy(tmp_path, body):
    p = tmp_path / "happy.vcf"; p.write_text(HDR + body); return p

def test_parse_vcf_counts_and_curve(tmp_path):
    # SNP TP: truth BD=TP, query BD=TP, qual 20
    # SNP FP: truth BD=N, query BD=FP, qual 5
    # SNP FN: truth BD=FN, query BD=N
    f = _happy(tmp_path,
        "chr1\t10\t.\tA\tT\t20\tPASS\t.\tGT:BD\t0/1:TP\t0/1:TP\n"
        "chr1\t20\t.\tG\tC\t5\tPASS\t.\tGT:BD\t./.:N\t0/1:FP\n"
        "chr1\t30\t.\tA\tG\t.\tPASS\t.\tGT:BD\t0/1:FN\t./.:N\n"
        "chr1\t100\t.\tA\t" + "A" + "C"*60 + "\t40\tPASS\t.\tGT:BD\t1/1:TP\t1/1:TP\n"
        "chr1\t200\t.\tA\t" + "A" + "C"*70 + "\t.\tPASS\t.\tGT:BD\t./.:N\t1/1:FP\n")
    counts, curve = ph.parse_vcf(str(f), sv_threshold=50)
    c = {r["size_class"]: r for r in counts}
    assert c["SNP"] == {"size_class": "SNP", "tp_query": 1, "tp_truth": 1, "fp": 1, "fn": 1}
    snp0 = next(r for r in curve if r["size_class"] == "SNP" and r["min_qual"] == 0.0)
    assert abs(snp0["precision"] - 0.5) < 1e-9   # 1 TP / (1 TP + 1 FP)
    assert abs(snp0["recall"] - 0.5) < 1e-9       # 1 query-TP / 2 truth
    snp20 = next(r for r in curve if r["size_class"] == "SNP" and r["min_qual"] == 20.0)
    assert abs(snp20["precision"] - 1.0) < 1e-9   # FP (qual 5) dropped
    assert c["SV"] == {"size_class": "SV", "tp_query": 1, "tp_truth": 1, "fp": 1, "fn": 0}
    sv0 = next(r for r in curve if r["size_class"] == "SV" and r["min_qual"] == 0.0)
    assert abs(sv0["precision"] - 0.5) < 1e-9   # missing-QUAL FP (-> 0.0) counted at t=0
    assert abs(sv0["recall"] - 1.0) < 1e-9
    sv40 = next(r for r in curve if r["size_class"] == "SV" and r["min_qual"] == 40.0)
    assert abs(sv40["precision"] - 1.0) < 1e-9  # missing-QUAL FP dropped at t=40

def test_qq_preferred_when_qual_missing(tmp_path):
    # hap.py leaves QUAL '.' and carries ROC quality in the QQ FORMAT field
    f = _happy(tmp_path,
        "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT:BD:QQ\t0/1:TP:30\t0/1:TP:30\n"
        "chr1\t20\t.\tG\tC\t.\tPASS\t.\tGT:BD:QQ\t./.:N:.\t0/1:FP:5\n")
    _counts, curve = ph.parse_vcf(str(f))
    quals = sorted({r["min_qual"] for r in curve if r["size_class"] == "SNP"})
    assert 30.0 in quals and 5.0 in quals            # QQ used, not collapsed to 0.0
    snp30 = next(r for r in curve if r["size_class"] == "SNP" and r["min_qual"] == 30.0)
    assert abs(snp30["precision"] - 1.0) < 1e-9       # FP (QQ=5) dropped at t=30
