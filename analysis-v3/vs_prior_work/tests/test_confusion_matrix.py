import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import matplotlib; matplotlib.use("Agg")
import confusion_matrix as cm

SUMMARY_HDR = ("##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
    '##FORMAT=<ID=BD,Number=1,Type=String,Description="decision">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n")
VCF_HDR = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")


def test_load_vcfdist_labels_from_summary(tmp_path):
    (tmp_path / "x.summary.vcf").write_text(SUMMARY_HDR +
        "chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT:BD\t0/1:TP\t0/1:TP\n"   # TP (query BD)
        "chr1\t20\t.\tG\tT\t30\tPASS\t.\tGT:BD\t./.:.\t0/1:FP\n"    # FP (query only)
        "chr1\t30\t.\tA\tG\t30\tPASS\t.\tGT:BD\t0/1:FN\t./.:.\n")   # FN (truth only)
    m = cm.load_vcfdist(str(tmp_path / "x."))
    assert m[("chr1", "10", "A", "C")] == "TP"
    assert m[("chr1", "20", "G", "T")] == "FP"
    assert m[("chr1", "30", "A", "G")] == "FN"


def test_load_vcfeval_labels(tmp_path):
    (tmp_path / "tp.vcf").write_text(VCF_HDR + "chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT\t0/1\n")
    (tmp_path / "fp.vcf").write_text(VCF_HDR + "chr1\t20\t.\tG\tT\t30\tPASS\t.\tGT\t0/1\n")
    (tmp_path / "fn.vcf").write_text(VCF_HDR + "chr1\t30\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n")
    m = cm.load_vcfeval(str(tmp_path))
    assert m == {("chr1", "10", "A", "C"): "TP",
                 ("chr1", "20", "G", "T"): "FP",
                 ("chr1", "30", "A", "G"): "FN"}


def test_build_cm_cross_tabulates_and_marks_absent():
    # v3: TP@k1, FP@k2, FN@k3 ; other: TP@k1, TP@k2 (disagrees), k3 absent -> N
    ref = {("chr1", "10", "A", "C"): "TP",
           ("chr1", "20", "G", "T"): "FP",
           ("chr1", "30", "A", "G"): "FN"}
    oth = {("chr1", "10", "A", "C"): "TP",
           ("chr1", "20", "G", "T"): "TP"}
    snp = cm.build_cm(ref, oth)["SNP"]
    assert snp["TP"]["TP"] == 1     # agree TP
    assert snp["FP"]["TP"] == 1     # v3 FP vs other TP (discordant)
    assert snp["FN"]["N"] == 1      # v3 FN, variant absent from other


def test_build_cm_size_class_split():
    # an SV insertion (>=50bp) keyed correctly under SV
    ref = {("chr1", "10", "A", "A" + "C" * 60): "TP"}
    oth = {("chr1", "10", "A", "A" + "C" * 60): "FP"}
    counts = cm.build_cm(ref, oth, sv_threshold=50)
    assert counts["SV"]["TP"]["FP"] == 1
    assert "SNP" not in counts
