import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import parse_vcfeval as pe

HDR = ("##fileformat=VCFv4.2\n"
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")

def _vcf(tmp_path, name, body):
    p = tmp_path / name; p.write_text(HDR + body); return p

def test_parse_dir_counts_and_curve(tmp_path):
    # query TPs: one SNP qual 20 (GT 0/1 -> 1 allele), one SV qual 30
    _vcf(tmp_path, "tp.vcf",
         "chr1\t10\t.\tA\tT\t20\tPASS\t.\tGT\t0/1\n"
         "chr1\t50\t.\tA\t" + "A" + "C"*60 + "\t30\tPASS\t.\tGT\t1/1\n")
    # query FP: one SNP qual 5
    _vcf(tmp_path, "fp.vcf", "chr1\t20\t.\tG\tC\t5\tPASS\t.\tGT\t0/1\n")
    # truth-side recovered (tp-baseline): the same SNP + SV
    _vcf(tmp_path, "tp-baseline.vcf",
         "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n"
         "chr1\t50\t.\tA\t" + "A" + "C"*60 + "\t.\tPASS\t.\tGT\t1/1\n")
    # truth-side missed (fn): one SNP
    _vcf(tmp_path, "fn.vcf", "chr1\t99\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n")

    counts, curve = pe.parse_dir(str(tmp_path), sv_threshold=50)
    c = {r["size_class"]: r for r in counts}
    assert c["SNP"]["fp"] == 1
    assert c["SNP"]["tp_query"] == 1          # 0/1 -> 1 allele
    assert c["SNP"]["fn"] == 1
    assert c["SNP"]["tp_truth"] == 2   # tp-baseline SNP 1/1 -> weight 2
    assert c["SNP"]["fn"] == 1         # fn SNP 0/1 -> weight 1 (distinct -> catches key swap)
    assert c["SV"]["tp_query"] == 2           # 1/1 -> 2 alleles
    # curve: SNP truth total = tp-baseline(1) + fn(1) = 2
    snp0 = next(r for r in curve if r["size_class"] == "SNP" and r["min_qual"] == 0.0)
    assert abs(snp0["recall"] - 0.5) < 1e-9   # 1 query-TP / 2 truth
    assert abs(snp0["precision"] - 0.5) < 1e-9  # 1 TP / (1 TP + 1 FP)
