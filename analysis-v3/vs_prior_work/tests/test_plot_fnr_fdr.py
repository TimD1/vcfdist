import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import matplotlib; matplotlib.use("Agg")
import plot_fnr_fdr as pf

COUNTS = ("tool\tdataset\tsize_class\ttp_query\ttp_truth\tfp\tfn\n"
          "vcfdist-v3\thprc\tSNP\t100\t100\t2\t5\n"
          "vcfeval\thprc\tSNP\t100\t100\t4\t8\n")

def test_fnr_fdr_math():
    assert abs(pf.fnr({"fn": 5, "tp_truth": 95}) - 0.05) < 1e-9
    assert abs(pf.fdr({"fp": 2, "tp_query": 98}) - 0.02) < 1e-9
    assert pf.fdr({"fp": 0, "tp_query": 0}) == 0.0

def test_main_writes_pdfs(tmp_path):
    c = tmp_path / "counts.tsv"; c.write_text(COUNTS)
    fnr_pdf, fdr_pdf = tmp_path / "fnr.pdf", tmp_path / "fdr.pdf"
    pf.main(str(c), str(fnr_pdf), str(fdr_pdf))
    assert fnr_pdf.exists() and fnr_pdf.stat().st_size > 0
    assert fdr_pdf.exists() and fdr_pdf.stat().st_size > 0
