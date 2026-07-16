import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import matplotlib; matplotlib.use("Agg")
import plot_pr_curves as pc

CURVE = ("tool\tdataset\tsize_class\tmin_qual\tprecision\trecall\n"
         "vcfdist-v3\thprc\tSNP\t0\t0.97\t0.95\n"
         "vcfdist-v3\thprc\tSNP\t10\t0.99\t0.90\n"
         "vcfeval\thprc\tSNP\t0\t0.96\t0.94\n")

def test_load_curve_sorted_by_recall():
    import tempfile, os as _os
    p = tempfile.mktemp(suffix=".tsv"); open(p, "w").write(CURVE)
    d = pc.load_curve(p)
    pts = d[("vcfdist-v3", "hprc", "SNP")]
    assert pts == [(0.90, 0.99), (0.95, 0.97)]   # ascending recall, independent expected
    _os.remove(p)

def test_main_writes_pdf(tmp_path):
    c = tmp_path / "pr.tsv"; c.write_text(CURVE)
    out = tmp_path / "pr.pdf"
    pc.main(str(c), str(out))
    assert out.exists() and out.stat().st_size > 0
