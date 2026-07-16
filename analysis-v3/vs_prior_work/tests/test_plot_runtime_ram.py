import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import matplotlib; matplotlib.use("Agg")
import plot_runtime_ram as pr

RUNTIME = ("tool\tdataset\twall_seconds\tmax_rss_kb\n"
           "vcfdist-v3\thprc\t120.0\t2048000\n"
           "vcfeval\thprc\t90.0\t1024000\n")

def test_load_runtime():
    import tempfile, os as _os
    p = tempfile.mktemp(suffix=".tsv"); open(p, "w").write(RUNTIME)
    d = pr.load_runtime(p)
    assert d[("vcfdist-v3", "hprc")]["wall_seconds"] == 120.0
    assert d[("vcfeval", "hprc")]["max_rss_kb"] == 1024000
    _os.remove(p)

def test_main_writes_pdf(tmp_path):
    p = tmp_path / "rt.tsv"; p.write_text(RUNTIME)
    out = tmp_path / "rt.pdf"
    pr.main(str(p), str(out))
    assert out.exists() and out.stat().st_size > 0
