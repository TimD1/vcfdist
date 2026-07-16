import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import aggregate as ag

TIME_LOG = """\tCommand being timed: "vcfdist ..."
\tUser time (seconds): 100.0
\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02:03
\tMaximum resident set size (kbytes): 2048000
"""

def test_parse_time_log(tmp_path):
    p = tmp_path / "t.log"; p.write_text(TIME_LOG)
    wall, rss = ag.parse_time_log(str(p))
    assert wall == 3723.0        # 1h 2m 3s
    assert rss == 2048000

def test_parse_time_log_mm_ss(tmp_path):
    p = tmp_path / "t.log"
    p.write_text("\tElapsed (wall clock) time (h:mm:ss or m:ss): 2:30.5\n"
                 "\tMaximum resident set size (kbytes): 1024\n")
    wall, rss = ag.parse_time_log(str(p))
    assert wall == 150.5
    assert rss == 1024

def test_write_tsv_roundtrip(tmp_path):
    p = tmp_path / "o.tsv"
    ag.write_tsv(str(p), [{"a": 1, "b": 2}], ["a", "b"])
    assert p.read_text() == "a\tb\n1\t2\n"

def test_main_vcfdist_end_to_end(tmp_path):
    summary = tmp_path / "s.tsv"; summary.write_text(
        "VAR_TYPE\tTHRESHOLD\tMIN_QUAL\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE\n"
        "SNP\tNONE\t0\t100\t101\t5\t3\t0.97\t0.95\t0.96\t14.0\n")
    curve = tmp_path / "c.tsv"; curve.write_text(
        "VAR_TYPE\tMIN_QUAL\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE\tTRUTH_TOTAL\tTRUTH_TP\tTRUTH_FN\tQUERY_TOTAL\tQUERY_TP\tQUERY_FP\n"
        "SNP\t0\t0.97\t0.95\t0.96\t14.0\t105\t100\t5\t104\t101\t3\n")
    tlog = tmp_path / "t.log"; tlog.write_text(
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:30.0\n"
        "\tMaximum resident set size (kbytes): 1000\n")
    co, cu, ro = tmp_path / "counts.tsv", tmp_path / "curve.tsv", tmp_path / "rt.tsv"
    ag.main(["--tool", "vcfdist-v3", "--dataset", "hprc", "--summary", str(summary),
             "--curve", str(curve), "--time-log", str(tlog),
             "--counts-out", str(co), "--curve-out", str(cu), "--runtime-out", str(ro)])
    cl = co.read_text().splitlines()
    assert cl[0] == "tool\tdataset\tsize_class\ttp_query\ttp_truth\tfp\tfn"
    assert cl[1] == "vcfdist-v3\thprc\tSNP\t101\t100\t3\t5"
    rl = ro.read_text().splitlines()
    assert rl[0] == "tool\tdataset\twall_seconds\tmax_rss_kb"
    assert rl[1] == "vcfdist-v3\thprc\t30.0\t1000"
    assert cu.read_text().splitlines()[0] == "tool\tdataset\tsize_class\tmin_qual\tprecision\trecall"
