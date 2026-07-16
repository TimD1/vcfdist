import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import parse_vcfdist as pv

SUMMARY = """VAR_TYPE\tTHRESHOLD\tMIN_QUAL\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE
SNP\tNONE\t0\t100\t101\t5\t3\t0.97\t0.95\t0.96\t14.0
SNP\tBEST\t6\t99\t100\t6\t2\t0.98\t0.94\t0.96\t14.0
INDEL\tNONE\t0\t40\t41\t4\t2\t0.95\t0.90\t0.92\t11.0
INDEL\tBEST\t3\t39\t40\t5\t1\t0.97\t0.88\t0.92\t11.0
SV\tNONE\t0\t10\t10\t2\t1\t0.90\t0.83\t0.86\t8.0
SV\tBEST\t4\t9\t9\t3\t0\t1.00\t0.75\t0.85\t8.0
ALL\tNONE\t0\t150\t152\t11\t6\t0.96\t0.93\t0.94\t12.0
ALL\tBEST\t5\t147\t149\t14\t3\t0.98\t0.91\t0.94\t12.0
"""

CURVE = """VAR_TYPE\tMIN_QUAL\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE\tTRUTH_TOTAL\tTRUTH_TP\tTRUTH_FN\tQUERY_TOTAL\tQUERY_TP\tQUERY_FP
SNP\t0\t0.970\t0.950\t0.96\t14.0\t105\t100\t5\t104\t101\t3
SNP\t1\t0.980\t0.940\t0.96\t14.0\t105\t99\t6\t101\t100\t1
SV\t0\t0.900\t0.830\t0.86\t8.0\t12\t10\t2\t11\t10\t1
INDEL\t0\t0.950\t0.900\t0.92\t11.0\t44\t40\t4\t42\t41\t2
ALL\t0\t0.960\t0.930\t0.94\t12.0\t150\t150\t11\t152\t152\t6
"""

def test_parse_counts(tmp_path):
    f = tmp_path / "s.tsv"; f.write_text(SUMMARY)
    rows = {r["size_class"]: r for r in pv.parse_counts(str(f))}
    assert set(rows) == {"SNP", "INDEL", "SV"}       # ALL excluded
    assert rows["SNP"] == {"size_class": "SNP", "tp_query": 101,
                           "tp_truth": 100, "fp": 3, "fn": 5}

def test_parse_curve(tmp_path):
    f = tmp_path / "c.tsv"; f.write_text(CURVE)
    rows = pv.parse_curve(str(f))
    snp0 = next(r for r in rows if r["size_class"] == "SNP" and r["min_qual"] == 0.0)
    assert abs(snp0["precision"] - 0.970) < 1e-9
    assert abs(snp0["recall"] - 0.950) < 1e-9
    classes = {r["size_class"] for r in rows}
    assert "ALL" not in classes          # ALL row excluded
    assert classes == {"SNP", "INDEL", "SV"}
