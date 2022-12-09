import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

names = ["A", "B", "C", "D"]
max_qual = 51

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,6))

# QUAL	SUB_DE	INS_DE	DEL_DE	INS_ED	DEL_ED	DE	ED	SCORE	SUB_DE_PCT	INDEL_DE_PCT	INDEL_ED_PCT
truth = "C"
for query in names:
    sub_des = []
    ins_des = []
    del_des = []
    des = []
    sub_des = []
    ins_eds = []
    del_eds = []
    eds = []
    quals = []
    with open(f"1_results/{query}-query_{truth}-truth_distance.tsv") as csv:
        next(csv)
        label = f"VCFdist Query={query} Truth={truth}"
        for line in csv:
            qual, sub_de, ins_de, del_de, ins_ed, del_ed, de, ed, score, sub_de_pct, indel_de_pct, indel_ed_pct = line.split("\t")
            quals.append(int(qual))

            sub_des.append(int(sub_de))
            ins_des.append(int(ins_de))
            del_des.append(int(del_de))
            des.append(int(de))

            ins_eds.append(int(ins_ed))
            del_eds.append(int(del_ed))
            eds.append(int(ed))

    ax1.plot(quals[:max_qual], sub_des[:max_qual], marker='.', color=f"C0", label=label)
    ax1.plot(quals[:max_qual], ins_des[:max_qual], marker='.', color=f"C1", label=label)
    ax1.plot(quals[:max_qual], del_des[:max_qual], marker='.', color=f"C2", label=label)
    ax1.plot(quals[:max_qual], des[:max_qual], marker='.', color=f"C3", label=label)

    ax2.plot(quals[:max_qual], sub_des[:max_qual], marker='.', color=f"C0", label=label)
    ax2.plot(quals[:max_qual], ins_eds[:max_qual], marker='.', color=f"C1", label=label)
    ax2.plot(quals[:max_qual], del_eds[:max_qual], marker='.', color=f"C2", label=label)
    ax2.plot(quals[:max_qual], eds[:max_qual], marker='.', color=f"C3", label=label)
 
ax1.set_xlabel("Quality Threshold")
ax1.set_ylabel("Distinct Edits (DE)")
ax1.set_yscale('log')
ax1.legend(["SUB", "INS", "DEL", "ALL"], loc="upper left")

ax2.set_xlabel("Quality Threshold")
ax2.set_ylabel("Edit Distance (ED)")
ax2.set_yscale('log')
ax2.legend(["SUB", "INS", "DEL", "ALL"], loc="upper left")

plt.tight_layout()
plt.savefig('img/4_sw_dist_curve.png', dpi=300)
