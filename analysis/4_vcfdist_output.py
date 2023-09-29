import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

names = ["O"]
max_qual = 31

fig, ax = plt.subplots(3, 2, figsize=(8,9))

# EDIT DISTANCE
sub_des = []
ins_des = []
del_des = []
des = []
sub_des = []
ins_eds = []
del_eds = []
eds = []
quals = []
with open(f"/x/vcfdist/data/pfda-v2/nist_vcfdist/K4GT3_HG002_O.distance.tsv") as csv:
    next(csv)
    label = f"O"
    for line in csv:
        qual, sub_de, ins_de, del_de, sub_ed, ins_ed, del_ed, de, ed, score, qscore = line.split("\t")
        quals.append(int(qual))

        sub_des.append(int(sub_de))
        ins_des.append(int(ins_de))
        del_des.append(int(del_de))
        des.append(int(de))

        ins_eds.append(int(ins_ed))
        del_eds.append(int(del_ed))
        eds.append(int(ed))

ax[1][0].plot(quals[:max_qual], sub_des[:max_qual], marker='.', linestyle='', color=f"C0", label=label)
ax[2][0].plot(quals[:max_qual], ins_des[:max_qual], marker='.', linestyle='', color=f"C1", label=label)
ax[2][0].plot(quals[:max_qual], del_des[:max_qual], marker='.', linestyle='', color=f"C2", label=label)

ax[1][1].plot(quals[:max_qual], sub_des[:max_qual], marker='.', linestyle='', color=f"C0", label=label)
ax[2][1].plot(quals[:max_qual], ins_eds[:max_qual], marker='.', linestyle='', color=f"C1", label=label)
ax[2][1].plot(quals[:max_qual], del_eds[:max_qual], marker='.', linestyle='', color=f"C2", label=label)
 
for i in range(2):
    ax[1][i].set_xlabel("Quality Threshold")
    ax[1][i].set_ylabel("Distinct Edits (DE)")
    ax[1][i].set_yscale('log')
    ax[1][i].legend(["SNP"], loc="upper left") 
for i in range(2):
    ax[2][i].set_xlabel("Quality Threshold")
    ax[2][i].set_ylabel("Edit Distance (ED)")
    ax[2][i].set_yscale('log')
    ax[2][i].legend(["INS", "DEL" ], loc="upper left")

# PRECISION RECALL
with open(f"/x/vcfdist/data/pfda-v2/nist_vcfdist/K4GT3_HG002_O.precision-recall.tsv") as csv:
    label = "O"
    indel_recall = []
    indel_prec = []
    snp_recall = []
    snp_prec = []
    next(csv)

    for line in csv:
        typ, qual, prec, recall, f1, f1q, truth_tot, truth_tp, truth_pp, truth_fn, \
            query_tot, query_tp, query_pp, query_fp = line.split('\t')
        if typ == "INDEL":
            indel_recall.append(float(truth_tp) / float(truth_tot))
            if int(truth_tp) + int(query_pp) + int(query_fp) == 0:
                indel_prec.append(1)
            else:
                indel_prec.append(float(truth_tp) / 
                        (float(truth_tp) + float(query_pp) + float(query_fp)))
        elif typ == "SNP":
            snp_recall.append(float(truth_tp) / float(truth_tot))
            if int(truth_tp) + int(query_pp) + int(query_fp) == 0:
                snp_prec.append(1)
            else:
                snp_prec.append(float(truth_tp) / 
                        (float(truth_tp) + float(query_pp) + float(query_fp)))

ax[0][0].set_title("SNPs")
ax[0][0].set_xlabel("Recall")
ax[0][0].set_ylabel("Precision")
ax[0][0].set_xlim(0.995, 1)
ax[0][0].set_xticks(np.arange(0.995, 1.0001, 0.001))
ax[0][0].set_ylim(0.995, 1)
ax[0][0].set_yticks(np.arange(0.995, 1.0001, 0.001))
ax[0][0].plot(snp_recall, snp_prec, linestyle='', 
        marker='.', color=f"C0", label=label)

ax[0][1].set_title("INDELs")
ax[0][1].set_xlabel("Recall")
ax[0][1].set_ylabel("Precision")
ax[0][1].set_xlim(0.98, 1)
ax[0][1].set_xticks(np.arange(0.98, 1.001, 0.005))
ax[0][1].set_ylim(0.98, 1)
ax[0][1].set_yticks(np.arange(0.98, 1.001, 0.005))
ax[0][1].plot(indel_recall, indel_prec, linestyle='', 
        marker='.', color=f"C0", label=label)

plt.tight_layout()
plt.savefig('img/4_vcfdist_plot.pdf', format="pdf")
