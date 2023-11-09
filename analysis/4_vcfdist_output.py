import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['text.usetex'] = True
plt.rcParams.update({"figure.facecolor": (0,0,0,0)})

def mm2in(x):
    return x / 25.4

names = ["O"]
max_qual = 31

fig, ax = plt.subplots(3, 2, figsize=(mm2in(80),mm2in(90)))

# source data headers
out_5bi = open("tsv/5bi.tsv", "w")
out_5bii = open("tsv/5bii.tsv", "w")
out_5biii = open("tsv/5biii.tsv", "w")
print("DATASET\tSUBMISSION\tVARIANT_TYPE\tQUALITY_SCORE\tPRECISION\tRECALL", file=out_5bi)
print("DATASET\tSUBMISSION\tVARIANT_TYPE\tQUALITY_SCORE\tEDIT_DISTANCE\tDISTINCT_EDITS", file=out_5bii)
print("DATASET\tSUBMISSION\tVARIANT_TYPE\tQUALITY_SCORE\tEDIT_DISTANCE\tDISTINCT_EDITS", file=out_5biii)

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
        print(f"NIST\tK4GT3\tSUB\t{qual}\t{sub_ed}\t{sub_de}", file=out_5bii)
        print(f"NIST\tK4GT3\tINS\t{qual}\t{ins_ed}\t{ins_de}", file=out_5biii)
        print(f"NIST\tK4GT3\tDEL\t{qual}\t{del_ed}\t{del_de}", file=out_5biii)

ax[1][0].plot(quals[:max_qual], sub_des[:max_qual], marker='.', markersize=1, linestyle='', color=f"C0", label=label)
ax[2][0].plot(quals[:max_qual], ins_des[:max_qual], marker='.', markersize=1, linestyle='', color=f"C1", label=label)
ax[2][0].plot(quals[:max_qual], del_des[:max_qual], marker='.', markersize=1, linestyle='', color=f"C2", label=label)

ax[1][1].plot(quals[:max_qual], sub_des[:max_qual], marker='.', markersize=1, linestyle='', color=f"C0", label=label)
ax[2][1].plot(quals[:max_qual], ins_eds[:max_qual], marker='.', markersize=1, linestyle='', color=f"C1", label=label)
ax[2][1].plot(quals[:max_qual], del_eds[:max_qual], marker='.', markersize=1, linestyle='', color=f"C2", label=label)
 
for i in range(2):
    ax[1][i].set_xlabel(r"\textbf{Quality Threshold}", fontsize=7)
    if i % 2:
        ax[1][i].set_ylabel(r"\textbf{Edit Distance}", fontsize=7)
    else:
        ax[1][i].set_ylabel(r"\textbf{Distinct Edits}", fontsize=7)
    ax[1][i].set_yscale('log')
    ax[1][i].yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    ax[1][i].tick_params(axis='both', which='both', labelsize=5)
    ax[1][i].legend(["substitution"], markerscale=3, loc="upper left", fontsize=5) 
for i in range(2):
    ax[2][i].set_xlabel(r"\textbf{Quality Threshold}", fontsize=7)
    if i % 2:
        ax[2][i].set_ylabel(r"\textbf{Edit Distance}", fontsize=7)
    else:
        ax[2][i].set_ylabel(r"\textbf{Distinct Edits}", fontsize=7)
    ax[2][i].set_yscale('log')
    ax[2][i].yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    ax[2][i].tick_params(axis='both', which='both', labelsize=5)
    ax[2][i].legend(["insertion", "deletion" ], markerscale=3, loc="upper left", fontsize=5)

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
            print(f"NIST\tK4GT3\tINDEL\t{qual}\t{indel_prec[-1]}\t{indel_recall[-1]}", file=out_5bi)
        elif typ == "SNP":
            snp_recall.append(float(truth_tp) / float(truth_tot))
            if int(truth_tp) + int(query_pp) + int(query_fp) == 0:
                snp_prec.append(1)
            else:
                snp_prec.append(float(truth_tp) / 
                        (float(truth_tp) + float(query_pp) + float(query_fp)))
            print(f"NIST\tK4GT3\tSNP\t{qual}\t{snp_prec[-1]}\t{snp_recall[-1]}", file=out_5bi)

ax[0][0].set_title(r"\textbf{SNPs}", fontsize=7)
ax[0][0].set_xlabel(r"\textbf{Recall}", fontsize=7)
ax[0][0].set_ylabel(r"\textbf{Precision}", fontsize=7)
ax[0][0].set_xlim(0.995, 1)
ax[0][0].set_xticks(np.arange(0.995, 1.0001, 0.0025))
ax[0][0].set_ylim(0.995, 1)
ax[0][0].set_yticks(np.arange(0.995, 1.0001, 0.001))
ax[0][0].tick_params(axis='both', which='both', labelsize=5)
ax[0][0].plot(snp_recall, snp_prec, linestyle='', 
        marker='.', color=f"C0", label=label, markersize=1)

ax[0][1].set_title(r"\textbf{INDELs}", fontsize=7)
ax[0][1].set_xlabel(r"\textbf{Recall}", fontsize=7)
ax[0][1].set_ylabel(r"\textbf{Precision}", fontsize=7)
ax[0][1].set_xlim(0.98, 1)
ax[0][1].set_xticks(np.arange(0.98, 1.001, 0.01))
ax[0][1].set_ylim(0.98, 1)
ax[0][1].set_yticks(np.arange(0.98, 1.001, 0.005))
ax[0][1].tick_params(axis='both', which='both', labelsize=5)
ax[0][1].plot(indel_recall, indel_prec, linestyle='', 
        marker='.', color=f"C0", label=label, markersize=1)

plt.tight_layout()
out_5bi.close()
out_5bii.close()
out_5biii.close()
plt.savefig('img/4_vcfdist_plot.pdf', format="pdf")
