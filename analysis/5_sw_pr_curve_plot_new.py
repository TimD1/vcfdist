import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# names = ["mm2-short", "mm2-ont", "mm2-hifi", "mm2-asm5", "mm2-asm10"]
names = ["A", "B", "C", "D"]

fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(20,8))
fig.suptitle("Varying Truth/Query Smith-Waterman Parameters")

idx = 0
for query in names:
    for truth in "C":
        with open(f"5_results/{query}-query_{truth}-truth.roc.all.csv") as csv:
            label = f"VCFeval Query={query} Truth={truth}"
            indel_recall = []
            indel_prec = []
            snp_recall = []
            snp_prec = []
            for line in csv:
                if line[:21] == "INDEL,*,*,PASS,*,QUAL":
                    indel_recall.append(float(line.split(',')[7]))
                    indel_prec.append(float(line.split(',')[8]))
                    indel_f1 = line.split(',')[10]
                elif line[:19] == "SNP,*,*,PASS,*,QUAL":
                    recall = float(line.split(',')[7])
                    snp_recall.append(recall)
                    prec = float(line.split(',')[8])
                    if prec == 0: prec = 1
                    snp_prec.append(prec)
                    snp_f1 = line.split(',')[10]
            ax1.plot(snp_recall, snp_prec, linestyle='', 
                    marker='.', color=f"C{idx}", label=label)
            ax2.plot(indel_recall, indel_prec, linestyle='', 
                    marker='.', color=f"C{idx}", label=label)

        with open(f"5_results/{query}-query_{truth}-truth_precision-recall.tsv") as csv:
            label = f"VCFdist Query={query} Truth={truth}"
            indel_recall = []
            indel_prec = []
            snp_recall = []
            snp_prec = []

            # precision/recall with no partial credit
            indel_recall2 = []
            indel_prec2 = []
            snp_recall2 = []
            snp_prec2 = []
            next(csv)

            for line in csv:
                typ, qual, prec, recall, f1, truth_tot, truth_tp, truth_pp, truth_fn, \
                    query_tot, query_tp, query_pp, query_fp = line.split('\t')
                if typ == "INDEL":
                    indel_recall.append(float(recall))
                    indel_prec.append(float(prec))
                    indel_recall2.append(float(truth_tp) / float(truth_tot))
                    if int(truth_tp) + int(query_pp) + int(query_fp) == 0:
                        indel_prec2.append(1)
                    else:
                        indel_prec2.append(float(truth_tp) / 
                                (float(truth_tp) + float(query_pp) + float(query_fp)))
                elif typ == "SNP":
                    snp_recall.append(float(recall))
                    snp_prec.append(float(prec))
                    snp_recall2.append(float(truth_tp) / float(truth_tot))
                    if int(truth_tp) + int(query_pp) + int(query_fp) == 0:
                        snp_prec2.append(1)
                    else:
                        snp_prec2.append(float(truth_tp) / 
                                (float(truth_tp) + float(query_pp) + float(query_fp)))
            ax1.plot(snp_recall, snp_prec, linestyle='', 
                    marker='+', color=f"C{idx}", label=label)
            ax2.plot(indel_recall, indel_prec, linestyle='', 
                    marker='+', color=f"C{idx}", label=label)
            ax1.plot(snp_recall2, snp_prec2, linestyle='', 
                    marker='x', color=f"C{idx}", label=label)
            ax2.plot(indel_recall2, indel_prec2, linestyle='', 
                    marker='x', color=f"C{idx}", label=label)
        ax1.set_title("SNPs")
        ax2.set_title("INDELs")
    idx += 1
 
ax1.set_xlabel("Recall")
ax1.set_ylabel("Precision")
ax1.set_xlim(0.9, 1)
ax1.set_xticks(np.arange(0.9, 1.01, 0.01))
ax1.set_ylim(0.995, 1.0001)
ax1.set_yticks(np.arange(0.995, 1.001, 0.001))
ax1.set_yticklabels([f"{x:.3f}" for x in np.arange(0.995, 1.001, 0.001)])
ax1.legend(loc="lower left")

ax2.set_xlabel("Recall")
ax2.set_ylabel("Precision")
ax2.set_xlim(0, 1)
ax2.set_xticks(np.arange(0, 1.1, 0.1))
ax2.set_ylim(0.6, 1.001)
ax2.set_yticks(np.arange(0.6, 1.05, 0.05))

plt.tight_layout()
plt.savefig('img/5_sw_pr_curve_new.png')
