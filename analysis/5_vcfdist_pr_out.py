import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# names = ["mm2-short", "mm2-ont", "mm2-hifi", "mm2-asm5", "mm2-asm10"]
# names = ["A", "B", "C", "D"]
names = ["O"]

fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(8,3))

for query in names:
    with open(f"/x/vcfdist/data/pfda-v2/nist_vcfdist/K4GT3_HG002_O.precision-recall.tsv") as csv:
        label = f"VCFdist {query}"
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
            typ, qual, prec, recall, f1, f1q, truth_tot, truth_tp, truth_pp, truth_fn, \
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
        # ax1.plot(snp_recall, snp_prec, linestyle='', 
        #         marker='+', color=f"C0", label=label)
        # ax2.plot(indel_recall, indel_prec, linestyle='', 
        #         marker='+', color=f"C0", label=label)
        ax1.plot(snp_recall2, snp_prec2, linestyle='', 
                marker='.', color=f"C0", label=label)
        ax2.plot(indel_recall2, indel_prec2, linestyle='', 
                marker='.', color=f"C0", label=label)
    ax1.set_title("SNPs")
    ax2.set_title("INDELs")
 
ax1.set_xlabel("Recall")
ax1.set_ylabel("Precision")
ax1.set_xlim(0.98, 1)
ax1.set_xticks(np.arange(0.98, 1.001, 0.005))
ax1.set_ylim(0.98, 1)
ax1.set_yticks(np.arange(0.98, 1.001, 0.005))
# ax1.legend(loc="lower left")

ax2.set_xlabel("Recall")
ax2.set_ylabel("Precision")
ax2.set_xlim(0.98, 1)
ax2.set_xticks(np.arange(0.98, 1.001, 0.005))
ax2.set_ylim(0.98, 1)
ax2.set_yticks(np.arange(0.98, 1.001, 0.005))

plt.tight_layout()
plt.savefig('img/5_vcfdist_pr_out.png', dpi=300)
