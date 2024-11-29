import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({"figure.facecolor": (0,0,0,0)})

fig, ax = plt.subplots(1, 2, figsize=(10,8))

with open(f"results/dev.precision-recall.tsv") as csv:
    indel_recall = []
    indel_prec = []
    snp_recall = []
    snp_prec = []
    next(csv) # skip header

    for line in csv:
        typ, qual, prec, recall, f1, f1q, truth_tot, truth_tp, truth_fn, \
            query_tot, query_tp, query_fp = line.split('\t')
        if typ == "INDEL":
            indel_recall.append(float(recall))
            indel_prec.append(float(prec))
        elif typ == "SNP":
            snp_recall.append(float(recall))
            snp_prec.append(float(prec))
    ax[0].plot(snp_recall, snp_prec, linestyle='', marker='.')
    ax[1].plot(indel_recall, indel_prec, linestyle='', marker='.')
     
    # SNP plot
    ax[0].set_title("SNPs")
    ax[0].set_xlabel("Recall", fontsize=15)
    ax[0].set_ylabel("Precision", fontsize=15)
    ax[0].set_xlim((0.99,1))
    ax[0].set_ylim((0.99,1))

    # INDEL plot
    ax[1].set_title("INDELs")
    ax[1].set_xlabel("Recall", fontsize=15)
    ax[1].set_ylabel("Precision", fontsize=15)
    ax[1].set_xlim((0.99,1))
    ax[1].set_ylim((0.99,1))

plt.tight_layout()
plt.savefig('results/pr_plot.png', dpi=200)
