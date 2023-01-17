import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

datasets = ["nist", "cmrg"]
names = [
    ["original", "bwa", "mm2-ont", "mm2-pb", "pbmm2"],
    ["original", "A", "B", "C", "D"]
]
colors = [
    ["k", "C7", "C2", "C3", "C6"],
    ["k", "C0", "C1", "C2", "C3"]
]

fig, ax = plt.subplots(2, 4, figsize=(20,8))

for di, dataset in enumerate(datasets):
    for ei, eval_list in enumerate(names):
        for ni, name in enumerate(eval_list):
            filename = "O" if name == "original" else name
            with open(f"/x/vcfdist/data/pfda-v2/{dataset}_vcfeval/K4GT3_HG002_{filename}.roc.all.csv") as csv:
                label = f"{name}"
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
                ax[di][ei*2].plot(snp_recall, snp_prec, linestyle='', 
                        marker='.', color=colors[ei][ni], label=label)
                ax[di][ei*2+1].plot(indel_recall, indel_prec, linestyle='', 
                        marker='.', color=colors[ei][ni], label=label)
     
        # SNP plot
        ax[di][ei*2].set_xlabel("Recall")
        ax[di][ei*2].set_ylabel("Precision")
        if dataset == "nist":
            ax[di][ei*2].set_xlim(0.99, 1)
            ax[di][ei*2].set_xticks(np.arange(0.99, 1.001, 0.002))
            ax[di][ei*2].set_ylim(0.999, 1.00001)
            ax[di][ei*2].set_yticks(np.arange(0.999, 1.0001, 0.0001))
            ax[di][ei*2].set_yticklabels([f"{x:.4f}" for x in 
                np.arange(0.999, 1.00001, 0.0001)])
        if dataset == "cmrg":
            ax[di][ei*2].set_xlim(0.9, 1)
            ax[di][ei*2].set_xticks(np.arange(0.9, 1.01, 0.01))
            ax[di][ei*2].set_ylim(0.995, 1.0001)
            ax[di][ei*2].set_yticks(np.arange(0.995, 1.001, 0.001))
            ax[di][ei*2].set_yticklabels([f"{x:.3f}" for x in 
                np.arange(0.995, 1.001, 0.001)])
        ax[di][ei*2].legend(loc="lower left")
        ax[di][ei*2].set_title(f"{dataset.upper()} SNPs")

        # INDEL plot
        ax[di][ei*2+1].set_xlabel("Recall")
        ax[di][ei*2+1].set_ylabel("Precision")
        if dataset == "nist":
            ax[di][ei*2+1].set_xlim(0.99, 1)
            ax[di][ei*2+1].set_xticks(np.arange(0.99, 1.001, 0.002))
            ax[di][ei*2+1].set_ylim(0.995, 1.0001)
            ax[di][ei*2+1].set_yticks(np.arange(0.995, 1.001, 0.001))
        if dataset == "cmrg":
            ax[di][ei*2+1].set_xlim(0.9, 1)
            ax[di][ei*2+1].set_xticks(np.arange(0.9, 1.01, 0.01))
            ax[di][ei*2+1].set_ylim(0.9, 1.001)
            ax[di][ei*2+1].set_yticks(np.arange(0.9, 1.01, 0.01))
        ax[di][ei*2+1].legend(loc="lower left")
        ax[di][ei*2+1].set_title(f"{dataset.upper()} INDELs")

plt.tight_layout()
plt.savefig('img/7_vcfeval_pr_plot.png')
