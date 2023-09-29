import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({"figure.facecolor": (0,0,0,0)})

datasets = ["nist", "cmrg"]
names = [
    ["original", "bwa", "mm2-ont", "mm2-pb", "pbmm2"],
    ["original", "A", "B", "C", "D"],
    ["originalstd", "Astd", "Bstd", "Cstd", "Dstd"]
]
colors = [
    ["k", "C7", "C2", "C3", "C6"],
    ["k", "C0", "C1", "C2", "C3"],
    ["k", "C0", "C1", "C2", "C3"]
]
sub_id = "K4GT3"

fig, ax = plt.subplots(3, 4, figsize=(20,12))

for di, dataset in enumerate(datasets):
    for ei, eval_list in enumerate(names):
        for ni, name in enumerate(eval_list):
            filename = name.replace("original", "O")
            try:
                with open(f"/x/vcfdist/data/pfda-v2/{dataset}_vcfeval/{sub_id}_HG002_{filename}.roc.all.csv") as csv:
                    label = f"{name}".replace("std","")
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
                    ax[ei][di*2].plot(snp_recall, snp_prec, linestyle='', 
                            marker='.', color=colors[ei][ni], label=label)
                    ax[ei][di*2+1].plot(indel_recall, indel_prec, linestyle='', 
                            marker='.', color=colors[ei][ni], label=label)
            except FileNotFoundError:
                print(f"'/x/vcfdist/data/pfda-v2/{dataset}_vcfeval/{sub_id}_HG002_{filename}.roc.all.csv' not found.")
                continue
     
        # SNP plot
        ax[ei][di*2].set_xlabel("Recall", fontsize=15)
        ax[ei][di*2].set_ylabel("Precision", fontsize=15)
        if dataset == "nist":
            ax[ei][di*2].set_xlim(0.99, 1)
            ax[ei][di*2].set_xticks(np.arange(0.99, 1.001, 0.002))
            ax[ei][di*2].set_ylim(0.999, 1.00001)
            ax[ei][di*2].set_yticks(np.arange(0.999, 1.0001, 0.0002))
            ax[ei][di*2].set_yticklabels([f"{x:.4f}" for x in 
                np.arange(0.999, 1.00001, 0.0002)])
        if dataset == "cmrg":
            ax[ei][di*2].set_xlim(0.9, 1)
            ax[ei][di*2].set_xticks(np.arange(0.9, 1.01, 0.02))
            ax[ei][di*2].set_ylim(0.95, 1.0001)
            ax[ei][di*2].set_yticks(np.arange(0.95, 1.001, 0.01))
            ax[ei][di*2].set_yticklabels([f"{x:.3f}" for x in 
                np.arange(0.95, 1.001, 0.01)])
        if di == 0:
            ax[ei][di*2].legend(loc="lower left", fontsize=15, markerscale=3)

        # INDEL plot
        ax[ei][di*2+1].set_xlabel("Recall", fontsize=15)
        ax[ei][di*2+1].set_ylabel("Precision", fontsize=15)
        if dataset == "nist":
            ax[ei][di*2+1].set_xlim(0.99, 1)
            ax[ei][di*2+1].set_xticks(np.arange(0.99, 1.001, 0.002))
            ax[ei][di*2+1].set_ylim(0.995, 1.0001)
            ax[ei][di*2+1].set_yticks(np.arange(0.995, 1.001, 0.001))
        if dataset == "cmrg":
            ax[ei][di*2+1].set_xlim(0.9, 1)
            ax[ei][di*2+1].set_xticks(np.arange(0.9, 1.01, 0.02))
            ax[ei][di*2+1].set_ylim(0.9, 1.001)
            ax[ei][di*2+1].set_yticks(np.arange(0.9, 1.01, 0.02))

plt.tight_layout()
plt.savefig('img/7_vcfeval_pr_plot.pdf', format="pdf")
