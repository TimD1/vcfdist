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

fig, ax = plt.subplots(3, 4, figsize=(mm2in(180), mm2in(108)))

for di, dataset in enumerate(datasets):
    for ei, eval_list in enumerate(names):
        out3 = open(f"tsv/3{chr(ord('a')+ei)}.tsv", "w")
        print("DATASET\tSUBMISSION\tREPRESENTATION\tVARIANT_TYPE\tPRECISION\tRECALL", file=out3)
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
                            print(f"{dataset.upper()}\t{sub_id}\t{label.upper()}\tINDEL\t{indel_prec[-1]}\t{indel_recall[-1]}", file=out3)
                        elif line[:19] == "SNP,*,*,PASS,*,QUAL":
                            recall = float(line.split(',')[7])
                            snp_recall.append(recall)
                            prec = float(line.split(',')[8])
                            if prec == 0: prec = 1
                            snp_prec.append(prec)
                            snp_f1 = line.split(',')[10]
                            print(f"{dataset.upper()}\t{sub_id}\t{label.upper()}\tSNP\t{snp_prec[-1]}\t{snp_recall[-1]}", file=out3)
                    ax[ei][di*2].plot(snp_recall, snp_prec, linestyle='', 
                            marker='.', markersize=1, color=colors[ei][ni], label=label)
                    ax[ei][di*2+1].plot(indel_recall, indel_prec, linestyle='', 
                            marker='.', markersize=1, color=colors[ei][ni], label=label)
            except FileNotFoundError:
                print(f"'/x/vcfdist/data/pfda-v2/{dataset}_vcfeval/{sub_id}_HG002_{filename}.roc.all.csv' not found.")
                continue
     
        # SNP plot
        ax[ei][0].set_ylabel(r"\textbf{Precision}", fontsize=7)
        ax[2][di*2].set_xlabel(r"\textbf{Recall}", fontsize=7)
        ax[2][di*2+1].set_xlabel(r"\textbf{Recall}", fontsize=7)
        if dataset == "nist":
            ax[ei][di*2].set_xlim(0.99, 1)
            ax[ei][di*2].set_xticks(np.arange(0.99, 1.001, 0.005))
            ax[ei][di*2].set_xticklabels([f"{x:.3f}" for x in 
                np.arange(0.99, 1.001, 0.005)], fontsize=5)
            ax[ei][di*2].set_ylim(0.999, 1.00002)
            ax[ei][di*2].set_yticks(np.arange(0.999, 1.00002, 0.0002))
            ax[ei][di*2].set_yticklabels([f"{x:.4f}" for x in 
                np.arange(0.999, 1.00002, 0.0002)], fontsize=5)
        if dataset == "cmrg":
            ax[ei][di*2].set_xlim(0.9, 1)
            ax[ei][di*2].set_xticks(np.arange(0.9, 1.01, 0.02))
            ax[ei][di*2].set_xticklabels([f"{x:.2f}" for x in 
                np.arange(0.9, 1.01, 0.02)], fontsize=5)
            ax[ei][di*2].set_ylim(0.95, 1.001)
            ax[ei][di*2].set_yticks(np.arange(0.95, 1.001, 0.01))
            ax[ei][di*2].set_yticklabels([f"{x:.3f}" for x in 
                np.arange(0.95, 1.001, 0.01)], fontsize=5)
        if di == 0:
            ax[ei][di*2].legend(loc="lower left", fontsize=5, markerscale=4)

        # INDEL plot
        # ax[ei][0].set_xlabel("Recall", fontsize=7)
        # ax[ei][di*2+1].set_ylabel("Precision", fontsize=7)
        if dataset == "nist":
            ax[ei][di*2+1].set_xlim(0.99, 1)
            ax[ei][di*2+1].set_xticks(np.arange(0.99, 1.001, 0.005))
            ax[ei][di*2+1].set_xticklabels([f"{x:.3f}" for x in 
                np.arange(0.99, 1.001, 0.005)], fontsize=5)
            ax[ei][di*2+1].set_ylim(0.995, 1.0001)
            ax[ei][di*2+1].set_yticks(np.arange(0.995, 1.0001, 0.001))
            ax[ei][di*2+1].set_yticklabels([f"{x:.3f}" for x in 
                np.arange(0.995, 1.0001, 0.001)], fontsize=5)
        if dataset == "cmrg":
            ax[ei][di*2+1].set_xlim(0.9, 1)
            ax[ei][di*2+1].set_xticks(np.arange(0.9, 1.01, 0.02))
            ax[ei][di*2+1].set_xticklabels([f"{x:.2f}" for x in 
                np.arange(0.9, 1.01, 0.02)], fontsize=5)
            ax[ei][di*2+1].set_ylim(0.9, 1.002)
            ax[ei][di*2+1].set_yticks(np.arange(0.9, 1.002, 0.02))
            ax[ei][di*2+1].set_yticklabels([f"{x:.2f}" for x in 
                np.arange(0.9, 1.01, 0.02)], fontsize=5)
        out3.close()

plt.tight_layout()
plt.savefig('img/7_vcfeval_pr_plot.pdf', format="pdf")
