import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

datasets = ["nist", "cmrg"]
names = ["original", "A", "B", "C", "D"]
evals = ["", "standard"]
colors = ["k", "C0", "C1", "C2", "C3"]

fig, ax = plt.subplots(2, 4, figsize=(20,8))

for di, dataset in enumerate(datasets):
    for ei, evaluation in enumerate(evals):
        for ni, name in enumerate(names):
            filename = "O" if name == "original" else name
            evalstring = "std" if evaluation == "standard" else ""
            with open(f"/x/vcfdist/data/pfda-v2/{dataset}_vcfdist/K4GT3_HG002_{filename}{evalstring}.precision-recall.tsv") as csv:
                label = f"{name}"
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
                ax[di][ei*2].plot(snp_recall, snp_prec, linestyle='', 
                        marker='.', color=colors[ni], label=label)
                ax[di][ei*2+1].plot(indel_recall, indel_prec, linestyle='', 
                        marker='.', color=colors[ni], label=label)
     
        # SNP plot
        ax[di][ei*2].set_xlabel("Recall")
        ax[di][ei*2].set_ylabel("Precision")
        if dataset == "nist":
            ax[di][ei*2].set_xlim(0.995, 1)
            ax[di][ei*2].set_xticks(np.arange(0.995, 1.001, 0.001))
            ax[di][ei*2].set_ylim(0.998, 1.00001)
            ax[di][ei*2].set_yticks(np.arange(0.998, 1.0001, 0.0004))
            ax[di][ei*2].set_yticklabels([f"{x:.4f}" for x in 
                np.arange(0.998, 1.00001, 0.0004)])
        if dataset == "cmrg":
            ax[di][ei*2].set_xlim(0.9, 1)
            ax[di][ei*2].set_xticks(np.arange(0.9, 1.01, 0.01))
            ax[di][ei*2].set_ylim(0.99, 1.0001)
            ax[di][ei*2].set_yticks(np.arange(0.99, 1.001, 0.001))
            ax[di][ei*2].set_yticklabels([f"{x:.3f}" for x in 
                np.arange(0.99, 1.001, 0.001)])
        ax[di][ei*2].legend(loc="lower left")
        ax[di][ei*2].set_title(f"{dataset.upper()} SNPs {evaluation}")

        # INDEL plot
        ax[di][ei*2+1].set_xlabel("Recall")
        ax[di][ei*2+1].set_ylabel("Precision")
        if dataset == "nist":
            ax[di][ei*2+1].set_xlim(0.98, 1)
            ax[di][ei*2+1].set_xticks(np.arange(0.98, 1.001, 0.004))
            ax[di][ei*2+1].set_ylim(0.98, 1.001)
            ax[di][ei*2+1].set_yticks(np.arange(0.98, 1.001, 0.004))
        if dataset == "cmrg":
            ax[di][ei*2+1].set_xlim(0.9, 1)
            ax[di][ei*2+1].set_xticks(np.arange(0.9, 1.01, 0.01))
            ax[di][ei*2+1].set_ylim(0.9, 1.001)
            ax[di][ei*2+1].set_yticks(np.arange(0.9, 1.01, 0.01))
        ax[di][ei*2+1].legend(loc="lower left")
        ax[di][ei*2+1].set_title(f"{dataset.upper()} INDELs {evaluation}")

plt.tight_layout()
plt.savefig('img/8_vcfdist_pr_plot.png')
