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
names = ["original", "A", "B", "C", "D"]
evals = ["vcfeval", "vcfdist_no_partial", "vcfdist", "vcfdist_standard"]
colors = ["k", "C0", "C1", "C2", "C3"]

fig, ax = plt.subplots(4, 4, figsize=(mm2in(180), mm2in(144)))
sub_id = "K4GT3"

for di, dataset in enumerate(datasets):
    for ei, evaluation in enumerate(evals):
        out4 = open(f"tsv/4{chr(ord('a')+ei)}.tsv", "w")
        print("DATASET\tSUBMISSION\tREPRESENTATION\tVARIANT_TYPE\tPRECISION\tRECALL", file=out4)

        for ni, name in enumerate(names):
            filename = "O" if name == "original" else name
            evalstring = "std" if evaluation == "vcfdist_standard" else ""
            if evaluation == "vcfeval":
                with open(f"/x/vcfdist/data/pfda-v2/{dataset}_vcfeval/{sub_id}_HG002_{filename}std.roc.all.csv") as csv:
                    indel_recall = []
                    indel_prec = []
                    snp_recall = []
                    snp_prec = []
                    for line in csv:
                        if line[:21] == "INDEL,*,*,PASS,*,QUAL":
                            indel_recall.append(float(line.split(',')[7]))
                            indel_prec.append(float(line.split(',')[8]))
                            indel_f1 = line.split(',')[10]
                            print(f"{dataset.upper()}\t{sub_id}\t{name.upper()}\tINDEL\t{indel_prec[-1]}\t{indel_recall[-1]}", file=out4)
                        elif line[:19] == "SNP,*,*,PASS,*,QUAL":
                            recall = float(line.split(',')[7])
                            snp_recall.append(recall)
                            prec = float(line.split(',')[8])
                            if prec == 0: prec = 1
                            snp_prec.append(prec)
                            print(f"{dataset.upper()}\t{sub_id}\t{name.upper()}\tSNP\t{indel_prec[-1]}\t{indel_recall[-1]}", file=out4)
                            snp_f1 = line.split(',')[10]
                    ax[ei][di*2].plot(snp_recall, snp_prec, linestyle='', 
                            marker='.', markersize=1, color=colors[ni], label=name)
                    ax[ei][di*2+1].plot(indel_recall, indel_prec, linestyle='', 
                            marker='.', markersize=1, color=colors[ni], label=name)
            else:
                with open(f"/x/vcfdist/data/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}{evalstring}.precision-recall.tsv") as csv:
                    label = f"{name}"
                    indel_recall = []
                    indel_prec = []
                    snp_recall = []
                    snp_prec = []
                    next(csv)

                    for line in csv:
                        typ, qual, prec, recall, f1, f1q, truth_tot, truth_tp, truth_pp, truth_fn, \
                            query_tot, query_tp, query_pp, query_fp = line.split('\t')
                        if evaluation == "vcfdist_no_partial":
                            if typ == "INDEL":
                                indel_recall.append(float(truth_tp) / float(truth_tot))
                                if int(truth_tp) + int(query_fp) == 0:
                                    indel_prec.append(1)
                                else:
                                    indel_prec.append(float(truth_tp) / 
                                            (float(truth_tp) + float(query_fp)))
                                print(f"{dataset.upper()}\t{sub_id}\t{name.upper()}\tINDEL\t{indel_prec[-1]}\t{indel_recall[-1]}", file=out4)
                            elif typ == "SNP":
                                snp_recall.append(float(truth_tp) / float(truth_tot))
                                if int(truth_tp) + int(query_fp) == 0:
                                    snp_prec.append(1)
                                else:
                                    snp_prec.append(float(truth_tp) / 
                                            (float(truth_tp) + float(query_fp)))
                                print(f"{dataset.upper()}\t{sub_id}\t{name.upper()}\tSNP\t{snp_prec[-1]}\t{snp_recall[-1]}", file=out4)
                        else: # normal vcfdist prec/recall calc, already done
                            if typ == "INDEL":
                                indel_recall.append(float(recall))
                                indel_prec.append(float(prec))
                                print(f"{dataset.upper()}\t{sub_id}\t{name.upper()}\tINDEL\t{indel_prec[-1]}\t{indel_recall[-1]}", file=out4)
                            elif typ == "SNP":
                                snp_recall.append(float(recall))
                                snp_prec.append(float(prec))
                                print(f"{dataset.upper()}\t{sub_id}\t{name.upper()}\tSNP\t{snp_prec[-1]}\t{snp_recall[-1]}", file=out4)
                    ax[ei][di*2].plot(snp_recall, snp_prec, linestyle='', 
                            marker='.', markersize=1, color=colors[ni], label=label)
                    ax[ei][di*2+1].plot(indel_recall, indel_prec, linestyle='', 
                            marker='.', markersize=1, color=colors[ni], label=label)
     
        # SNP plot
        ax[ei][0].set_ylabel(r"\textbf{Precision}", fontsize=7)
        ax[3][di*2].set_xlabel(r"\textbf{Recall}", fontsize=7)
        ax[3][di*2+1].set_xlabel(r"\textbf{Recall}", fontsize=7)
        if dataset == "nist":
            ax[ei][di*2].set_xlim(0.99, 1)
            ax[ei][di*2].set_xticks(np.arange(0.99, 1.0002, 0.005))
            ax[ei][di*2].set_xticklabels([f"{x:.3f}" for x in 
                np.arange(0.99, 1.0002, 0.005)], fontsize=5)
            ax[ei][di*2].set_ylim(0.998, 1.00004)
            ax[ei][di*2].set_yticks(np.arange(0.998, 1.00004, 0.0004))
            ax[ei][di*2].set_yticklabels([f"{x:.4f}" for x in 
                np.arange(0.998, 1.00004, 0.0004)], fontsize=5)
        if dataset == "cmrg":
            ax[ei][di*2].set_xlim(0.9, 1)
            ax[ei][di*2].set_xticks(np.arange(0.9, 1.002, 0.02))
            ax[ei][di*2].set_xticklabels([f"{x:.2f}" for x in 
                np.arange(0.9, 1.002, 0.02)], fontsize=5)
            ax[ei][di*2].set_ylim(0.95, 1.001)
            ax[ei][di*2].set_yticks(np.arange(0.95, 1.001, 0.01))
            ax[ei][di*2].set_yticklabels([f"{x:.2f}" for x in 
                np.arange(0.95, 1.001, 0.01)], fontsize=5)
        if di == 0:
            ax[ei][di*2].legend(loc="lower left", fontsize=5, markerscale=4)

        # INDEL plot
        # ax[ei][di*2+1].set_xlabel("Recall", fontsize=15)
        # ax[ei][di*2+1].set_ylabel("Precision", fontsize=15)
        if dataset == "nist":
            ax[ei][di*2+1].set_xlim(0.95, 1)
            ax[ei][di*2+1].set_xticks(np.arange(0.95, 1.001, 0.01))
            ax[ei][di*2+1].set_xticklabels([f"{x:.2f}" for x in 
                np.arange(0.95, 1.001, 0.01)], fontsize=5)
            ax[ei][di*2+1].set_ylim(0.98, 1.0004)
            ax[ei][di*2+1].set_yticks(np.arange(0.98, 1.0004, 0.004))
            ax[ei][di*2+1].set_yticklabels([f"{x:.3f}" for x in 
                np.arange(0.98, 1.0004, 0.004)], fontsize=5)
        if dataset == "cmrg":
            ax[ei][di*2+1].set_xlim(0.9, 1)
            ax[ei][di*2+1].set_xticks(np.arange(0.9, 1.002, 0.02))
            ax[ei][di*2+1].set_xticklabels([f"{x:.2f}" for x in 
                np.arange(0.9, 1.002, 0.02)], fontsize=5)
            ax[ei][di*2+1].set_ylim(0.9, 1.002)
            ax[ei][di*2+1].set_yticks(np.arange(0.9, 1.002, 0.02))
            ax[ei][di*2+1].set_yticklabels([f"{x:.2f}" for x in 
                np.arange(0.9, 1.002, 0.02)], fontsize=5)
        out4.close()

plt.tight_layout()
plt.savefig('img/8_vcfdist_pr_plot.pdf', format="pdf")
