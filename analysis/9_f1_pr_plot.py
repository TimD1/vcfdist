import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import scipy.stats
import matplotlib as mpl
import statistics
import bisect
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['text.usetex'] = True
plt.rcParams.update({"figure.facecolor": (0,0,0,0)})

data = "/home/timdunn/vcfdist/data"
datasets = ["nist", "cmrg"]
tools = ["vcfeval", "vcfdist"]
names = ["original", "A", "B", "C", "D"]
colors = ["k", "C0", "C1", "C2", "C3"]

# sub_ids = [
#     "0GOOR", "23O09", "4GKR1", "61YRJ", "8H0ZB", "B1S5A", "C6JUX", "EIUT6", "H9OJ3", "IA789",
#     "KXBR8", "MECML", "NWQ6Y", "R9CXN", "SEX9X", "UYMUW", "W91C1", "XC97E", "YBE9U", "YUI27",
#     "4HL0B", "7NJ5C", "9DGOR", "BARQS", "CN36L", "ES3XW", "HB8P3", "ISKOS", "JIPSI", "KFA0T",
#     "PGXA4", "RU88N", "TG5TE", "VES2R", "WGQ43", "XV7ZN", "YGOTK", "13678", "32LOW", "60Z59", 
#     "JBBK0", "K4GT3", "7RR8Z", "ASJT6", "BSODP", "CZA1Y", "0O7FL", "2OT9Q", "FFFGB", "HF8CT", "Y8QR8", "YJN61", "LR1FD", "MT57N",
#     "J04SL", "K33QJ", "KPXT2", "M9KLP", "NFT0L", "QUE7Q", "S7K7S", "TZMTX", "W607K", "WX8VK"
# ]
# markers = ['.']*len(sub_ids)

sub_ids = [
    "K4GT3", # Google Health
    "W91C1", # Sentieon
    "IA789"  # Roche
]
markers = ['s', 'o', '^']

fig, ax = plt.subplots(2, 4, figsize=(20,8))

for di, dataset in enumerate(datasets):
    for ti, tool in enumerate(tools):

        # VCFEVAL
        if tool == "vcfeval":
            all_snp_f1_x = []
            all_snp_f1_y = []
            all_indel_f1_x = []
            all_indel_f1_y = []
            all_snp_f1_meds = []
            all_indel_f1_meds = []
            for si, sub_id in enumerate(sub_ids):
                print(f"\nvcfeval {dataset} {sub_id}: ", end="")
                snp_f1_list = []
                indel_f1_list = []
                for ni, name in enumerate(names):
                    filename = "O" if name == "original" else name
                    try:
                        with open(f"{data}/pfda-v2/{dataset}_vcfeval/{sub_id}_HG002_{filename}.summary.csv") as csv:
                            label = f"{name}"
                            lines = [line.rstrip() for line in csv]
                            snp_line = lines[4]
                            indel_line = lines[2]
                            try:
                                snp_f1 = float(snp_line.split(",")[13])
                            except ValueError:
                                snp_f1 = 0
                            try:
                                indel_f1 = float(indel_line.split(",")[13])
                            except ValueError:
                                indel_f1 = 0
                            snp_f1q = -10 * np.log10(1-snp_f1)
                            indel_f1q = -10 * np.log10(1-indel_f1)
                            snp_f1_list.append(snp_f1q)
                            indel_f1_list.append(indel_f1q)
                            print(f"{filename}={snp_f1q:6.4f},{indel_f1q:6.4f}\t", end="")
                    except FileNotFoundError:
                        print(f"File '{data}/pfda-v2/{dataset}_vcfeval/{sub_id}_HG002_{filename}.summary.csv' not found.")
                        continue
                if len(snp_f1_list) != len(names): continue
                snp_f1_mean = sum(snp_f1_list) / len(snp_f1_list)
                indel_f1_mean = sum(indel_f1_list) / len(indel_f1_list)
                for ni, name in enumerate(names):
                    ax[ti][di*2].plot(snp_f1_mean, snp_f1_list[ni], linestyle='', 
                            marker=markers[si], color=colors[ni], label=label)
                    ax[ti][di*2+1].plot(indel_f1_mean, indel_f1_list[ni], linestyle='', 
                            marker=markers[si], color=colors[ni], label=label)
                all_snp_f1_y.extend(snp_f1_list)
                all_snp_f1_x.extend([snp_f1_mean]*len(names))
                all_snp_f1_meds.append(statistics.median(snp_f1_list))
                all_indel_f1_y.extend(indel_f1_list)
                all_indel_f1_x.extend([indel_f1_mean]*len(names))
                all_indel_f1_meds.append(statistics.median(indel_f1_list))

            m, b, r, p, std_err = scipy.stats.linregress(all_snp_f1_x, all_snp_f1_y)
            ax[ti][di*2].set_xlabel("Avg SNP F1 Q-score", fontsize=15)
            ax[ti][di*2].set_ylabel("SNP F1 Q-score", fontsize=15)
            ax[ti][di*2].text(0.05,0.9,"$R^2$" +f": {r*r:.6f}", fontsize=15, transform=ax[ti][di*2].transAxes)

            m, b, r, p, std_err = scipy.stats.linregress(all_indel_f1_x, all_indel_f1_y)
            ax[ti][di*2+1].set_xlabel("Avg INDEL F1 Q-score", fontsize=15)
            ax[ti][di*2+1].set_ylabel("INDEL F1 Q-score", fontsize=15)
            ax[ti][di*2+1].text(0.05,0.9,"$R^2$" +f": {r*r:.6f}", fontsize=15, transform=ax[ti][di*2+1].transAxes)

            # second pass, calculate average max rank change
            all_snp_f1_meds.sort()
            all_indel_f1_meds.sort()
            snp_rank_change = []
            indel_rank_change = []
            for si, sub_id in enumerate(sub_ids):
                snp_f1_list = []
                indel_f1_list = []
                for ni, name in enumerate(names):
                    filename = "O" if name == "original" else name
                    try:
                        with open(f"{data}/pfda-v2/{dataset}_vcfeval/{sub_id}_HG002_{filename}.summary.csv") as csv:
                            label = f"{name}"
                            lines = [line.rstrip() for line in csv]
                            snp_line = lines[4]
                            indel_line = lines[2]
                            try:
                                snp_f1 = float(snp_line.split(",")[13])
                            except ValueError:
                                snp_f1 = 0
                            try:
                                indel_f1 = float(indel_line.split(",")[13])
                            except ValueError:
                                indel_f1 = 0
                            snp_f1q = -10 * np.log10(1-snp_f1)
                            indel_f1q = -10 * np.log10(1-indel_f1)
                            snp_f1_list.append(snp_f1q)
                            indel_f1_list.append(indel_f1q)
                    except FileNotFoundError:
                        break
                if len(snp_f1_list) != len(names): continue
                # remove this median from list before doing rank calcs
                snp_f1_meds = all_snp_f1_meds[:]
                snp_f1_med = statistics.median(snp_f1_list)
                snp_f1_meds.remove(snp_f1_med)
                snp_min_rank = bisect.bisect_left(snp_f1_meds, min(snp_f1_list))
                snp_max_rank = bisect.bisect_left(snp_f1_meds, max(snp_f1_list))
                snp_rank_change.append(snp_max_rank - snp_min_rank)
                indel_f1_meds = all_indel_f1_meds[:]
                indel_f1_med = statistics.median(indel_f1_list)
                indel_f1_meds.remove(indel_f1_med)
                indel_min_rank = bisect.bisect_left(indel_f1_meds, min(indel_f1_list))
                indel_max_rank = bisect.bisect_left(indel_f1_meds, max(indel_f1_list))
                indel_rank_change.append(indel_max_rank - indel_min_rank)
            snp_amrc = np.mean(snp_rank_change)
            indel_amrc = np.mean(indel_rank_change)
            ax[ti][di*2].text(0.05,0.8,f"AMRC: {snp_amrc:.2f}", fontsize=15, transform=ax[ti][di*2].transAxes)
            ax[ti][di*2+1].text(0.05,0.8,f"AMRC: {indel_amrc:.2f}", fontsize=15, transform=ax[ti][di*2+1].transAxes)

        # VCFDIST
        elif tool == "vcfdist":
            all_snp_f1_x = []
            all_snp_f1_y = []
            all_indel_f1_x = []
            all_indel_f1_y = []
            all_snp_f1_meds = []
            all_indel_f1_meds = []
            for si, sub_id in enumerate(sub_ids):
                print(f"\nvcfdist {dataset} {sub_id}: ", end="")
                snp_f1_list = []
                indel_f1_list = []
                for ni, name in enumerate(names):
                    filename = "O" if name == "original" else name
                    try:
                        with open(f"{data}/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}.precision-recall-summary.tsv") as tsv:
                            label = f"{name}"
                            lines = [line.rstrip() for line in tsv]
                            snp_line = lines[4]
                            indel_line = lines[2]
                            snp_f1q = float(snp_line.split("\t")[9])
                            indel_f1q = float(indel_line.split("\t")[9])
                            snp_f1_list.append(snp_f1q)
                            indel_f1_list.append(indel_f1q)
                            print(f"{filename}={snp_f1q:6.4f},{indel_f1q:6.4f}\t", end="")
                    except FileNotFoundError:
                        print(f"File '{data}/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}.precision-recall-summary.tsv' not found.")
                        continue
                if len(snp_f1_list) != len(names): continue
                snp_f1_mean = sum(snp_f1_list) / len(snp_f1_list)
                indel_f1_mean = sum(indel_f1_list) / len(indel_f1_list)
                for ni, name in enumerate(names):
                    ax[ti][di*2].plot(snp_f1_mean, snp_f1_list[ni], linestyle='', 
                            marker=markers[si], color=colors[ni], label=label)
                    ax[ti][di*2+1].plot(indel_f1_mean, indel_f1_list[ni], linestyle='', 
                            marker=markers[si], color=colors[ni], label=label)
                all_snp_f1_y.extend(snp_f1_list)
                all_snp_f1_x.extend([snp_f1_mean]*len(names))
                all_snp_f1_meds.append(statistics.median(snp_f1_list))
                all_indel_f1_y.extend(indel_f1_list)
                all_indel_f1_x.extend([indel_f1_mean]*len(names))
                all_indel_f1_meds.append(statistics.median(indel_f1_list))

            m, b, r, p, std_err = scipy.stats.linregress(all_snp_f1_x, all_snp_f1_y)
            ax[ti][di*2].set_xlabel("Avg SNP F1 Q-score", fontsize=15)
            ax[ti][di*2].set_ylabel("SNP F1 Q-score", fontsize=15)
            ax[ti][di*2].text(0.05,0.9,"$R^2$" +f": {r*r:.6f}", fontsize=15, transform=ax[ti][di*2].transAxes)

            m, b, r, p, std_err = scipy.stats.linregress(all_indel_f1_x, all_indel_f1_y)
            ax[ti][di*2+1].set_xlabel("Avg INDEL F1 Q-score", fontsize=15)
            ax[ti][di*2+1].set_ylabel("INDEL F1 Q-score", fontsize=15)
            ax[ti][di*2+1].text(0.05,0.9,"$R^2$" +f": {r*r:.6f}", fontsize=15, transform=ax[ti][di*2+1].transAxes)

            # second pass, calculate average max rank change
            all_snp_f1_meds.sort()
            all_indel_f1_meds.sort()
            snp_rank_change = []
            indel_rank_change = []
            for si, sub_id in enumerate(sub_ids):
                snp_f1_list = []
                indel_f1_list = []
                for ni, name in enumerate(names):
                    filename = "O" if name == "original" else name
                    try:
                        with open(f"{data}/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}.precision-recall-summary.tsv") as tsv:
                            label = f"{name}"
                            lines = [line.rstrip() for line in tsv]
                            snp_line = lines[4]
                            indel_line = lines[2]
                            snp_f1q = float(snp_line.split("\t")[9])
                            indel_f1q = float(indel_line.split("\t")[9])
                            snp_f1_list.append(snp_f1q)
                            indel_f1_list.append(indel_f1q)
                    except FileNotFoundError:
                        break
                if len(snp_f1_list) != len(names): continue
                # remove this median from list before doing rank calcs
                snp_f1_meds = all_snp_f1_meds[:]
                snp_f1_med = statistics.median(snp_f1_list)
                snp_f1_meds.remove(snp_f1_med)
                snp_min_rank = bisect.bisect_left(snp_f1_meds, min(snp_f1_list))
                snp_max_rank = bisect.bisect_left(snp_f1_meds, max(snp_f1_list))
                snp_rank_change.append(snp_max_rank - snp_min_rank)
                indel_f1_meds = all_indel_f1_meds[:]
                indel_f1_med = statistics.median(indel_f1_list)
                indel_f1_meds.remove(indel_f1_med)
                indel_min_rank = bisect.bisect_left(indel_f1_meds, min(indel_f1_list))
                indel_max_rank = bisect.bisect_left(indel_f1_meds, max(indel_f1_list))
                indel_rank_change.append(indel_max_rank - indel_min_rank)
            snp_amrc = np.mean(snp_rank_change)
            indel_amrc = np.mean(indel_rank_change)
            ax[ti][di*2].text(0.05,0.8,f"AMRC: {snp_amrc:.2f}", fontsize=15, transform=ax[ti][di*2].transAxes)
            ax[ti][di*2+1].text(0.05,0.8,f"AMRC: {indel_amrc:.2f}", fontsize=15, transform=ax[ti][di*2+1].transAxes)

square = mlines.Line2D([], [], linestyle="", color='k', marker='s',
                                  markersize=12, label='Google Health')
circle = mlines.Line2D([], [], linestyle="", color='k', marker='o',
                                  markersize=12, label='Sentieon')
triangle = mlines.Line2D([], [], linestyle="", color='k', marker='^',
                                  markersize=12, label='Roche')
ax[0][0].legend(handles=[square, circle, triangle], loc="center right", fontsize=15)

O = mlines.Line2D([], [], linestyle="", color='k', marker='.',
                                  markersize=12, label='original')
A = mlines.Line2D([], [], linestyle="", color='C0', marker='.',
                                  markersize=12, label='A')
B = mlines.Line2D([], [], linestyle="", color='C1', marker='.',
                                  markersize=12, label='B')
C = mlines.Line2D([], [], linestyle="", color='C2', marker='.',
                                  markersize=12, label='C')
D = mlines.Line2D([], [], linestyle="", color='C3', marker='.',
                                  markersize=12, label='D')
ax2 = ax[0][0].twinx()

# ax2.legend(handles=[O,A,B,C,D], loc="center right", fontsize=15)
ax2.legend(handles=[O,A,B,C,D], loc="center left", fontsize=15)


plt.tight_layout()
plt.savefig('img/9_f1_pr_plot.pdf', format="pdf")
