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
plt.rcParams.update({"figure.facecolor": (0,0,0,0)})
mpl.rcParams['text.usetex'] = True

data = "/home/timdunn/vcfdist/data"
datasets = ["nist", "cmrg"]
names = ["original", "A", "B", "C", "D"]
colors = ["k", "C0", "C1", "C2", "C3"]

sub_ids = [
    "0GOOR", "23O09", "4GKR1", "61YRJ", "8H0ZB", "B1S5A", "C6JUX", "EIUT6", "H9OJ3", "IA789",
    "KXBR8", "MECML", "NWQ6Y", "R9CXN", "SEX9X", "UYMUW", "W91C1", "XC97E", "YBE9U", "YUI27",
    "4HL0B", "7NJ5C", "9DGOR", "BARQS", "CN36L", "ES3XW", "HB8P3", "ISKOS", "JIPSI", "KFA0T",
    "PGXA4", "RU88N", "TG5TE", "VES2R", "WGQ43", "XV7ZN", "YGOTK", "13678", "32LOW", "60Z59", 
    "JBBK0", "K4GT3", "7RR8Z", "ASJT6", "BSODP", "CZA1Y", "0O7FL", "2OT9Q", "FFFGB", "HF8CT", "Y8QR8", "YJN61", "LR1FD", "MT57N",
    "J04SL", "K33QJ", "KPXT2", "M9KLP", "NFT0L", "QUE7Q", "S7K7S", "TZMTX", "W607K", "WX8VK"
]
markers = ['.']*len(sub_ids)

# sub_ids = [
#     "K4GT3", # Google Health
#     "W91C1", # Sentieon
#     "IA789"  # Roche
# ]
# markers = ['s', 'o', '^']

fig, ax = plt.subplots(1, 6, figsize=(20,4))

for di, dataset in enumerate(datasets):

    # VCFDIST ED/DE SUMMARY METRICS
    all_ed_f1_x = []
    all_ed_f1_y = []
    all_ed_f1_meds = []
    all_de_f1_x = []
    all_de_f1_y = []
    all_de_f1_meds = []
    all_aln_f1_x = []
    all_aln_f1_y = []
    all_aln_f1_meds = []
    for si, sub_id in enumerate(sub_ids):
        print(f"\n{dataset} {sub_id}: ", end="")
        ed_f1_list = []
        de_f1_list = []
        aln_f1_list = []
        for ni, name in enumerate(names):
            filename = "O" if name == "original" else name
            try:
                with open(f"{data}/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}.distance-summary.tsv") as tsv:
                    label = f"{name}"
                    lines = [line.rstrip() for line in tsv]
                    line = lines[2]
                    ed_f1q = float(line.split("\t")[4])
                    de_f1q = float(line.split("\t")[5])
                    aln_f1q = float(line.split("\t")[6])
                    ed_f1_list.append(ed_f1q)
                    de_f1_list.append(de_f1q)
                    aln_f1_list.append(aln_f1q)
                    print(f"{filename}={ed_f1q:.4f}\t", end="")
            except FileNotFoundError:
                # print(f"File '{data}/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}.distance-summary.tsv' not found.")
                break
            except IndexError:
                print(f"File '{data}/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}.distance-summary.tsv' empty.")
                break
        if len(ed_f1_list) != 5: continue
        ed_f1_mean = np.mean(ed_f1_list)
        de_f1_mean = np.mean(de_f1_list)
        aln_f1_mean = np.mean(aln_f1_list)
        for ni, name in enumerate(names):
            ax[di*3].plot(ed_f1_mean, ed_f1_list[ni], linestyle='', 
                    marker=markers[si], color=colors[ni], label=label)
            ax[di*3+1].plot(de_f1_mean, de_f1_list[ni], linestyle='', 
                    marker=markers[si], color=colors[ni], label=label)
            ax[di*3+2].plot(aln_f1_mean, aln_f1_list[ni], linestyle='', 
                    marker=markers[si], color=colors[ni], label=label)
        all_ed_f1_y.extend(ed_f1_list)
        all_ed_f1_x.extend([ed_f1_mean]*5)
        all_ed_f1_meds.append(statistics.median(ed_f1_list))
        all_de_f1_y.extend(de_f1_list)
        all_de_f1_x.extend([de_f1_mean]*5)
        all_de_f1_meds.append(statistics.median(de_f1_list))
        all_aln_f1_y.extend(aln_f1_list)
        all_aln_f1_x.extend([aln_f1_mean]*5)
        all_aln_f1_meds.append(statistics.median(aln_f1_list))

    print(" ")
    m, b, r, p, std_err = scipy.stats.linregress(all_ed_f1_x, all_ed_f1_y)
    ax[di*3].set_xlabel("Avg ED Q-score", fontsize=15)
    ax[di*3].set_ylabel("ED Q-score", fontsize=15)
    ax[di*3].text(0.05,0.9,"$R^2$" +f": {r*r:.6f}", fontsize=15, 
            transform=ax[di*3].transAxes)

    m, b, r, p, std_err = scipy.stats.linregress(all_de_f1_x, all_de_f1_y)
    ax[di*3+1].set_xlabel("Avg DE Q-score", fontsize=15)
    ax[di*3+1].set_ylabel("DE Q-score", fontsize=15)
    ax[di*3+1].text(0.05,0.9,"$R^2$" +f": {r*r:.6f}", fontsize=15, 
            transform=ax[di*3+1].transAxes)

    m, b, r, p, std_err = scipy.stats.linregress(all_aln_f1_x, all_aln_f1_y)
    ax[di*3+2].set_xlabel("Avg ALN Q-score", fontsize=15)
    ax[di*3+2].set_ylabel("ALN Q-score", fontsize=15)
    ax[di*3+2].text(0.05,0.9,"$R^2$" +f": {r*r:.6f}", fontsize=15, 
            transform=ax[di*3+2].transAxes)

    # second pass, calculate average max rank change
    all_ed_f1_meds.sort()
    all_de_f1_meds.sort()
    all_aln_f1_meds.sort()
    ed_rank_change = []
    de_rank_change = []
    aln_rank_change = []
    for si, sub_id in enumerate(sub_ids):
        ed_f1_list = []
        de_f1_list = []
        aln_f1_list = []
        for ni, name in enumerate(names):
            filename = "O" if name == "original" else name
            try:
                with open(f"{data}/pfda-v2/{dataset}_vcfdist/{sub_id}_HG002_{filename}.distance-summary.tsv") as tsv:
                    label = f"{name}"
                    lines = [line.rstrip() for line in tsv]
                    line = lines[2]
                    ed_f1q = float(line.split("\t")[4])
                    de_f1q = float(line.split("\t")[5])
                    aln_f1q = float(line.split("\t")[6])
                    ed_f1_list.append(ed_f1q)
                    de_f1_list.append(de_f1q)
                    aln_f1_list.append(aln_f1q)
            except FileNotFoundError:
                break
            except IndexError:
                break
        if len(ed_f1_list) != 5: continue
        # remove this median from list before doing rank calcs
        ed_f1_meds = all_ed_f1_meds[:]
        ed_f1_med = statistics.median(ed_f1_list)
        ed_f1_meds.remove(ed_f1_med)
        ed_min_rank = bisect.bisect_left(ed_f1_meds, min(ed_f1_list))
        ed_max_rank = bisect.bisect_left(ed_f1_meds, max(ed_f1_list))
        ed_rank_change.append(ed_max_rank - ed_min_rank)
        de_f1_meds = all_de_f1_meds[:]
        de_f1_med = statistics.median(de_f1_list)
        de_f1_meds.remove(de_f1_med)
        de_min_rank = bisect.bisect_left(de_f1_meds, min(de_f1_list))
        de_max_rank = bisect.bisect_left(de_f1_meds, max(de_f1_list))
        de_rank_change.append(de_max_rank - de_min_rank)
        aln_f1_meds = all_aln_f1_meds[:]
        aln_f1_med = statistics.median(aln_f1_list)
        aln_f1_meds.remove(aln_f1_med)
        aln_min_rank = bisect.bisect_left(aln_f1_meds, min(aln_f1_list))
        aln_max_rank = bisect.bisect_left(aln_f1_meds, max(aln_f1_list))
        aln_rank_change.append(aln_max_rank - aln_min_rank)
    ed_amrc = np.mean(ed_rank_change)
    de_amrc = np.mean(de_rank_change)
    aln_amrc = np.mean(aln_rank_change)
    ax[di*3].text(0.05,0.8,f"AMRC: {ed_amrc:.2f}", fontsize=15, transform=ax[di*3].transAxes)
    ax[di*3+1].text(0.05,0.8,f"AMRC: {de_amrc:.2f}", fontsize=15, transform=ax[di*3+1].transAxes)
    ax[di*3+2].text(0.05,0.8,f"AMRC: {aln_amrc:.2f}", fontsize=15, transform=ax[di*3+2].transAxes)

plt.tight_layout()
plt.savefig('img/9_f1_ed_plot.pdf', format="pdf")
