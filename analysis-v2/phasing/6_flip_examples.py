import numpy as np
import matplotlib.pyplot as plt
from random import sample

NONE = 0
FLIP = 1

data_fn = "flips.tsv"
cat_strs = ["NONE", "FLIP"]
datasets = ["hprc", "pav", "giab-tr"]
colors = {"hprc": "Greens", "pav": "Reds", "giab-tr": "Blues"}
names = {"t2t-q100": "Q100-dipcall", "hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR"}

for ds in datasets:
    vcfdist_results = np.zeros((2,2))
    whatshap_results = np.zeros((2,2))
    with open(data_fn, "r") as data:
        next(data) # skip header
        for line in data:
            dataset, vcfdist, whatshap, truth = line.strip().split()
            if dataset == ds:
                vcfdist_results[int(vcfdist == "TRUE")][int(truth == "TRUE")] += 1
                whatshap_results[int(whatshap == "TRUE")][int(truth == "TRUE")] += 1

    # whatshap
    fig, ax = plt.subplots(figsize=(1.55,1.55))
    ax.matshow(np.log(whatshap_results + 0.1), cmap=colors[ds])
    plt.title(f"{names[ds]}\nmanual examination\nof UNKNOWNs", fontsize=7)
    plt.ylabel("WhatsHap", fontsize=7)
    ax.set_yticks(list(range(2)))
    ax.set_yticklabels(cat_strs, fontsize=5)
    plt.xlabel("Truth", fontsize=7)
    ax.set_xticks(list(range(2)))
    ax.set_xticklabels(cat_strs, fontsize=5)
    for (i,j), z in np.ndenumerate(whatshap_results):
        ax.text(j, i, f"{int(z)}", ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'), fontsize=5)
    plt.tight_layout()
    plt.savefig(f"img/{ds}_wh_ex_flip_cm.pdf", format='pdf')

    # vcfdist
    fig, ax = plt.subplots(figsize=(1.33,1.33))
    ax.matshow(np.log(vcfdist_results + 0.1), cmap=colors[ds])
    plt.title(f"                      ", fontsize=7)
    plt.ylabel("vcfdist", fontsize=7)
    ax.set_yticks(list(range(2)))
    ax.set_yticklabels(cat_strs, fontsize=5)
    plt.xlabel("Truth", fontsize=7)
    ax.set_xticks(list(range(2)))
    ax.set_xticklabels(cat_strs, fontsize=5)
    for (i,j), z in np.ndenumerate(vcfdist_results):
        ax.text(j, i, f"{int(z)}", ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'), fontsize=5)
    plt.tight_layout()
    plt.savefig(f"img/{ds}_vd_ex_flip_cm.pdf", format='pdf')
