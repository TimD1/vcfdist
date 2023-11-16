from collections import defaultdict
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

datasets = ["hprc", "pav", "giab-tr"]
query_data, truth_data = {}, {}
for ds in datasets:
    query_data[ds] = defaultdict(int)
    truth_data[ds] = defaultdict(int)
keys = set()
results = np.zeros((len(datasets),2,2), dtype=float)
cmaps = ["Greens", "Reds", "Blues"]
fig, ax = plt.subplots(1, len(datasets), figsize=(len(datasets)*1.5, 1.5))

QUERY = 0
TRUTH = 1
SINGLE = 0
MULTI = 1

TP = 0
PP = 1
FP = 2
FN = 3

for ds_idx, ds in enumerate(datasets):

    results_fn = f"/home/timdunn/vcfdist/analysis-v2/multi_match/evals/vcfdist/{ds}.summary.vcf"
    with open(results_fn, 'r') as results_file:
        line_ct = 0
        for line in results_file:
            if line[0] == "#": # header
                continue
            contig, pos, var_id, ref, alt, qual, filt, info, fmt, truth, query = line.strip().split('\t')
            tGT, tBD, tBC, tBK, tQQ, tSC, tSG, tPB, tPS, tPF = truth.split(":")
            qGT, qBD, qBC, qBK, qQQ, qSC, qSG, qPB, qPS, qPF = query.split(":")

            # flip qGT if necessary, to compare to correct truth variant
            if qGT == "0|1":
                if (int(qPS) + int(qPF)) % 2:
                    qGT = "1|0"
            elif qGT == "1|0":
                if (int(qPS) + int(qPF)) % 2:
                    qGT = "0|1"

            # parse truth variant
            if tBC == ".":
                pass
            elif float(tBC) >= 0.7:
                truth_data[ds][f"{contig}:{tGT}:{tSC}:{tSG}"] += 1
                keys.add(f"{contig}:{tGT}:{tSC}:{tSG}")
                # print(f"Truth TP, {contig}:{tGT}:{tSC}:{tSG}")

            # parse query variant
            if qBC == ".":
                pass
            elif float(qBC) >= 0.7:
                query_data[ds][f"{contig}:{qGT}:{qSC}:{qSG}"] += 1
                keys.add(f"{contig}:{qGT}:{qSC}:{qSG}")
                # print(f"Query TP, {contig}:{qGT}:{qSC}:{qSG}")

            line_ct += 1
            # if line_ct > 1000: break

    # plot results
    for key in keys:
        # not in either
        if truth_data[ds][key] == 0 and query_data[ds][key] == 0:
            pass
        # PP/TP, should be present in both
        elif truth_data[ds][key] > 1 and query_data[ds][key] > 1: # MULTI_MULTI
            results[ds_idx][MULTI][MULTI] += truth_data[ds][key] + query_data[ds][key]
        elif truth_data[ds][key] > 1 and query_data[ds][key] == 1:
            results[ds_idx][SINGLE][MULTI] += truth_data[ds][key] + query_data[ds][key]
        elif truth_data[ds][key] == 1 and query_data[ds][key] > 1:
            results[ds_idx][MULTI][SINGLE] += truth_data[ds][key] + query_data[ds][key]
        elif truth_data[ds][key] == 1 and query_data[ds][key] == 1: # SINGLE-SINGLE
            results[ds_idx][SINGLE][SINGLE] += truth_data[ds][key] + query_data[ds][key]

    props = [100*results[ds] / results[ds].sum(axis=None) for ds in range(len(datasets))]
    ax[ds_idx].matshow(results[ds_idx], cmap=cmaps[ds_idx])
    ax[ds_idx].set_title(f"{ds.upper()}", fontsize=7)
    ax[ds_idx].set_yticks([0, 1])
    ax[ds_idx].set_xticks([0, 1])
    ax[ds_idx].set_yticklabels(["SINGLE", "MULTI"], fontsize=5)
    ax[ds_idx].set_xticklabels(["SINGLE", "MULTI"], fontsize=5)
    ax[ds_idx].set_xlabel("Query", fontsize=7)
    ax[ds_idx].set_ylabel("Truth", fontsize=7)

    for (i,j), z in np.ndenumerate(props[ds_idx]):
        ax[ds_idx].text(j, i, f"{z:.2f}%", ha='center', va='center', fontsize=5)
            # bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'))

plt.tight_layout()
plt.savefig("single_multi_frac.pdf", format="pdf")
