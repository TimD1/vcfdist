from collections import defaultdict
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

query_data, truth_data = {}, {}
BKs = ["TP", "PP", "FP", "FN"]
for BK in BKs:
    query_data[BK] = defaultdict(int)
    truth_data[BK] = defaultdict(int)
keys = set()

QUERY = 0
TRUTH = 1
SINGLE = 0
MULTI = 1

TP = 0
PP = 1
FP = 2
FN = 3
results = np.zeros((4,2,2), dtype=float)

results_fn = f"/home/timdunn/vcfdist/analysis-v2/multi_match/evals/vcfdist/giab-tr.summary.vcf"
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
        elif float(tBC) == 0:
            truth_data["FN"][f"{contig}:{tGT}:{tSC}:{tSG}"] += 1
            keys.add(f"{contig}:{tGT}:{tSC}:{tSG}")
            # print(f"Truth FN, {contig}:{tGT}:{tSC}:{tSG}")
        elif float(tBC) < 1:
            truth_data["PP"][f"{contig}:{tGT}:{tSC}:{tSG}"] += 1
            keys.add(f"{contig}:{tGT}:{tSC}:{tSG}")
            # print(f"Truth PP, {contig}:{tGT}:{tSC}:{tSG}")
        elif float(tBC) == 1:
            truth_data["TP"][f"{contig}:{tGT}:{tSC}:{tSG}"] += 1
            keys.add(f"{contig}:{tGT}:{tSC}:{tSG}")
            # print(f"Truth TP, {contig}:{tGT}:{tSC}:{tSG}")
        else:
            print("ERROR: unexpected tBC =", tBC)

        # parse query variant
        if qBC == ".":
            pass
        elif float(qBC) == 0:
            query_data["FP"][f"{contig}:{qGT}:{qSC}:{qSG}"] += 1
            keys.add(f"{contig}:{qGT}:{qSC}:{qSG}")
            # print(f"Query FP, {contig}:{qGT}:{qSC}:{qSG}")
        elif float(qBC) < 1:
            query_data["PP"][f"{contig}:{qGT}:{qSC}:{qSG}"] += 1
            keys.add(f"{contig}:{qGT}:{qSC}:{qSG}")
            # print(f"Query PP, {contig}:{qGT}:{qSC}:{qSG}")
        elif float(qBC) == 1:
            query_data["TP"][f"{contig}:{qGT}:{qSC}:{qSG}"] += 1
            keys.add(f"{contig}:{qGT}:{qSC}:{qSG}")
            # print(f"Query TP, {contig}:{qGT}:{qSC}:{qSG}")
        else:
            print("ERROR: unexpected qBC =", qBC)

        line_ct += 1
        # if line_ct > 1000: break

# plot results
for key in keys:
    for bk, BK in enumerate(BKs):
        # not in either
        if truth_data[BK][key] == 0 and query_data[BK][key] == 0:
            pass
        # PP/TP, should be present in both
        elif truth_data[BK][key] > 1 and query_data[BK][key] > 1:
            results[bk][MULTI][MULTI] += truth_data[BK][key] + query_data[BK][key]
        elif truth_data[BK][key] > 1 and query_data[BK][key] == 1:
            results[bk][SINGLE][MULTI] += truth_data[BK][key] + query_data[BK][key]
        elif truth_data[BK][key] == 1 and query_data[BK][key] > 1:
            results[bk][MULTI][SINGLE] += truth_data[BK][key] + query_data[BK][key]
        elif truth_data[BK][key] == 1 and query_data[BK][key] == 1:
            results[bk][SINGLE][SINGLE] += truth_data[BK][key] + query_data[BK][key]
        else: # present in one
            if BK == "FP":
                results[bk][SINGLE][SINGLE] += 1
                if truth_data[BK][key] > 0:
                    print("ERROR: FP key present in truth")
                    print(BK, key)

            elif BK == "FN":
                results[bk][SINGLE][SINGLE] += 1
                if query_data[BK][key] > 0:
                    print("ERROR: FN key present in query")
                    print(BK, key)
            else:
                print("ERROR: PP/TP key not present for truth and query")
                print(BK, key)

fracs = [100*results[x].sum(axis=None) / results.sum(axis=None) for x in range(4)]
props = [100*results[x] / results[x].sum(axis=None) for x in range(4)]
# results = 100 * results / results.sum(axis=None)
# print(results)
cmaps = ["Greens", "Oranges", "Reds", "Blues"]
fig, ax = plt.subplots(1, 5, figsize=(15,4))
for x in range(4):
    ax[x].matshow(results[x], cmap=cmaps[x])
    ax[x].set_title(f"{BKs[x]} ({fracs[x]:.3f}%)")
    ax[x].set_yticks([0, 1])
    ax[x].set_xticks([0, 1])
    ax[x].set_yticklabels(["SINGLE", "MULTI"])
    ax[x].set_xticklabels(["SINGLE", "MULTI"])
    ax[x].set_xlabel("Query")
    ax[x].set_ylabel("Truth")

    for (i,j), z in np.ndenumerate(props[x]):
        ax[x].text(j, i, f"{z:.3f}%", ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'))

ax[4].set_title("Variant Calls")
ax[4].pie(fracs, labels=BKs, autopct='%1.1f%%', explode=(0, 0.6, 0.3, 0),
        colors=['green', 'orange', 'red', 'blue'])
plt.tight_layout()
plt.savefig("results.png")
