import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_right

vcf_fn = "/home/timdunn/vcfdist/data/t2t-q100-v0.9/split/all.vcf"

bins = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1_000, 2000,5000, 10_000,20000, 50000, 
        100_000, 200000, 500000, 1_000_000, 2000000, 5000000, 10_000_000, 20000000, 50000000, 
        100_000_000, 200000000, 500000000, 1_000_000_000]
labels = [1, 10, 100, "1K", "10K", "100K", "1M", "10M", "100M", "1G"]
all_bins = [-x for x in bins[::-1]] + [0] + bins
all_labels = [f"-{l}" for l in labels[::-1]] + [0] + labels

counts = [0]*len(all_bins)
with open(vcf_fn, "r") as vcf:
    for line in vcf:
        if line[0] == "#": continue
        fields = line.split()
        ctg, pos, uid, ref, alt = fields[:5]
        varlen = len(alt) - len(ref)
        idx = bisect_right(all_bins, varlen)
        counts[idx] += 1

fig, ax = plt.subplots(1,1, figsize=(15,5))
for x in range(56):
    plt.axvline(x, color='k', linestyle=':')

ax.bar(np.arange(-0.5, len(all_bins)-0.5), counts, width=1)
ax.set_xlim(0,len(all_bins))
ax.set_xticks([0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 28, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56])
ax.set_xticklabels(all_labels)
ax.set_yscale('log')
ax.set_ylim(0.5, 2000000)
ax.set_xlabel("Variant Size")
ax.set_ylabel("Counts")
plt.savefig(f"img/var-lengths.png")
