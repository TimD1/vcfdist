import matplotlib.pyplot as plt
import numpy as np

gaps = [10, 20, 50, 100, 200] # 500, 1000 didn't finish
data = "/home/timdunn/vcfdist/analysis-v2/clustering"

bins = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1_000, 2000,5000, 10_000,20000, 50000, 
        100_000, 200000, 500000, 1_000_000, 2000000, 5000000, 10_000_000, 20000000, 50000000, 
        100_000_000, 200000000, 500000000, 1_000_000_000]
for gap in gaps:
    print(f"parsing gap {gap} superclusters")
    counts = [0]*len(bins)
    bases, bases2 = 0, 0
    with open(f"{data}/gaps/{gap}.superclusters.tsv", "r") as bedfile:
        next(bedfile)
        for line in bedfile:
            fields = line.split()
            ctg, start, stop, size = fields[0], int(fields[1]), int(fields[2]), int(fields[3])

            for bin_idx, bin_val in enumerate(bins):
                if size < bin_val:
                    counts[bin_idx] += 1
                    bases += size
                    bases2 += size*size
                    break
    fig, ax = plt.subplots(1,1)
    ax.bar(np.arange(-0.5, len(bins)-0.5), counts, width=1)
    ax.set_xlim(0,len(bins))
    ax.set_xticks(range(0, len(bins), 3))
    ax.set_xticklabels([1, 10, 100, "1K", "10K", "100K", "1M", "10M", "100M", "1G"])
    ax.set_yscale('log')
    ax.set_ylim(0.5, 20000000)
    ax.set_xlabel("Region Size")
    ax.set_ylabel("Counts")
    ax.text(0.5, 0.8, f"$\sum$bases: {bases:,}", transform=ax.transAxes)
    ax.text(0.5, 0.75, f"$\sum$bases$^2$: {bases2:,}", transform=ax.transAxes)
    ax.text(0.5, 0.7, f"regions: {sum(counts):,}", transform=ax.transAxes)
    plt.savefig(f"img/gaps-{gap}.png")
