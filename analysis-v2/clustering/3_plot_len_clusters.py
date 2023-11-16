import matplotlib.pyplot as plt
import numpy as np
import statistics

# lens = [10, 50, 100, 500, 1000] # 5000, 10000 didn't finish
lens = [1000]
data = "/home/timdunn/vcfdist/analysis-v2/clustering"

bins = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1_000, 2000,5000, 10_000,20000, 50000, 
        100_000, 200000, 500000, 1_000_000, 2000000, 5000000, 10_000_000, 20000000, 50000000, 
        100_000_000, 200000000, 500000000, 1_000_000_000]
for length in lens:
    print(f"Parsing WFA variant max length {length} superclusters")
    counts = [0]*len(bins)
    bases, bases2 = 0, 0
    regions = []
    with open(f"{data}/sweep_varlens/pav.wfa.len{length}.superclusters.tsv", "r") as bedfile:
        next(bedfile)
        for line in bedfile:
            fields = line.split()
            ctg, sc, start, stop, size = fields[0], int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4])

            for bin_idx, bin_val in enumerate(bins):
                if size < bin_val:
                    counts[bin_idx] += 1
                    regions.append(size)
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
    ax.text(0.5, 0.65, f"median: {statistics.median(regions)}", transform=ax.transAxes)
    ax.text(0.5, 0.6, f"mean: {statistics.mean(regions)}", transform=ax.transAxes)
    ax.text(0.5, 0.55, f"max: {max(regions)}", transform=ax.transAxes)
    plt.savefig(f"img/lens-{length}.png")
