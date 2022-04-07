import matplotlib.pyplot as plt
import numpy as np

bed_filename = "../r10.4_chr20/ont-case-study/input/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
inside, outside = [], []

# read high-conf bed file
with open(bed_filename, 'r') as bed_file:
    last_contig = ""
    last_start = last_stop = 0
    for line in bed_file:

        # parse bed line
        contig, start, stop = line.split()
        start = int(start)
        stop = int(stop)

        inside.append(stop-start)
        if contig == last_contig: # continue
            outside.append(start-last_stop)
        else: # start
            outside.append(start)

        # setup for next iter
        last_start = start
        last_stop = stop
        last_contig = contig

plt.ylabel('Regions (count)')
plt.xlabel('Region size (bases)')
plt.yscale('log')
plt.xscale('log')
plt.hist(inside, bins=np.logspace(0, np.log10(max(inside)), 50), alpha=0.6)
plt.hist(outside, bins=np.logspace(0, np.log10(max(outside)), 50), alpha=0.6)
plt.tight_layout()
plt.legend(['inside high-conf', 'outside high-conf'])
plt.savefig('img/bed_regions.png')
