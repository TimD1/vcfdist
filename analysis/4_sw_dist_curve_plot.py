import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

names = ["A", "B", "C", "D"]
quals = list(range(11))

plt.figure(figsize=(15,10))

truth = "C"
for calls in names:
    dists = []
    for q in quals:
        with open(f"4_results/{calls}-calls_C-truth_{q}_distance.tsv") as csv:
            label = f"VCFdist Calls={calls} Truth={truth}"
            for line in csv:
                dist = int(line)
                dists.append(dist)
                break
    plt.plot(quals, dists, linestyle='', marker='.', label=label)
 
plt.xlabel("Quality Threshold")
plt.ylabel("Edit Distance")
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig('img/4_sw_dist_curve.png')
