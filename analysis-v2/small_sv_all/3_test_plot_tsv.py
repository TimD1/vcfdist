from collections import defaultdict
import json
import matplotlib.pyplot as plt
import numpy as np

vcfs = ["t2t-q100", "hprc", "pav", "giab-tr"]
variant_types = ["snv", "indel", "sv"]
vcf_types = ["snv", "indel", "sv", "small", "large", "all"]

# initialize counts
counts = {}
partial = {}
for vcf in vcfs:
    counts[vcf] = {"snv": defaultdict(int), "indel": defaultdict(int), "sv": defaultdict(int)}
    partial[vcf] = {"snv": defaultdict(int), "indel": defaultdict(int), "sv": defaultdict(int)}
colors = ["yellow", "blue", "red", "green", "orange", "purple"]

# count variants
for vcf in vcfs:
    for typ in vcf_types:
        print(f"{vcf} {typ}")
        results_fn = f"vcfdist/{vcf}.{typ}.query.tsv"
        with open(results_fn, 'r') as results_file:
            next(results_file) # skip header
            line_ct = 0
            for line in results_file:
                contig, pos, hap, ref, alt, qual, var_type, err_type, credit, clust, sc, loc = line.strip().split('\t')
                if err_type == "TP":
                    if len(ref) == len(alt) == 1 and var_type == "SNP": # snp
                        counts[vcf]["snv"][typ] += 1
                    elif len(ref) == 0 and var_type == "INS": # insertion
                        if len(alt) >= 50:
                            counts[vcf]["sv"][typ] += 1
                        else:
                            counts[vcf]["indel"][typ] += 1
                    elif len(alt) == 0 and var_type == "DEL": # deletion
                        if len(ref) >= 50:
                            counts[vcf]["sv"][typ] += 1
                        else:
                            counts[vcf]["indel"][typ] += 1
                    else:
                        print("ERROR: unexpected variant type")
                if err_type == "PP":
                    if len(ref) == len(alt) == 1 and var_type == "SNP": # snp
                        partial[vcf]["snv"][typ] += 1
                    elif len(ref) == 0 and var_type == "INS": # insertion
                        if len(alt) >= 50:
                            partial[vcf]["sv"][typ] += 1
                        else:
                            partial[vcf]["indel"][typ] += 1
                    elif len(alt) == 0 and var_type == "DEL": # deletion
                        if len(ref) >= 50:
                            partial[vcf]["sv"][typ] += 1
                        else:
                            partial[vcf]["indel"][typ] += 1
                    else:
                        print("ERROR: unexpected variant type")
                line_ct += 1
                # if line_ct > 1000: break
print(json.dumps(counts, indent=4))

fig, ax = plt.subplots(1, 3, figsize=(15,5))
indices = np.arange(len(vcfs))
width = 0.1
for var_type_idx, var_type in enumerate(variant_types):
    for vcf_type_idx, vcf_type in enumerate(vcf_types):
        ax[var_type_idx].bar(indices-2*width+vcf_type_idx*width, 
                [counts[vcf][var_type][vcf_type] for vcf in vcfs], width, color=colors[vcf_type_idx])
        ax[var_type_idx].bar(indices-2*width+vcf_type_idx*width, 
                [partial[vcf][var_type][vcf_type] for vcf in vcfs], width, color=colors[vcf_type_idx],
                bottom=[counts[vcf][var_type][vcf_type] for vcf in vcfs], alpha=0.5)
    ax[var_type_idx].set_title(var_type)
    ax[var_type_idx].set_ylabel("Counts")
    ax[var_type_idx].set_xlabel("VCFs")
    ax[var_type_idx].set_xticks(indices)
    ax[var_type_idx].set_xticklabels(vcfs)
ax[0].set_ylim(0, 6_000_000)
ax[1].set_ylim(0, 1_400_000)
ax[2].set_ylim(0, 35_000)
fig.tight_layout()
ax[2].legend(vcf_types, loc=(0.6,0.6))
leg = ax[2].get_legend()
for i in range(len(colors)):
    leg.legendHandles[i].set_color(colors[i])
plt.savefig("counts.png", dpi=300)
