from collections import defaultdict
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

vcfs = ["t2t-q100", "hprc", "pav", "giab-tr"]
variant_types = ["snv", "indel", "sv"]
vcf_types = ["snv", "indel", "sv", "small", "large", "all"]
regions = [
    "summary",
    # "alldifficultregions",
    # "AllHomopolymers_ge7bp_imperfectge11bp_slop5",
    # "AllTandemRepeats",
    # "AllTandemRepeatsandHomopolymers_slop5",
    # "MHC",
    # "satellites_slop5",
    # "segdups",
    # "alllowmapandsegdupregions"
]

for region in regions:
# initialize fp_counts
    print(region)
    fp_counts = {}
    tp_counts = {}
    for vcf in vcfs:
        fp_counts[vcf] = {"snv": defaultdict(int), "indel": defaultdict(int), "sv": defaultdict(int)}
        tp_counts[vcf] = {"snv": defaultdict(int), "indel": defaultdict(int), "sv": defaultdict(int)}
    colors = ["yellow", "blue", "red", "green", "orange", "purple"]

# count variants
    for vcf in vcfs:
        print(f"    {vcf}")
        for typ in vcf_types:
            results_fn = f"vcfdist/{vcf}.{typ}.{region}.vcf"
            with open(results_fn, 'r') as results_file:
                line_ct = 0
                for line in results_file:
                    if line[0] == "#": # header
                        continue
                    contig, pos, var_id, ref, alt, qual, filt, info, fmt, truth, query = line.strip().split('\t')
                    gt, decision, credit, ga4gh_cat, qual2, sc, swap = query.split(":")
                    if credit == '.': continue
                    haps = sum([0 if x == '.' else int(x) for x in gt.split('|')])
                    if decision == "FP":
                        if len(ref) == len(alt) == 1: # snp
                            fp_counts[vcf]["snv"][typ] += haps
                        elif len(ref) == 1 and len(alt) > 1: # insertion
                            if len(alt) > 50:
                                fp_counts[vcf]["sv"][typ] += haps
                            else:
                                fp_counts[vcf]["indel"][typ] += haps
                        elif len(alt) == 1 and len(ref) > 1: # deletion
                            if len(ref) > 50:
                                fp_counts[vcf]["sv"][typ] += haps
                            else:
                                fp_counts[vcf]["indel"][typ] += haps
                        else:
                            print("ERROR: unexpected variant type")
                    else: # TP
                        if len(ref) == len(alt) == 1: # snp
                            tp_counts[vcf]["snv"][typ] += haps
                        elif len(ref) == 1 and len(alt) > 1: # insertion
                            if len(alt) > 50:
                                tp_counts[vcf]["sv"][typ] += haps
                            else:
                                tp_counts[vcf]["indel"][typ] += haps
                        elif len(alt) == 1 and len(ref) > 1: # deletion
                            if len(ref) > 50:
                                tp_counts[vcf]["sv"][typ] += haps
                            else:
                                tp_counts[vcf]["indel"][typ] += haps
                        else:
                            print("ERROR: unexpected variant type")
                    line_ct += 1
                    # if line_ct > 1000: break
    print(json.dumps(fp_counts, indent=4))

    fig, ax = plt.subplots(1, 3, figsize=(7,2.5))
    indices = np.arange(len(vcfs)-1)
    width = 0.1
    thirty_less_yquals = [0, 3.01, 6.99, 10, 13.01, 16.99, 20, 23.01, 26.99, 30] 
    ylabels = ["0.1%", "0.2%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%", "100%"]
    for var_type_idx, var_type in enumerate(variant_types):
        for vcf_type_idx, vcf_type in enumerate(vcf_types):

            # skip certain plots (few variants due to complex not getting filtered)
            if var_type == "snv" and vcf_type in ["indel", "large", "sv"]: continue
            if var_type == "indel" and vcf_type in ["snv", "sv"]: continue
            if var_type == "sv" and vcf_type in ["snv", "indel", "small"]: continue

            fp = [fp_counts[vcf][var_type][vcf_type]/max(1, tp_counts[vcf][var_type][vcf_type]) for vcf in vcfs[1:]]
            fpq = [0 if not frac else -10*np.log10(frac) for frac in fp]

            ax[var_type_idx].bar(indices-2*width+vcf_type_idx*width, 
                    [30-fpqual for fpqual in fpq], width, color=colors[vcf_type_idx])
        for yqual in thirty_less_yquals:
            ax[var_type_idx].axhline(y=yqual, color='k', alpha=0.5, linestyle=':', ms=0.5, zorder=-1)
        ax[var_type_idx].set_title(f"{var_type.upper()} evaluation", fontsize=7)
        ax[var_type_idx].set_xlabel("VCFs", fontsize=7)
        ax[var_type_idx].set_xticks(indices)
        ax[var_type_idx].set_xticklabels([x.upper() for x in vcfs[1:]], fontsize=5)
        ax[var_type_idx].set_yticks(thirty_less_yquals)
        ax[var_type_idx].set_yticklabels(ylabels, fontsize=5)
    patches = [mpatches.Patch(color=c, label=f"{l.upper()} variants") 
            for c,l in zip(colors, vcf_types)]
    ax[0].legend(handles=patches, loc=(0.05,0.6), fontsize=5)
    ax[0].set_ylabel("False Positive Rate", fontsize=7)
    # plt.suptitle(f"{region}")
    plt.tight_layout()
    plt.savefig(f"img/counts-{region}-fp.pdf", format='pdf')
