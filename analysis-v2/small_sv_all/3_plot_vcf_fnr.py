from collections import defaultdict
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

vcfs = ["t2t-q100", "hprc", "pav", "giab-tr"]
names = {"t2t-q100": "Q100-dipcall", "hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR"}
variant_types = ["snp", "indel", "sv"]
vcf_types = ["snp", "indel", "sv", "small", "large", "all"]
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

colors = ["#9e9ac8", "#6aadd5", "#fa6949", "#fc8c3b", "#959595", "#73c375"]
# colors = ["#6950a2", "#2070b4", "#ca171c", "#d64701", "#505050", "#228a44"]

for region in regions:
# initialize counts
    print(region)
    counts = {}
    for vcf in vcfs:
        counts[vcf] = {"snp": defaultdict(int), "indel": defaultdict(int), "sv": defaultdict(int)}

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
                    gt, decision, credit, ga4gh_cat, qual2, sc, sg, ps, pb, bs, fe = truth.split(":")
                    if credit == '.': continue
                    haps = sum([0 if x == '.' else int(x) for x in gt.split('|')])
                    if float(credit) >= 0.7: # TP
                        if len(ref) == len(alt) == 1: # snp
                            counts[vcf]["snp"][typ] += haps
                        elif len(ref) == 1 and len(alt) > 1: # insertion
                            if len(alt) > 50:
                                counts[vcf]["sv"][typ] += haps
                            else:
                                counts[vcf]["indel"][typ] += haps
                        elif len(alt) == 1 and len(ref) > 1: # deletion
                            if len(ref) > 50:
                                counts[vcf]["sv"][typ] += haps
                            else:
                                counts[vcf]["indel"][typ] += haps
                        else:
                            print("ERROR: unexpected variant type")
                    line_ct += 1
                    # if line_ct > 1000: break
    print(json.dumps(counts, indent=4))

    # count separate vs joint FNs
    sep_fns = {"snp": 0, "indel": 0, "sv": 0}
    jnt_fns = {"snp": 0, "indel": 0, "sv": 0}
    for vcf in vcfs[1:]:
        sep_fns["snp"] += counts["t2t-q100"]["snp"]["small"] - counts[vcf]["snp"]["small"]
        sep_fns["indel"] += counts["t2t-q100"]["indel"]["small"] - counts[vcf]["indel"]["small"]
        sep_fns["sv"] += counts["t2t-q100"]["sv"]["sv"] - counts[vcf]["sv"]["sv"]
        jnt_fns["snp"] += counts["t2t-q100"]["snp"]["all"] - counts[vcf]["snp"]["all"]
        jnt_fns["indel"] += counts["t2t-q100"]["indel"]["all"] - counts[vcf]["indel"]["all"]
        jnt_fns["sv"] += counts["t2t-q100"]["sv"]["all"] - counts[vcf]["sv"]["all"]
    print(sep_fns)
    print(jnt_fns)

    fig, ax = plt.subplots(1, 3, figsize=(7,2.5))
    indices = np.arange(len(vcfs)-1)
    width = 0.1
    yquals = [0, 3.01, 6.99, 10, 13.01, 16.99, 20, 23.01, 26.99, 30, # 33.01, 36.99, 40, 43.01, 46.99, 50
            ]
    ylabels = ["0%", "50%", "80%", "90%", "95%", "98%", "99%", "99.5%", "99.8%", "99.9%", # "99.95%", "99.98%", "99.99%", "99.995%", "99.998%", "99.999%"
            ]
    yquals2 = [0, 3.01, 6.99, 10, 13.01, 16.99, 20, 23.01, 26.99, 30] 
    ylabels2 = ["0.1%", "0.2%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%", "100%"]
    for var_type_idx, var_type in enumerate(variant_types):
        for vcf_type_idx, vcf_type in enumerate(vcf_types):

            # skip certain plots (few variants due to complex not getting filtered)
            if var_type == "snp" and vcf_type in ["indel", "large", "sv"]: continue
            if var_type == "indel" and vcf_type in ["snp", "sv"]: continue
            if var_type == "sv" and vcf_type in ["snp", "indel", "small"]: continue

            fn_fracs = [1 - counts[vcf][var_type][vcf_type] / 
                    max(1, counts["t2t-q100"][var_type][vcf_type]) for vcf in vcfs[1:]]
            fn_qscores = [0 if not frac else -10*np.log10(frac) for frac in fn_fracs]

            ax[var_type_idx].bar(indices-2.5*width+vcf_type_idx*width, 
                    [30-x for x in fn_qscores], width, color=colors[vcf_type_idx])
        for yqual in yquals:
            ax[var_type_idx].axhline(y=yqual, color='k', alpha=0.5, linestyle=':', ms=0.5, zorder=-1)
        ax[var_type_idx].set_title(f"{var_type.upper()} evaluation", fontsize=7)
        ax[var_type_idx].set_xlabel("VCFs", fontsize=7)
        ax[var_type_idx].set_xticks(indices)
        ax[var_type_idx].set_xticklabels([names[x] for x in vcfs[1:]], fontsize=5)
        ax[var_type_idx].set_yticks(yquals2)
        ax[var_type_idx].set_yticklabels(ylabels2, fontsize=5)
        ax[var_type_idx].set_ylim(0,30)
    patches = [mpatches.Patch(color=c, label=f"{l.upper()} variants")for c,l in zip(colors, vcf_types)]
    ax[0].legend(handles=patches, loc=('upper left'), fontsize=5)
    ax[0].set_ylabel("False Negative Rate", fontsize=7)
    # plt.suptitle(f"{region}")
    plt.tight_layout()
    plt.savefig(f"./img/{region}-fnr.pdf", format='pdf')
