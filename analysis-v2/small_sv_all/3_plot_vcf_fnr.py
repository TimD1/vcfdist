from collections import defaultdict
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

vcfs = ["t2t-q100", "hprc", "pav", "giab-tr"]
names = {"t2t-q100": "Q100-dipcall", "hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR"}
variant_types = ["snp", "indel", "sv"]
vcf_types = ["sv", "small", "all"]
regions = ["summary"]
colors = ["#fa6949", "#9e9ac8", "#73c375"]

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
                    gt, decision, credit, rd, qd, ga4gh_cat, qual2, sc, sg, ps, pb, bs, fe = truth.split(":")
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

    # save truth results
    tp_counts_json = json.dumps(counts, indent=4)
    with open(f"{region}-truth-tp.json", "w") as json_file:
        json_file.write(tp_counts_json)

    # FALSE NEGATIVE RATE PLOT
    fig, ax = plt.subplots(1, 3, figsize=(7,2.5))
    indices = np.arange(len(vcfs)-1)
    width = 0.25
    yquals = [0, 3.01, 6.99, 10, 13.01, 16.99, 20, 23.01, 26.99, 30]
    ylabels = ["0.1%", "0.2%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%", "100%"]
    for var_type_idx, var_type in enumerate(variant_types):
        for vcf_type_idx, vcf_type in enumerate(vcf_types):

            # skip certain plots (few variants due to complex not getting filtered)
            if var_type == "snp" and vcf_type in ["indel", "large", "sv"]: continue
            if var_type == "indel" and vcf_type in ["snp", "sv"]: continue
            if var_type == "sv" and vcf_type in ["snp", "indel", "small"]: continue

            fn_fracs = [1 - counts[vcf][var_type][vcf_type] / 
                    max(1, counts["t2t-q100"][var_type][vcf_type]) for vcf in vcfs[1:]]
            fn_qscores = [0 if not frac else -10*np.log10(frac) for frac in fn_fracs]

            ax[var_type_idx].bar(indices+(width/2 if vcf_type == "all" else -width/2), 
                    [30-x for x in fn_qscores], width, color=colors[vcf_type_idx])
        for yqual in yquals:
            ax[var_type_idx].axhline(y=yqual, color='k', alpha=0.5, linestyle=':', ms=0.5, zorder=-1)
        ax[var_type_idx].set_title(f"{var_type.upper()} evaluation", fontsize=7)
        ax[var_type_idx].set_xticks(indices)
        ax[var_type_idx].set_xticklabels([names[x] for x in vcfs[1:]], fontsize=5)
        ax[var_type_idx].set_yticks(yquals)
        ax[var_type_idx].set_yticklabels(ylabels, fontsize=5)
        ax[var_type_idx].set_ylim(0,30)
    patches = [mpatches.Patch(color=c, label=f"{l.upper()} variants")for c,l in zip(colors, vcf_types)]
    ax[0].legend(handles=patches, loc=('upper left'), fontsize=5)
    ax[0].set_ylabel("False Negative Rate", fontsize=7)
    plt.tight_layout()
    plt.savefig(f"./img/{region}-fnr.pdf", format='pdf')
