from collections import defaultdict
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

datasets = ["hprc", "pav", "giab-tr"]
names = {"hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR", "t2t-q100": "Q100-dipcall"}
tools = ["vcfdist", "vcfeval", "truvari", "truvari-wfa", "truvari-mafft"]
colors = {"vcfdist": "#73c375", "vcfeval": "#fa6949", "truvari": "#6aadd5", "truvari-wfa": "#9e9ac8", "truvari-mafft": "#959595"}

sizes = ["snp", "indel", "sv"]
SZ_SNP   = 0
SZ_INDEL = 1
SZ_SV    = 2
SZ_DIMS  = 3

categories = ["tp-baseline", "tp", "fn", "fp"]
T_TP = 0
Q_TP = 1
T_FN = 2
Q_FP = 3
CATS = 4 # categories

counts = np.zeros((len(sizes), len(tools), len(datasets), CATS), dtype=int)

def get_size_idx(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return SZ_SNP
    else:
        size_diff = abs(len(ref) - len(alt))
        if size_diff == 0:
            return -1
        elif size_diff < 50:
            return SZ_INDEL
        else:
            return SZ_SV

def count(gt):
    # count non-zero alleles, ignoring phasing
    return sum([int(x) != 0 for x in gt.replace('/', '|').split("|")])

# count variants
for ds_idx, ds in enumerate(datasets):
    print(f"    {ds}")
    for tool_idx, tool in enumerate(tools):
        print(f"        {tool}")
        tool_type = tool[:7]

        if tool_type == "vcfdist":
            results_fn = f"{tool}/{ds}.precision-recall-summary.tsv"
            with open(results_fn, 'r') as results_file:
                for line in results_file:
                    try:
                        line = next(results_file)
                        var_size, minq, t_tp, q_tp, t_fn, q_fp, prec, recall, f1, f1q = \
                                line.strip().split()
                        var_size_idx = sizes.index(var_size.lower())
                        counts[var_size_idx][tool_idx][ds_idx][T_TP] = int(t_tp)
                        counts[var_size_idx][tool_idx][ds_idx][Q_TP] = int(q_tp)
                        counts[var_size_idx][tool_idx][ds_idx][T_FN] = int(t_fn)
                        counts[var_size_idx][tool_idx][ds_idx][Q_FP] = int(q_fp)
                    except StopIteration:
                        pass

        elif tool_type == "vcfeval":
            for cat_idx, cat in enumerate(categories):
                vcf = open(f"{tool}/{ds}/{cat}.vcf", 'r')
                for record in vcf:
                    if record[0] == "#": continue
                    fields = record.strip().split('\t')
                    ref = fields[3]
                    alt = fields[4]
                    if alt == "*": continue
                    size_idx = get_size_idx(ref, alt)
                    if size_idx < 0: continue # INV

                    fmt = fields[8].split(':')
                    gt_ct = 1
                    for idx, fmt_field in enumerate(fmt):
                        if fmt_field == "GT":
                            gt = fields[9].split(':')[idx]
                            gt_ct = count(gt)
                    counts[size_idx][tool_idx][ds_idx][cat_idx] += gt_ct

        elif tool == "truvari":
            for callset in ["query", "truth"]:
                vcf = open(f"{tool}/{ds}/result_{callset}.vcf", 'r')
                for record in vcf:
                    if record[0] == "#": continue
                    fields = record.strip().split('\t')
                    ref = fields[3]
                    alt = fields[4]
                    if alt == "*": continue
                    size_idx = get_size_idx(ref, alt)
                    if size_idx < 0: continue # INV

                    fmt = fields[8].split(':')
                    gt_ct = 1
                    for idx, fmt_field in enumerate(fmt):
                        if fmt_field == "GT":
                            gt = fields[9].split(':')[idx]
                            gt_ct = count(gt)

                    for idx, fmt_field in enumerate(fmt):
                        if fmt_field == "BD":
                            bd = fields[9].split(':')[idx]
                            if bd == "TP" and callset == "query":
                                counts[size_idx][tool_idx][ds_idx][Q_TP] += gt_ct
                            elif bd == "TP" and callset == "truth":
                                counts[size_idx][tool_idx][ds_idx][T_TP] += gt_ct
                            elif bd == "FN":
                                counts[size_idx][tool_idx][ds_idx][T_FN] += gt_ct
                            elif bd == "FP":
                                counts[size_idx][tool_idx][ds_idx][Q_FP] += gt_ct
                            else:
                                print(f"ERROR: bd = {bd}, fields = {fields}")

        else: # truvari refine wfa/mafft
            for output in ["phab_bench", "no_phab_bench"]:
                for callset in ["query", "truth"]:
                    vcf = open(f"{tool}/{ds}/{output}/result_{callset}.vcf", 'r')
                    for record in vcf:
                        if record[0] == "#": continue
                        fields = record.strip().split('\t')
                        ref = fields[3]
                        alt = fields[4]
                        if alt == "*": continue
                        size_idx = get_size_idx(ref, alt)
                        if size_idx < 0: continue # INV

                        fmt = fields[8].split(':')
                        gt_ct = 1
                        for idx, fmt_field in enumerate(fmt):
                            if fmt_field == "GT":
                                gt = fields[9].split(':')[idx]
                                gt_ct = count(gt)

                        for idx, fmt_field in enumerate(fmt):
                            if fmt_field == "BD":
                                bd = fields[9].split(':')[idx]
                                if bd == "TP" and callset == "query":
                                    counts[size_idx][tool_idx][ds_idx][Q_TP] += gt_ct
                                elif bd == "TP" and callset == "truth":
                                    counts[size_idx][tool_idx][ds_idx][T_TP] += gt_ct
                                elif bd == "FN":
                                    counts[size_idx][tool_idx][ds_idx][T_FN] += gt_ct
                                elif bd == "FP":
                                    counts[size_idx][tool_idx][ds_idx][Q_FP] += gt_ct
                                else:
                                    print(f"ERROR: bd = {bd}, fields = {fields}")

for size_idx, size in enumerate(sizes):
    print(f"{size.upper()}:")
    for tool_idx, tool in enumerate(tools):
        print(f"  {tool}")
        for ds_idx, ds in enumerate(datasets):
            print(f"    {ds}: \t{counts[size_idx][tool_idx][ds_idx][Q_TP]} Q_TP\t{counts[size_idx][tool_idx][ds_idx][Q_FP]} Q_FP\t{counts[size_idx][tool_idx][ds_idx][T_TP]} T_TP\t{counts[size_idx][tool_idx][ds_idx][T_FN]} T_FN")

fig, ax = plt.subplots(1, 3, figsize=(7,2.5))
indices = np.arange(len(datasets))
width = 0.15
yquals = [0, 3.01, 6.99, 10, 13.01, 16.99, 20, 23.01, 26.99, 30, # 33.01, 36.99, 40, 43.01, 46.99, 50
        ]
ylabels = ["0%", "50%", "80%", "90%", "95%", "98%", "99%", "99.5%", "99.8%", "99.9%", # "99.95%", "99.98%", "99.99%", "99.995%", "99.998%", "99.999%"
        ]
yquals2 = [0, 3.01, 6.99, 10, 13.01, 16.99, 20, 23.01, 26.99, 30] 
ylabels2 = ["0.1%", "0.2%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%", "100%"]
for var_size_idx, var_size in enumerate(sizes):
    for tool_idx, t in enumerate(tools):

        fnr_fracs = [1 - counts[var_size_idx][tool_idx][ds_idx][T_TP] / 
                max(1, counts[var_size_idx][tool_idx][ds_idx][T_TP] + \
                        counts[var_size_idx][tool_idx][ds_idx][T_FN]) for ds_idx in range(len(datasets))]
        fnr_qscores = [30 if frac == 0 else -10*np.log10(frac) for frac in fnr_fracs]

        ax[var_size_idx].bar(indices - 2*width + tool_idx*width, 
                [30-x for x in fnr_qscores], width, color=colors[t])
    for yqual in yquals:
        ax[var_size_idx].axhline(y=yqual, color='k', alpha=0.5, linestyle=':', ms=0.5, zorder=-1)
    ax[var_size_idx].set_title(f"{var_size.upper()} evaluation", fontsize=7)
    ax[var_size_idx].set_xticks(range(len(datasets)))
    ax[var_size_idx].tick_params(axis='x')
    ax[var_size_idx].set_xticklabels([names[ds] for ds in datasets], fontsize=5)
    ax[var_size_idx].set_yticks(yquals2)
    ax[var_size_idx].set_yticklabels(ylabels2, fontsize=5)
    ax[var_size_idx].set_ylim(0,30)
patches = [mpatches.Patch(color=colors[t], label=t.upper()) for t in tools]
ax[0].legend(handles=patches, loc='upper left', fontsize=5)
ax[0].set_ylabel("False Negative Rate", fontsize=7)
plt.tight_layout()
plt.savefig(f"./img/fnr.pdf", format='pdf')

fig, ax = plt.subplots(1, 3, figsize=(7,2.5))
for var_size_idx, var_size in enumerate(sizes):
    for tool_idx, t in enumerate(tools):

        fpr_fracs = [counts[var_size_idx][tool_idx][ds_idx][Q_FP] / 
                max(1, counts[var_size_idx][tool_idx][ds_idx][Q_FP] + \
                        counts[var_size_idx][tool_idx][ds_idx][Q_TP]) for ds_idx in range(len(datasets))]
        fpr_qscores = [50 if not frac else -10*np.log10(frac) for frac in fpr_fracs]

        ax[var_size_idx].bar(indices - 2*width + tool_idx*width, 
                [30-x for x in fpr_qscores], width, color=colors[t])
    for yqual in yquals2:
        ax[var_size_idx].axhline(y=yqual, color='k', alpha=0.5, linestyle=':', ms=0.5, zorder=-1)
    ax[var_size_idx].set_title(f"{var_size.upper()} evaluation", fontsize=7)
    ax[var_size_idx].set_xticks(range(len(datasets)))
    ax[var_size_idx].tick_params(axis='x')
    ax[var_size_idx].set_xticklabels([names[ds] for ds in datasets], fontsize=5)
    ax[var_size_idx].set_yticks(yquals2)
    ax[var_size_idx].set_yticklabels(ylabels2, fontsize=5)
    ax[var_size_idx].set_ylim(0,30)
patches = [mpatches.Patch(color=colors[t], label=t.upper()) for t in tools]
ax[0].legend(handles=patches, loc='upper left', fontsize=5)
ax[0].set_ylabel("False Positive Rate", fontsize=7)
plt.tight_layout()
plt.savefig(f"./img/fpr.pdf", format='pdf')
