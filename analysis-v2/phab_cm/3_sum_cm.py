import vcf
import numpy as np
import matplotlib.pyplot as plt

datasets = ["hprc", "pav", "giab-tr"]

SIZE     = 0
SZ_SNP   = 0
SZ_INDEL = 1
SZ_SV    = 2
SZ_DIMS  = 3
sizes = ["SNP", "INDEL", "SV"]

NONE = 0
FP   = 1
TP   = 2
DIMS = 3
cats = ["None", "FP", "TP"]

tv_min_seq_pct = 0.7
tv_min_size_pct = 0.7
tv_min_ovlp_pct = 0.0

counts = np.zeros((SZ_DIMS, DIMS-1, 2*(DIMS-1)))
skips = np.zeros((2,2,2), dtype=int)

def record_string(record):
    return f"{record.CHROM}:{record.POS} {record.REF} {record.ALT[0]} {record.genotype('QUERY')['GT']}"

def get_size(record):
    ref = record.REF
    alt = record.ALT[0]
    if len(ref) == 1 and len(alt) == 1:
        return SZ_SNP
    else:
        size_diff = abs(len(ref) - len(alt))
        if size_diff == 0:
            print("ERROR: size 0 INDEL")
            exit(1)
        elif size_diff < 50:
            return SZ_INDEL
        else:
            return SZ_SV


def get_vd_type(vd_rec):
    credit = float(vd_rec.genotype('QUERY')['BC'])
    if credit >= 0 and credit < 0.7:
        return FP
    elif credit >= 0.7 and credit <= 1:
        return TP
    else:
        print("ERROR: credit out of range")
        exit(1)


def get_tv_type(tv_rec):
    decision = tv_rec.genotype('QUERY')['BD']
    seqsim = tv_rec.INFO['PctSeqSimilarity']
    sizesim = tv_rec.INFO['PctSizeSimilarity']
    ovlp = tv_rec.INFO['PctRecOverlap']
    if seqsim == None or sizesim == None or ovlp == None:
        return FP
    if decision == "TP":
        return TP
    if decision == "FP":
        return FP


def get_ve_type(ve_rec):
    decision = ve_rec.genotype('QUERY')['BD']
    if decision == "TP":
        return TP
    if decision == "FP":
        return FP
    if decision == "N":
        return NONE


def count_gts(record):
    if record.genotype('QUERY')['GT'] == "1|1" or \
            record.genotype('QUERY')['GT'] == "1/1":
        return 2
    elif record.genotype('QUERY')['GT'] == '0|1' or \
           record.genotype('QUERY')['GT'] == '0/1' or \
           record.genotype('QUERY')['GT'] == '1|0' or \
           record.genotype('QUERY')['GT'] == '1/0':
        return 1
    return 0

def count(gt):
    # count non-zero alleles, ignoring phasing
    return sum([int(x) != 0 for x in gt.replace('/', '|').split("|")])


def match(record1, record2):
    return  record1.CHROM == record2.CHROM and \
            record1.POS == record2.POS and \
            record1.REF == record2.REF and \
            record1.ALT[0] == record2.ALT[0] # always one ALT (split multiallelics)

def hap1(var):
    return var.replace('/', '|').replace('1|1', '1|0')
def hap2(var):
    return var.replace('/', '|').replace('1|1', '0|1')



am = pp = am_diff = snp = 0
for ds in datasets:
    truvari_records = {}
    vcfdist_records = {}
    vcfeval_records = {}

    # parse summary VCFs
    vcfdist_vcf = vcf.Reader(open(f"./phab-vcfdist/{ds}.summary.vcf", "r"))
    vcfeval_vcf = vcf.Reader(open(f"./phab-vcfeval/{ds}/output.vcf", "r"))
    truvari_vcf = vcf.Reader(open(f"./phab-truvari/{ds}/result_query.vcf", "r"))
    records = set()

    for record in vcfdist_vcf:
        var = record_string(record)
        if record.genotype('QUERY')['BD'] == None: continue
        vcfdist_records[var] = (get_size(record), get_vd_type(record), 
                record.genotype('QUERY')['BK'])
        records.add(var)

    for record in vcfeval_vcf:
        var = record_string(record)
        if record.ALT[0] == "*": continue
        if record.ALT[0] == None: continue
        if record.genotype('QUERY')['BD'] == None: continue
        if record.genotype('QUERY')['BD'] == 'N': continue
        if count_gts(record) == 2:
            vcfeval_records[hap1(var)] = (get_size(record), 
                    get_ve_type(record), record.genotype('QUERY')['BK'])
            records.add(hap1(var))
            vcfeval_records[hap2(var)] = (get_size(record), 
                    get_ve_type(record), record.genotype('QUERY')['BK'])
            records.add(hap2(var))
        else:
            vcfeval_records[var.replace('/', '|')] = (get_size(record), get_ve_type(record), 
                    record.genotype('QUERY')['BK'])
            records.add(var.replace('/', '|'))

    for record in truvari_vcf:
        var = record_string(record)
        if record.ALT[0] == "*": continue
        if count_gts(record) == 2:
            truvari_records[hap1(var)] = (get_size(record), 
                    get_tv_type(record), record.genotype('QUERY')['BK'])
            truvari_records[hap2(var)] = (get_size(record), 
                    get_tv_type(record), record.genotype('QUERY')['BK'])
            records.add(hap1(var))
            records.add(hap2(var))
        else:
            truvari_records[var.replace('/', '|')] = (get_size(record), 
                    get_tv_type(record), record.genotype('QUERY')['BK'])
            records.add(var.replace('/', '|'))

    for var in records:
        size = None
        vd_type, ve_type, tv_type = NONE, NONE, NONE
        gts = 0
        if var in vcfdist_records: (size, vd_type, vd_bk) = vcfdist_records[var]
        if var in vcfeval_records: (size, ve_type, ve_bk) = vcfeval_records[var]
        if var in truvari_records: (size, tv_type, tv_bk) = truvari_records[var]
        skips[0 if vd_type else 1][0 if ve_type else 1][0 if tv_type else 1] += 1
        if vd_type and ve_type and tv_type:
            i = vd_type-1
            j = (ve_type-1)*2 + tv_type-1
            counts[size][i][j] += 1
            if i == 0 and j == 1:
                if tv_bk == "am": # truvari allele match
                    am += 1
                    print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")
                else:
                    # print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")
                    pass
            if i == 0 and j == 2:
                # vcfeval considers two INSs at same location
                # print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")
                pass
            if i == 0 and j == 3:
                pass
                # print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")
            if i == 1 and j == 0:
                pass
                # print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")
            if i == 1 and j == 1:
                if vd_bk == "lm":
                    pp += 1
                elif vd_bk == "gm" and ve_bk == "am" and tv_bk == "am":
                    am_diff += 1
                else:
                    # print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")
                    pass
            if i == 1 and j == 2:
                pass
                # print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")

        elif size == SZ_SNP and vd_type and ve_type and not tv_type: # Truvari skips SNPs
            counts[size][vd_type-1][(ve_type-1)*2] += 1
            snp += 1
        elif ve_type and tv_type and not vd_type:
            # print(f"{ds} {var} VD={cats[vd_type]}:{vd_bk} VE={cats[ve_type]}:{ve_bk} TV={cats[tv_type]}:{tv_bk}")
            pass


sum_counts = counts[SZ_INDEL][:][:] + counts[SZ_SV][:][:]
print("[truvari-only TP] more lenient allele match:", am)
print("[vcfeval-only FP] inexact match:", pp)
print("[vcfeval-only FP] allele match stringency diff:", am_diff)

# print("         all keep:", skips[0][0][0])
# print("         all skip:", skips[1][1][1])
# print("vcfdist-only skip:", skips[1][0][0])
# print("vcfeval-only skip:", skips[0][1][0])
# print("truvari-only skip:", skips[0][0][1], f" of which {snp} are SNPs")
# print("vcfdist-only keep:", skips[0][1][1])
# print("vcfeval-only keep:", skips[1][0][1])
# print("truvari-only keep:", skips[1][1][0])

fig, ax = plt.subplots(figsize=(3,6))
ax.matshow(np.log(sum_counts + 0.1), cmap='Greys')
plt.title(f"INDEL + SV Confusion Matrix", fontsize=7)
plt.ylabel("vcfdist", fontsize=7)
ax.set_yticks(list(range(DIMS-1)))
ax.set_yticklabels(cats[1:], fontsize=5)
plt.xlabel("vcfeval / Truvari", fontsize=7)
ax.set_xticks(list(range((DIMS-1)*2)))
ax.set_xticklabels([f"FP / {c}" for c in cats[1:]] + [f"TP / {c}" for c in cats[1:]], fontsize=5)
for (i,j), z in np.ndenumerate(sum_counts):
    ax.text(j, i, f"{int(z)}", ha='center', va='center',
        bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'), fontsize=5)
plt.tight_layout()
plt.savefig(f"img/sum_cm.pdf", format='pdf')

fig, ax = plt.subplots(figsize=(3,6))
ax.matshow(np.log(sum_counts + 0.1), cmap='Greys')
plt.title(f"INDEL + SV Confusion Matrix", fontsize=7)
plt.ylabel("vcfdist", fontsize=7)
ax.set_yticks(list(range(DIMS-1)))
ax.set_yticklabels(cats[1:], fontsize=5)
plt.xlabel("vcfeval / Truvari", fontsize=7)
ax.set_xticks(list(range((DIMS-1)*2)))
ax.set_xticklabels([f"FP / {c}" for c in cats[1:]] + [f"TP / {c}" for c in cats[1:]], fontsize=5)
for (i,j), z in np.ndenumerate(sum_counts):
    text = ""
    if i == 0 and j == 0:
        text = f"all\ncall FP"
    elif i == 0 and j == 1:
        text = f"a) only\ntruvari\ncalls TP"
    elif i == 0 and j == 2:
        text = f"b) only\nvcfeval\ncalls TP"
    elif i == 0 and j == 3:
        text = f"c) only\nvcfdist\ncalls FP"
    elif i == 1 and j == 0:
        text = f"d) only\nvcfdist\ncalls TP"
    elif i == 1 and j == 1:
        text = f"e) only\nvcfeval\ncalls FP"
    elif i == 1 and j == 2:
        text = f"f) only\ntruvari\ncalls FP"
    else:
        text = f"all\ncall TP"
    ax.text(j, i, f"{text}", ha='center', va='center',
        bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'), fontsize=5)
plt.tight_layout()
plt.savefig(f"img/label_cm.pdf", format='pdf')
