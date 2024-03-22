import vcf
import numpy as np
import matplotlib.pyplot as plt

# chr20, with phab
phab = True
ds = "pav"
colors = {"hprc": "Greens", "pav": "Reds", "giab-tr": "Blues"}
names = {"t2t-q100": "Q100-dipcall", "hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR"}
if phab: phab_prefix = "phab-"
else: phab_prefix = ""
truvari_prefix = f"./{phab_prefix}truvari/{ds}/result_"
vcfdist_prefix = f"./{phab_prefix}vcfdist/{ds}."
vcfeval_prefix = f"./{phab_prefix}vcfeval/{ds}/"

SIZE     = 0
SZ_SNP   = 0
SZ_INDEL = 1
SZ_SV    = 2
SZ_DIMS  = 3
sizes = ["SNP", "INDEL", "SV"]

VE_NONE = 0
VE_FP   = 1
VE_TP   = 2
VE_DIMS = 3
ve_cats = ["None", "FP", "TP"]

VD_NONE = 0
VD_FP   = 1
VD_TP   = 2
VD_DIMS = 3
vd_cats = ["None", "FP", "TP"]

TV_NONE = 0
TV_FP   = 1
TV_TP   = 2
TV_DIMS = 3
tv_cats = ["None", "FP", "TP"]

tv_min_seq_pct = 0.7
tv_min_size_pct = 0.7
tv_min_ovlp_pct = 0.0

counts = np.zeros((SZ_DIMS, VD_DIMS-1, 2*(TV_DIMS-1)))

truvari_records = {}
vcfdist_records = {}
vcfeval_records = {}

def record_string(record):
    return f"{record.CHROM}:{record.POS} {record.REF} {record.ALT[0]}"

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
    credit = float(vd_rec.genotype(name)['BC'])
    if credit >= 0 and credit < 0.7:
        return VD_FP
    elif credit >= 0.7 and credit <= 1:
        return VD_TP
    else:
        print("ERROR: credit out of range")
        exit(1)


def get_tv_type(tv_rec):
    decision = tv_rec.genotype(name)['BD']
    seqsim = tv_rec.INFO['PctSeqSimilarity']
    sizesim = tv_rec.INFO['PctSizeSimilarity']
    ovlp = tv_rec.INFO['PctRecOverlap']
    if seqsim == None or sizesim == None or ovlp == None:
        return TV_FP
    if decision == "TP":
        return TV_TP
    if decision == "FP":
        return TV_FP


def get_ve_type(ve_rec):
    decision = ve_rec.genotype(name)['BD']
    if decision == "TP":
        return VE_TP
    if decision == "FP":
        return VE_FP
    if decision == "N":
        return VE_NONE


def count_gts(record):
    if record.genotype(name)['GT'] == "1|1" or \
            record.genotype(name)['GT'] == "1/1":
        return 2
    elif record.genotype(name)['GT'] == '0|1' or \
           record.genotype(name)['GT'] == '0/1' or \
           record.genotype(name)['GT'] == '1|0' or \
           record.genotype(name)['GT'] == '1/0':
        return 1
    return 0


def match(record1, record2):
    return  record1.CHROM == record2.CHROM and \
            record1.POS == record2.POS and \
            record1.REF == record2.REF and \
            record1.ALT[0] == record2.ALT[0] # always one ALT (split multiallelics)


# parse vcfdist and truvari summary VCFs
vcfdist_vcf = vcf.Reader(open(f"{vcfdist_prefix}summary.vcf", "r"))
vcfeval_vcf = vcf.Reader(open(f"{vcfeval_prefix}output.vcf", "r"))
truvari_vcf = vcf.Reader(open(f"{truvari_prefix}query.vcf", "r"))
records = set()
name = "QUERY"

for record in vcfdist_vcf:
    var = record_string(record)
    if record.genotype(name)['BD'] == None: continue
    if var not in vcfdist_records.keys():
        vcfdist_records[var] = (get_size(record), get_vd_type(record), 1)
    else:
        vcfdist_records[var] = (get_size(record), get_vd_type(record), 2)
    records.add(var)

for record in vcfeval_vcf:
    var = record_string(record)
    if record.ALT[0] == "*": continue
    if record.ALT[0] == None: continue
    if record.genotype(name)['BD'] == None: continue
    vcfeval_records[var] = (get_size(record), get_ve_type(record), count_gts(record))
    records.add(var)

for record in truvari_vcf:
    var = record_string(record)
    if record.ALT[0] == "*": continue
    truvari_records[var] = (get_size(record), get_tv_type(record), count_gts(record))
    records.add(var)

for var in records:
    size = None
    vd_type = VD_NONE
    ve_type = VE_NONE
    tv_type = TV_NONE
    gts = 0
    if var in vcfdist_records: (size, vd_type, gts) = vcfdist_records[var]
    if var in vcfeval_records: (size, ve_type, gts) = vcfeval_records[var]
    if var in truvari_records: (size, tv_type, gts) = truvari_records[var]
    if vd_type and ve_type and tv_type:
        counts[size][vd_type-1][(ve_type-1)*2 + tv_type-1] += gts
    elif size == SZ_SNP and vd_type and ve_type and not tv_type:
        counts[size][vd_type-1][(ve_type-1)*2] += gts

for size_idx in range(SZ_DIMS):
    fig, ax = plt.subplots(figsize=(2,4))
    ax.matshow(np.log(counts[size_idx] + 0.1), cmap=colors[ds])
    plt.title(f"{names[ds]}\n{sizes[size_idx]} Confusion Matrix", fontsize=7)
    plt.ylabel("vcfdist", fontsize=7)
    ax.set_yticks(list(range(VD_DIMS-1)))
    ax.set_yticklabels(vd_cats[1:], fontsize=5)
    plt.xlabel("vcfeval / Truvari", fontsize=7)
    ax.set_xticks(list(range((TV_DIMS-1)*2)))
    ax.set_xticklabels([f"FP / {c}" for c in tv_cats[1:]] + [f"TP / {c}" for c in tv_cats[1:]], fontsize=5)
    for (i,j), z in np.ndenumerate(counts[size_idx]):
        ax.text(j, i, f"{int(z)}", ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'), fontsize=5)
    plt.tight_layout()
    plt.savefig(f"img/{ds}_{sizes[size_idx].lower()}_cm.pdf", format='pdf')
