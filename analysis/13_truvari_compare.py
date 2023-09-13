import vcf
import numpy as np
import matplotlib.pyplot as plt

# full genome, no phab
# truvari_prefix = "../out/giabtr/truvari-init/result_"
# vcfdist_prefix = "../out/giabtr/vcfdist-keep/"
# chr20, with phab
truvari_prefix = "../out/giabtr/truvari-norm-phab-chr20/result_"
vcfdist_prefix = "../out/giabtr/vcfdist-norm-phab-chr20/"
do_print = True

SIZE             = 0
SZ_SNP           = 0
SZ_INDEL_1_10    = 1
SZ_INDEL_10_50   = 2
SZ_INDEL_50_500  = 3
SZ_INDEL_500PLUS = 4
SZ_DIMS = 5
sizes = ["SNP", "INDEL 1-10", "INDEL 10-50", "INDEL 50-500", "INDEL 500+"]

VCF_DIST     = 1
VD_NONE      = 0
VD_FPN       = 1
VD_PP_0_25   = 2
VD_PP_25_50  = 3
VD_PP_50_75  = 4
VD_PP_75_100 = 5
VD_TP        = 6
VD_DIMS = 7
vd_cats = ["None", "FP/FN", "PP (0, .25)", "PP (.25, .5)", "PP(.5, .75)", "PP(.75, 1)", "TP"]

TRU_VARI    = 2
TV_NONE     = 0
TV_FPN_ANY  = 1
TV_FPN_SEQ  = 2
TV_FPN_SIZE = 3
TV_FPN_OVLP = 4
TV_TP       = 5
TV_DIMS = 6
tv_cats = ["None", "FP/FN any", "FP/FN seq", "FP/FN size", "FP/FN overlap", "TP"]

tv_min_seq_pct = 0.7
tv_min_size_pct = 0.7
tv_min_ovlp_pct = 0.0

counts = np.zeros((SZ_DIMS, VD_DIMS, TV_DIMS))

def get_size(ref : str, alt : str):
    if len(ref) == 1 and len(alt) == 1:
        return SZ_SNP
    else:
        size_diff = abs(len(ref) - len(alt))
        if size_diff == 0:
            print("ERROR: size 0 INDEL")
            exit(1)
        if size_diff <= 10:
            return SZ_INDEL_1_10
        elif size_diff <= 50:
            return SZ_INDEL_10_50
        elif size_diff <= 500:
            return SZ_INDEL_50_500
        else:
            return SZ_INDEL_500PLUS

def get_vd_type(credit : float):
    if credit == 0:
        return VD_FPN
    elif credit > 0 and credit <= 0.25:
        return VD_PP_0_25
    elif credit > 0 and credit <= 0.50:
        return VD_PP_25_50
    elif credit > 0 and credit <= 0.75:
        return VD_PP_50_75
    elif credit > 0 and credit < 1:
        return VD_PP_75_100
    elif credit == 1:
        return VD_TP
    else:
        print("ERROR: credit out of range")
        exit(1)


def get_tv_type(tv_type : str):
    if tv_type == "TP":
        return TV_TP
    elif tv_type == "FP" or tv_type == "FN":
        return TV_FPN_ANY


# for callset in ["query", "truth"]:
for callset in ["query"]:

    # parse vcfdist and truvari summary VCFs
    vcfdist_vcf = vcf.Reader(open(f"{vcfdist_prefix}summary.vcf", "r"))
    truvari_vcf = vcf.Reader(open(f"{truvari_prefix}{callset}.vcf", "r"))
    name = callset.upper()
    print(name)
    tv_used, vd_used = True, True
    this_ctg = "chr1"
    print(this_ctg)

    while True:

        # get next valid records for each
        if tv_used: tv_rec = next(truvari_vcf, None)
        if vd_used: vd_rec = next(vcfdist_vcf, None)
        while vd_rec != None and vd_rec.genotype(name)['BD'] == None: # skip other callset
            vd_rec = next(vcfdist_vcf, None)
        while tv_rec != None and tv_rec.ALT[0] == "*": # skip nested var
            tv_rec = next(truvari_vcf, None)
        tv_used, vd_used = False, False

        if do_print:
            print("============================================================")
            print("Truvari:", tv_rec)
            print("vcfdist:", vd_rec)

        # we've finished iterating through both VCFs
        if tv_rec == None and vd_rec == None: break

        # we've finished this contig for both VCFs
        if (tv_rec == None or tv_rec.CHROM != this_ctg) and \
                (vd_rec == None or vd_rec.CHROM != this_ctg):
            if tv_rec == None:
                this_ctg = vd_rec.CHROM
                print(this_ctg)
            elif vd_rec == None:
                this_ctg = tv_rec.CHROM
                print(this_ctg)
            elif tv_rec.CHROM == vd_rec.CHROM:
                this_ctg = tv_rec.CHROM
                print(this_ctg)
            else:
                print("ERROR: different contigs up next")

        # if we've finished only one VCF, set high position
        tv_pos = 2_000_000_000 if tv_rec == None else (
                1_000_000_000 if tv_rec.CHROM != this_ctg else tv_rec.POS)
        vd_pos = 2_000_000_000 if vd_rec == None else (
                1_000_000_000 if vd_rec.CHROM != this_ctg else vd_rec.POS)

        if tv_pos == vd_pos: # position match
            if do_print: print(f"TV VD {name} {tv_rec.genotype(name)['GT']} {vd_rec.genotype(name)['GT']} {tv_rec.CHROM}:{tv_rec.POS}\t{tv_rec.REF}\t{tv_rec.ALT[0]}")
            if tv_rec.REF == vd_rec.REF and tv_rec.ALT[0] == vd_rec.ALT[0]: # full match
                tv_used, vd_used = True, True
                size = get_size(tv_rec.REF, tv_rec.ALT[0])
                vd_type = get_vd_type(float(vd_rec.genotype(name)['BC']))
                tv_type = get_tv_type(tv_rec.genotype(name)['BD'])
                counts[size][vd_type][tv_type] += 1
                if vd_type == VD_TP and tv_type == TV_FPN_ANY:
                    if do_print: print("vcfdist TP, Truvari FP")
                if vd_type == VD_FPN and tv_type == TV_TP:
                    if do_print: print("vcfdist FP, Truvari TP")

                # skip if info unavailable
                if tv_rec.INFO['PctSeqSimilarity'] == None or \
                        tv_rec.INFO['PctSizeSimilarity'] == None or \
                        tv_rec.INFO['PctRecOverlap'] == None:
                    counts[size][vd_type][TV_FPN_OVLP] += 1
                    continue

                # count truvari filter fails
                if tv_rec.INFO['PctSeqSimilarity'] < tv_min_seq_pct:
                    counts[size][vd_type][TV_FPN_SEQ] += 1
                if tv_rec.INFO['PctSizeSimilarity'] < tv_min_size_pct:
                    counts[size][vd_type][TV_FPN_SIZE] += 1
                if tv_rec.INFO['PctRecOverlap'] < tv_min_ovlp_pct:
                    counts[size][vd_type][TV_FPN_OVLP] += 1

            # vcfdist SNP adjacent to INS
            elif len(vd_rec.REF) == 1 and len(vd_rec.ALT[0]) == 1:
                vd_used = True
                size = get_size(vd_rec.REF, vd_rec.ALT[0])
                vd_type = get_vd_type(float(vd_rec.genotype(name)['BC']))
                tv_type = TV_NONE
                counts[size][vd_type][tv_type] += 1

            else: # position but not allele match, likely different phasing order

                # get next record for each
                tv_rec_next = next(truvari_vcf, None)
                vd_rec_next = next(vcfdist_vcf, None)
                while vd_rec_next != None and vd_rec_next.genotype(name)['BD'] == None: # skip other callset
                    vd_rec_next = next(vcfdist_vcf, None)
                while tv_rec_next != None and tv_rec_next.ALT[0] == "*": # skip nested var
                    tv_rec_next = next(truvari_vcf, None)

                # test if swapping causes matches
                if tv_rec.REF == vd_rec_next.REF and \
                        tv_rec.ALT[0] == vd_rec_next.ALT[0] and \
                        tv_rec_next.REF == vd_rec.REF and \
                        tv_rec_next.ALT[0] == vd_rec.ALT[0]:
                    # 1: tv_rec, vd_rec_next
                    size1 = get_size(tv_rec.REF, tv_rec.ALT[0])
                    vd_type1 = get_vd_type(float(vd_rec_next.genotype(name)['BC']))
                    tv_type1 = get_tv_type(tv_rec.genotype(name)['BD'])
                    counts[size1][vd_type1][tv_type1] += 1
                    # 2: tv_rec_next, vd_rec
                    size2 = get_size(tv_rec_next.REF, tv_rec_next.ALT[0])
                    vd_type2 = get_vd_type(float(vd_rec.genotype(name)['BC']))
                    tv_type2 = get_tv_type(tv_rec_next.genotype(name)['BD'])
                    counts[size2][vd_type2][tv_type2] += 1
                    tv_used, vd_used = True, True
                else:
                    if do_print: print("ERROR: failed to match")

                    # discard current two, pretend we haven't looked at next two
                    counts[size][VD_NONE][TV_NONE] += 2
                    tv_rec = tv_rec_next
                    vd_rec = vd_rec_next


        elif vd_pos < tv_pos: # vcfdist only
            vd_used = True
            size = get_size(vd_rec.REF, vd_rec.ALT[0])
            if do_print: 
                print(f"   VD {name} {vd_rec.genotype(name)['GT']} {vd_rec.CHROM}:{vd_rec.POS}\t{vd_rec.REF}\t{vd_rec.ALT[0]}\t")
                if size != SZ_SNP:
                    print("WARN: vcfdist only, not SNP")
            vd_type = get_vd_type(float(vd_rec.genotype(name)['BC']))
            tv_type = TV_NONE
            counts[size][vd_type][tv_type] += 1

        elif tv_pos < vd_pos: # truvari only
            if do_print:
                print(f"TV    {name} {tv_rec.genotype(name)['GT']} {tv_rec.CHROM}:{tv_rec.POS}\t{tv_rec.REF}\t{tv_rec.ALT[0]}")
                print("WARN: Truvari only")
            tv_used = True
            size = get_size(tv_rec.REF, tv_rec.ALT[0])
            vd_type = VD_NONE
            tv_type = get_tv_type(tv_rec.genotype(name)['BD'])
            counts[size][vd_type][tv_type] += 1

            # skip if info unavailable
            if tv_rec.INFO['PctSeqSimilarity'] == None or \
                    tv_rec.INFO['PctSizeSimilarity'] == None or \
                    tv_rec.INFO['PctRecOverlap'] == None:
                counts[size][vd_type][TV_FPN_OVLP] += 1
                continue

            # count truvari filter fails
            if tv_rec.INFO['PctSeqSimilarity'] < tv_min_seq_pct:
                counts[size][vd_type][TV_FPN_SEQ] += 1
            if tv_rec.INFO['PctSizeSimilarity'] < tv_min_size_pct:
                counts[size][vd_type][TV_FPN_SIZE] += 1
            if tv_rec.INFO['PctRecOverlap'] < tv_min_ovlp_pct:
                counts[size][vd_type][TV_FPN_OVLP] += 1
 
    for size_idx in range(SZ_DIMS):
        fig, ax = plt.subplots(figsize=(12,12))
        ax.matshow(np.log(counts[size_idx] + 0.1), cmap="Blues")
        plt.title(f"{name} {sizes[size_idx]} Confusion Matrix")
        plt.ylabel("vcdist")
        ax.set_yticks(list(range(VD_DIMS)))
        ax.set_yticklabels(vd_cats)
        plt.xlabel("TruVari")
        ax.set_xticks(list(range(TV_DIMS)))
        ax.set_xticklabels(tv_cats)
        for (i,j), z in np.ndenumerate(counts[size_idx]):
            ax.text(j, i, f"{int(z)}", ha='center', va='center',
                bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'))
        plt.savefig(f"img/13_{callset}_{size_idx}_cm.png", dpi=200)
