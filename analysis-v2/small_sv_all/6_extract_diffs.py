import sys

datasets = ["hprc", "pav", "giab-tr"]
vcf_types = ["small", "sv"]

def RED(text):
    return "\x1b[0;31;49m" + text + "\033[0m"
def YLW(text):
    return "\x1b[0;33;49m" + text + "\033[0m"
def GRN(text):
    return "\x1b[0;32;49m" + text + "\033[0m"
# def RED(text):
#     return "\x1b[7;49;31m" + text + "\033[0m"
# def YLW(text):
#     return "\x1b[7;49;33m" + text + "\033[0m"
# def GRN(text):
#     return "\x1b[7;49;32m" + text + "\033[0m"

all_sc_dict = {}
for ds in datasets:
    for vcf in vcf_types:
        print(f"DATASET: {ds}, SUBSET: {vcf}")
        # subset, not substitution mutation
        sub_file = open(f"vcfdist/{ds}.{vcf}.summary.vcf", 'r')
        all_file = open(f"vcfdist/{ds}.all.summary.vcf", 'r')
        sub_vars = [l for l in sub_file.readlines() if l[0] != "#"]
        all_vars = [l for l in all_file.readlines() if l[0] != "#"]
        out_file = open(f"diffs/{ds}.{vcf}_vs_all.txt", 'w')

        # get all differences in subset (per supercluster)
        sub_scs = []
        all_scs = []
        sc_ctgs = []
        svi = 0
        for avi, all_line in enumerate(all_vars):
            # extract fields
            all_ctg, all_pos, all_name, all_ref, all_alt, all_qual, all_filt, \
                    all_info, all_fmt, all_truth, all_query = \
                    all_line.strip().split('\t')
            if svi >= len(sub_vars): break
            sub_line = sub_vars[svi]
            sub_ctg, sub_pos, sub_name, sub_ref, sub_alt, sub_qual, sub_filt, \
                    sub_info, sub_fmt, sub_truth, sub_query = \
                    sub_line.strip().split('\t')

            # extract tags
            all_qGT, all_qBD, all_qBC, all_qBK, all_qQQ, all_qSC, all_qSP = all_query.split(":")
            all_tGT, all_tBD, all_tBC, all_tBK, all_tQQ, all_tSC, all_tSP = all_truth.split(":")
            all_SC = all_tSC if all_qSC == "." else all_qSC

            sub_qGT, sub_qBD, sub_qBC, sub_qBK, sub_qQQ, sub_qSC, sub_qSP = sub_query.split(":")
            sub_tGT, sub_tBD, sub_tBC, sub_tBK, sub_tQQ, sub_tSC, sub_tSP = sub_truth.split(":")
            sub_SC = sub_tSC if sub_qSC == "." else sub_qSC

            # same variant
            if sub_ctg == all_ctg and sub_pos == all_pos and \
                    sub_ref == all_ref and sub_alt == all_alt:

                # create many-to-one mapping from sub_scs to all_scs (larger due to more variants)
                if sub_SC != "." and all_SC != ".":
                    all_sc_dict[sub_ctg+":"+sub_SC] = all_SC

                # skip exact matches
                if sub_qBD == all_qBD and sub_qBC == all_qBC:
                    svi += 1
                    continue

                # different classification
                if not len(all_scs) or int(all_SC) != all_scs[-1]:
                    sub_scs.append(int(sub_SC))
                    all_scs.append(int(all_SC))
                    sc_ctgs.append(sub_ctg)
                    print(sub_ctg, sub_SC, all_SC)
                svi += 1

            else: # different variant
                if sub_ctg == all_ctg and int(sub_pos) < int(all_pos):
                    # this only occurs with variant missing from 'all' file due to overlap
                    svi += 1

        # print superclusters
        svi = 0
        avi = 0
        for sc_ctg, (sub_sc, all_sc) in zip(sc_ctgs, zip(sub_scs, all_scs)):
            print(f"ctg: {sc_ctg}\t{vcf.upper()} supercluster: {sub_sc}\tALL supercluster: {all_sc}\t{RED('FP/FN')} {YLW('PP')} {GRN('TP')}", file=out_file)
            print("\tCONTIG\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY", file=out_file)

            # get info for all vars
            all_line = all_vars[avi]
            all_ctg, all_pos, all_name, all_ref, all_alt, all_qual, all_filt, \
                    all_info, all_fmt, all_truth, all_query = \
                    all_line.strip().split('\t')
            all_qGT, all_qBD, all_qBC, all_qBK, all_qQQ, all_qSC, all_qSP = all_query.split(":")
            all_tGT, all_tBD, all_tBC, all_tBK, all_tQQ, all_tSC, all_tSP = all_truth.split(":")
            all_SC = all_tSC if all_qSC == "." else all_qSC
            all_BC = all_tBC if all_qBC == "." else all_qBC
            all_BC = max(float(0 if all_tBC == "." else all_tBC), \
                         float(0 if all_qBC == "." else all_qBC))

            # find supercluster start
            while int(all_SC) < all_sc or all_ctg != sc_ctg:
                avi += 1
                all_line = all_vars[avi]
                all_ctg, all_pos, all_name, all_ref, all_alt, all_qual, all_filt, \
                        all_info, all_fmt, all_truth, all_query = \
                        all_line.strip().split('\t')
                all_qGT, all_qBD, all_qBC, all_qBK, all_qQQ, all_qSC, all_qSP = all_query.split(":")
                all_tGT, all_tBD, all_tBC, all_tBK, all_tQQ, all_tSC, all_tSP = all_truth.split(":")
                all_SC = all_tSC if all_qSC == "." else all_qSC
                all_BC = max(float(0 if all_tBC == "." else all_tBC), \
                             float(0 if all_qBC == "." else all_qBC))

            # get info for sub vars
            sub_line = sub_vars[svi]
            sub_ctg, sub_pos, sub_name, sub_ref, sub_alt, sub_qual, sub_filt, \
                    sub_info, sub_fmt, sub_truth, sub_query = \
                    sub_line.strip().split('\t')
            sub_qGT, sub_qBD, sub_qBC, sub_qBK, sub_qQQ, sub_qSC, sub_qSP = sub_query.split(":")
            sub_tGT, sub_tBD, sub_tBC, sub_tBK, sub_tQQ, sub_tSC, sub_tSP = sub_truth.split(":")
            sub_SC = sub_tSC if sub_qSC == "." else sub_qSC
            sub_BC = max(float(0 if sub_tBC == "." else sub_tBC), \
                         float(0 if sub_qBC == "." else sub_qBC))

            # find supercluster start
            while int(sub_SC) < sub_sc or sub_ctg != sc_ctg:
                svi += 1
                sub_line = sub_vars[svi]
                sub_ctg, sub_pos, sub_name, sub_ref, sub_alt, sub_qual, sub_filt, \
                        sub_info, sub_fmt, sub_truth, sub_query = \
                        sub_line.strip().split('\t')
                sub_qGT, sub_qBD, sub_qBC, sub_qBK, sub_qQQ, sub_qSC, sub_qSP = sub_query.split(":")
                sub_tGT, sub_tBD, sub_tBC, sub_tBK, sub_tQQ, sub_tSC, sub_tSP = sub_truth.split(":")
                sub_SC = sub_tSC if sub_qSC == "." else sub_qSC
                sub_BC = max(float(0 if sub_tBC == "." else sub_tBC), \
                             float(0 if sub_qBC == "." else sub_qBC))

            # print supercluster
            while sub_ctg+":"+sub_SC in all_sc_dict and all_sc_dict[sub_ctg+":"+sub_SC] == all_SC:
                if float(sub_BC) == 0:
                    print(f"{vcf.upper().rjust(6)}: {RED(sub_line)}", end="", file=out_file)
                elif float(sub_BC) == 1:
                    print(f"{vcf.upper().rjust(6)}: {GRN(sub_line)}", end="", file=out_file)
                else:
                    print(f"{vcf.upper().rjust(6)}: {YLW(sub_line)}", end="", file=out_file)
                svi += 1
                sub_line = sub_vars[svi]
                sub_ctg, sub_pos, sub_name, sub_ref, sub_alt, sub_qual, sub_filt, \
                        sub_info, sub_fmt, sub_truth, sub_query = \
                        sub_line.strip().split('\t')
                sub_qGT, sub_qBD, sub_qBC, sub_qBK, sub_qQQ, sub_qSC, sub_qSP = sub_query.split(":")
                sub_tGT, sub_tBD, sub_tBC, sub_tBK, sub_tQQ, sub_tSC, sub_tSP = sub_truth.split(":")
                sub_SC = sub_tSC if sub_qSC == "." else sub_qSC
                sub_BC = max(float(0 if sub_tBC == "." else sub_tBC), \
                             float(0 if sub_qBC == "." else sub_qBC))

            # print supercluster
            print(" ", file=out_file)
            while int(all_SC) == all_sc and all_ctg == sc_ctg:
                if float(all_BC) == 0:
                    print(f"   ALL: {RED(all_line)}", end="", file=out_file)
                elif float(all_BC) == 1:
                    print(f"   ALL: {GRN(all_line)}", end="", file=out_file)
                else:
                    print(f"   ALL: {YLW(all_line)}", end="", file=out_file)
                avi += 1
                all_line = all_vars[avi]
                all_ctg, all_pos, all_name, all_ref, all_alt, all_qual, all_filt, \
                        all_info, all_fmt, all_truth, all_query = \
                        all_line.strip().split('\t')
                all_qGT, all_qBD, all_qBC, all_qBK, all_qQQ, all_qSC, all_qSP = all_query.split(":")
                all_tGT, all_tBD, all_tBC, all_tBK, all_tQQ, all_tSC, all_tSP = all_truth.split(":")
                all_SC = all_tSC if all_qSC == "." else all_qSC
                all_BC = max(float(0 if all_tBC == "." else all_tBC), \
                             float(0 if all_qBC == "." else all_qBC))
            print("\n", file=out_file)
