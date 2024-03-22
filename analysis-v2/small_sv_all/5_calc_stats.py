import json
import copy

datasets = ["t2t-q100", "hprc", "pav", "giab-tr"]
names = {"t2t-q100": "Q100-dipcall", "hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR"}
variant_types = ["snp", "indel", "sv"]
categories = ["small", "small", "sv"]
subsets = ["sv", "small", "all"]
regions = ["summary"]

for region in regions:
    print(f"{region} region:")

    # load in FP and TP
    fp_json= None
    with open(f"{region}-query-fp.json", "r") as fp_json_file:
        fp_json = json.load(fp_json_file)
    tp_json= None
    with open(f"{region}-truth-tp.json", "r") as tp_json_file:
        tp_json = json.load(tp_json_file)

    # calculate fn from tp
    fn_json = copy.deepcopy(tp_json)
    for ds in datasets:
        for vt in variant_types:
            if vt in tp_json[ds] and vt in tp_json["t2t-q100"]:
                for ss in subsets:
                    if ss in tp_json[ds][vt] and ss in tp_json["t2t-q100"][vt]:
                        fn_json[ds][vt][ss] = \
                                tp_json["t2t-q100"][vt][ss] - tp_json[ds][vt][ss]

    for ds in datasets:
        print(f"  {ds} dataset:")
        for vt, cat in zip(variant_types, categories):
            if cat in fp_json[ds][vt] and cat in fn_json[ds][vt]:
                cat_errs = fp_json[ds][vt][cat] + fn_json[ds][vt][cat]
                all_errs = fp_json[ds][vt]["all"] + fn_json[ds][vt]["all"]
                red = 1 - (all_errs / cat_errs)
                print(f"    {vt} variants:\t{cat} errors: {cat_errs}\tall errors: {all_errs}\t{red*100:.2f}% reduction")

