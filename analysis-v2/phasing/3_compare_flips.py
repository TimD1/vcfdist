import numpy as np
import matplotlib.pyplot as plt
from random import sample

VD_NONE = 0
VD_UNKNOWN_NONE = 1
VD_UNKNOWN_FLIP = 2
VD_FLIP = 3
VD_CATS = 4

WH_NONE = 0
WH_FLIP = 1

PHASE_ORIG = 0
PHASE_SWAP = 1

wh_cat_strs = ["NONE", "FLIP"]
vd_cat_strs = ["NONE", "UNKNOWN / NONE", "UNKNOWN / FLIP", "FLIP"]
colors = {"hprc": "Greens", "pav": "Reds", "giab-tr": "Blues"}
names = {"t2t-q100": "Q100-dipcall", "hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR"}
datasets = ["hprc", "pav", "giab-tr"]
tests = { "hprc": ["chr9:81725221-81725225", "chr7:127706261-127706267", "chr14:39094226-39094447", "chr20:21598710-21598773", "chr8:95596931-95597391", "chr1:11080562-11080570", "chr7:73123787-73123789", "chr1:93507055-93507114", "chr11:124299970-124299974", "chr17:1060281-1067362", "chr3:180295301-180295321", "chr6:103637864-103637922", "chr3:41721332-41721366", "chr4:100361368-100361376", "chr7:15097103-15097107", "chr17:16092645-16092703"],
        "giab-tr": [ "chr3:184752974-184754943", "chr3:195477243-195508035", "chr12:29643690-29644112", "chr2:242003448-242012851", "chr13:114042597-114045525", "chr4:75571670-75571681", "chr10:130788761-130792518", "chr5:83361-85442", "chr3:195775381-195790683", "chr18:79012005-79020176", "chr16:85397299-85419607", "chr7:61145443-61145458", "chrX:268020-273177", "chr7:35044332-35046217", "chr3:76256352-76256358", "chr4:184785599-184785623"],
        "pav": [ "chr4:71156554-71156589", "chr16:1113777-1114085", "chr4:122180024-122180048", "chr19:58427287-58427509", "chr1:975056-979561", "chr4:148492402-148494852", "chr4:167498959-167499009", "chr4:108123431-108123435", "chr4:118299711-118299766", "chr4:150250951-150251882", "chr6:157309761-157315305", "chr2:169465347-169468449", "chr11:17288879-17288964", "chr4:72529977-72530137", "chr18:79012005-79020176", "chr5:648572-651181"]
        }

class SC():

    def __init__(self, line):
        fields = line.strip().split()
        self.ctg = fields[0]
        self.sc = int(fields[1])
        self.start = int(fields[2])
        self.stop = int(fields[3])
        self.orig_ed = int(fields[9])
        self.swap_ed = int(fields[10])
        self.phase = int(fields[11])
        self.flip = int(fields[15])

    def __str__(self):
        return f"{self.ctg}:{self.start}-{self.stop} sc={self.sc} ed={self.orig_ed},{self.swap_ed} phase={self.phase} {'flip' if self.flip else ''}"

    def __repr__(self):
        return self.__str__()


class WHFlip():
    def __init__(self, ctg, start):
        self.ctg = ctg
        self.start = start

    def __str__(self):
        return f"{self.ctg}:{self.start}"


for ds in datasets:
    results = np.zeros((VD_CATS,2))
    unknowns = []

    # get whatshap flips
    whatshap_bed = open(f"whatshap/{ds}.switches.bed", "r")
    wh_flips = []
    prev_start, prev_stop = 0, 0
    prev_ctg = ""
    for line in whatshap_bed:
        fields = line.strip().split()
        ctg = fields[0]
        start = int(fields[1])
        stop = int(fields[2])
        if start == prev_stop and ctg == prev_ctg:
            wh_flips.append( WHFlip(ctg, start) )
        prev_start, prev_stop = start, stop
        prev_ctg = ctg

    print("whatshap flips:", len(wh_flips))

    # get all vcfdist superclusters
    vcfdist_scs = open(f"vcfdist/{ds}.superclusters.tsv", "r")
    next(vcfdist_scs) # skip header
    vd_flip_scs = []
    vd_scs = []
    for line in vcfdist_scs:
        vd_scs.append(SC(line))
        if int(line.strip().split()[-1]):
            vd_flip_scs.append(SC(line))
    print("vcfdist flips:", len(vd_flip_scs))

    # map whatshap flips to vcfdist categories
    for wh_flip in wh_flips:
        for vd_sc in vd_scs:

            if wh_flip.ctg == vd_sc.ctg and \
                    wh_flip.start >= vd_sc.start and \
                    wh_flip.start <= vd_sc.stop:

                # pos = f"{vd_sc.ctg}:{vd_sc.start}-{vd_sc.stop}"
                # if pos in tests[ds]:
                #     print(pos, "flip")

                category = -1
                if vd_sc.phase == PHASE_ORIG:
                    if vd_sc.orig_ed == 0:
                        category = VD_NONE
                    elif vd_sc.swap_ed == 0:
                        category = VD_FLIP
                    elif vd_sc.swap_ed >= vd_sc.orig_ed:
                        category = VD_UNKNOWN_NONE
                        unknowns.append(vd_sc)
                    elif 1 - float(vd_sc.swap_ed) / vd_sc.orig_ed >= 0.7:
                        category = VD_UNKNOWN_FLIP
                        unknowns.append(vd_sc)
                    else:
                        category = VD_UNKNOWN_NONE
                        unknowns.append(vd_sc)

                else: # PHASE_SWAP
                    if vd_sc.swap_ed == 0:
                        category = VD_NONE
                    elif vd_sc.orig_ed == 0:
                        category = VD_FLIP
                    elif vd_sc.orig_ed >= vd_sc.swap_ed:
                        category = VD_UNKNOWN_NONE
                        unknowns.append(vd_sc)
                    elif 1 - float(vd_sc.orig_ed) / vd_sc.swap_ed >= 0.7:
                        category = VD_UNKNOWN_FLIP
                        unknowns.append(vd_sc)
                    else:
                        category = VD_UNKNOWN_NONE
                        unknowns.append(vd_sc)

                results[category][WH_FLIP] += 1
                # if category != VD_SC_ED_RATIO_PASS:
                #     print(vd_sc)

                break # found corresponding supercluster

    for vd_sc in vd_flip_scs:
        found_match = False
        for wh_flip in wh_flips:
            if wh_flip.ctg == vd_sc.ctg and \
                    wh_flip.start >= vd_sc.start and \
                    wh_flip.start <= vd_sc.stop:
                found_match = True
        if not found_match:
            pos = f"{vd_sc.ctg}:{vd_sc.start}-{vd_sc.stop}"
            if pos in tests[ds]:
                print(pos, "none")
            category = -1
            if vd_sc.phase == PHASE_ORIG:
                if vd_sc.orig_ed == 0:
                    category = VD_NONE
                elif vd_sc.swap_ed == 0:
                    category = VD_FLIP
                elif vd_sc.swap_ed >= vd_sc.orig_ed:
                    category = VD_UNKNOWN_NONE
                    unknowns.append(vd_sc)
                elif 1 - float(vd_sc.swap_ed) / vd_sc.orig_ed >= 0.7:
                    category = VD_UNKNOWN_FLIP
                    unknowns.append(vd_sc)
                else:
                    category = VD_UNKNOWN_NONE
                    unknowns.append(vd_sc)

            else: # PHASE_SWAP
                if vd_sc.swap_ed == 0:
                    category = VD_NONE
                elif vd_sc.orig_ed == 0:
                    category = VD_FLIP
                elif vd_sc.orig_ed >= vd_sc.swap_ed:
                    category = VD_UNKNOWN_NONE
                    unknowns.append(vd_sc)
                elif 1 - float(vd_sc.orig_ed) / vd_sc.swap_ed >= 0.7:
                    category = VD_UNKNOWN_FLIP
                    unknowns.append(vd_sc)
                else:
                    category = VD_UNKNOWN_NONE
                    unknowns.append(vd_sc)

            results[category][WH_NONE] += 1
            # print(vd_sc)

    fig, ax = plt.subplots(figsize=(4,2))
    ax.matshow(np.log(results + 0.1), cmap=colors[ds])
    plt.title(f"{names[ds]}\nphasing flip errors", fontsize=7)
    plt.ylabel("Truth / vcfdist", fontsize=7)
    ax.set_yticks(list(range(VD_CATS)))
    ax.set_yticklabels(vd_cat_strs, fontsize=5)
    plt.xlabel("WhatsHap", fontsize=7)
    ax.set_xticks(list(range(2)))
    ax.set_xticklabels(wh_cat_strs, fontsize=5)
    for (i,j), z in np.ndenumerate(results):
        ax.text(j, i, f"{int(z)}", ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'), fontsize=5)
    plt.tight_layout()
    plt.savefig(f"img/{ds}_flip_cm.pdf", format='pdf')

    # x = sample(unknowns, 16)
    # print("Random samples:")
    # for y in x:
    #     print(y)
