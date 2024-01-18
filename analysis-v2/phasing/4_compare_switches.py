import numpy as np
import matplotlib.pyplot as plt

WH_NONE = 0
WH_SWITCH = 1
WH_CATS = 2

VD_NONE = 0
VD_SWITCH = 1
VD_FLIP = 2
VD_CATS = 3

wh_cat_strs = ["NONE", "SWITCH"]
vd_cat_strs = ["NONE", "SWITCH", "FLIP"]
colors = {"hprc": "Greens", "pav": "Reds", "giab-tr": "Blues"}
names = {"t2t-q100": "Q100-dipcall", "hprc": "hifiasm-dipcall", "pav": "Q100-PAV", "giab-tr": "hifiasm-GIAB-TR"}
datasets = ["hprc", "pav", "giab-tr"]

class VDSwitch():

    def __init__(self, line):
        fields = line.strip().split()
        self.ctg = fields[0]
        self.start = int(fields[1])
        self.stop = int(fields[2])
        if fields[3] == "FLIP_BEG":
            self.type = VD_FLIP
        elif fields[3] == "FLIP_END":
            self.type = VD_FLIP
        elif fields[3] == "SWITCH_ERR":
            self.type = VD_SWITCH
        else:
            print(f"ERROR: unexpected vcfdist switch type '{fields[3]}'")

    def __str__(self):
        return f"{self.ctg}:{self.start}-{self.stop} {vd_cat_strs[self.type]}"


class WHSwitch():
    def __init__(self, line):
        fields = line.strip().split()
        self.ctg = fields[0]
        self.start = int(fields[1])
        self.stop = int(fields[2])

    def __str__(self):
        return f"{self.ctg}:{self.start}-{self.stop}"


for ds in datasets:
    results = np.zeros((VD_CATS, WH_CATS))

    # get whatshap flips
    whatshap_bed = open(f"whatshap/{ds}.switches.bed", "r")
    wh_switches = []
    prev_start, prev_stop = 0, 0
    prev_ctg = ""
    for line in whatshap_bed:
        wh_switches.append( WHSwitch(line) )
    wh_matches = [False] * len(wh_switches)
    print("whatshap switches:", len(wh_switches))

    # get vcfdist flip superclusters
    vcfdist_tsv = open(f"vcfdist/{ds}.switchflips.tsv", "r")
    next(vcfdist_tsv) # skip header
    vd_switches = []
    for line in vcfdist_tsv:
        vd_switches.append( VDSwitch(line) )
    print("vcfdist switches:", len(vd_switches))
    vd_matches = [False] * len(vd_switches)

    def overlap(wh_sw, vd_sw):
        if wh_sw.ctg != vd_sw.ctg:
            return False
        if wh_sw.stop < vd_sw.start:
            return False
        if wh_sw.start > vd_sw.stop:
            return False
        return True

    # find all pairs of matches
    for vdi, vd_sw in enumerate(vd_switches):
        for whi, wh_sw in enumerate(wh_switches):
            if overlap(wh_sw, vd_sw):
                results[vd_sw.type][WH_SWITCH] += 1
                vd_matches[vdi] = True
                wh_matches[whi] = True
                break

    for whi, wh_sw in enumerate(wh_switches):
        if not wh_matches[whi]:
            results[VD_NONE][WH_SWITCH] += 1

    for vdi, vd_sw in enumerate(vd_switches):
        if not vd_matches[vdi]:
            results[vd_sw.type][WH_NONE] += 1


    fig, ax = plt.subplots(figsize=(3,1.7))
    ax.matshow(np.log(results + 0.1), cmap=colors[ds])
    plt.title(f"{names[ds]}\nphasing switch errors", fontsize=7)
    plt.ylabel("vcfdist", fontsize=7)
    ax.set_yticks(list(range(VD_CATS)))
    ax.set_yticklabels(vd_cat_strs, fontsize=5)
    plt.xlabel("WhatsHap", fontsize=7)
    ax.set_xticks(list(range(WH_CATS)))
    ax.set_xticklabels(wh_cat_strs, fontsize=5)
    for (i,j), z in np.ndenumerate(results):
        ax.text(j, i, f"{int(z)}", ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.3'), fontsize=5)
    plt.tight_layout()
    plt.savefig(f"img/{ds}_switch_cm.pdf", format='pdf')
