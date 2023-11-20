import numpy as np
import matplotlib.pyplot as plt

VD_NONE = 0
VD_PROB_NONE = 1
VD_PROB_FLIP = 2
VD_FLIP = 3
VD_CATS = 4

WH_NONE = 0
WH_FLIP = 1

PHASE_ORIG = 0
PHASE_SWAP = 1

results = np.zeros((VD_CATS,2))
wh_cat_strs = ["NONE", "FLIP"]
vd_cat_strs = ["NONE", "LIKELY NONE", "LIKELY FLIP", "FLIP"]
colors = {"hprc": "Greens", "pav": "Reds", "giab-tr": "Blues"}
ds = "hprc"

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


class WHFlip():
    def __init__(self, ctg, start):
        self.ctg = ctg
        self.start = start

    def __str__(self):
        return f"{self.ctg}:{self.start}"


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

# # get vcfdist flip superclusters
# vcfdist_tsv = open(f"vcfdist/{ds}.switchflips.tsv", "r")
# next(vcfdist_tsv) # skip header
# vd_flip_scs = []
# for line in vcfdist_tsv:
#     fields = line.strip().split()
#     fliptype = fields[3]
#     sc = int(fields[4])
#     if fliptype == "FLIP_BEG":
#         vd_flip_scs.append(sc)
# print("vcfdist flips:", len(vd_flip_scs))

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

            category = -1
            if vd_sc.phase == PHASE_ORIG:
                if vd_sc.orig_ed == 0:
                    category = VD_NONE
                elif vd_sc.swap_ed == 0:
                    category = VD_FLIP
                elif vd_sc.swap_ed >= vd_sc.orig_ed:
                    category = VD_PROB_NONE
                elif 1 - float(vd_sc.swap_ed) / vd_sc.orig_ed >= 0.7:
                    category = VD_PROB_FLIP
                else:
                    category = VD_PROB_NONE

            else: # PHASE_SWAP
                if vd_sc.swap_ed == 0:
                    category = VD_NONE
                elif vd_sc.orig_ed == 0:
                    category = VD_FLIP
                elif vd_sc.orig_ed >= vd_sc.swap_ed:
                    category = VD_PROB_NONE
                elif 1 - float(vd_sc.orig_ed) / vd_sc.swap_ed >= 0.7:
                    category = VD_PROB_FLIP
                else:
                    category = VD_PROB_NONE

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
        category = -1
        if vd_sc.phase == PHASE_ORIG:
            if vd_sc.orig_ed == 0:
                category = VD_NONE
            elif vd_sc.swap_ed == 0:
                category = VD_FLIP
            elif vd_sc.swap_ed >= vd_sc.orig_ed:
                category = VD_PROB_NONE
            elif 1 - float(vd_sc.swap_ed) / vd_sc.orig_ed >= 0.7:
                category = VD_PROB_FLIP
            else:
                category = VD_PROB_NONE

        else: # PHASE_SWAP
            if vd_sc.swap_ed == 0:
                category = VD_NONE
            elif vd_sc.orig_ed == 0:
                category = VD_FLIP
            elif vd_sc.orig_ed >= vd_sc.swap_ed:
                category = VD_PROB_NONE
            elif 1 - float(vd_sc.orig_ed) / vd_sc.swap_ed >= 0.7:
                category = VD_PROB_FLIP
            else:
                category = VD_PROB_NONE

        results[category][WH_NONE] += 1
        # print(vd_sc)

fig, ax = plt.subplots(figsize=(4,2))
ax.matshow(np.log(results + 0.1), cmap=colors[ds])
plt.title(f"{ds.upper()} Phase Flip CM", fontsize=7)
plt.ylabel("vcfdist", fontsize=7)
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
