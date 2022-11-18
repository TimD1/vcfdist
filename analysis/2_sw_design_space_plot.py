import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['text.usetex'] = True

plt.figure(figsize=(12, 5))
xmax = 3
ymax = 1.5
label_pad = 0.03
plt.xlim(0,xmax)
plt.ylim(0,ymax)
x = list(np.arange(0, xmax+1, 0.1))

# G_open > 0
plt.annotate("$G_{\mathrm{open}} < 0$", (1.25,1.3), rotation=40, color="white")
plt.fill_between(x, x, ymax, alpha=0.8, facecolor='k')

# G_start + G_start > P_sub  (otherwise all SUBs->INDELs)
plt.axvspan(0, 0.5, alpha=0.1, color='k', zorder=-100)
plt.annotate("$2G_{\mathrm{start}} < P_{\mathrm{sub}}$", (0.46,1.13), 
        rotation=90, color="white")

# G_start + G_extend > P_sub (same DE, decreases ED)
plt.fill_between(x, [1-i for i in x], 0, facecolor='k', alpha=0.1)
plt.annotate("$G_{\mathrm{start}} + G_{\mathrm{extend}} < P_{\mathrm{sub}}$", 
        (label_pad,0.55), rotation=-40, color="white")

# G_extend + G_extend < P_sub (shouldn't create new SUB to reduce INDEL length)
plt.fill_between(x, 0.5, ymax, facecolor='k', alpha=0.1)
plt.annotate("$2G_{\mathrm{extend}} < P_{\mathrm{sub}}$", 
        (label_pad,0.55), ha="left", va="top", color="white")

# min ED
i = 0
plt.scatter(1,1)
plt.annotate("edit-dist", (1+label_pad,1), ha="left", va="center", color=f"C{i}")
i += 1

# minimap2
mm2_presets = {
    "mm2-short":    { "a": 0, "b": 10, "o1": 13, "e1":3, "o2":25,  "e2": 2 },
    "mm2-ont":  { "a": 0, "b": 6,  "o1": 5,  "e1":3, "o2":25,  "e2": 2 },
    "mm2-hifi": { "a": 0, "b": 10, "o1": 13, "e1":5, "o2":53,  "e2": 3 },
    "mm2-asm5":     { "a": 0, "b": 40, "o1": 79, "e1":7, "o2":163, "e2": 3 },
    "mm2-asm10":    { "a": 0, "b": 20, "o1": 33, "e1":5, "o2":83,  "e2": 3 },
}
for name, x in mm2_presets.items():
    start = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"], color=f"C{i}")
    end = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e2"]/x["b"], color=f"C{i}")
    plt.arrow((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"]-0.02, 0, x["e2"]/x["b"]-x["e1"]/x["b"]+0.04, color='k',
            length_includes_head=True, head_width=0.02, head_length = 0.016)
    plt.annotate(name, ((x["o1"] + x["e1"])/x["b"]+label_pad, x["e1"]/x["b"]), ha="left", va="center", color=f"C{i}")
    i += 1

# ngmlr
ngmlr_presets = {
    "ngmlr-ont":  { "a": 0, "b": 12, "o": 11, "e_max": 11, "e_min": 3, "e_dec": 0.15},
    "ngmlr-hifi": { "a": 0, "b": 12, "o": 5,  "e_max": 5,  "e_min": 4, "e_dec": 0.15}
}
for name, x in ngmlr_presets.items():
    for e in np.arange(x["e_max"], x["e_min"], -x["e_dec"]):
        plt.scatter(x["o"]/x["b"], e/x["b"], color=f"C{i}")
    plt.annotate(name, (x["o"]/x["b"]+label_pad, (x["e_max"]+x["e_min"])/(2*x["b"])), ha="left", va="center", color=f"C{i}")
    i += 1

# npore
points = [[1.7, 0.5], [1.7,0.7], [1.7, 1], [1.7, 0.8], [1.7, 1.1], [2.3, 1.8], [2.8, 1], [3.2, 1.1], [3.2, 2.1], [3.2, 1], [3.7, 1.9], [3.7, 2], [3.7, 1], [3.7, 0.6], [3.8, 2.4], [3.8, 1.7], [3.8, 0.8], [4.6, 2.5], [4.6, 2.3], [4.6, 1.5]]
for p in points:
    plt.scatter(p[0]/5.72, p[1]/5.72, color=f"C{i}")
plt.annotate("npore", (0.5, 0.25), ha="left", va="center", color=f"C{i}")
plt.scatter(6/5.72, 1/5.72, color=f"C{i}")
plt.annotate("npore", (6/5.72+label_pad, 1/5.72), ha="left", va="center", color=f"C{i}")
i += 1

# verkko
plt.scatter(0, 0, color=f"C{i}")
plt.annotate("verkko", (label_pad, 0), ha="left", va="center", color=f"C{i}")

# points
plt.scatter(0.3, 0.2, color=f"k")
plt.annotate("A", (0.3+label_pad, 0.2), ha="left", va="center", color="k")
plt.scatter(0.7, 0.2, color=f"k")
plt.annotate("B", (0.7+label_pad, 0.2), ha="left", va="center", color="k")
plt.scatter(1.7, 0.2, color=f"k")
plt.annotate("C", (1.7+label_pad, 0.2), ha="left", va="center", color="k")
plt.scatter(1, 0.8, color=f"k")
plt.annotate("D", (1+label_pad, 0.8), ha="left", va="center", color="k")


plt.title("Gap Penalty Design Space")
plt.xlabel(r'$\frac{G_{\mathrm{open}}+G_{\mathrm{extend}}}{P_\mathrm{sub}} = \frac{G_{\mathrm{start}}}{P_{\mathrm{sub}}}$')
plt.ylabel(r'$\frac{G_{\mathrm{extend}}}{P_{\mathrm{sub}}}$')
plt.savefig("img/2_sw_design_space.png")
