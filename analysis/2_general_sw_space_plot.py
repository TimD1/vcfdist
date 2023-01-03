import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['text.usetex'] = True

plt.figure(figsize=(12, 5))
xmin, xmax = 0, 3.2
ymin, ymax = -0.1, 1.5
label_pad = 0.03
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
x = list(np.arange(xmin, xmax+1, 0.1))

# G_open > 0
plt.annotate("$o < 0$", (1.32,1.37), rotation=35, color="white")
plt.fill_between(x, x, ymax, alpha=0.8, facecolor='k')

# vertical
plt.annotate("$2(o+e) < x$", (0.46,1.18), rotation=90, color="white")
plt.axvspan(0, 0.5, alpha=0.1, color='k', zorder=-100)
plt.annotate("$3x < o+e$", (3.01,1.22), rotation=90, color="black")
plt.axvspan(3, 3.2, alpha=0.1, color='k', zorder=-100)

# horizontal
plt.annotate("$x < 2e$", (label_pad,0.57), ha="left", va="top", color="white")
plt.fill_between(x, 0.5, ymax, facecolor='k', alpha=0.1)
plt.annotate("$e < 0$", (label_pad,-0.03), ha="left", va="top", color="white")
plt.fill_between(x, -0.2, 0, facecolor='k', alpha=0.8)

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
    plt.arrow(x["o"]/x["b"], x["e_max"]/x["b"], 0, x["e_min"]/x["b"]-x["e_max"]/x["b"], color='k',
            length_includes_head=True, head_width=0.02, head_length = 0.016)

# npore
points = [[1.7, 0.5], [1.7,0.7], [1.7, 1], [1.7, 0.8], [1.7, 1.1], [2.3, 1.8], [2.8, 1], [3.2, 1.1], [3.2, 2.1], [3.2, 1], [3.7, 1.9], [3.7, 2], [3.7, 1], [3.7, 0.6], [3.8, 2.4], [3.8, 1.7], [3.8, 0.8], [4.6, 2.5], [4.6, 2.3], [4.6, 1.5]]
for p in points:
    plt.scatter(p[0]/5.72, p[1]/5.72, color=f"C{i}")
plt.annotate("npore", (0.5, 0.25), ha="left", va="center", color=f"C{i}")
plt.scatter(6/5.72, 1/5.72, color=f"C{i}")
plt.scatter(6/5.72, 1/5.72, color=f"k", marker="*", s=5)
plt.annotate("npore", (6/5.72+label_pad, 1/5.72), ha="left", va="center", color=f"C{i}")
i += 1

# verkko
plt.scatter(0, 0, color=f"C{i}", zorder=100)
plt.annotate("HP-compression", (label_pad, 0), ha="left", va="bottom", color=f"C{i}")

# points
plt.scatter(0.3, 0.1, color=f"k")
plt.annotate("A", (0.3+label_pad, 0.1), ha="left", va="center", color="k")
plt.scatter(1, 0.333, color=f"k")
plt.annotate("B", (1+label_pad, 0.333), ha="left", va="center", color="k")
plt.scatter(1.6, 0.4, color=f"k")
plt.annotate("C", (1.6+label_pad, 0.4), ha="left", va="center", color="k")
plt.scatter(1.5, 0.05, color=f"k")
plt.annotate("D", (1.5+label_pad, 0.05), ha="left", va="center", color="k")

# arrows
plt.arrow(2.5, 0.9,  0.1,  0.05, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # right
plt.annotate("prefer\nsubstitutions", (2.6+label_pad, 0.95), ha="left", va="center")
plt.arrow(2.5, 0.9, -0.1, -0.05, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # left
plt.annotate("prefer\nINDELs", (2.4-label_pad, 0.85), ha="right", va="center")
plt.arrow(2.5, 0.9, -0.05,  0.1, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # up
plt.annotate("prefer\nseparate\ngaps", (2.45, 1+label_pad), ha="center", va="bottom")
plt.arrow(2.5, 0.9,  0.05, -0.1, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # down
plt.annotate("prefer\nmerged\ngaps", (2.55, 0.8-label_pad), ha="center", va="top")

# labels
plt.xlabel(r'\LARGE{$\frac{o+e}{x}$}')
plt.ylabel(r'\LARGE{$\frac{e}{x}$}')
plt.tight_layout()
plt.savefig("img/2_general_sw_space.png")
