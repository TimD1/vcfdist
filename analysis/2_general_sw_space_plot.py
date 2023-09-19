import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from numpy import sin, cos, pi, linspace
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
plt.axvline(0.5, alpha=0.8, color='k', zorder=-100)
# plt.annotate("$3x < o+e$", (3.01,1.22), rotation=90, color="black")
# plt.axvspan(3, 3.2, alpha=0.1, color='k', zorder=-100)

# horizontal
# plt.annotate("$x < 2e$", (label_pad,0.57), ha="left", va="top", color="white")
# plt.plot(x, [0.5]*len(x), color='k', alpha=0.8)
plt.annotate("$e < 0$", (label_pad,-0.03), ha="left", va="top", color="white")
plt.fill_between(x, -0.2, 0, facecolor='k', alpha=0.8)

# min ED
i = 0
plt.scatter(1,1, marker="X")
plt.annotate("edlib", (1+label_pad,1), ha="left", va="center", color=f"C{i}")
i += 1

# minimap2
mm2_presets = {
        r"mm2 \texttt{sr}":     { "a": 0, "b": 10, "o1": 13, "e1":3, "o2":25,  "e2": 2, "m":"^" },
        r"\noindent mm2 \texttt{map-ont}\\winnowmap":{ "a": 0, "b": 6,  "o1": 5,  "e1":3, "o2":25, "e2": 2, "m":"s" },
        r"mm2 \texttt{map-pb}": { "a": 0, "b": 10, "o1": 13, "e1":5, "o2":53, "e2": 3, "m":"s" },
        r"mm2 \texttt{asm5}":   { "a": 0, "b": 40, "o1": 79, "e1":7, "o2":163, "e2": 3, "m":"P" },
        r"mm2 \texttt{asm10}":  { "a": 0, "b": 20, "o1": 33, "e1":5, "o2":83, "e2": 3, "m":"P" },
    "pbmm2":                { "a": 0, "b": 7, "o1": 6, "e1":5, "o2":57, "e2": 3, "m":"s" },
}
for name, x in mm2_presets.items():
    start = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"], color=f"C{i}", marker=x["m"])
    end = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e2"]/x["b"], color=f"C{i}", marker=x["m"])
    plt.arrow((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"]-0.02, 0, x["e2"]/x["b"]-x["e1"]/x["b"]+0.04, color='k',
            length_includes_head=True, head_width=0.02, head_length = 0.016)
    plt.annotate(name, ((x["o1"] + x["e1"])/x["b"]+label_pad, x["e1"]/x["b"]), ha="left", va="center", color=f"C{i}")
    i += 1

# BWA, GRAF, DRAGEN
plt.scatter(1.6, 0.3, color=f"C{i}", marker="^", zorder=100)
plt.annotate(r"\noindent DRAGEN\\ BWA\\GRAF", (1.6+label_pad, 0.275), ha="right", va="center", color=f"C{i}")
i += 1

# ngmlr
ngmlr_presets = {
    r"ngmlr \texttt{ont}":  { "a": 0, "b": 12, "o": 11, "e_max": 11, "e_min": 3, "e_dec": 0.15},
    r"ngmlr \texttt{hifi}": { "a": 0, "b": 12, "o": 5,  "e_max": 5,  "e_min": 4, "e_dec": 0.15}
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
    plt.scatter(p[0]/5.72, p[1]/5.72, color=f"purple")
plt.annotate("npore", (0.5, 0.25), ha="left", va="center", color=f"purple")
plt.scatter(6/5.72, 1/5.72, marker="s", color=f"purple")
plt.annotate("npore", (6/5.72+label_pad, 1/5.72), ha="left", va="center", color=f"purple")
i += 1

# verkko
plt.scatter(0, 0, color=f"goldenrod", zorder=100)
plt.annotate("verkko", (label_pad, 0), ha="left", va="bottom", color=f"goldenrod")

# points
plt.scatter(0.4, 0.3, marker="*", color=f"k")
plt.annotate("A", (0.4+label_pad, 0.3), ha="left", va="center", color="k")
plt.scatter(1, 0.333, marker="*", color=f"k")
plt.annotate("B", (1+label_pad, 0.333), ha="left", va="center", color="k")
plt.scatter(1.6, 0.4, marker="*", color=f"k")
plt.annotate("C", (1.6+label_pad, 0.4), ha="left", va="center", color="k")
plt.scatter(2, 0.2, marker="*", color=f"k")
plt.annotate("D", (2+label_pad, 0.2), ha="left", va="center", color="k")

# arrows
plt.arrow(2.1, 1.1,  0.1,  0.05, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # right
plt.annotate("prefer\nsubstitutions", (2.2+label_pad, 1.15), ha="left", va="center")
plt.arrow(2.1, 1.1, -0.1, -0.05, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # left
plt.annotate("prefer\nINDELs", (2.0-label_pad, 1.05), ha="right", va="center")
plt.arrow(2.065, 1.175, -0.01,  0.01, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # up
plt.annotate("prefer\nseparate\ngaps", (2.05, 1.2+label_pad), ha="center", va="bottom")
plt.arrow(2.13, 1.05,  0, -0.02, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # down
plt.annotate("prefer\nmerged\ngaps", (2.15, 1-label_pad), ha="center", va="top")
arc_angles = linspace(0 * pi, pi/4, 20)
r = 0.2
arc_xs = r * cos(arc_angles)
arc_ys = r * sin(arc_angles)
plt.plot(arc_xs+1.93, arc_ys+1.03, color = 'k')

# legend
sra = mlines.Line2D([], [], color='k', marker='^', linestyle='None',
        markersize=7, label = 'Short Read Aligners')
lra = mlines.Line2D([], [], color='k', marker='s', linestyle='None',
        markersize=7, label = 'Long Read Aligners')
aa = mlines.Line2D([], [], color='k', marker='P', linestyle='None',
        markersize=7, label = 'Assembly Aligners')
sva = mlines.Line2D([], [], color='k', marker='o', linestyle='None',
        markersize=7, label = 'SV/CNV Aligners')
ed = mlines.Line2D([], [], color='k', marker='X', linestyle='None',
        markersize=7, label = 'Edit Distance Align')
points = mlines.Line2D([], [], color='k', marker='*', linestyle='None',
        markersize=9, label = 'Select Design Points')
leg1 = plt.legend(handles = [sra, lra, aa, sva, ed, points], loc = (0.8,0.2))

m = mlines.Line2D([], [], linestyle='None', label = r'$m$ : match penalty ($0$)')
x = mlines.Line2D([], [], linestyle='None', label = r'$x$ : mismatch penalty')
o = mlines.Line2D([], [], linestyle='None', label = r'$o$ : gap opening penalty')
e = mlines.Line2D([], [], linestyle='None', label = r'$e$ : gap extension penalty')
plt.legend(handles = [m, x, o, e], loc=(0.78, 0.6), frameon=False)
plt.gca().add_artist(leg1)

# labels
plt.xlabel(r'\LARGE{$(o+e) / x$}')
plt.ylabel(r'\LARGE{$e / x$}')
plt.tight_layout()
plt.savefig("img/2_general_sw_space.png", dpi=200)
