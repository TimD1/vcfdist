import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from numpy import sin, cos, pi, linspace
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['text.usetex'] = True

def mm2in(x):
    return x / 25.4

plt.figure(figsize=(mm2in(180), mm2in(75)))
xmin, xmax = 0, 3.2
ymin, ymax = -0.1, 1.5
label_pad = 0.03
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
x = list(np.arange(xmin, xmax+1, 0.1))

# G_open > 0
plt.annotate("$o < 0$", (1.31,1.36), rotation=35, color="white", fontsize=7)
plt.fill_between(x, x, ymax, alpha=0.8, facecolor='k')

# vertical
plt.annotate("$2(o+e) < x$", (0.45,1.08), rotation=90, color="white", fontsize=7)
plt.axvline(0.5, alpha=0.8, color='k', zorder=-100)
# plt.annotate("$3x < o+e$", (3.01,1.22), rotation=90, color="black")
# plt.axvspan(3, 3.2, alpha=0.1, color='k', zorder=-100)

# horizontal
# plt.annotate("$x < 2e$", (label_pad,0.57), ha="left", va="top", color="white")
# plt.plot(x, [0.5]*len(x), color='k', alpha=0.8)
plt.annotate("$e < 0$", (label_pad,-0.03), ha="left", va="top", color="white", fontsize=7)
plt.fill_between(x, -0.2, 0, facecolor='k', alpha=0.8)

# min ED
i = 0
plt.scatter(1,1, marker="X", s=7)
plt.annotate("edlib", (1+label_pad,1), ha="left", va="center", color=f"C{i}", fontsize=7)
i += 1

# minimap2
mm2_presets = {
    r"mm2 \texttt{sr}":     
        { "a": 0, "b": 10, "o1": 13, "e1":3, "o2":25,  "e2": 2, "m":"^", "l":True },
    r"\noindent mm2 \texttt{map-ont}\vspace{-2mm}\\winnowmap":
        { "a": 0, "b": 6,  "o1": 5,  "e1":3, "o2":25, "e2": 2, "m":"s", "l":False },
    r"mm2 \texttt{map-pb}": 
        { "a": 0, "b": 10, "o1": 13, "e1":5, "o2":53, "e2": 3, "m":"s", "l":True },
    r"mm2 \texttt{asm5}":   
        { "a": 0, "b": 40, "o1": 79, "e1":7, "o2":163, "e2": 3, "m":"P", "l":False },
    r"mm2 \texttt{asm10}":  
        { "a": 0, "b": 20, "o1": 33, "e1":5, "o2":83, "e2": 3, "m":"P", "l":True },
    r"mm2 \texttt{asm20}":  
        { "a": 0, "b": 10, "o1": 13, "e1":5, "o2":53, "e2": 3, "m":"P", "l":False },
    "pbmm2":                
        { "a": 0, "b": 7, "o1": 6, "e1":5, "o2":57, "e2": 3, "m":"s", "l":False },
}
for name, x in mm2_presets.items():
    start = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"], color=f"C{i}", marker=x["m"], s=7)
    end = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e2"]/x["b"], color=f"C{i}", marker=x["m"], s=7)
    plt.arrow((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"]-0.02, 0, x["e2"]/x["b"]-x["e1"]/x["b"]+0.04, color='k',
            length_includes_head=True, head_width=0.02, head_length = 0.016)
    plt.annotate(name, 
            ((x["o1"] + x["e1"])/x["b"]+(label_pad if x["l"] else -label_pad), 
                x["e1"]/x["b"]), ha="left" if x["l"] else "right", va="center", 
            color=f"C{i}", fontsize=7)
    i += 1

# BWA, GRAF, DRAGEN
plt.scatter(1.6, 0.3, color=f"C{i}", marker="^", zorder=100, s=7)
plt.annotate(r"BWA", (1.6-label_pad, 0.275), ha="right", va="center", color=f"C{i}", fontsize=7)
i += 1

# ngmlr
ngmlr_presets = {
    r"ngmlr \texttt{ont}":  { "a": 0, "b": 12, "o": 11, "e_max": 11, "e_min": 3, "e_dec": 0.15},
    r"ngmlr \texttt{hifi}": { "a": 0, "b": 12, "o": 5,  "e_max": 5,  "e_min": 4, "e_dec": 0.15}
}
for name, x in ngmlr_presets.items():
    for e in np.arange(x["e_max"], x["e_min"], -x["e_dec"]):
        plt.scatter(x["o"]/x["b"], e/x["b"], color=f"C{i}", s=7)
    plt.annotate(name, (x["o"]/x["b"]-label_pad, (x["e_max"]+x["e_min"])/(2*x["b"])), ha="right", va="center", color=f"C{i}", fontsize=7)
    i += 1
    plt.arrow(x["o"]/x["b"], x["e_max"]/x["b"], 0, x["e_min"]/x["b"]-x["e_max"]/x["b"], color='k',
            length_includes_head=True, head_width=0.02, head_length = 0.016)

# npore
points = [[1.7, 0.5], [1.7,0.7], [1.7, 1], [1.7, 0.8], [1.7, 1.1], [2.3, 1.8], [2.8, 1], [3.2, 1.1], [3.2, 2.1], [3.2, 1], [3.7, 1.9], [3.7, 2], [3.7, 1], [3.7, 0.6], [3.8, 2.4], [3.8, 1.7], [3.8, 0.8], [4.6, 2.5], [4.6, 2.3], [4.6, 1.5]]
for p in points:
    plt.scatter(p[0]/5.72, p[1]/5.72, color=f"purple", s=7)
plt.annotate("npore", (0.5, 0.25), ha="left", va="center", color=f"purple", fontsize=7)
plt.scatter(6/5.72, 1/5.72, marker="s", color=f"purple", s=7)
plt.annotate("npore", (6/5.72+label_pad, 1/5.72), ha="left", va="center", color=f"purple", fontsize=7)
i += 1

# verkko
plt.scatter(0, 0, color=f"goldenrod", zorder=100, s=7)
plt.annotate("verkko", (label_pad, 0), ha="left", va="bottom", color=f"goldenrod", fontsize=7)

# VCFs
vcf_scores = {
        r"\noindent Q100-dipcall\vspace{-2mm}\\hifiasm-dipcall":     
            { "a": 0, "b": 40, "o1": 79, "e1":7, "o2":163, "e2": 3},
        r"Q100-PAV":
            { "a": 0, "b": 12, "o1": 11, "e1":9, "o2":113, "e2": 3},
        r"hifiasm-GIAB-TR":
            { "a": 0, "b": 10, "o1": 11, "e1":9, "o2":113, "e2": 3},
}
for name, x in vcf_scores.items():
    start = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"], color=f"k", marker='*', s=7)
    end = plt.scatter((x["o1"] + x["e1"])/x["b"], x["e2"]/x["b"], color=f"k", marker='*', s=7)
    plt.arrow((x["o1"] + x["e1"])/x["b"], x["e1"]/x["b"]-0.02, 0, x["e2"]/x["b"]-x["e1"]/x["b"]+0.04, color='k',
            length_includes_head=True, head_width=0.02, head_length = 0.016)
    plt.annotate(name, ((x["o1"] + x["e1"])/x["b"]+label_pad, x["e1"]/x["b"]), ha="left", va="center", color='k', fontsize=7)
    i += 1

# arrows
plt.arrow(2.75, 1.1,  0.1,  0.05, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # right
plt.annotate("prefer\nsubstitutions", (2.85+label_pad, 1.15), ha="left", va="center", fontsize=7)
plt.arrow(2.75, 1.1, -0.1, -0.05, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # left
plt.annotate("prefer\nINDELs", (2.65-label_pad, 1.05), ha="right", va="center", fontsize=7)
plt.arrow(2.715, 1.175, -0.01,  0.01, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # up
plt.annotate("prefer\nseparate\ngaps", (2.7, 1.2+label_pad), ha="center", va="bottom", fontsize=7)
plt.arrow(2.78, 1.05,  0, -0.02, color='k', length_includes_head=True, head_width=0.02, head_length=0.016) # down
plt.annotate("prefer\nmerged\ngaps", (2.8, 1-label_pad), ha="center", va="top", fontsize=7)
arc_angles = linspace(0 * pi, pi/4, 20)
r = 0.2
arc_xs = r * cos(arc_angles)
arc_ys = r * sin(arc_angles)
plt.plot(arc_xs+2.58, arc_ys+1.03, color = 'k')

# legend
sra = mlines.Line2D([], [], color='k', marker='^', linestyle='None',
        markersize=6, label = 'Short Read Aligners')
lra = mlines.Line2D([], [], color='k', marker='s', linestyle='None',
        markersize=6, label = 'Long Read Aligners')
aa = mlines.Line2D([], [], color='k', marker='P', linestyle='None',
        markersize=6, label = 'Assembly Aligners')
sva = mlines.Line2D([], [], color='k', marker='o', linestyle='None',
        markersize=6, label = 'SV/CNV Aligners')
ed = mlines.Line2D([], [], color='k', marker='X', linestyle='None',
        markersize=6, label = 'Edit Distance Align')
points = mlines.Line2D([], [], color='k', marker='*', linestyle='None',
        markersize=8, label = 'Select Design Points')
leg1 = plt.legend(handles = [sra, lra, aa, sva, ed, points], loc = (0.8,0.1), fontsize=7)

m = mlines.Line2D([], [], linestyle='None', label = r'$m$ : match penalty ($0$)')
x = mlines.Line2D([], [], linestyle='None', label = r'$x$ : mismatch penalty')
o = mlines.Line2D([], [], linestyle='None', label = r'$o$ : gap opening penalty')
e = mlines.Line2D([], [], linestyle='None', label = r'$e$ : gap extension penalty')
plt.legend(handles = [m, x, o, e], loc=(0.5, 0.7), frameon=False, fontsize=7)
plt.gca().add_artist(leg1)

# labels
plt.xlabel(r'\textbf{Starting a Gap:} $(o+e) / x$', fontsize=7)
plt.ylabel(r'\textbf{Extending a Gap:} $e / x$', fontsize=7)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.tight_layout()
plt.savefig("sw_space.pdf", format="pdf")
