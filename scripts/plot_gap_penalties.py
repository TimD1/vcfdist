import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(10,10))
max_gap = 50

# affine: minimap2
mm2_presets = {
        "short":    { "a": 0, "b": 10, "o1": 14, "e1":4, "o2":26, "e2": 3 },
        "map-ont":  { "a": 0, "b": 6,  "o1": 6,  "e1":4, "o2":26, "e2": 3 },
        "map-hifi": { "a": 0, "b": 5,  "o1": 7,  "e1":3, "o2":27, "e2": 2 },
        "asm5":     { "a": 0, "b": 20, "o1": 40, "e1":4, "o2":82, "e2": 2 },
        "asm10":    { "a": 0, "b": 10, "o1": 17, "e1":3, "o2":42, "e2": 2 },
}
for pre, val in mm2_presets.items():
    plt.plot(range(max_gap), [0] + [min(val["o2"] + val["e2"]*i, 
        val["o1"] + val["e1"]*i)/val["b"] for i in range(1,max_gap)])

# affine: npore
plt.plot(range(max_gap), [0] + [(4 + i)/5 for i in range(1,max_gap)])

# linear: edit/indel dist
plt.plot(range(max_gap), range(max_gap))
plt.plot(range(max_gap), [y/2 for y in range(max_gap)])

# convex: ngmlr
ngmlr_presets = {
        "ont": { "a": 0, "b": 6, "o": 4, "e_max": 4, "e_min": 3.5, "e_dec": 0.15},
        "pb": {  "a": 0, "b": 7, "o": 7, "e_max": 7, "e_min": 3, "e_dec": 0.15},
        }
for pre, val in ngmlr_presets.items():

    e_val = val["e_max"]
    gap_scores = [0, (val["o"]+e_val)/val["b"]]
    for x in range(2, max_gap):
        e_val = max(val["e_min"], e_val - val["e_dec"])
        gap_scores.append(gap_scores[-1] + e_val/val["b"])
    plt.plot(range(max_gap), gap_scores)

# actual: measured
plt.plot(range(1,max_gap), 
        # from polyfit2 
        [-0.0003596*x**2+0.036792*x+1.5675427 for x in range(1,max_gap)])


# plot
plt.legend(
        [f"minimap2 {x}" for x in mm2_presets.keys()] + 
        ["npore"] +
        ["edit-dist", "indel-dist"] + 
        [f"ngmlr {x}" for x in ngmlr_presets.keys()] +
        ["measured HG002-GRCh38"]
        )
plt.xlabel('Gap Length')
plt.xlim(0,20)
plt.ylim(0,10)
plt.xticks(range(21))
plt.yticks(range(11))
plt.ylabel('Gap Penalty (relative to substitution)')
plt.tight_layout()
plt.savefig("img/gap_penalties.png")
