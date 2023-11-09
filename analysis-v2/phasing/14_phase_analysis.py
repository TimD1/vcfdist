from collections import defaultdict

prefix="../out/"
prefix="../data/pfda-v2/nist_vcfdist/K4GT3_HG002_C."

# get phase block phasings
pb_phases = defaultdict(list)
pb_fn = f"{prefix}phase-blocks.tsv"
pb_file = open(pb_fn, "r")
next(pb_file)
for line in pb_file:
    ctg, start, stop, size, scs, phase = line.strip().split()
    pb_phases[ctg].append(int(phase))
print(pb_phases)

# count incorrectly phased and total superclusters with 1 or 2 variants
sc_fn = f"{prefix}superclusters.tsv"
sc_file = open(sc_fn, "r")
next(sc_file)
sc_phase_errs = defaultdict(int)
sc_counts = defaultdict(int)
best_tot_ed = 0
orig_tot_ed = 0
swap_tot_ed = 0
real_tot_ed = 0
for line in sc_file:
    ctg, start, stop, size, q1v, q2v, t1v, t2v, orig_ed, swap_ed, phase, pb = line.strip().split()
    vc1 = min(int(q1v), int(q2v)) # variant counts
    vc2 = max(int(q1v), int(q2v))
    sc_counts[f"{vc1},{vc2}"] += 1
    pb, orig_ed, swap_ed = int(pb), int(orig_ed), int(swap_ed)

    # calculate total edit distances
    best_tot_ed += min(orig_ed, swap_ed)
    orig_tot_ed += orig_ed
    swap_tot_ed += swap_ed
    if pb_phases[ctg][pb] == 1:
        real_tot_ed += swap_ed
    else:
        real_tot_ed += orig_ed

    # calculate phasing errors
    if phase == "?":
        continue
    elif phase == "=":
        if pb_phases[ctg][pb] == 1:
            sc_phase_errs[f"{vc1},{vc2}"] += 1
    elif phase == "X":
        if pb_phases[ctg][pb] == 0:
            sc_phase_errs[f"{vc1},{vc2}"] += 1
    else:
        print("ERROR: unexpected phase")
print("sc phase errors:", sc_phase_errs)
print("sc counts:", sc_counts)
print("best_tot_ed:", best_tot_ed)
print("orig_tot_ed:", orig_tot_ed)
print("swap_tot_ed:", swap_tot_ed)
print("real_tot_ed:", real_tot_ed)
