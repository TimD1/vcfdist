import matplotlib.pyplot as plt
import numpy as np

def get_cluster_sizes(var_filename):

    cluster_vars = []  # variants per cluster
    cluster_sizes = [] # bases per cluster
    with open(var_filename, "r") as varfile:
        first_position = [0, 0]
        last_position = [0, 0]
        dependencies = set()
        cluster = [-1, -1]
        variants = [0, 0]
        prev_ctg = ""
        next(varfile) # skip header
        for variant in varfile:
            ctg, pos, hap, ref, alt, q, typ, errtype, credit, orig_gt, \
                    clust, supclust, loc = variant.split("\t")
            hap = int(hap)
            pos = int(pos)
            clust = int(clust)

            # new contig
            if ctg != prev_ctg:
                first_position = [0, 0]
                last_position = [0, 0]
                cluster = [-1, -1]
                variants = [0, 0]
                prev_ctg = ctg

            # same cluster
            if clust == cluster[hap]:
                dependencies.add(f"{ctg}_{hap}:{last_position[hap]}:{pos}")
                variants[hap] += 1
                last_position[hap] = pos + len(ref)

            else: # new cluster
                # save previous cluster stats
                if cluster[hap] != -1:
                    cluster_vars.append(variants[hap])
                    cluster_sizes.append(last_position[hap]-first_position[hap])

                # reset
                first_position[hap] = pos
                last_position[hap] = pos + len(ref)
                cluster[hap] = clust
                variants[hap] = 1

        cluster_vars.append(variants[hap])
        cluster_sizes.append(last_position[hap]-first_position[hap])

    return cluster_vars, cluster_sizes, dependencies

        
gap_filename = "3_results/gap_truth.tsv"
sw_filename = "3_results/sw_truth.tsv"

sw_vars, sw_bases, sw_deps =  get_cluster_sizes(sw_filename)
gap_vars, gap_bases, gap_deps =  get_cluster_sizes(gap_filename)

fig, (vars_ax, bases_ax) = plt.subplots(1, 2, figsize=(10,6))

max_vars = max(sw_vars + gap_vars)
vars_ax.hist(sw_vars, histtype="step", bins=range(1, max_vars+1))
vars_ax.hist(gap_vars, histtype="step", bins=range(1, max_vars+1))
vars_ax.set_xlabel("Variants in Cluster")
vars_ax.set_ylabel("Counts")
vars_ax.set_yscale("log")

max_bases = max(sw_bases + gap_bases)
bases_ax.hist(sw_bases, histtype="step", bins=np.linspace(1, max_bases, 50))
bases_ax.hist(gap_bases, histtype="step", bins=np.linspace(1, max_bases, 50))
bases_ax.set_xlabel("Bases in Cluster")
bases_ax.set_ylabel("Counts")
bases_ax.set_yscale("log")

print(f"\ntotal variants: {sum(gap_vars)}")
print(f"gap clusters: {len(gap_vars)}")
print(f"s-w clusters: {len(sw_vars)}")

print(f"\ndeps broken: {len(gap_deps-sw_deps)}")
print(f"deps added: {len(sw_deps-gap_deps)}")
# for x in sw_deps-gap_deps:
#     print(f"    {x}")

print(f"\ngap avg variants: {np.mean(gap_vars):.2f}")
print(f"s-w avg variants: {np.mean(sw_vars):.2f}")

print(f"\ngap avg bases: {np.mean(gap_bases):.2f}")
print(f"s-w avg bases: {np.mean(sw_bases):.2f}")

print(f"\ngap total bases: {np.sum([x*x for x in gap_bases])}")
print(f"s-w total bases: {np.sum([x*x for x in sw_bases])}")
print(f"reduction factor: {np.sum([x*x for x in gap_bases])/np.sum([x*x for x in sw_bases]):.2f}x\n")

plt.legend(["Smith-Waterman Clustering", "50 Base Gap Clustering"])
plt.savefig("img/3_cluster_sizes.png")
