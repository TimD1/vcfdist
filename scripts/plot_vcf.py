import pysam
import matplotlib.pyplot as plt
import numpy as np

vcf_filename = "../r10.4_chr20/ont-case-study/input/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
vcf_file = pysam.VariantFile(vcf_filename, 'r')

prev_contig = prev_start = prev_stop = None
curr_contig = curr_start = curr_stop = None
distances = []            # between each variant
variant_region_bases = [] # bases in region
variant_region_sizes = [] # variants in region
insertion_sizes = []        # INDEL sizes
deletion_sizes = []        # INDEL sizes
variant_qualities = []    # variant Q-scores
variant_types = [0]*9     # SUB, INS, DEL
variant_genotypes = [0]*6 # phasings
types = ["S", "I", "D", "S,S", "S,I", "S,D", "I,I", "I,D", "D,D"]
SUB = 0
INS = 1
DEL = 2
SUB_SUB = 3
SUB_INS = 4
SUB_DEL = 5
INS_INS = 6
INS_DEL = 7
DEL_DEL = 8
phasings = ["1|1", "1|0", "0|1", "0|0", "1|2", "2|1"]
_1_1 = 0
_1_0 = 1
_0_1 = 2
_0_0 = 3
_1_2 = 4
_2_1 = 5


def get_type(phasing, alleles):
    ''' Determine SUB/INS/DEL from phasing and alleles. '''

    if phasing == _0_0:
        print(f"ERROR, unexpected phasing: {phasings[phasing]} {alleles}")
        exit()

    if phasing == _1_2 or phasing == _2_1: # two ALTs
        if len(alleles) != 3:
            print(f"ERROR, expected 3 alleles: {phasings[phasing]} {alleles}")
            exit()
        indel_len1 = len(alleles[1]) - len(alleles[0])
        indel_len2 = len(alleles[2]) - len(alleles[0])
        if indel_len1 < 0 and indel_len2 < 0:
            return DEL_DEL
        elif indel_len1 > 0 and indel_len2 < 0 or \
             indel_len1 < 0 and indel_len2 > 0:
            return INS_DEL
        elif indel_len1 > 0 and indel_len2 > 0:
            return INS_INS
        elif indel_len1 == 0 and indel_len2 > 0 or \
             indel_len2 == 0 and indel_len1 > 0:
                 return SUB_INS
        elif indel_len1 == 0 and indel_len2 < 0 or \
             indel_len2 == 0 and indel_len1 < 0:
                 return SUB_DEL
        else:
            if len(alleles[0]) != 1:
                print(f"ERROR, expected SUB: {phasings[phasing]} {alleles}")
                exit()

            if (alleles[0] == alleles[1] or
                   alleles[0] == alleles[2] or
                   alleles[1] == alleles[2]):
                print(f"ERROR, expected different SUBs: {phasings[phasing]} {alleles}")
                exit()

            return SUB_SUB

    else: # single ALT

        if len(alleles) != 2:
            print(f"WARNING, expected two alleles: {phasings[phasing]} {alleles}")

        indel_len = len(alleles[1]) - len(alleles[0])
        if indel_len < 0:
            deletion_sizes.append(-indel_len)
            return DEL
        elif indel_len > 0:
            insertion_sizes.append(indel_len)
            return INS
        else:
            if len(alleles[0]) != 1:
                print(f"ERROR, expected SUB: {phasings[phasing]} {alleles}")
                exit()
            return SUB



def get_phasing(gt):
    if   gt[0] == 1 and gt[1] == 1: return _1_1
    elif gt[0] == 1 and gt[1] == 0: return _1_0
    elif gt[0] == 0 and gt[1] == 1: return _0_1
    elif gt[0] == 0 and gt[1] == 0: return _0_0
    elif gt[0] == 1 and gt[1] == 2: return _1_2
    elif gt[0] == 2 and gt[1] == 1: return _2_1
    else:
        print(f"ERROR: unexpected phasing ({gt[0]}, {gt[1]})")
        exit()



region_bases = 0
region_size = 0
for record in vcf_file.fetch():

    # get genotype
    for sample in record.samples:
        curr_gt = record.samples[sample]['GT']
        break

    # categorize phasing
    curr_phasing = get_phasing(curr_gt)
    variant_genotypes[curr_phasing] += 1

    # unpack variant info
    curr_contig = record.contig
    curr_start = record.start
    curr_stop = record.stop

    # save quality
    if not record.qual:
        variant_qualities.append(-1)
    else:
        variant_qualities.append(record.qual)

    # categorize type
    curr_type = get_type(curr_phasing, record.alleles)
    variant_types[curr_type] += 1

    # calculate region-based stats
    if prev_start is None: # first variant
        distances.append(curr_start)
        region_bases = curr_stop - curr_start
        region_size = 1

    elif prev_contig != curr_contig: # new contig
        variant_region_bases.append(region_bases)
        variant_region_sizes.append(region_size)
        distances.append(curr_start)
        region_bases = curr_stop - curr_start
        region_size = 1

    else: # same contig
        distances.append(curr_start - prev_stop)

        # variants within 50 in same region
        if curr_start - prev_stop < 50:
            region_bases += curr_stop - prev_stop
            region_size += 1
        else:
            variant_region_bases.append(region_bases)
            variant_region_sizes.append(region_size)
            region_bases = curr_stop - curr_start
            region_size = 1

    # set data for next variant
    prev_contig = curr_contig
    prev_start = curr_start
    prev_stop = curr_stop

variant_region_bases.append(region_bases)
variant_region_sizes.append(region_size)


plt.figure()
plt.hist(variant_qualities, bins=50)
plt.ylabel('Counts')
plt.xlabel('Variant Quality')
plt.yscale('log')
plt.savefig("img/variant_qualities.png")

plt.figure()
plt.hist(insertion_sizes, 
        bins=np.logspace(0, np.log10(max(insertion_sizes)), 50))
plt.ylabel('Counts')
plt.xlabel('Insertion Size')
plt.yscale('log')
plt.xscale('log')
plt.savefig("img/insertion_sizes.png")

plt.figure()
plt.hist(deletion_sizes, 
        bins=np.logspace(0, np.log10(max(deletion_sizes)), 50))
plt.ylabel('Counts')
plt.xlabel('Deletion Size')
plt.yscale('log')
plt.xscale('log')
plt.savefig("img/deletion_sizes.png")

plt.figure()
plt.hist(distances, 
        bins=np.logspace(0, np.log10(max(distances)), 50))
plt.ylabel('Counts')
plt.xlabel('Distance between variants')
plt.yscale('log')
plt.xscale('log')
plt.savefig("img/distances.png")

plt.figure()
plt.hist(variant_region_bases,
        bins=np.logspace(0, np.log10(max(variant_region_bases)), 50))
plt.ylabel('Counts')
plt.xlabel('Bases in Region')
plt.yscale('log')
plt.xscale('log')
plt.savefig("img/variant_region_bases.png")

plt.figure()
plt.hist(variant_region_sizes, 
        bins=np.logspace(0, np.log10(max(variant_region_sizes)), 50))
plt.xlabel('Variants in Region')
plt.ylabel('Counts')
plt.yscale('log')
plt.xscale('log')
plt.savefig("img/variant_region_sizes.png")

plt.figure()
plt.bar(range(9), variant_types)
plt.xticks(range(9), labels=types)
plt.ylabel('Counts')
plt.yscale('log')
plt.xlabel('Variant Type')
plt.savefig("img/variant_types.png")

plt.figure()
plt.bar(range(6), variant_genotypes)
plt.xticks(range(6), labels=phasings)
plt.yscale('log')
plt.ylabel('Counts')
plt.xlabel('Phasing')
plt.savefig("img/variant_genotypes.png")
