import argparse, subprocess
from time import perf_counter
from collections import defaultdict

from lib.vcf import *
from lib.cig import *

def main():

    print("> getting regions")
    get_bam_regions()
    print([f"{ctg}:{start}-{stop}" for (ctg, start ,stop) in cfg.args.regions])

    print("> extracting read data")
    try:
        bam = pysam.AlignmentFile(cfg.args.bam, "rb")
    except ValueError: 
        print(f"ERROR: BAM file '{cfg.args.bam}' contains no data.")
        exit(1)

    ctg1 = hap1 = pos1 = seq1 = ref1 = cig = None
    ctg2 = hap2 = pos2 = seq2 = ref2 = cig = None
    for ctg, beg, end in cfg.args.regions:
        for read in bam.fetch(ctg, beg, end):
            if read.query_name == "chr20-hap1":
                ctg1 = read.reference_name
                hap1 = 1
                pos1 = read.reference_start
                seq1 = read.query_alignment_sequence.upper()
                ref1 = read.get_reference_sequence().upper()
                cig1 = read.cigartuples
            elif read.query_name == "chr20-hap2":
                ctg2 = read.reference_name
                hap2 = 2
                pos2 = read.reference_start
                seq2 = read.query_alignment_sequence.upper()
                ref2 = read.get_reference_sequence().upper()
                cig2 = read.cigartuples
            else:
                print("ERROR: unexpected haplotype.")

    if ctg1 == None:
        print("ERROR: haplotypes unaligned.")
        exit(1)

    hap1_data = [(ctg1, hap1, pos1, seq1, ref1, cig1)]
    hap2_data = [(ctg2, hap2, pos2, seq2, ref2, cig2)]

    print('> generating standardized vcfs')
    filename = cfg.args.bam[:-4]
    vcf1 = gen_vcf(hap1_data, 1, filename)
    vcf2 = gen_vcf(hap2_data, 2, filename)

    print(f"> merging vcfs")
    out_fn = f"{filename}.vcf.gz"
    merge_vcfs(vcf1, vcf2, out_fn)
    subprocess.run(['tabix', '-f', '-p', 'vcf', out_fn])




def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("ref")
    parser.add_argument("contig")
    parser.add_argument("--contigs")
    parser.add_argument("--contig_beg")
    parser.add_argument("--contig_end")
    parser.add_argument("--min_qual", type=int, default=0)
    return parser


if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
