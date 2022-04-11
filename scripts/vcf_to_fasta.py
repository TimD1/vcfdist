import argparse, subprocess

from lib.vcf import *
from lib.cig import *

def main():

    print("> getting regions")
    get_vcf_regions()
    print(cfg.args.regions)

    print("> splitting vcf")
    filename = cfg.args.vcf[:-7]
    vcf1, vcf2 = split_vcf(cfg.args.vcf, vcf_out_pre=filename)

    print("> applying hap1")
    hap1_data = apply_vcf(vcf1, 1)

    print("> applying hap2")
    hap2_data = apply_vcf(vcf2, 2)

    print("> saving FASTA")
    fasta_filename = f"{filename}.fasta"
    first = True
    with open(fasta_filename, "w") as fasta:
        for (ctg1, hap1, seq1, ref1, cig1), (ctg2, hap2, seq2, ref2, cig2) in \
                zip(hap1_data, hap2_data):
            assert(ctg1 == ctg2)
            n = '\n' if not first else ''
            fasta.write(f"{n}>{ctg1}-hap{hap1}")
            fasta.write(f"\n{seq1}")
            fasta.write(f"\n>{ctg2}-hap{hap2}")
            fasta.write(f"\n{seq2}")
            first = False

    print("> saving SAM")
    sam_filename = f"{filename}_orig.sam"
    with open(sam_filename, "w") as sam:
        sam.write("@HD\tVN:1.6\tSO:coordinate")
        for ctg, start, stop in cfg.args.regions:
            sam.write(f"\n@SQ\tSN:{ctg}\tLN:{stop+1}")
        sam.write(f"\n@PG\tPN:vcfdist\tID:vcf_to_fasta\tVN:0.0.1\tCL:{' '.join(sys.argv)}")
        for (ctg1, hap1, seq1, ref1, cig1), (ctg2, hap2, seq2, ref2, cig2) in \
                zip(hap1_data, hap2_data):
            assert(ctg1 == ctg2)
            sam.write(f"\n{ctg1}-hap{hap1}\t0\t{ctg1}\t1\t60\t{collapse_cigar(cig1)}\t*\t0\t{ref_len(cig1)}\t{seq1}\t{'*'*len(seq1)}")
            sam.write(f"\n{ctg2}-hap{hap2}\t0\t{ctg2}\t1\t60\t{collapse_cigar(cig2)}\t*\t0\t{ref_len(cig2)}\t{seq2}\t{'*'*len(seq2)}")

    print("> converting to BAM")
    subprocess.Popen(["samtools", "view", "-b", sam_filename], 
            stdout=open(f"{filename}_orig.bam", "w"))




def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
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
