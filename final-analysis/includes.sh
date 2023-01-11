#!/bin/bash

data="/x/vcfdist/data"
# data="/home/timdunn/vcfdist/data"
parallel="/home/timdunn/parallel-20211222/src/parallel"
vcfdist="/home/timdunn/vcfdist/src/vcfdist"

# get IDs for each pFDA v2 submission VCF
# sub_ids=( $(ls $data/pfda-v2/submission_vcfs) )
# sub_ids=("0GOOR" "23O09" "4GKR1" "61YRJ" "8H0ZB" "B1S5A" "C6JUX" "EIUT6" "H9OJ3" "IA789") #mbit3
# sub_ids=("KXBR8" "MECML" "NWQ6Y" "R9CXN" "SEX9X" "UYMUW" "W91C1" "XC97E" "YBE9U" "YUI27") #mbit4
# sub_ids=("4HL0B" "7NJ5C" "9DGOR" "BARQS" "CN36L" "ES3XW" "HB8P3" "ISKOS" "JIPSI" "KFA0T") #mbit5
# sub_ids=("PGXA4" "RU88N" "TG5TE" "VES2R" "WGQ43" "XV7ZN" "YGOTK" "13678" "32LOW" "60Z59" "JBBK0" "K4GT3" "7RR8Z" "ASJT6" "BSODP" "CZA1Y" "0O7FL" "2OT9Q" "FFFGB" "HF8CT" "Y8QR8" "YJN61" "LR1FD" "MT57N") # mbit7
# sub_ids=("J04SL" "K33QJ" "KPXT2" "M9KLP" "NFT0L" "QUE7Q" "S7K7S" "TZMTX" "W607K" "WX8VK") # mbit8

sub_ids=("K4GT3")

# define shorthand IDs for reference sequences (used for directories)
ref_id="GRCh38"

# define full filenames for corresponding FASTAs
ref_name="GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
truth_names=(
    "nist-v4.2.1"
    "cmrg-v1.00"
)
truth_ids=(
    "nist"
    "cmrg"
)
truth_vcfs=(
    "HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz"
    "HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz"
)
truth_beds=(
    "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
    "HG002_GRCh38_CMRG_smallvar_v1.00.bed"
)

# shorthand IDs for each sequencing dataset
read_ids=(
    "pacbio-hifi"
    "nanopore-prom"
    "illumina-pairs1"
    "illumina-pairs2"
)

# full filenames for corresponding FASTQs
read_names=(
    "HG002_35x_PacBio_14kb-15kb.fastq.gz"
    "HG002_GM24385_1_2_3_Guppy_3.6.0_prom.fastq.gz"
    "HG002.novaseq.pcr-free.35x.R1.fastq.gz"
    "HG002.novaseq.pcr-free.35x.R2.fastq.gz"
)

# save read types for minimap2 alignment
read_types=(
    "map-pb"
    "map-ont"
    "sr"
    "sr"
)

# data points for evaluation
pts=(
    "A"
    "B"
    "C"
    "D"
)
x=(
    "10"
    "3"
    "5"
    "5"
)
o=(
    "1"
    "2"
    "6"
    "9"
)
e=(
    "3"
    "1"
    "2"
    "1"
)
# all representations for evaluation (include original)
reps=(
    "O"
    "A"
    "B"
    "C"
    "D"
)
