#!/bin/bash

data="/home/timdunn/vcfdist/data"
out="/home/timdunn/vcfdist/analysis-v2/phab_cm"
parallel="/home/timdunn/parallel-20211222/src/parallel"
vcfdist="/home/timdunn/vcfdist/src/vcfdist"
rtg="/home/timdunn/software/rtg-tools-3.12.1/rtg"
tabix="/home/timdunn/software/htslib-1.16/tabix"
bgzip="/home/timdunn/software/htslib-1.16/bgzip"
timer="/usr/bin/time"

# define reference FASTA
ref_name="GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"

# T2T information
truth_name="t2t-q100"
truth_version="v0.9"

# query VCF information
query_names=(
    # "hprc"
    # "pav"
    "giab-tr"
)
query_versions=(
    # "v1"
    # "v4"
    "v4.20"
)
