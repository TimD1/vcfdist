#!/bin/bash

data="/home/timdunn/vcfdist/data"
dir="/home/timdunn/vcfdist/analysis-v2/small_sv_all"
parallel="/home/timdunn/parallel-20211222/src/parallel"
vcfdist="/home/timdunn/vcfdist/src/vcfdist"
rtg="/home/timdunn/software/rtg-tools-3.12.1/rtg"
tabix="/home/timdunn/software/htslib-1.16/tabix"
bgzip="/home/timdunn/software/htslib-1.16/bgzip"
timer="/usr/bin/time"

# define reference FASTA
ref_name="${data}/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"

# filename information for each dataset
ds_names=(
    "t2t-q100" # must always be included as first VCF
    "hprc"
    "pav"
    "giab-tr"
)
ds_versions=(
    "v0.9"
    "v1"
    "v4"
    "v4.20"
)
splits=(
    "snp"
    "indel"
    "sv"
    "small"
    "large"
    "all"
)
bednames=(
    "alldifficultregions"
    "AllHomopolymers_ge7bp_imperfectge11bp_slop5"
    "AllTandemRepeats"
    "AllTandemRepeatsandHomopolymers_slop5"
    "MHC"
    "satellites_slop5"
    "segdups"
    "alllowmapandsegdupregions"
)
