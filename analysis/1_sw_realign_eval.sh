#!/bin/bash

calls='../test/calls.vcf.gz'
truth='../test/truth.vcf.gz'
ref='/x/gm24385/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
bed='../scripts/data/bench_20.bed'
sub_penalties=(
    # "10"
    # "6"
    # "10"
    # "40"
    # "20"
    "10"
    "10"
    "10"
    "10"
)
open_penalties=(
    # "13"
    # "5"
    # "13"
    # "79"
    # "33"
    "1"
    "5"
    "15"
    "2"
)
extend_penalties=(
    # "3"
    # "3"
    # "5"
    # "7"
    # "5"
    "2"
    "2"
    "2"
    "8"
)
names=(
    # "mm2-short"
    # "mm2-ont"
    # "mm2-hifi"
    # "mm2-asm5"
    # "mm2-asm10"
    "A"
    "B"
    "C"
    "D"
)

# # VCFdist realign
# parallel -j4 --joblog "vcfdist-realign.log" \
# "../src/vcfdist \
#     $calls $truth $ref -b $bed \
#     -r 1_results/{4}- \
#     -s {1} -o {2} -e {3} -x" ::: \
#     ${sub_penalties[@]} :::+ ${open_penalties[@]} :::+ \
#     ${extend_penalties[@]} :::+ ${names[@]}

# VCFdist compare
parallel -j8 --joblog "vcfdist-compare.log" \
"../src/vcfdist \
    1_results/{1}-calls.vcf \
    1_results/{2}-truth.vcf \
    $ref \
    -b $bed \
    --keep-query --keep-truth \
    -r 1_results/{1}-calls_{2}-truth_" ::: ${names[@]} ::: ${names[@]}

# # hap.py evaluation
# source ~/software/happy/venv2/bin/activate
# parallel -j25 \
#     "python ~/software/happy/install/bin/hap.py \
#     1_results/{1}-truth.vcf \
#     1_results/{2}-calls.vcf \
#     -r $ref \
#     --engine-vcfeval-template /x/gm24385/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.sdf \
#     -T $bed \
#     --roc QUAL \
#     --threads 1 \
#     --write-counts \
#     --engine vcfeval \
#     -o 1_results/{2}-calls_{1}-truth" ::: \
#     ${names[@]} ::: ${names[@]}

# # unzip results
# parallel -j25 \
#     "gunzip -f 1_results/{1}-calls_{2}-truth.*.gz" ::: \
#     ${names[@]} ::: ${names[@]}
