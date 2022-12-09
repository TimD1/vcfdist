#!/bin/bash

query='../test/query.vcf.gz'
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
    "5"
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
    "6"
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
mkdir -p 5_results

# parallel -j8 --joblog "vcfdist-realign.log" \
# "../src/vcfdist \
#     $query \
#     $truth \
#     $ref \
#     -b $bed \
#     -s {2} -o {3} -e {4} \
#     --exit-after-realign \
#     -r 5_results/{1}_ > 5_results/{1}.out" ::: \
#         ${names[@]} :::+ ${sub_penalties[@]} :::+ \
#         ${open_penalties[@]} :::+ ${extend_penalties[@]} 

# VCFdist compare
parallel -j8 --joblog "vcfdist-compare.log" \
"../src/vcfdist \
    5_results/{1}_query.vcf \
    5_results/{2}_truth.vcf \
    $ref \
    -b $bed \
    -p 1 \
    --keep-query --keep-truth \
    --new-prec-calc \
    -r 5_results/{1}-query_{2}-truth_ > 5_results/{1}-query_{2}-truth.out" ::: \
        ${names[@]} ::: ${names[@]}

# # hap.py evaluation
# source ~/software/happy/venv2/bin/activate
# parallel -j25 \
#     "python ~/software/happy/install/bin/hap.py \
#     5_results/{1}_truth.vcf \
#     5_results/{2}_query.vcf \
#     -r $ref \
#     --engine-vcfeval-template /x/gm24385/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.sdf \
#     -T $bed \
#     --roc QUAL \
#     --threads 1 \
#     --write-counts \
#     --engine vcfeval \
#     -o 5_results/{2}-query_{1}-truth" ::: \
#     ${names[@]} ::: ${names[@]}

# # unzip results
# parallel -j25 \
#     "gunzip -f 5_results/{1}-query_{2}-truth.*.gz" ::: \
#     ${names[@]} ::: ${names[@]}
