#!/bin/bash

calls='../test/calls.vcf.gz'
truth='../test/truth.vcf.gz'
ref='/x/gm24385/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
bed='../scripts/data/bench_20.bed'
names=(
    "A"
    "B"
    "C"
    "D"
)

# VCFdist compare
parallel -j8 --joblog "vcfdist-compare.log" \
"../src/vcfdist \
    1_results/{1}-calls.vcf \
    1_results/C-truth.vcf \
    $ref \
    -b $bed \
    --keep-query --keep-truth \
    --simple-cluster \
    -q {2} \
    -r 4_results/{1}-calls_C-truth_{2}_" ::: ${names[@]} ::: $(seq 0 50)
