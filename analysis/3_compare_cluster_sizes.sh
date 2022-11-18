#!/bin/bash

time ../src/vcfdist \
    ../test/calls.vcf.gz \
    ../test/truth.vcf.gz \
    /x/gm24385/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b ../scripts/data/bench_20.bed \
    --keep-calls --keep-truth \
    -r 3_results/sw_

time ../src/vcfdist \
    ../test/calls.vcf.gz \
    ../test/truth.vcf.gz \
    /x/gm24385/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b ../scripts/data/bench_20.bed \
    --simple-cluster \
    --keep-calls --keep-truth \
    -r 3_results/gap_
