#!/bin/bash

time ../src/vcfdist \
    ../test/query.vcf.gz \
    ../test/truth.vcf.gz \
    /x/gm24385/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b ../scripts/data/bench_20.bed \
    --keep-query --keep-truth \
    -r 3_results/sw_

time ../src/vcfdist \
    ../test/query.vcf.gz \
    ../test/truth.vcf.gz \
    /x/gm24385/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b ../scripts/data/bench_20.bed \
    --simple-cluster \
    --keep-query --keep-truth \
    -r 3_results/gap_
