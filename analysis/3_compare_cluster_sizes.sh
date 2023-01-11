#!/bin/bash

data="/x/vcfdist/data"

time ../src/vcfdist \
    $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz \
    $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz \
    $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    --keep-query --keep-truth \
    -r 3_results/sw_

time ../src/vcfdist \
    $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz \
    $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz \
    $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    --simple-cluster \
    --keep-query --keep-truth \
    -r 3_results/gap_
