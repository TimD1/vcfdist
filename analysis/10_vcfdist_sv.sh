#!/bin/bash

source ../pipeline/includes.sh

# # Zook GIAB TR dataset (0:56, 3:40, 7:30)
# $timer -v ../src/vcfdist \
#     $data/zook/verkko.vcf.gz \
#     $data/giab-tr-v4.20/GIABTR.HG002.benchmark.vcf.gz \
#     $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
#     -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
#     -l 10000 \
#     -t -q \
#     -p ../out/giabtr/vcfdist-keep/

# Zook GIAB TR dataset
$timer -v ../src/vcfdist \
    $data/zook/verkko.vcf.gz \
    $data/giab-tr-v4.20/GIABTR.HG002.benchmark.vcf.gz \
    $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.small.bed \
    -l 1000 \
    -p ../out/small
