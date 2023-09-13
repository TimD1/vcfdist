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
    -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
    -l 10000 \
    -p ../out/giabtr/vcfdist-chr20

# Benchmark on Truvari phab variants (chr20)
# file="../out/giabtr/phab-chr20/output.vcf.gz"
# for sample in `bcftools query -l $file`; do
#     bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
# done
# $timer -v ../src/vcfdist \
#     ../out/giabtr/phab-chr20/norm.QUERY.vcf.gz \
#     ../out/giabtr/phab-chr20/norm.TRUTH.vcf.gz \
#     $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
#     -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#     -l 10000 \
#     -t -q \
#     -p ../out/giabtr/vcfdist-norm-phab-chr20/
