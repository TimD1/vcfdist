#!/bin/bash

source ~/truvari/venv3.10/bin/activate
source ../pipeline/includes.sh

# bcftools annotate -x FMT/AD GRCh38_HG2-verrkoV1.1-V0.7_dipcall-z2k.vcf.gz > verkko0.vcf.gz
# bcftools norm -m-any --check-ref -w -f GCA... verkko0.vcf.gz -o verkko-norm.vcf.gz
# bcftools sort verkko-norm.vcf.gz -o verkko.vcf.gz
# tabix -p vcf verkko.vcf.gz

# # TruVari without harmonization
# rm -r ../out/giabtr/truvari-initial
# $timer -v truvari bench \
#     -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.vcf.gz \
#     -c $data/zook/verkko.vcf.gz \
#     -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
#     --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
#     --no-ref a \
#     --sizemin 1 --sizefilt 1 --sizemax 10000 \
#     --pick multi \
#     --typeignore \
#     --dup-to-ins \
#     -o ../out/giabtr/truvari-initial/

# updated Truvari benchmarking
rm -r ../out/giabtr/truvari-phab
$timer -v truvari phab \
    -r $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
    -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.vcf.gz \
    -c $data/zook/verkko.vcf.gz \
    --bSample HG002 \
    --cSample syndip \
    -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -o ../out/giabtr/truvari-phab/output.vcf.gz

rm -r ../out/giabtr/truvari-final
$timer -v truvari bench \
    -b ../out/giabtr/truvari-phab/output.vcf.gz \
    -c ../out/giabtr/truvari-phab/output.vcf.gz \
    --bSample HG002 \
    --cSample syndip \
    -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
    --no-ref a \
    --sizemin 1 --sizefilt 1 --sizemax 10000 \
    --pick multi \
    --typeignore \
    --dup-to-ins \
    -o ../out/giabtr/truvari-final/

deactivate
