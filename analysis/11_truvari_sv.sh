#!/bin/bash

source ~/truvari/venv3.10/bin/activate
source ../pipeline/includes.sh

# bcftools annotate -x FMT/AD GRCh38_HG2-verrkoV1.1-V0.7_dipcall-z2k.vcf.gz > verkko0.vcf.gz
# bcftools norm -m-any --check-ref -w -f GCA... verkko0.vcf.gz -o verkko-norm.vcf.gz
# bcftools sort verkko-norm.vcf.gz -o verkko.vcf.gz
# tabix -p vcf verkko.vcf.gz

# # TruVari without harmonization
# rm -r ../out/giabtr/truvari-init
# truvari bench \
#     -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.vcf.gz \
#     -c $data/zook/verkko.vcf.gz \
#     -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
#     --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.lregions.bed \
#     --no-ref a \
#     --sizemin 1 --sizefilt 1 --sizemax 10000 \
#     --pick multi \
#     --typeignore \
#     --dup-to-ins \
#     -o ../out/giabtr/truvari-init/

# laytr tru2ga \
#     -i ../out/giabtr/truvari-init/ \
    # -o ../out/giabtr/truvari-init/result

# Truvari benchmarking chr20
giab_dir="../out/giabtr"
phab_dir="$giab_dir/phab-chr20"
truvari_phab_dir="$giab_dir/truvari-norm-phab-chr20"

mkdir -p $phab_dir
$timer -v truvari phab \
    -r $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
    -b $data/giab-tr-v4.20/GIABTR.HG002.benchmark.vcf.gz \
    -c $data/zook/verkko.vcf.gz \
    --bSamples HG002\
    --cSamples syndip\
    --align wfa \
    -t 56 \
    -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -o $phab_dir/output.vcf.gz

# file="../out/giabtr/phab-chr20/output.vcf.gz"
# for sample in `bcftools query -l $file`; do
#     bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
#     tabix -p vcf ${file/.vcf*/.$sample.vcf.gz}
# done
samples=(
    "QUERY"
    "TRUTH"
)
# for sample in ${samples[@]}; do
#     bcftools norm \
#         -a \
#         -m-any \
#         --check-ref w \
#         -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
#         $phab_dir/output.${sample}.vcf.gz \
#         -o $phab_dir/norm.${sample}.vcf.gz
#     tabix -p vcf $phab_dir/norm.${sample}.vcf.gz
# done

# $timer -v truvari bench \
#     -b $phab_dir/norm.TRUTH.vcf.gz \
#     -c $phab_dir/norm.QUERY.vcf.gz \
#     -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
#     --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#     --no-ref a \
#     --sizemin 1 --sizefilt 1 --sizemax 10000 \
#     --pick multi \
#     --typeignore \
#     --dup-to-ins \
#     -o $truvari_phab_dir

# laytr tru2ga \
#     -i $truvari_phab_dir/ \
#     -o $truvari_phab_dir/result
# gunzip $truvari_phab_dir/result*.gz

# deactivate
