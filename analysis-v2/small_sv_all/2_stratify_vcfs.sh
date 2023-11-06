#!/bin/bash

source globals.sh

$parallel \
    "bgzip ./vcfdist/{1}.{2}.summary.vcf" ::: \
    ${ds_names[@]} ::: ${splits[@]}

$parallel \
    "tabix -p vcf ./vcfdist/{1}.{2}.summary.vcf.gz" ::: \
    ${ds_names[@]} ::: ${splits[@]}

$parallel \
    "bcftools view \
        -R $data/genome-stratifications-v3.3/GRCh38_{3}.bed.gz \
        -Ou \
        ./vcfdist/{1}.{2}.summary.vcf.gz \
        -o ./vcfdist/{1}.{2}.{3}.vcf" ::: \
        ${ds_names[@]} ::: ${splits[@]} ::: ${bednames[@]}

$parallel \
    "gunzip ./vcfdist/{1}.{2}.summary.vcf.gz" ::: \
    ${ds_names[@]} ::: ${splits[@]}
