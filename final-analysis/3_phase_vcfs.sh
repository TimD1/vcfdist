#!/bin/bash

source includes.sh

source ~/miniconda3/etc/profile.d/conda.sh
conda activate whatshap-env

mkdir -p $data/pfda-v2/phased_vcfs
for ref_id in "${!ref_ids[@]}"; do
    mkdir -p $data/pfda-v2/phased_vcfs/${ref_ids[ref_id]}
done

# TODO: fix so it'll run with multiple references
phased_bams=""
for ref_id in "${ref_ids[@]}"; do
    for rd_id in "${read_ids[@]}"; do
        phased_bams="$phased_bams $data/pfda-v2/phased_bams/${ref_id}/${rd_id}.bam"
    done
done
phased_bams="${phased_bams:1}"

$parallel -j5 --joblog 3_phase_vcfs.log \
    "whatshap phase \
        ${data}/pfda-v2/submission_vcfs/{3}/{3}_HG002.vcf.gz \
        ${phased_bams} \
        --output ${data}/pfda-v2/phased_vcfs/{2}/{3}_HG002_O.vcf.gz \
        --reference ${data}/refs/{1} \
        --ignore-read-groups \
        --indels; \
    tabix -f -p vcf ${data}/pfda-v2/phased_vcfs/{2}/{3}_HG002_O.vcf.gz" ::: \
    ${ref_names[@]} :::+ ${ref_ids[@]} ::: ${sub_ids[@]}

conda deactivate
