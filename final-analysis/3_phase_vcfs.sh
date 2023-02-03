#!/bin/bash

source includes.sh

source ~/miniconda3/etc/profile.d/conda.sh
conda activate whatshap-env

mkdir -p $data/pfda-v2/phased_vcfs/$ref_id

phased_bams=""
for rd_id in "${read_ids[@]}"; do
    phased_bams="$phased_bams $data/pfda-v2/phased_bams/${ref_id}/${rd_id}.bam"
done
phased_bams="${phased_bams:1}"

$parallel -j5 --joblog 3_phase_vcfs.log \
    "whatshap phase \
        ${data}/pfda-v2/submission_vcfs/{1}/{1}_HG002.vcf.gz \
        ${phased_bams} \
        --output ${data}/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_O.vcf.gz \
        --reference ${data}/refs/$ref_name \
        --ignore-read-groups \
        --indels; \
    tabix -f -p vcf ${data}/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_O.vcf.gz" ::: ${sub_ids[@]}

conda deactivate
