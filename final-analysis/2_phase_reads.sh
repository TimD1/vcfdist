#!/bin/bash

source includes.sh

source ~/miniconda3/etc/profile.d/conda.sh
conda activate whatshap-env

for ref_id in "${!ref_ids[@]}"; do
    mkdir -p $data/pfda-v2/phased_bams/${ref_ids[ref_id]}
done

parallel -j4 --joblog 2_phase_reads.log \
    "whatshap haplotag \
        --output $data/pfda-v2/phased_bams/{5}/{2}.bam \
        --reference $data/refs/{4} \
        --ignore-read-groups \
        --skip-missing-contigs \
        --output-threads=10 \
        $data/nist-v4.2.1/{5}/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz \
        $data/pfda-v2/aligned_bams/{5}/{2}.bam;
    samtools index -@ 12 $data/pfda-v2/phased_bams/{5}/{2}.bam" ::: \
    ${read_names[@]} :::+ ${read_ids[@]} :::+ ${read_types[@]} ::: \
    ${ref_names[@]} :::+ ${ref_ids[@]}

conda deactivate
