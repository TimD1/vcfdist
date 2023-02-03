#!/bin/bash

source includes.sh

for ref_id in "${!ref_ids[@]}"; do
    mkdir -p $data/pfda-v2/aligned_bams/${ref_ids[ref_id]}
done

parallel -j4 --joblog 1_align_reads.log \
    "minimap2 \
        -ax {3} \
        --eqx \
        -t 10 \
        $data/refs/{4} \
        $data/pfda-v2/input_fastqs/{1} |
    samtools view \
        -@ 10 -b -h -F 2304 |
    samtools sort \
        -@ 10 > \
        $data/pfda-v2/aligned_bams/{5}/{2}.bam;
    samtools index -@ 10 $data/pfda-v2/aligned_bams/{5}/{2}.bam" ::: \
    ${read_names[@]} :::+ ${read_ids[@]} :::+ ${read_types[@]} ::: \
    ${ref_names[@]} :::+ ${ref_ids[@]}

# parallel -j4 --joblog 1_align_reads.log \
#     "samtools sort \
#         -@ 10 \
#         $data/pfda-v2/aligned_bams/{5}/{2}.bam > \
#         $data/pfda-v2/aligned_bams/{5}/{2}s.bam;
#     samtools index -@10 $data/pfda-v2/aligned_bams/{5}/{2}s.bam" ::: \
#     ${read_names[@]} :::+ ${read_ids[@]} :::+ ${read_types[@]} ::: \
#     ${ref_names[@]} :::+ ${ref_ids[@]}
