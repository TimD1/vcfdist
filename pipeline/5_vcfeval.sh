#!/bin/bash

source includes.sh

source ~/software/happy/venv2/bin/activate
export HGREF="$data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"

for truth in "${truth_ids[@]}"; do
    mkdir -p $data/pfda-v2/${truth}_vcfeval
done

# # rtgtools vcfeval timing evaluation
# $parallel -j1 --joblog 5_vcfeval.log \
#     "$timer -v -o 5_{4}_{2}.txt ~/software/rtg-tools-3.12.1/rtg RTG_MEM=8G vcfeval \
#         -b $data/{3}/$ref_id/{5} \
#         -c $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
#         -t $data/refs/$ref_name.sdf \
#         -e $data/{3}/$ref_id/{6} \
#         --all-records \
#         -o $data/pfda-v2/{4}_vcfeval/{1}_HG002_{2}" ::: \
#     ${sub_ids[@]} ::: ${reps[@]} ::: \
#     ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}

# hap.py vcfeval timing evaluation
# "$timer -v -o 5_{4}_{2}.txt python ~/software/happy/install/bin/hap.py \
$parallel -j5 --joblog 5_vcfeval.log \
    "$timer -v -o 5_{4}_{2}.txt taskset 1 python ~/software/happy/install/bin/hap.py \
        $data/{3}/$ref_id/{5} \
        $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
        -r $data/refs/$ref_name \
        -T $data/{3}/$ref_id/{6} \
        --threads 1 \
        --roc QUAL \
        --write-counts \
        --engine vcfeval \
        -o $data/pfda-v2/{4}_vcfeval/{1}_HG002_{2}" ::: \
    ${sub_ids[@]} ::: ${reps[@]} ::: \
    ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}

# # standardized representation (vcfdist output)
# $parallel -j2 --joblog 5_vcfevalstd.log \
#     "python ~/software/happy/install/bin/hap.py \
#         $data/{3}/$ref_id/{5} \
#         $data/pfda-v2/{4}_vcfdist/{1}_HG002_{2}std.query.vcf \
#         -r $data/refs/$ref_name \
#         -T $data/{3}/$ref_id/{6} \
#         --threads 5 \
#         --roc QUAL \
#         --write-counts \
#         --engine vcfeval \
#         -o $data/pfda-v2/{4}_vcfeval/{1}_HG002_{2}std" ::: \
#     ${sub_ids[@]} ::: ${reps[@]} ::: \
#     ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}
# $parallel -j5 \
#     "gunzip -f $data/pfda-v2/{3}_vcfeval/{1}_HG002_{2}std.*.gz" ::: \
#     ${sub_ids[@]} ::: ${reps[@]} ::: ${truth_ids[@]}

# $parallel -j5 --joblog 5_vcfeval.log \
#     "python ~/software/happy/install/bin/hap.py \
#         $data/{3}/$ref_id/{5} \
#         $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
#         -r $data/refs/$ref_name \
#         -T $data/{3}/$ref_id/{6} \
#         --roc QUAL \
#         --write-counts \
#         --engine vcfeval \
#         -o $data/pfda-v2/{4}_vcfeval/{1}_HG002_{2}" ::: \
#     ${sub_ids[@]} ::: ${reps[@]} ::: \
#     ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}
# $parallel -j5 \
#     "gunzip -f $data/pfda-v2/{3}_vcfeval/{1}_HG002_{2}.*.gz" ::: \
#     ${sub_ids[@]} ::: ${reps[@]} ::: ${truth_ids[@]}

deactivate
