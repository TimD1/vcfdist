#!/bin/bash

source includes.sh

for truth in "${truth_ids[@]}"; do
    mkdir -p $data/pfda-v2/${truth}_vcfdist
done

# # without standardization, use original representations
# $parallel -j5 --joblog 6_vcfdist.log \
# "$vcfdist  \
#     $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
#     $data/{3}/$ref_id/{5} \
#     $data/refs/$ref_name \
#     --keep-query --keep-truth \
#     -b $data/{3}/$ref_id/{6} \
#     -p $data/pfda-v2/{4}_vcfdist/{1}_HG002_{2}." ::: \
#     ${sub_ids[@]} ::: ${reps[@]} ::: \
#     ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}

# # standardize representations
# $parallel -j5 --joblog 6_vcfdiststd.log \
# "$vcfdist  \
#     $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
#     $data/{3}/$ref_id/{5} \
#     $data/refs/$ref_name \
#     -b $data/{3}/$ref_id/{6} \
#     -p $data/pfda-v2/{4}_vcfdist/{1}_HG002_{2}std." ::: \
#     ${sub_ids[@]} ::: ${reps[@]} ::: \
#     ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}

# standard reps
$parallel -j1 --joblog 6_vcfdist.log \
"$timer -v -o vcfdist_{1}_{4}_{2}.log $vcfdist  \
    $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
    $data/{3}/$ref_id/{5} \
    $data/refs/$ref_name \
    -b $data/{3}/$ref_id/{6} \
    -p $data/pfda-v2/{4}_vcfdist/{1}_HG002_{2}." ::: \
    ${sub_ids[@]} ::: ${reps[@]} ::: \
    ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}

