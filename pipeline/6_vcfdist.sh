#!/bin/bash

source includes.sh

for truth in "${truth_ids[@]}"; do
    mkdir -p $data/pfda-v2/${truth}_vcfdist
done

$parallel -j1 --joblog 6_vcfdist.log \
"$timer -v -o vcfdist_{1}_{4}_{2}.log $vcfdist  \
    $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
    $data/{3}/{5} \
    $data/refs/$ref_name \
    -b $data/{3}/{6} \
    -l 500 \
    -p $data/pfda-v2/{4}_vcfdist/{1}_HG002_{2}." ::: \
    ${sub_ids[@]} ::: ${reps[@]} ::: \
    ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}
    # -t -q \
