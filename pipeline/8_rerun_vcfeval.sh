#!/bin/bash

source includes.sh

source ~/software/happy/venv2/bin/activate
export HGREF="$data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"

sub_ids=(
    "61YRJ"
    "61YRJ"
    "NWQ6Y"
    "R9CXN"
    "R9CXN"
    "UYMUW"
    "UYMUW"
    "CN36L"
    "13678"
    "60Z59"
    "JBBK0"
    "CZA1Y"
    "FFFGB"
    "HF8CT"
    "HF8CT"
    "MT57N"
    "M9KLP"
    "TZMTX"
    "TZMTX"
    "WX8VK"
)
reps=(
    "O"
    "A"
    "A"
    "B"
    "D"
    "C"
    "D"
    "O"
    "O"
    "A"
    "D"
    "O"
    "A"
    "A"
    "D"
    "D"
    "B"
    "O"
    "A"
    "O"
)

$parallel -j1 --joblog 8_vcfeval_rerun.log \
    "python ~/software/happy/install/bin/hap.py \
        $data/{3}/$ref_id/{5} \
        $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
        -r $data/refs/$ref_name \
        -T $data/{3}/$ref_id/{6} \
        --roc QUAL \
        --write-counts \
        --engine vcfeval \
        -o $data/pfda-v2/{4}_vcfeval/{1}_HG002_{2}" ::: \
    ${sub_ids[@]} :::+ ${reps[@]} ::: \
    ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}

$parallel -j5 \
    "gunzip -f $data/pfda-v2/{3}_vcfeval/{1}_HG002_{2}.*.gz" ::: \
    ${sub_ids[@]} :::+ ${reps[@]} ::: ${truth_ids[@]}

deactivate
