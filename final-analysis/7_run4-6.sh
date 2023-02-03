#!/bin/bash

echo "./4_vcf_reps.sh"
source includes.sh

# realign using selected design pts
$parallel -j5 --joblog "4_vcf_reps.log" \
"$vcfdist  \
    $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_O.vcf.gz \
    $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz \
    $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b all.bed \
    -x {3} -o {4} -e {5} \
    --realign-only \
    -p $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}." ::: \
    ${sub_ids[@]} ::: ${pts[@]} :::+ ${x[@]} :::+ ${o[@]} :::+ ${e[@]}

# remove other outputs
$parallel -j5 \
    "rm $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.truth.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j5 \
    "rm $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.orig_truth.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j5 \
    "rm $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.orig_query.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}

# rename vcfdist output, compress, and index
$parallel -j5 \
    "mv $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.query.vcf \
        $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j5 \
    "$bgzip -f $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j5 \
    "$tabix -f -p vcf $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.vcf.gz" ::: \
        ${sub_ids[@]} ::: ${pts[@]}

echo "./5_vcfeval.sh"
source ~/software/happy/venv2/bin/activate
export HGREF="$data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"

for truth in "${truth_ids[@]}"; do
    mkdir -p $data/pfda-v2/${truth}_vcfeval
done

$parallel -j5 --joblog 5_vcfeval.log \
    "python ~/software/happy/install/bin/hap.py \
        $data/{3}/$ref_id/{5} \
        $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
        -r $data/refs/$ref_name \
        -T $data/{3}/$ref_id/{6} \
        --roc QUAL \
        --write-counts \
        --engine vcfeval \
        -o $data/pfda-v2/{4}_vcfeval/{1}_HG002_{2}" ::: \
    ${sub_ids[@]} ::: ${reps[@]} ::: \
    ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}

$parallel -j5 \
    "gunzip -f $data/pfda-v2/{3}_vcfeval/{1}_HG002_{2}.*.gz" ::: \
    ${sub_ids[@]} ::: ${reps[@]} ::: ${truth_ids[@]}

deactivate

echo "./6_vcfdist.sh"
for truth in "${truth_ids[@]}"; do
    mkdir -p $data/pfda-v2/${truth}_vcfdist
done

# standard reps
$parallel -j5 --joblog 6_vcfdist.log \
"$vcfdist  \
    $data/pfda-v2/phased_vcfs/$ref_id/{1}_HG002_{2}.vcf.gz \
    $data/{3}/$ref_id/{5} \
    $data/refs/$ref_name \
    -b $data/{3}/$ref_id/{6} \
    -p $data/pfda-v2/{4}_vcfdist/{1}_HG002_{2}." ::: \
    ${sub_ids[@]} ::: ${reps[@]} ::: \
    ${truth_names[@]} :::+ ${truth_ids[@]} :::+ ${truth_vcfs[@]} :::+ ${truth_beds[@]}
