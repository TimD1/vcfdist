#!/bin/bash

source includes.sh

# realign using A,B,C,D
$parallel -j10 --joblog "4_vcf_reps.log" \
"$vcfdist  \
    $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_O.vcf.gz \
    $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz \
    $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -b $data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    -x {3} -o {4} -e {5} \
    -v 1 \
    --realign-only \
    -p $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}." ::: \
    ${sub_ids[@]} ::: ${pts[@]} :::+ ${x[@]} :::+ ${o[@]} :::+ ${e[@]}

# remove other outputs
$parallel -j10 \
    "rm $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.truth.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j10 \
    "rm $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.orig_truth.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j10 \
    "rm $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.orig_query.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}

# rename vcfdist output, compress, and index
$parallel -j10 \
    "mv $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.query.vcf \
        $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j10 \
    "bgzip $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.vcf" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
$parallel -j10 \
    "tabix -p vcf $data/pfda-v2/phased_vcfs/GRCh38/{1}_HG002_{2}.vcf.gz" ::: \
        ${sub_ids[@]} ::: ${pts[@]}
