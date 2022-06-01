#!/bin/bash

main_dir="/home/timdunn/vcfdist/scripts/data"
out_dir="/home/timdunn/vcfdist/out"
ont_dir="/home/timdunn/vcfdist/r10.4_chr20/ont-case-study/input/data"
pfda_dir="/home/timdunn/vcfdist/pfda-v2/submission_vcfs"

chr=20

callvcfs=(
    # "$main_dir/truth_orig.vcf.gz"
    # "$main_dir/truth_npore.vcf.gz"
    # "$pfda_dir/2OT9Q/2OT9Q_HG002.vcf.gz"
    "$out_dir/orig_truth.vcf.gz"
    "$out_dir/truth.vcf.gz"
)
callvcfnames=(
    # 'truth_orig'
    # 'truth_npore'
    # '2OT9Q'
    'out_truth_orig'
    'out_truth_min_ed'
)

# callvcfnames=( $(ls $pfda_dir) )
# callvcfs=()
# for callvcf in ${callvcfnames[@]}; do
#     callvcfs+=($pfda_dir/${callvcf}/${callvcf}_HG002.vcf.gz)
# done

truthvcfs=(
    # "$ont_dir/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    "$main_dir/truth_orig.vcf.gz"
    # "$main_dir/truth_npore.vcf.gz"
)

truthvcfnames=(
    # 'truth421'
    'truth_orig'
    # 'truth_npore'
)

fullbeds=(
    "$main_dir/all.bed"
)
bednames=(
    'all'
)

# create list of bed file names for just eval chrs
beds=()
for bed in ${bednames[@]}; do
    beds+=($main_dir/${bed}_${chr}.bed)
done

# create BED region files for just eval chrs
for i in ${!fullbeds[@]}; do
    grep "chr${chr}" ${fullbeds[i]} > ${beds[i]}
done

evalbeds=(
    "$ont_dir/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
)
evalbednames=(
    'eval421'
)

ref_fasta="GRCh38_no_alt.chr20.fa"
source ~/software/happy/venv2/bin/activate
mkdir -p results/

parallel -j25 \
    "python ~/software/happy/install/bin/hap.py \
        {3} \
        {1} \
        -r $ont_dir/$ref_fasta \
        -T {7} \
        -R {5} \
        --roc QUAL \
        --write-counts \
        --engine vcfeval \
        -o results/{2}-{6}-{4}-{8}" ::: \
            ${callvcfs[@]} :::+ ${callvcfnames[@]} ::: \
            ${truthvcfs[@]} :::+ ${truthvcfnames[@]} ::: \
            ${beds[@]} :::+ ${bednames[@]} ::: \
            ${evalbeds[@]} :::+ ${evalbednames[@]}

parallel -j25 \
    "gunzip -f results/{2}-{6}-{4}-{8}.*.gz" ::: \
        ${callvcfs[@]} :::+ ${callvcfnames[@]} ::: \
        ${truthvcfs[@]} :::+ ${truthvcfnames[@]} ::: \
        ${beds[@]} :::+ ${bednames[@]} ::: \
        ${evalbeds[@]} :::+ ${evalbednames[@]}

deactivate
