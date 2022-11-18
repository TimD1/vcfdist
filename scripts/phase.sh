#!/bin/bash

pfda_dir="/home/timdunn/vcfdist/pfda-v2/submission_vcfs"
sub_ids=( $(ls $pfda_dir) )

mkdir -p phase
for id in ${sub_ids[@]}; do
    cp $pfda_dir/$id/${id}_HG002.vcf.gz ./phase
    gunzip ./phase/${id}_HG002.vcf.gz
done
