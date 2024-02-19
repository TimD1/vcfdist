#!/bin/bash

source globals.sh

# vcfeval evaluation
mkdir -p vcfeval
for i in "${!query_names[@]}"; do
    echo "vcfeval: evaluating '${query_names[i]}'"
    rm -r vcfeval/${query_names[i]}
    $timer -v $rtg vcfeval \
        -b ./vcfs/${truth_name}.vcf.gz \
        -c ./vcfs/${query_names[i]}.vcf.gz \
        -t $data/refs/${ref_name}.sdf \
        --threads 64 \
        --ref-overlap \
        -m ga4gh \
        --vcf-score-field=QUAL \
        -o vcfeval/${query_names[i]} \
        2> vcfeval/${query_names[i]}.log
    gunzip vcfeval/${query_names[i]}/*.vcf.gz
done
