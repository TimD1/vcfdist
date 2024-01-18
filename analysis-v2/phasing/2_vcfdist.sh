#!/bin/bash

source globals.sh

mkdir -p vcfdist
for q in "${!query_names[@]}"; do
    echo "vcfdist: evaluating '${query_names[q]}'"
    $timer -v ../../src/vcfdist \
        ./vcfs/${query_names[q]}.vcf \
        ./vcfs/${truth_name}.vcf \
        $data/refs/$ref_name \
        -l 1000 \
        -p ./vcfdist/${query_names[q]}. \
        2> ./vcfdist/${query_names[q]}.log
done
