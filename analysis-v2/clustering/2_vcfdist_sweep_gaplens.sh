#!/bin/bash

source globals.sh

dists=( 
        # "10"
        # "50"
        # "100"
        "500"
)

# vcfdist evaluation
mkdir -p $out/sweep_gaplens
for q in "${!query_names[@]}"; do
    for d in "${dists[@]}"; do
        echo "vcfdist: evaluating len 1000 '${query_names[q]}' gap $d"
        $timer -v ../../src/vcfdist \
            $data/${query_names[q]}-${query_versions[q]}/split/${query_names[q]}.all.vcf.gz \
            $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
            $data/refs/$ref_name \
            --bed $data/${truth_name}-${truth_version}/split/bench.bed \
            --simple-cluster \
            -g $d \
            -l 1000 \
            -p $out/sweep_gaplens/${query_names[q]}.gap${d}.len1000. \
            2> $out/sweep_gaplens/${query_names[q]}.gap${d}.len1000.log
    done
done
