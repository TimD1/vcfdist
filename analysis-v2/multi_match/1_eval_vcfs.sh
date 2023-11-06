#!/bin/bash

source globals.sh

        # --bed $data/cmrg-v1.00/HG002_GRCh38_CMRG_smallvar_v1.00.bed \
        # --max-threads 1 \
# vcfdist evaluation
mkdir -p $out/evals/vcfdist
for i in "${!query_names[@]}"; do
    echo "vcfdist: evaluating '${query_names[i]}'"
    $timer -v ../../src/vcfdist \
        $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
        $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
        $data/refs/$ref_name \
        --bed $data/${truth_name}-${truth_version}/split/bench.bed \
        --keep-query --keep-truth \
        -l 1000 \
        -p $out/evals/vcfdist/${query_names[i]}. \
        2> $out/evals/vcfdist/${query_names[i]}.log
done
