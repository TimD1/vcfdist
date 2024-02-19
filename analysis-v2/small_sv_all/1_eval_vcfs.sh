#!/bin/bash

source globals.sh

# vcfdist evaluation
mkdir -p vcfdist
for cat in "${splits[@]}"; do
    for i in "${!ds_names[@]}"; do
        echo "vcfdist: evaluating '${ds_names[i]}' ${cat}"
        $timer -v $vcfdist \
            $data/${ds_names[i]}-${ds_versions[i]}/split/${ds_names[i]}.${cat}.vcf.gz \
            $data/t2t-q100-v0.9/split/t2t-q100.${cat}.vcf.gz \
            $ref_name \
            --bed $data/t2t-q100-v0.9/GRCh38_HG2-T2TQ100-V0.9_dipcall-z2k.benchmark.bed \
            -l 1000 \
            -p $dir/vcfdist/${ds_names[i]}.${cat}. \
            2> $dir/vcfdist/${ds_names[i]}.${cat}.log
    done
done
