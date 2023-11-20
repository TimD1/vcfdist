#!/bin/bash

source globals.sh

# create VCFs with only benchmarking regions
echo "filtering $truth_name"
bcftools filter \
    $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
    -R $data/${truth_name}-${truth_version}/split/bench.bed \
    -o ./vcfs/${truth_name}.vcf

for i in ${!query_names[@]}; do
    echo "filtering ${query_names[i]}"
    bcftools filter \
        $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
        -R $data/${truth_name}-${truth_version}/split/bench.bed \
        -o ./vcfs/${query_names[i]}.vcf
done
