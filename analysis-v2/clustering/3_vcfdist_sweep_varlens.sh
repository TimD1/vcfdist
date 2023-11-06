#!/bin/bash

source globals.sh

lens=( 
        # "10"
        # "50"
        # "100"
        # "500"
        "1000"
        # "5000"
        # "10000"
)

mkdir -p $out/sweep_varlens

# # wfa-based superclustering
# for q in "${!query_names[@]}"; do
#     for l in "${lens[@]}"; do
#         echo "vcfdist: evaluating wfa '${query_names[q]}' len $l"
#         $timer -v ../../src/vcfdist \
#             $data/${query_names[q]}-${query_versions[q]}/split/${query_names[q]}.all.vcf.gz \
#             $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#             $data/refs/$ref_name \
#             --bed $data/${truth_name}-${truth_version}/split/bench.bed \
#             -l $l \
#             -p $out/sweep_varlens/${query_names[q]}.wfa.len${l}. \
#             2> $out/sweep_varlens/${query_names[q]}.wfa.len${l}.log
#     done
# done

# wfa-based superclustering
for q in "${!query_names[@]}"; do
    for l in "${lens[@]}"; do
        echo "vcfdist: evaluating wfa '${query_names[q]}' len $l"
        $timer -v ../../src/vcfdist \
            $data/${query_names[q]}-${query_versions[q]}/split/${query_names[q]}.all.vcf.gz \
            $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
            $data/refs/$ref_name \
            --bed $data/${truth_name}-${truth_version}/split/bench.bed \
            -l $l \
            -p $out/sweep_varlens/${query_names[q]}.wfa.len${l}ra. \
            2> $out/sweep_varlens/${query_names[q]}.wfa.len${l}ra.log
    done
done

# # gap-based superclustering
# for q in "${!query_names[@]}"; do
#     for l in "${lens[@]}"; do
#         echo "vcfdist: evaluating gap 100 '${query_names[q]}' len $l"
#         $timer -v ../../src/vcfdist \
#             $data/${query_names[q]}-${query_versions[q]}/split/${query_names[q]}.all.vcf.gz \
#             $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#             $data/refs/$ref_name \
#             --bed $data/${truth_name}-${truth_version}/split/bench.bed \
#             --simple-cluster \
#             -g 100 \
#             -l $l \
#             -p $out/sweep_varlens/${query_names[q]}.gap100.len${l}. \
#             2> $out/sweep_varlens/${query_names[q]}.gap100.len${l}.log
#     done
# done
