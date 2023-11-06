#!/bin/bash

source globals.sh

# # vcfdist evaluation
# mkdir -p $out/evals/vcfdist
# for i in "${!query_names[@]}"; do
#     echo "vcfdist: evaluating '${query_names[i]}'"
#     $timer -v ../../src/vcfdist \
#         $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         $data/refs/$ref_name \
#         --bed $data/${truth_name}-${truth_version}/split/bench.bed \
#         --keep-query --keep-truth \
#         -l 1000 \
#         -p $out/evals/vcfdist/${query_names[i]}. \
#         2> $out/evals/vcfdist/${query_names[i]}.log
# done

# # vcfeval evaluation
# mkdir -p $out/evals/vcfeval
# for i in "${!query_names[@]}"; do
#     echo "vcfeval: evaluating '${query_names[i]}'"
#     $timer -v $rtg vcfeval \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         -t $data/refs/${ref_name}.sdf \
#         -e $data/${truth_name}-${truth_version}/split/bench.bed \
#         --threads 64 \
#         --all-records \
#         -o $out/evals/vcfeval/${query_names[i]} \
#         2> $out/evals/vcfeval/${query_names[i]}.log
# done

source ~/truvari/venv3.10/bin/activate

# # truvari evaluation
# mkdir -p $out/evals/truvari
# for i in "${!query_names[@]}"; do
#     echo "truvari: evaluating '${query_names[i]}'"
#     truvari bench \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         --bSample HG002 \
#         --cSample HG002 \
#         --includebed $data/${truth_name}-${truth_version}/split/bench.bed \
#         -f $data/refs/$ref_name \
#         -s 1 -S 1 \
#         --sizemax 1000 \
#         --no-ref a \
#         --pick multi \
#         -o $out/evals/truvari/${query_names[i]} \
#         2> $out/evals/truvari/${query_names[i]}.log
# done

# wfa refine
for i in "${!query_names[@]}"; do
    if [[ $i -eq  1 ]]; then
        continue
    fi
    echo "truvari: wfa refine '${query_names[i]}'"
    rm -rf $out/evals/truvari/${query_names[i]}_wfa
    cp -r $out/evals/truvari/${query_names[i]}_\
          $out/evals/truvari/${query_names[i]}_wfa
    truvari refine \
        -f $data/refs/$ref_name \
        -t 64 \
        --use-original-vcfs \
        --align wfa \
        $out/evals/truvari/${query_names[i]}_wfa \
        2> $out/evals/truvari/${query_names[i]}_wfa.log
done

# for rep in "${reps[@]}"; do
#     for i in "${!ds_names[@]}"; do
#         echo "truvari: mafft refine $rep '${ds_names[i]}'"
#         rm -rf $out/evals/$rep/truvari/${ds_names[i]}_mafft
#         cp -r $out/evals/$rep/truvari/${ds_names[i]} \
#               $out/evals/$rep/truvari/${ds_names[i]}_mafft
#         truvari refine \
#             -f $data/refs/$ref_name \
#             -t 64 \
#             --use-original-vcfs \
#             --align mafft \
#             $out/evals/$rep/truvari/${ds_names[i]}_mafft
#     done
# done

deactivate
