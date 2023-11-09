#!/bin/bash

source globals.sh

# # vcfdist evaluation
# mkdir -p $out/evals/vcfdist_tr
# for i in "${!query_names[@]}"; do
#     echo "vcfdist: evaluating '${query_names[i]}'"
#     $timer -v ../../src/vcfdist \
#         $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         $data/refs/$ref_name \
#         --bed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
#         -l 1000 \
#         -p $out/evals/vcfdist_tr/${query_names[i]}. \
#         2> $out/evals/vcfdist_tr/${query_names[i]}.log
# done

# # vcfeval evaluation
# mkdir -p $out/evals/vcfeval_tr
# for i in "${!query_names[@]}"; do
#     echo "vcfeval: evaluating '${query_names[i]}'"
#     $timer -v $rtg vcfeval \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         -t $data/refs/${ref_name}.sdf \
#         -e $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
#         --threads 64 \
#         --all-records \
#         -o $out/evals/vcfeval_tr/${query_names[i]} \
#         2> $out/evals/vcfeval_tr/${query_names[i]}.log
# done

source ~/truvari/venv3.10/bin/activate

# # truvari evaluation
# mkdir -p $out/evals/truvari_tr
# for i in "${!query_names[@]}"; do
#     echo "truvari: evaluating '${query_names[i]}'"
#     truvari bench \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         --bSample HG002 \
#         --cSample HG002 \
#         --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
#         -f $data/refs/$ref_name \
#         -s 1 -S 1 \
#         --sizemax 1000 \
#         --no-ref a \
#         --pick multi \
#         -o $out/evals/truvari_tr/${query_names[i]} \
#         2> $out/evals/truvari_tr/${query_names[i]}.log
# done

# # wfa refine
# for i in "${!query_names[@]}"; do
#     echo "truvari: wfa refine '${query_names[i]}'"
#     rm -rf $out/evals/truvari_tr/${query_names[i]}_wfa
#     cp -r $out/evals/truvari_tr/${query_names[i]} \
#           $out/evals/truvari_tr/${query_names[i]}_wfa
#     truvari refine \
#         -f $data/refs/$ref_name \
#         -t 64 \
#         --use-original-vcfs \
#         --align wfa \
#         $out/evals/truvari_tr/${query_names[i]}_wfa \
#         2> $out/evals/truvari_tr/${query_names[i]}_wfa.log
# done

# mafft refine
for i in "${!query_names[@]}"; do
    echo "truvari: mafft refine '${query_names[i]}'"
    rm -rf $out/evals/truvari_tr/${query_names[i]}_mafft
    cp -r $out/evals/truvari_tr/${query_names[i]} \
          $out/evals/truvari_tr/${query_names[i]}_mafft
    truvari refine \
        -f $data/refs/$ref_name \
        -t 64 \
        --use-original-vcfs \
        --align mafft \
        $out/evals/truvari_tr/${query_names[i]}_mafft \
        2> $out/evals/truvari_tr/${query_names[i]}_mafft.log
done

deactivate
