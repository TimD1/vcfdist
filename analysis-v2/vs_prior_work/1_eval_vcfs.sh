#!/bin/bash

source globals.sh

# EVALUATE WHOLE-GENOME BENCHMARKING REGIONS

# # vcfdist evaluation
# mkdir -p ./vcfdist
# for i in "${!query_names[@]}"; do
#     echo "vcfdist: evaluating '${query_names[i]}'"
#     $timer -v ../../src/vcfdist \
#         $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
#         $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
#         $data/refs/$ref_name \
#         --bed $data/${truth_name}-${truth_version}/split/bench.bed \
#         -t 64 -r 128 \
#         -l 1000 \
#         -p ./vcfdist/${query_names[i]}. \
#         2> ./vcfdist/${query_names[i]}.log
# done

# # vcfeval evaluation
# mkdir -p ./vcfeval
# for i in "${!query_names[@]}"; do
#     echo "vcfeval: evaluating '${query_names[i]}'"
#     rm -r ./vcfeval/${query_names[i]}
#     $timer -v $rtg vcfeval \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
#         -t $data/refs/${ref_name}.sdf \
#         -e $data/${truth_name}-${truth_version}/split/bench.bed \
#         --threads 64 \
#         --all-records \
#         --ref-overlap \
#         -o ./vcfeval/${query_names[i]} \
#         2> ./vcfeval/${query_names[i]}.log
#     gunzip ./vcfeval/${query_names[i]}/*.vcf.gz
# done

source ~/software/Truvari-4.2.1/venv3.10/bin/activate

# # truvari evaluation
# mkdir -p ./truvari
# for i in "${!query_names[@]}"; do
#     echo "truvari: evaluating '${query_names[i]}'"
#     rm -rf ./truvari/${query_names[i]}
#     $timer -v truvari bench \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
#         --bSample HG002 \
#         --cSample HG002 \
#         --includebed $data/${truth_name}-${truth_version}/split/bench.bed \
#         --no-ref a \
#         --sizemin 1 --sizefilt 1 --sizemax 1000 \
#         --pick single \
#         --typeignore \
#         --dup-to-ins \
#         -o ./truvari/${query_names[i]} \
#         2> ./truvari/${query_names[i]}.log
#     truvari ga4gh \
#         -i ./truvari/${query_names[i]}/ \
#         -o ./truvari/${query_names[i]}/result
#     gunzip ./truvari/${query_names[i]}/result*.gz
# done

# wfa/mafft refine
for aln in "${truvari_refine_options[@]}"; do
    for i in "${!query_names[@]}"; do
        echo "truvari: ${aln} refine '${query_names[i]}'"
        mkdir -p ./truvari-${aln}
        rm -rf ./truvari-${aln}/${query_names[i]}
        cp -r ./truvari/${query_names[i]} \
              ./truvari-${aln}/${query_names[i]}

        # refine difficult regions
        $timer -v truvari refine \
            -f $data/refs/$ref_name \
            -t 64 \
            --regions ./truvari-${aln}/${query_names[i]}/candidate.refine.bed \
            --use-original-vcfs \
            --use-region-coords \
            --recount \
            --align ${aln} \
            ./truvari-${aln}/${query_names[i]} \
            2> ./truvari-${aln}/${query_names[i]}.log
        truvari ga4gh \
            -i ./truvari-${aln}/${query_names[i]}/phab_bench/ \
            -o ./truvari-${aln}/${query_names[i]}/phab_bench/result
        gunzip ./truvari-${aln}/${query_names[i]}/phab_bench/result*.gz
        rm -f ./truvari-${aln}/${query_names[i]}/*.gz
        rm -f ./truvari-${aln}/${query_names[i]}/*.vcf
        rm -f ./truvari-${aln}/${query_names[i]}/*.tbi

        # normal bench in unrefined regions
        echo "bedtools: exclude phab '${query_names[i]}'"
        bedtools subtract \
            -a $data/${truth_name}-${truth_version}/split/bench.bed \
            -b ./truvari-${aln}/${query_names[i]}/candidate.refine.bed \
            > ./truvari-${aln}/${query_names[i]}/no-refine.bed
        echo "truvari: bench excluded phab '${query_names[i]}'"
        rm -rf ./truvari-${aln}/${query_names[i]}/no_phab_bench
        $timer -v truvari bench \
            -b $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
            -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
            --bSample HG002 \
            --cSample HG002 \
            --includebed ./truvari-${aln}/${query_names[i]}/no-refine.bed \
            --no-ref a \
            --sizemin 1 --sizefilt 1 --sizemax 1000 \
            --pick single \
            --typeignore \
            --dup-to-ins \
            -o ./truvari-${aln}/${query_names[i]}/no_phab_bench \
            2> ./truvari-${aln}/${query_names[i]}/no_phab_bench.log
        truvari ga4gh \
            -i ./truvari-${aln}/${query_names[i]}/no_phab_bench/ \
            -o ./truvari-${aln}/${query_names[i]}/no_phab_bench/result
        gunzip ./truvari-${aln}/${query_names[i]}/no_phab_bench/result*.gz
    done
done

deactivate
