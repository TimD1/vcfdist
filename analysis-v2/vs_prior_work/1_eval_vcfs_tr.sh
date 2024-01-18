#!/bin/bash

source globals.sh

# EVALUATE TANDEM REPEAT REGIONS

# vcfdist evaluation
mkdir -p ./vcfdist_tr
for i in "${!query_names[@]}"; do
    echo "vcfdist: evaluating '${query_names[i]}'"
    $timer -v ../../src/vcfdist \
        $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
        $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
        $data/refs/$ref_name \
        --bed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
        -t 64 -r 128 \
        -l 1000 \
        -p ./vcfdist_tr/${query_names[i]}. \
        2> ./vcfdist_tr/${query_names[i]}.log
done

# vcfeval evaluation
mkdir -p ./vcfeval_tr
for i in "${!query_names[@]}"; do
    echo "vcfeval: evaluating '${query_names[i]}'"
    $timer -v $rtg vcfeval \
        -b $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
        -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
        -t $data/refs/${ref_name}.sdf \
        -e $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
        --threads 64 \
        --all-records \
        --ref-overlap \
        -o ./vcfeval_tr/${query_names[i]} \
        2> ./vcfeval_tr/${query_names[i]}.log
    gunzip ./vcfeval_tr/${query_names[i]}/*.vcf.gz
done

source ~/software/Truvari-4.2.0-dev/venv3.10/bin/activate

# truvari evaluation
mkdir -p ./truvari_tr
for i in "${!query_names[@]}"; do
    echo "truvari: evaluating '${query_names[i]}'"
    rm -rf ./truvari_tr/${query_names[i]}
    $timer -v truvari bench \
        -b $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
        -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
        -f $data/refs/$ref_name \
        --bSample HG002 \
        --cSample HG002 \
        --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
        --no-ref a \
        --sizemin 1 --sizefilt 1 --sizemax 1000 \
        --pick single \
        --typeignore \
        --dup-to-ins \
        -o ./truvari_tr/${query_names[i]} \
        2> ./truvari_tr/${query_names[i]}.log
    laytr tru2ga \
        -i ./truvari_tr/${query_names[i]}/ \
        -o ./truvari_tr/${query_names[i]}/result
    gunzip ./truvari_tr/${query_names[i]}/result*.gz
done

# wfa/mafft refine
for aln in "${truvari_refine_options[@]}"; do
    for i in "${!query_names[@]}"; do
        echo "truvari: ${aln} refine '${query_names[i]}'"
        mkdir -p ./truvari_tr-${aln}
        rm -rf ./truvari_tr-${aln}/${query_names[i]}
        cp -r ./truvari_tr/${query_names[i]} \
              ./truvari_tr-${aln}/${query_names[i]}
        $timer -v truvari refine \
            -f $data/refs/$ref_name \
            -t 64 \
            --regions ./truvari_tr-${aln}/${query_names[i]}/candidate.refine.bed \
            --use-original-vcfs \
            --recount \
            --align ${aln} \
            ./truvari_tr-${aln}/${query_names[i]} \
            2> ./truvari_tr-${aln}/${query_names[i]}.log
        echo "laytr: ${aln} refine '${query_names[i]}'"
        laytr tru2ga \
            -i ./truvari_tr-${aln}/${query_names[i]}/phab_bench/ \
            -o ./truvari_tr-${aln}/${query_names[i]}/phab_bench/result
        gunzip ./truvari_tr-${aln}/${query_names[i]}/phab_bench/result*.gz
        rm -f ./truvari_tr-${aln}/${query_names[i]}/*.gz
        rm -f ./truvari_tr-${aln}/${query_names[i]}/*.vcf
        rm -f ./truvari_tr-${aln}/${query_names[i]}/*.tbi
        echo "bedtools: exclude phab '${query_names[i]}' $aln"
        bedtools subtract \
            -a $data/giab-tr-v4.20/GIABTR.HG002.benchmark.regions.bed \
            -b ./truvari_tr-${aln}/${query_names[i]}/candidate.refine.bed \
            > ./truvari_tr-${aln}/${query_names[i]}/no-refine.bed
        echo "truvari: bench excluded phab '${query_names[i]}' $aln"
        rm -rf ./truvari_tr-${aln}/${query_names[i]}/no_phab_bench
        $timer -v truvari bench \
            -b $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
            -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
            -f $data/refs/$ref_name \
            --bSample HG002 \
            --cSample HG002 \
            --includebed ./truvari_tr-${aln}/${query_names[i]}/no-refine.bed \
            --no-ref a \
            --sizemin 1 --sizefilt 1 --sizemax 1000 \
            --pick single \
            --typeignore \
            --dup-to-ins \
            -o ./truvari_tr-${aln}/${query_names[i]}/no_phab_bench \
            2> ./truvari_tr-${aln}/${query_names[i]}/no_phab_bench.log
        echo "laytr: ${aln} refine '${query_names[i]}' $aln"
        laytr tru2ga \
            -i ./truvari_tr-${aln}/${query_names[i]}/no_phab_bench/ \
            -o ./truvari_tr-${aln}/${query_names[i]}/no_phab_bench/result
        gunzip ./truvari_tr-${aln}/${query_names[i]}/no_phab_bench/result*.gz
    done
done

deactivate
