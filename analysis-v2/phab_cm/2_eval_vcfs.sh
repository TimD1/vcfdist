#!/bin/bash

source globals.sh

# vcfdist evaluation
mkdir -p $out/vcfdist
for i in "${!query_names[@]}"; do
    echo "vcfdist: evaluating '${query_names[i]}'"
    $timer -v /home/timdunn/vcfdist/src/vcfdist \
        $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
        $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
        $data/refs/$ref_name \
        --bed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
        -t 64 -r 128 \
        -l 1000 \
        -p $out/vcfdist/${query_names[i]}. \
        2> $out/vcfdist/${query_names[i]}.log
done

# vcfeval evaluation
mkdir -p $out/vcfeval
for i in "${!query_names[@]}"; do
    echo "vcfeval: evaluating '${query_names[i]}'"
    rm -r $out/vcfeval/${query_names[i]}
    $timer -v $rtg vcfeval \
        -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
        -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
        -t $data/refs/${ref_name}.sdf \
        -e $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
        --threads 64 \
        --ref-overlap \
        --vcf-score-field=QUAL \
        -o $out/vcfeval/${query_names[i]} \
        2> $out/vcfeval/${query_names[i]}.log
    gunzip $out/vcfeval/${query_names[i]}/*.vcf.gz
done

source ~/software/Truvari-4.1.0/venv3.10/bin/activate

# # truvari evaluation
# mkdir -p $out/truvari
# for i in "${!query_names[@]}"; do
#     echo "truvari: evaluating '${query_names[i]}'"
#     rm -r $out/truvari/${query_names[i]}
#     $timer -v truvari bench \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         -f $data/refs/$ref_name \
#         --bSample HG002 \
#         --cSample HG002 \
#         --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.minus1.bed \
#         --no-ref a \
#         --sizemin 1 --sizefilt 1 --sizemax 1000 \
#         --pick single \
#         --typeignore \
#         --dup-to-ins \
#         -o $out/truvari/${query_names[i]} \
#         2> $out/truvari/${query_names[i]}.log
#     laytr tru2ga \
#         -i $out/truvari/${query_names[i]}/ \
#         -o $out/truvari/${query_names[i]}/result
#     gunzip $out/truvari/${query_names[i]}/result*.gz
# done

# # wfa/mafft refine
# for aln in "${truvari_refine_options[@]}"; do
#     for i in "${!query_names[@]}"; do
#         echo "\ntruvari: ${aln} refine '${query_names[i]}'"
#         mkdir -p $out/truvari-${aln}
#         rm -rf $out/truvari-${aln}/${query_names[i]}
#         cp -r $out/truvari/${query_names[i]} \
#               $out/truvari-${aln}/${query_names[i]}
#         $timer -v truvari refine \
#             -f $data/refs/$ref_name \
#             -t 64 \
#             --use-original-vcfs \
#             --recount \
#             --align ${aln} \
#             $out/truvari-${aln}/${query_names[i]} \
#             2> $out/truvari-${aln}/${query_names[i]}.log
#         echo "laytr: ${aln} refine '${query_names[i]}'"
#         laytr tru2ga \
#             -i $out/truvari-${aln}/${query_names[i]}/phab_bench/ \
#             -o $out/truvari-${aln}/${query_names[i]}/phab_bench/result
#         gunzip $out/truvari-${aln}/${query_names[i]}/phab_bench/result*.gz
#         rm -f $out/truvari-${aln}/${query_names[i]}/*.gz
#         rm -f $out/truvari-${aln}/${query_names[i]}/*.vcf
#         rm -f $out/truvari-${aln}/${query_names[i]}/*.tbi
#         echo "bedtools: exclude phab '${query_names[i]}'"
#         bedtools subtract \
#             -a $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.minus1.bed \
#             -b $out/truvari-${aln}/${query_names[i]}/candidate.refine.bed \
#             > $out/truvari-${aln}/${query_names[i]}/no-refine.bed
#         echo "truvari: bench excluded phab '${query_names[i]}'"
#         rm -rf $out/truvari-${aln}/${query_names[i]}/no_phab_bench
#         truvari bench \
#             -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#             -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#             -f $data/refs/$ref_name \
#             --bSample HG002 \
#             --cSample HG002 \
#             --includebed $out/truvari-${aln}/${query_names[i]}/no-refine.bed \
#             --no-ref a \
#             --sizemin 1 --sizefilt 1 --sizemax 1000 \
#             --pick single \
#             --typeignore \
#             --dup-to-ins \
#             -o $out/truvari-${aln}/${query_names[i]}/no_phab_bench \
#             2> $out/truvari-${aln}/${query_names[i]}/no_phab_bench.log
#         echo "laytr: ${aln} refine '${query_names[i]}'"
#         laytr tru2ga \
#             -i $out/truvari-${aln}/${query_names[i]}/no_phab_bench/ \
#             -o $out/truvari-${aln}/${query_names[i]}/no_phab_bench/result
#         gunzip $out/truvari-${aln}/${query_names[i]}/no_phab_bench/result*.gz
#     done
# done

deactivate
