#!/bin/bash

source globals.sh

# # truvari phab VCF normalization v4.2.1
# source ~/software/Truvari-4.2.1/venv3.10/bin/activate
# mkdir -p $out/phab
# for i in "${!query_names[@]}"; do
#     echo "phab: mafft for '${query_names[i]}'"
#     truvari phab \
#         -r $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.most.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.most.vcf.gz \
#         --bSamples HG002 \
#         --cSamples HG002 \
#         -f $data/refs/$ref_name \
#         -t 64 \
#         -o $out/phab/${query_names[i]}.vcf.gz \
#         2> $out/phab/${query_names[i]}.log
# done
# deactivate

# # split phab normalized VCF into truth/query
# in_samples=(
#     "HG002"
#     "p:HG002"
# )
# out_samples=(
#     "TRUTH"
#     "QUERY"
# )
# for s in "${out_samples[@]}"; do 
#     echo $s > "${s}.txt" 
# done
# for i in "${!query_names[@]}"; do
#     for s in "${!in_samples[@]}"; do
#         bcftools view -c1 \
#                 -s ${in_samples[s]} \
#                 $out/phab/${query_names[i]}.vcf.gz \
#             | bcftools reheader -s ${out_samples[s]}.txt \
#             | sed "s/\//|/g" \
#             | bcftools norm -a -m-any \
#             | bcftools view -Oz \
#             &> $out/phab/${query_names[i]}.${out_samples[s]}.vcf.gz
#         tabix -f -p vcf $out/phab/${query_names[i]}.${out_samples[s]}.vcf.gz
#     done
# done

# vcfdist evaluation
mkdir -p $out/phab-vcfdist
for i in "${!query_names[@]}"; do
    echo "phab-vcfdist: evaluating '${query_names[i]}'"
    $timer -v /home/timdunn/vcfdist/src/vcfdist \
        $out/phab/${query_names[i]}.QUERY.vcf.gz \
        $out/phab/${query_names[i]}.TRUTH.vcf.gz \
        $data/refs/$ref_name \
        --bed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
        -l 1000 \
        -p $out/phab-vcfdist/${query_names[i]}. \
        2> $out/phab-vcfdist/${query_names[i]}.log
done

# vcfeval evaluation
mkdir -p $out/phab-vcfeval
for i in "${!query_names[@]}"; do
    echo "phab-vcfeval: evaluating '${query_names[i]}'"
    rm -r $out/phab-vcfeval/${query_names[i]}
    $timer -v $rtg vcfeval \
        -b $out/phab/${query_names[i]}.TRUTH.vcf.gz \
        -c $out/phab/${query_names[i]}.QUERY.vcf.gz \
        -t $data/refs/${ref_name}.sdf \
        --bed-regions $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
        --evaluation-regions $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
        --threads 64 \
        -m ga4gh \
        --ref-overlap \
        --all-records \
        -o $out/phab-vcfeval/${query_names[i]} \
        2> $out/phab-vcfeval/${query_names[i]}.log
    gunzip $out/phab-vcfeval/${query_names[i]}/output.vcf.gz
done

# Truvari evaluation
mkdir -p $out/phab-truvari
source ~/software/Truvari-4.2.1/venv3.10/bin/activate
for i in "${!query_names[@]}"; do
    echo "phab-truvari: evaluating '${query_names[i]}'"
    rm -r $out/phab-truvari/${query_names[i]}
    $timer -v truvari bench \
        -b $out/phab/${query_names[i]}.TRUTH.vcf.gz \
        -c $out/phab/${query_names[i]}.QUERY.vcf.gz \
        --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
        --no-ref a \
        --sizemin 1 --sizefilt 1 --sizemax 1000 \
        --pick single \
        --typeignore \
        --dup-to-ins \
        -o $out/phab-truvari/${query_names[i]} \
        2> $out/phab-truvari/${query_names[i]}.log
    truvari ga4gh \
        -i $out/phab-truvari/${query_names[i]}/ \
        -o $out/phab-truvari/${query_names[i]}/result
    gunzip $out/phab-truvari/${query_names[i]}/result*.gz
done
deactivate
