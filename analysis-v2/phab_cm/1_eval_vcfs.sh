#!/bin/bash

source globals.sh

# # truvari phab VCF normalization v4.1.0 (doesn't work)
# mkdir -p $out/phab-mafft
# mkdir -p $out/phab-wfa
# source ~/truvari/venv3.10/bin/activate
# for i in "${!query_names[@]}"; do
#     echo "truvari v4.1: phab wfa for '${query_names[i]}'"
#     truvari phab \
#         -r $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         --bSamples HG002 \
#         --cSamples HG002 \
#         -f $data/refs/$ref_name \
#         --align wfa \
#         -t 64 \
#         -o $out/phab-wfa/${query_names[i]}.vcf.gz \
#         2> $out/phab-wfa/${query_names[i]}.log
#     echo "truvari v4.1: phab mafft for '${query_names[i]}'"
#     truvari phab \
#         -r $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
#         --bSamples HG002 \
#         --cSamples HG002 \
#         -f $data/refs/$ref_name \
#         --align mafft \
#         -t 64 \
#         -o $out/phab-mafft/${query_names[i]}.vcf.gz \
#         2> $out/phab-mafft/${query_names[i]}.log
# done
# deactivate

# # truvari phab VCF normalization v4.0
# source ~/software/Truvari-v4.0.0/venv3.10/bin/activate
# mkdir -p $out/phab
# for i in "${!query_names[@]}"; do
#     echo "truvari v4.0: phab mafft for '${query_names[i]}'"
#     truvari phab \
#         -r $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#         -b $data/${truth_name}-${truth_version}/split/${truth_name}.all.vcf.gz \
#         -c $data/${query_names[i]}-${query_versions[i]}/split/${query_names[i]}.all.vcf.gz \
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
#     in_file="$out/phab/${query_names[i]}.vcf.gz"
#     for s in "${!in_samples[@]}"; do
#         bcftools view -c1 -s ${in_samples[s]} $in_file | bcftools reheader -s ${out_samples[s]}.txt | 
#         sed "s/\//|/g" | bcftools view -Oz &> $out/phab/${query_names[i]}.${out_samples[s]}.vcf.gz
#         tabix -p vcf $out/phab/${query_names[i]}.${out_samples[s]}.vcf.gz
#     done
# done

# # vcfdist evaluation
# mkdir -p $out/vcfdist
# for i in "${!query_names[@]}"; do
#     echo "vcfdist: evaluating '${query_names[i]}'"
#     $timer -v /home/timdunn/vcfdist/src/vcfdist \
#         $out/phab/${query_names[i]}.QUERY.vcf.gz \
#         $out/phab/${query_names[i]}.TRUTH.vcf.gz \
#         $data/refs/$ref_name \
#         --bed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#         -l 1000 \
#         -p $out/vcfdist/${query_names[i]}. \
#         2> $out/vcfdist/${query_names[i]}.log
# done

# # vcfeval evaluation
# mkdir -p $out/vcfeval
# for i in "${!query_names[@]}"; do
#     echo "vcfeval: evaluating '${query_names[i]}'"
#     $timer -v $rtg vcfeval \
#         -b $out/phab/${query_names[i]}.TRUTH.vcf.gz \
#         -c $out/phab/${query_names[i]}.QUERY.vcf.gz \
#         -t $data/refs/${ref_name}.sdf \
#         -e $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
#         --threads 64 \
#         --all-records \
#         -o $out/vcfeval/${query_names[i]} \
#         2> $out/vcfeval/${query_names[i]}.log
# done

# Truvari evaluation
mkdir -p $out/truvari
source ~/software/Truvari-v4.0.0/venv3.10/bin/activate
for i in "${!query_names[@]}"; do
    echo "Truvari: evaluating '${query_names[i]}'"
    $timer -v truvari bench \
        -b $out/phab/${query_names[i]}.TRUTH.vcf.gz \
        -c $out/phab/${query_names[i]}.QUERY.vcf.gz \
        -f $data/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
        --includebed $data/giab-tr-v4.20/GIABTR.HG002.benchmark.chr20.bed \
        --no-ref a \
        --sizemin 1 --sizefilt 1 --sizemax 1000 \
        --pick multi \
        --typeignore \
        --dup-to-ins \
        -o $out/truvari/${query_names[i]}

    laytr tru2ga \
        -i $out/truvari/${query_names[i]}/ \
        -o $out/truvari/${query_names[i]}/result
    gunzip $out/truvari/${query_names[i]}/result*.gz
done
deactivate
