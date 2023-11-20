#!/bin/bash

source globals.sh

for i in ${!query_names[@]}; do
    /home/timdunn/truvari/venv3.10/bin/whatshap compare \
        --names ${query_names[i]},${truth_name} \
        ./vcfs/${query_names[i]}.vcf \
        ./vcfs/${truth_name}.vcf \
        --tsv-pairwise whatshap/${query_names[i]}.summary.tsv \
        --switch-error-bed whatshap/${query_names[i]}.switches.bed
done

for i in ${!query_names[@]}; do
    echo "${query_names[i]} switches: `tail -n +2 whatshap/${query_names[i]}.summary.tsv | cut -f 12 | cut -f1 -d/ | paste -sd+ - | bc`"
    echo "${query_names[i]}    flips: `tail -n +2 whatshap/${query_names[i]}.summary.tsv | cut -f 12 | cut -f2 -d/ | paste -sd+ - | bc`"
done
