#!/bin/bash

Y='\033[0;33m'
W='\033[0m'
source lib/venv3/bin/activate

# files="orig mm2-map-ont mm2-asm5 mm2-asm10 mm2-asm20"
files="npore"

for file in $files; do

    echo -e "\n${Y}> truth_${file}${W}"
    python3.7 bam_to_vcf.py \
        data/truth_${file}.bam \
        data/ref.fasta \
        chr20

    # echo -e "\n${Y}> calls_${file}${W}"
    # python3.7 bam_to_vcf.py \
    #     data/calls_${file}.bam \
    #     data/ref.fasta \
    #     chr20

done
deactivate
