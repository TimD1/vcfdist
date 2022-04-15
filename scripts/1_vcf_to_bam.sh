#!/bin/bash

ref="data/ref.fasta"
fastas="data/truth data/calls"
Y='\033[0;33m'
W='\033[0m'

source lib/venv3/bin/activate

for fasta in $fastas; do
    echo -e "\n\n${Y}[$fasta VCF to orig BAM]${W}"
    python3.7 vcf_to_sam.py \
        ${fasta}.vcf.gz \
        data/ref.fasta \
        chr20 |
    samtools view \
        -@ `nproc` \
        -b |
    samtools sort \
        -@ `nproc` |
    samtools calmd \
        -@ `nproc` \
        -b \
        -Q \
        - \
        $ref \
    > ${fasta}_orig.bam
    samtools index ${fasta}_orig.bam
done
deactivate
