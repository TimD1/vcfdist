#!/bin/bash

Y='\033[0;33m'
W='\033[0m'

source lib/venv3/bin/activate
echo -e "\n\n${Y}[truth VCF to FASTA]${W}"
python3.7 vcf_to_fasta.py \
    data/truth.vcf.gz \
    data/ref.fasta \
    chr20
echo -e "\n\n${Y}[calls VCF to FASTA]${W}"
python3.7 vcf_to_fasta.py \
    data/calls.vcf.gz \
    data/ref.fasta \
    chr20
deactivate
