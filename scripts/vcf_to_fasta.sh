#!/bin/bash

G='\033[0;32m'
Y='\033[0;33m'
W='\033[0m'

source lib/venv3/bin/activate
echo -e "\n\n${G}[truth VCF to FASTA]${W}"
python3.7 vcf_to_fasta.py input/truth.vcf.gz input/ref.fasta chr20
echo -e "\n\n${G}[calls VCF to FASTA]${W}"
python3.7 vcf_to_fasta.py input/calls.vcf.gz input/ref.fasta chr20
deactivate
