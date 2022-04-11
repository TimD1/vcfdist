#!/bin/bash

G='\033[0;32m'
Y='\033[0;33m'
W='\033[0m'
source lib/venv3/bin/activate
echo -e "\n\n${G}[BAM to VCF]${W}"
python3.7 bam_to_vcf.py input/truth_orig.bam input/ref.fasta chr20
deactivate
