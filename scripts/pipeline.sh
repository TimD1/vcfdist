#!/bin/bash

G='\033[0;32m'
W='\033[0m'

echo -e "\n\n${G}[converting truth/call VCFs to BAMs]${W}"
./1_vcf_to_bam.sh

echo -e "\n\n${G}[converting truth/call VCFs to FASTAs]${W}"
./1a_vcf_to_fasta.sh

echo -e "\n\n${G}[(re)aligning FASTAs to BAMs]${W}"
./1b_fasta_to_bam.sh

echo -e "\n\n${G}[converting BAMs to VCFs]${W}"
./2_bam_to_vcf.sh
