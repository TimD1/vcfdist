#!/bin/bash

/home/timdunn/truvari/venv3.10/bin/whatshap compare \
    --names pav,t2t-q100 \
    ./vcfs/pav.vcf \
    ./vcfs/t2t-q100.vcf \
    --tsv-pairwise pav.tsv \
    --switch-error-bed switches.bed
