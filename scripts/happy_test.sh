#!/bin/bash
source ~/software/happy/venv2/bin/activate

python ~/software/happy/install/bin/hap.py \
    ../test/truth.vcf.gz \
    ../test/calls.vcf.gz \
    -r /x/gm24385/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    -T ../scripts/data/bench_20.bed \
    --roc QUAL \
    --write-counts \
    --engine vcfeval \
    -o ../out/vcfeval

gunzip -f ../out/vcfeval.*.gz
deactivate
