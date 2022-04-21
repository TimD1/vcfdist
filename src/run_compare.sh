#!/bin/bash

./compare \
    ../scripts/data/calls.vcf.gz \
    ../scripts/data/truth.vcf.gz \
    ../scripts/data/ref.fasta \
    -b ../scripts/data/all.bed
