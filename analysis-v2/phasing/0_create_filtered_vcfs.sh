#!/bin/bash

echo "filtering pav"
bcftools filter \
    ~/vcfdist/data/pav-v4/split/pav.all.vcf.gz \
    -R ~/vcfdist/data/t2t-q100-v0.9/split/bench.bed \
    -o ./vcfs/pav.vcf

echo "filtering t2t-q100"
bcftools filter \
    ~/vcfdist/data/t2t-q100-v0.9/split/t2t-q100.all.vcf.gz \
    -R ~/vcfdist/data/t2t-q100-v0.9/split/bench.bed \
    -o ./vcfs/t2t-q100.vcf
