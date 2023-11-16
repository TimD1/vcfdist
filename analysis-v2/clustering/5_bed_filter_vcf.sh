#!/bin/bash

bcftools view \
    -i 'TYPE=="SNP" || (ILEN < 1000 && ILEN > -1000)'\
    -R /home/timdunn/vcfdist/data/t2t-q100-v0.9/GRCh38_HG2-T2TQ100-V0.9_dipcall-z2k.benchmark.bed \
    -o ./t2t_t2t-bed.vcf \
    /home/timdunn/vcfdist/data/t2t-q100-v0.9/split/t2t-q100.all.vcf.gz

bcftools view \
    -i 'TYPE=="SNP" || (ILEN < 1000 && ILEN > -1000)'\
    -R /home/timdunn/vcfdist/data/t2t-q100-v0.9/GRCh38_HG2-T2TQ100-V0.9_dipcall-z2k.benchmark.bed \
    -o ./pav_t2t-bed.tmp.vcf \
    /home/timdunn/vcfdist/data/pav-v4/split/pav.all.vcf.gz

cat ./pav_t2t-bed.tmp.vcf | grep -v "INV"> pav_t2t-bed.vcf # remove large inversions
