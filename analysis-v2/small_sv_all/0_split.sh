#!/bin/bash
data="/home/timdunn/vcfdist/data"
src="pav"
version="v4"
base="${data}/${src}-${version}/split/${src}"

# split VCF into sub-VCFs of SNPs, INDELs, SVs
echo "$src snp"
bcftools view \
    -i 'TYPE=="SNP"'\
    -Oz \
    -o ${base}.snp.vcf.gz \
    ${base}.all.vcf.gz
tabix -p vcf ${base}.snp.vcf.gz

echo "$src indel"
bcftools view \
    -i 'ILEN < 50 && ILEN > -50'\
    -Oz \
    -o ${base}.indel.vcf.gz \
    ${base}.all.vcf.gz
tabix -p vcf ${base}.indel.vcf.gz

echo "$src sv"
bcftools view \
    -i 'ILEN <= -50 || ILEN >= 50'\
    -Oz \
    -o ${base}.sv.vcf.gz \
    ${base}.all.vcf.gz
tabix -p vcf ${base}.sv.vcf.gz

echo "$src small"
bcftools view \
    -i 'TYPE=="SNP" || (ILEN < 50 && ILEN > -50)'\
    -Oz \
    -o ${base}.small.vcf.gz \
    ${base}.all.vcf.gz
tabix -p vcf ${base}.small.vcf.gz

echo "$src large"
bcftools view \
    -i 'ILEN <= -1 || ILEN >= 1'\
    -Oz \
    -o ${base}.large.vcf.gz \
    ${base}.all.vcf.gz
tabix -p vcf ${base}.large.vcf.gz
