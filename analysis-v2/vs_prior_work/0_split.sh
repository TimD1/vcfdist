#!/bin/bash
data="/home/timdunn/vcfdist/data"
src="pav"
version="v4"
base="${data}/${src}-${version}/split/${src}"

# filter SVs > 1000bp and INV
echo "$src most"
bcftools view \
    -i 'TYPE=="SNP" || (ILEN < 1000 && ILEN > -1000)'\
    ${base}.all.vcf.gz |
    grep -v "INV" > ${base}.most.vcf
bgzip -f ${base}.most.vcf
tabix -p vcf ${base}.most.vcf.gz
