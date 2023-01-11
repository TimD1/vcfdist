#!/bin/bash

source includes.sh

orig_truth="$data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf"
truth="$data/nist-v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_phased.vcf.gz"

# decompress truth
bgzip -d ${orig_truth}.gz

# incorrect FORMAT tags in MHC (on chr6)
sed -i 's/PATMAT:\.:\./PATMAT/g' ${orig_truth}
sed -i 's/HOMVAR:\.:\./HOMVAR/g' ${orig_truth}

# remove PS tags
bcftools annotate \
    -x FORMAT/PS ${orig_truth} | bgzip -c > ${truth}

# index truth
tabix -p vcf ${truth}

# recompress truth
bgzip ${orig_truth}
