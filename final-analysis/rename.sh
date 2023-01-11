#!/bin/bash

source includes.sh

# $parallel -j20 \
#     "mv $data/pfda-v2/phased_vcfs/{1}/{2}_HG002.vcf.gz $data/pfda-v2/phased_vcfs/{1}/{2}_HG002_O.vcf.gz" ::: \
#             ${ref_ids[@]} ::: ${sub_ids[@]}

$parallel -j20 \
    "mv $data/pfda-v2/phased_vcfs/{1}/{2}_HG002.vcf.gz.tbi $data/pfda-v2/phased_vcfs/{1}/{2}_HG002_O.vcf.gz.tbi" ::: \
            ${ref_ids[@]} ::: ${sub_ids[@]}
