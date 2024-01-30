#!/bin/bash
echo -e "Running vcfdist..."
vcfdist \
    query.vcf \
    nist-v4.2.1_chr1_5Mb.vcf.gz \
    GRCh38_chr1_5Mb.fa \
    -b nist-v4.2.1_chr1_5Mb.bed \
    -p results/ \
    -v 0
echo -e "done!"

echo -en "Plotting results..."
python3 pr_plot.py
echo -e "done!"
