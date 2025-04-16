We currently output an intermediate VCF in GA4GH compatible format, meaning the results can be stratified and analyzed by `hap.py`'s quantification helper script `qfy.py`.
In order to use `qfy.py` please install <a href="https://github.com/Illumina/hap.py">`hap.py`</a>.
`tabix` and `bgzip` should already be included as part of HTSlib.

```bash
> ./vcfdist \                                 # run vcfdist
    query.vcf.gz \
    truth.vcf.gz \
    reference.fasta \
    -b analysis-regions.bed \
    -p output-prefix/
> bgzip output-prefix/summary.vcf             # compress summary VCF
> tabix -p vcf output-prefix/summary.vcf.gz   # index summary VCF
> export HGREF=/path/to/reference.fasta       # set reference path
> source /path/to/happy/venv2/bin/activate    # activate hap.py virtualenv
> python /path/to/happy/install/bin/qfy.py \  # run quantification script
    -t ga4gh \
    --write-vcf \
    --write-counts \
    --stratification strat.tsv \
    --roc QUAL \
    --o results/qfy-output-prefix \
    output-prefix/summary.vcf.gz
```
Ensure that `strat.tsv` contains one stratification region per line; each line consists of a region name and BED file name separated by a tab.
GIAB stratification regions for GRCh38 can be found <a href="https://github.com/genome-in-a-bottle/genome-stratifications/tree/master/GRCh38">here</a>.