# vcfdist: benchmarking phased small variant calls

<!-- [![DOI](https://zenodo.org/badge/365294513.svg)](https://zenodo.org/badge/latestdoi/365294513) -->

## Introduction
vcfdist is a distance-based small variant calling evaluator that:
- standardizes query and truth VCF variants to a consistent representation
- discovers long-range variant representation dependencies using a novel clustering algorithm
- requires local phasing of both input VCFs and enforces correct local phasing of variants
- gives partial credit to variant calls which are mostly correct

This results in more stable SNP and INDEL precision-recall curves than vcfeval, particularly for complex variants. vcfdist also reports alignment distance based metrics for evaluation which are entirely independent of variant representation, providing greater insight into variant calling performance.

This project is currently under active development. We welcome the submission of any feedback, issues, or suggestions for improvement!


## Citation
Please cite the following pre-print if you use vcfdist:

<details>
<summary>
<a href="https://www.biorxiv.org/content/early/2023/03/12/2023.03.10.532078"><b>[bioRxiv]</b> vcfdist: Accurately benchmarking phased small variant calls in human genomes</a>
</summary>

<pre>
@article {dunn2023vcfdist,
author = {Dunn, Tim and Narayanasamy, Satish},
title = {vcfdist: Accurately benchmarking phased small variant calls in human genomes},
elocation-id = {2023.03.10.532078},
year = {2023},
doi = {10.1101/2023.03.10.532078},
publisher = {Cold Spring Harbor Laboratory},
URL = {https://www.biorxiv.org/content/early/2023/03/12/2023.03.10.532078},
eprint = {https://biorxiv.org/content/early/2023/03/12/2023.03.10.532078.full.pdf},
journal = {bioRxiv}
}
</pre>
</details>


## Contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Variant Stratification](#variant-stratification)
* [Acknowledgements](#acknowledgements)
* [Limitations](#limitations)


## Installation

vcfdist's only dependencies are GCC v8+ and HTSlib. Please note that on Mac, `g++` is aliased to `clang`, which is currently not supported. If you don't have HTSlib already, please set it up as follows:
```bash
> wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
> tar -xvf htslib-1.17.tar.bz2
> cd htslib-1.17
> ./configure --prefix=/usr/local
> make
> sudo make install
> export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```
If you do already have HTSlib installed elsewhere, make sure you've added it to your `LD_LIBRARY_PATH`. At this point, installation is as simple as cloning the repository and building the executable:

```bash
> git clone https://github.com/timd1/vcfdist
> cd vcfdist/src
> make
> ./vcfdist -h
```


## Usage

Here's an example usage of vcfdist:

```bash
> ./vcfdist \
    query.vcf.gz \
    truth.vcf.gz \
    reference.fasta \
    -b analysis-regions.bed \
    -p output-prefix/
```
Please see additional options documented <a href="./src/README.md">here</a>, or run `./vcfdist --help`.

The output TSV files are documented <a href="./docs/outputs.md">here</a>.


## Variant Stratification

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
Ensure that `strat.tsv` contains one stratification region per line. 
Each line must contain a region name and BED file name, separated by a tab (`\t`).
GIAB stratification regions for GRCh38 can be found <a href="https://github.com/genome-in-a-bottle/genome-stratifications/tree/master/GRCh38">here</a>.


## Acknowledgements
Datasets used in the evaluation of the accompanying paper are listed <a href="./data/README.md">here</a>.

## Limitations
The current version of vcfdist is not designed to support:
- overlapping or unphased variants
- haploid or polyploid chromosomes
- structural variant evaluations
