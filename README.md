# vcfdist: benchmarking phased variant calls
![build](https://github.com/timd1/vcfdist/actions/workflows/build.yml/badge.svg)
<!-- ![Github All Releases](https://img.shields.io/github/downloads/timd1/vcfdist/total.svg) -->

## Overview

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Variant Stratification](#variant-stratification)
* [Acknowledgements](#acknowledgements)
* [Limitations](#limitations)
* [License](#license)


## Introduction
vcfdist is a distance-based **variant calling evaluation tool** that:
- simultaneously evaluates **SNPS, small INDELs, complex, tandem repeat, and structural variants**
- gives **partial credit** to variant calls which are mostly correct
- **standardizes** query and truth VCF **variants** to a consistent representation
- **requires local phasing** of both input VCFs and enforces correct local phasing of variants
- **discovers** long-range variant **representation dependencies** using a novel clustering algorithm
- works on **monoploid and diploid** VCF contigs

This results in more stable SNP and INDEL precision-recall curves than previous work, particularly for complex variants. vcfdist also reports alignment distance based metrics for evaluation which are entirely independent of variant representation, providing greater insight into variant calling performance.

This project is currently under active development. We welcome the submission of any feedback, issues, or suggestions for improvement!


### Citation
Please cite the following pre-print if you use vcfdist:

<details>
<summary>
<a href="https://www.biorxiv.org/content/10.1101/2023.03.10.532078v2" target="_blank"><b>[bioRxiv]</b> vcfdist: Accurately benchmarking phased small variant calls in human genomes</a>
</summary>

<pre>
@article {dunn2023vcfdist,
    author = {Dunn, Tim and Narayanasamy, Satish},
    title = {vcfdist: Accurately benchmarking phased small variant calls in human genomes},
    elocation-id = {2023.03.10.532078},
    year = {2023},
    doi = {10.1101/2023.03.10.532078},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/10.1101/2023.03.10.532078v2},
    eprint = {https://biorxiv.org/content/10.1101/2023.03.10.532078.full.pdf},
    journal = {bioRxiv}
}
</pre>
</details>


## Installation

vcfdist is developed for Linux and its only dependencies are GCC v8+ and HTSlib. Please note that on Mac, `g++` is aliased to `clang`, which is currently not supported. If you don't have HTSlib already, please set it up as follows:
```bash
> wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
> tar -xvf htslib-1.17.tar.bz2
> cd htslib-1.17
> ./configure --prefix=/usr/local
> make
> sudo make install
> export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```
If you do already have HTSlib installed elsewhere, make sure you've added it to your `LD_LIBRARY_PATH`. At this point, installation is as simple as cloning the repository and building the executable. It should compile in less than one minute.

```bash
> git clone https://github.com/timd1/vcfdist
> cd vcfdist/src
> make
> ./vcfdist --version
vcfdist v2.0.0
```


## Usage

The <a href="./demo">`demo`</a> directory contains a <a href="./demo/demo.sh">demo script</a> (shown below) and all required inputs. It operates on the first 5 million bases on `chr1`, and should run in about 3 seconds.
```bash
../src/vcfdist \
    query.vcf \
    nist-v4.2.1_chr1_5Mb.vcf.gz \
    GRCh38_chr1_5Mb.fa \
    -b nist-v4.2.1_chr1_5Mb.bed \
    -p results/ \
    -v 0
```

You can expect to see <a href="./demo/output.txt">this output</a>.

To include more details on intermediate results, run it again at higher verbosity by removing the `-v 0` flag.
Please note that your results may not be identical, since vcfdist is under active development and handling of edge-cases may differ between versions.

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
Ensure that `strat.tsv` contains one stratification region per line; each line consists of a region name and BED file name separated by a tab.
GIAB stratification regions for GRCh38 can be found <a href="https://github.com/genome-in-a-bottle/genome-stratifications/tree/master/GRCh38">here</a>.


## Acknowledgements
Datasets used in the evaluation of the accompanying paper are listed <a href="./data/README.md">here</a>.

## Limitations
The current version of vcfdist is not designed to support:
- overlapping or unphased variants
- polyploid chromosomes

## License
This project is covered under the <a href="LICENSE">GNU GPL v3</a> license.
