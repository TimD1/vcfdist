# vcfdist: Accurately benchmarking complex small variant calls

<!-- [![DOI](https://zenodo.org/badge/365294513.svg)](https://zenodo.org/badge/latestdoi/365294513) -->

## Introduction
vcfdist is a distance-based small variant calling evaluator that:
- standardizes query and truth VCF variants to a consistent representation
- discovers long-range variant representation dependencies using a novel clustering algorithm
- requires local phasing of both input VCFs and enforces correct local phasing of variants
- gives partial credit to variant calls which are mostly correct

This results in more stable SNP and INDEL precision-recall curves than vcfeval, particularly for complex variants. vcfeval also reports alignment distance based metrics for evaluation which are entirely independent of variant representation, providing greater insight into variant calling performance.

<div align="center">
<img src="img/f1_overview.png" width="900p" alt="vcfdist vs vcfeval precision-recall curve stability">
</div>


## Citation
Please cite the following pre-print if you use vcfdist:

<details>
<summary>
<a href=""><b>[bioRxiv]</b> vcfdist: Accurately benchmarking small complex variant calls in human genomes</a>
</summary>

<pre>
@article {dunn-vcfdist,
author = {Dunn, Tim and Narayanasamy, Satish},
title = {vcfdist: vcfdist: Accurately benchmarking small complex variant calls in human genomes},
elocation-id = {},
year = {2023},
doi = {},
publisher = {},
URL = {},
eprint = {},
journal = {bioRxiv}
}
</pre>
</details>

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Acknowledgements](#acknowledgements)

## Installation

vcfdist's only dependencies are a C++20 compliant compiler (`g++` will do the trick) and HTSlib, so if you don't have it already, please download from <a href="http://www.htslib.org/">here</a>. Then, add HTSlib to your library path:
```bash
> export LD_LIBRARY_PATH=/path/to/htslib-1.16:$LD_LIBRARY_PATH
```
At this point, installation is as simple as cloning the repository and building the executable:

```bash
> git clone https://github.com/timd1/vcfdist
> cd vcfdist/src
> make
> ./vcfdist --version
vcfdist 1.0.0
```


## Usage

Here's an example usage vcfdist:

```bash
./vcfdist \
    query.vcf.gz \
    truth.vcf.gz \
    reference.fasta \
    -b analysis-regions.bed \
    -p results/output-prefix
```
Please see additional options documented <a href="./img">here</a>, or run `./vcfdist --help`.


## Acknowledgements
Datasets used in the evaluation of the accompanying paper are listed <a href="./data/README.md">here</a>.