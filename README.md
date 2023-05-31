# vcfdist: benchmarking phased small variant calls
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
vcfdist is a distance-based small variant calling evaluator that:
- gives partial credit to variant calls which are mostly correct
- standardizes query and truth VCF variants to a consistent representation
- **requires local phasing of both input VCFs** and enforces correct local phasing of variants
- discovers long-range variant representation dependencies using a novel clustering algorithm
- works on monoploid and diploid VCF contigs

This results in more stable SNP and INDEL precision-recall curves than vcfeval, particularly for complex variants. vcfdist also reports alignment distance based metrics for evaluation which are entirely independent of variant representation, providing greater insight into variant calling performance.

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
> git clone --branch v1.3.1 https://github.com/timd1/vcfdist
> cd vcfdist/src
> make
> ./vcfdist --version
vcfdist v1.3.1
```


## Usage

The `demo` directory contains all input files required to run `vcfdist`. This demonstration operates on the first 5 million bases on `chr1`, and should run in about 3 seconds.
```bash
./vcfdist \
    ../demo/query.vcf \
    ../demo/nist-v4.2.1_chr1_5Mb.vcf.gz \
    ../demo/GRCh38_chr1_5Mb.fa \
    -b ../demo/nist-v4.2.1_chr1_5Mb.bed \
    -p ../demo/results/ \
    -v 0
```

You can expect to see the following output:

```
PRECISION-RECALL SUMMARY
 
TYPE   MIN_QUAL        TRUTH_TP        QUERY_TP        TRUTH_FN        QUERY_FP        PREC            RECALL          F1_SCORE        F1_QSCORE
SNP    Q >= 0          8220            8220            5               6               0.999149        0.999271        0.999210        31.023401
SNP    Q >= 0          8220            8220            5               6               0.999149        0.999271        0.999210        31.023401
 
TYPE   MIN_QUAL        TRUTH_TP        QUERY_TP        TRUTH_FN        QUERY_FP        PREC            RECALL          F1_SCORE        F1_QSCORE
INDEL  Q >= 0          932             932             3               3               0.996793        0.996261        0.996527        24.592749
INDEL  Q >= 0          932             932             3               3               0.996793        0.996261        0.996527        24.592749
 
 
ALIGNMENT DISTANCE SUMMARY
 
TYPE   MIN_QUAL        EDIT_DIST       DISTINCT_EDITS  ED_QSCORE       DE_QSCORE       ALN_QSCORE
ALL    Q >= 0          26              16              26.566509       27.579178       27.154125
ALL    Q >= 0          26              16              26.566509       27.579178       27.154125
ALL    Q >= 61         11793           9163            0.000000        0.000000        0.000000
 
TYPE   MIN_QUAL        EDIT_DIST       DISTINCT_EDITS  ED_QSCORE       DE_QSCORE
SNP    Q >= 0          10              10              29.151360       29.151360
SNP    Q >= 0          10              10              29.151360       29.151360
SNP    Q >= 61         8225            8225            0.000000        0.000000
 
TYPE   MIN_QUAL        EDIT_DIST       DISTINCT_EDITS  ED_QSCORE       DE_QSCORE
INDEL  Q >= 0          16              6               23.483049       21.940516
INDEL  Q >= 0          16              6               23.483049       21.940516
INDEL  Q >= 61         3568            938             0.000000        0.000000
```

To include more details on intermediate results, run it again at higher verbosity by removing the `-v 0` flag.
Please note that your results may not be identical to the example shown, since vcfdist is under active development and handling of certain edge cases may differ between versions.

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
- polyploid chromosomes

## License
This project is covered under the <a href="LICENSE">GNU GPL v3</a> license.
