# vcfdist: benchmarking phased variant calls
![build](https://github.com/timd1/vcfdist/actions/workflows/build.yml/badge.svg)
[![DOI](https://zenodo.org/badge/472945373.svg)](https://zenodo.org/badge/latestdoi/472945373)

## Overview

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Wiki](#wiki)
* [License](#license)


## Introduction
vcfdist is a distance-based **germline variant calling evaluation tool** that:
- simultaneously evaluates **SNPS, INDELs, complex, tandem repeat, and structural variants**
- **requires local phasing** information for truth and query variants
- **discovers** long-range variant **representation dependencies**
- works on **monoploid and diploid** contigs
- can **identify** variant calls which are **partially correct**
- can **standardize** query and truth VCF **variants** to a consistent representation
- can report **alignment distance** based **metrics**

This results in more stable and accurate SNP, INDEL, and SV precision-recall curves than previous work, particularly when complex variants are involved.

This project is currently under active development. We welcome the submission of any feedback, issues, or suggestions for improvement! Check out the [wiki](https://github.com/TimD1/vcfdist/wiki) for more information.


### Citation
Please cite the following works if you use vcfdist:

<details>
<summary>
<a href="https://doi.org/10.1038/s41467-023-43876-x" target="_blank"><b>[Nature Comms]</b> vcfdist: Accurately benchmarking phased small variant calls in human genomes</a>
</summary>

<pre>
@article{dunn2023vcfdist,
  author={Dunn, Tim and Narayanasamy, Satish},
  title={vcfdist: Accurately benchmarking phased small variant calls in human genomes},
  journal={Nature Communications},
  year={2023},
  volume={14},
  number={1},
  pages={8149},
  issn={2041-1723},
  doi={10.1038/s41467-023-43876-x},
  URL={https://doi.org/10.1038/s41467-023-43876-x}
}
</pre>
</details>

<details>
<summary>
<a href="https://doi.org/10.1101/2024.01.23.575922" target="_blank"><b>[bioRxiv]</b> Jointly benchmarking small and structural variant calls with vcfdist</a>
</summary>

<pre>
@article{dunn2024vcfdist,
  author={Dunn, Tim and Zook, Justin M and Holt, James M and Narayanasamy, Satish},
  title={Jointly benchmarking small and structural variant calls with vcfdist},
  journal={bioRxiv},
  year={2024},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2024.01.23.575922},
  URL={https://doi.org/10.1101/2024.01.23.575922}
}
</pre>
</details>


## Installation

### Option 1: GitHub Source
vcfdist is developed for Linux and its only dependencies are GCC v8+ and HTSlib. If you don't have HTSlib already, please set it up as follows:
```bash
> wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
> tar -xvf htslib-1.17.tar.bz2
> cd htslib-1.17
> ./configure --prefix=/usr/local
> make
> sudo make install
> export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```
If you do already have HTSlib installed elsewhere, make sure you've added it to your `LD_LIBRARY_PATH`, and that the HTSlib headers are included during compilation. At this point, installation is as simple as cloning the repository and building the executable. It should compile in less than one minute.

```bash
> git clone https://github.com/timd1/vcfdist
> cd vcfdist/src
> make
> sudo make install
> vcfdist --version
vcfdist v2.3.2
```

### Option 2: Docker Image
A pre-built Docker Hub image can be downloaded from <a href="https://hub.docker.com/r/timd1/vcfdist">here</a> using:
```bash
sudo docker pull timd1/vcfdist
sudo docker run -it timd1/vcfdist:latest vcfdist --help
```

## Usage

The <a href="./demo">`demo`</a> directory contains a <a href="./demo/demo.sh">demo script</a> (shown below) and all required inputs. It operates on the first 5 million bases on `chr1`, and should run in about 3 seconds.
```bash
vcfdist \
    query.vcf \
    nist-v4.2.1_chr1_5Mb.vcf.gz \
    GRCh38_chr1_5Mb.fa \
    -b nist-v4.2.1_chr1_5Mb.bed \
    -p results/ \
    -v 0
```

<!--- You can expect to see <a href="./demo/output.txt">this output</a>. -->
You can expect to see the following output:
```
PRECISION-RECALL SUMMARY

TYPE   THRESHOLD     TRUTH_TP  QUERY_TP  TRUTH_FN  QUERY_FP  PREC     RECALL   F1_SCORE  F1_QSCORE
SNP    NONE Q >= 0   8222      8222      1         2         0.9997   0.9998   0.9998    37.3885
SNP    BEST Q >= 0   8222      8222      1         2         0.9997   0.9998   0.9998    37.3885

INDEL  NONE Q >= 0   876       876       51        12        0.9864   0.9449   0.9652    14.5953
INDEL  BEST Q >= 0   876       876       51        12        0.9864   0.9449   0.9652    14.5953

SV     NONE Q >= 0   0         0         0         0         1.0000   1.0000   1.0000    100.000
SV     BEST Q >= 0   0         0         0         0         1.0000   1.0000   1.0000    100.000

ALL    NONE Q >= 0   9098      9098      52        14        0.9984   0.9943   0.9963    24.4200
ALL    BEST Q >= 0   9098      9098      52        14        0.9984   0.9943   0.9963    24.4200
```

To include more details on intermediate results, run it again at higher verbosity by removing the `-v 0` flag.

## Wiki

The [vcfdist wiki](https://github.com/TimD1/vcfdist/wiki) is currently a work-in-progress, but has helpful information on [command-line parameters](https://github.com/TimD1/vcfdist/wiki/02-Parameters-and-Usage) and [output documentation](https://github.com/TimD1/vcfdist/wiki/09-Outputs). If something isn't covered yet, just start a [discussion](https://github.com/TimD1/vcfdist/discussions) or file an [issue](https://github.com/TimD1/vcfdist/issues) and I'd be happy to answer.

## License
This project is covered under the <a href="LICENSE">GNU GPL v3</a> license.
