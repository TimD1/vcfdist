### Data
This directory includes the following subdirectories containing the inputs used to evaluate VCFdist. If you would like to reproduce our analyses using our scripts from the `pipeline` directory, please create the directory structure shown below, using publicly available data.

```
pfda-v2/
nist-v4.2.1/
cmrg-v1.00/
refs/
```

#### Precision FDA Challenge v2 Submission Inputs
```
pfda-v2/                                              [1] NIST pFDA TCV2 Submission Data
    submission_vcfs/                                  (2) pFDA TCV2 Submission VCFs
        <id>/
            <id>_HG002.vcf.gz
    input_fastqs/                                     (3) pFDA TCV2 Read FASTQs
        HG002*.fastq.gz
```
1. [NIST pFDA Truth Challenge V2 Submission Data](https://data.nist.gov/od/id/mds2-2336)


#### NIST Genome In A Bottle Benchmarks v4.2.1
```
nist-v4.2.1/                                          [1]
    GRCh38/
        HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
        HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf.gz [2]
```
1. [NIST/GIAB Benchmarking BEDs](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/)
2. [Phased GRCh38 Truth VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/)


#### Challenging Medically Relevent Genes Benchmark v1.00
```
cmrg-v1.00/                                           [1]
    GRCh38/
        HG002_GRCh38_CMRG_smallvar_v1.00_draft.bed
        HG002_GRCh38_CMRG_smallvar_v1.00_draft.vcf.gz
```
1. [CMRG Truth VCFs and Benchmarking BEDs](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/)

#### References
```
refs/
    GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta [1]
```
1. [GRCh38 Reference](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38)
