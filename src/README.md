
### Usage
```
Usage: vcfdist <query.vcf> <truth.vcf> <ref.fasta> [options]

Required:
  <STRING>  query.vcf   phased VCF file containing variant calls to evaluate 
  <STRING>  truth.vcf   phased VCF file containing ground truth variant calls 
  <STRING>  ref.fasta   FASTA file containing draft reference sequence 

Options:
  -b, --bed <STRING>
    BED file containing regions to evaluate
  
  -p, --prefix <STRING> [./]
    prefix for output files (directory needs a trailing slash)
  
  -v, --verbosity <INTEGER> [1]
    printing verbosity (0: succinct, 1: default, 2:verbose)
  
  -r, --realign-only
    standardize truth and query variant representations, then exit
  
  -q, --keep-query
    do not realign query variants, keep original representation
  
  -t, --keep-truth
    do not realign truth variants, keep original representation

  -x, --mismatch-penalty <INTEGER> [3]
    Smith-Waterman mismatch (substitution) penalty

  -o, --gap-open-penalty <INTEGER> [2]
    Smith-Waterman gap opening penalty

  -e, --gap-extend-penalty <INTEGER> [1]
    Smith-Waterman gap extension penalty

  --min-qual <INTEGER> [0]
    minimum variant quality, lower qualities ignored

  --max-qual <INTEGER> [60]
    maximum variant quality, higher qualities kept but thresholded

  -s, --smallest-variant <INTEGER> [1]
    minimum variant size, smaller variants ignored (SNPs are size 1)

  -l, --largest-variant <INTEGER> [5000]
    maximum variant size, larger variants ignored

  -i, --max-iterations <INTEGER> [4]
    maximum iterations for expanding/merging clusters

  -g, --supercluster-gap <INTEGER> [50]
    minimum base gap between independent superclusters

  -h, --help
    show this help message

  -a, --advanced
    show advanced options

  -c, --citation
    please cite vcfdist if used in your analyses

  -v, --version
    print vcfdist version (v1.2.3)
