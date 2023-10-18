### Usage
```
Usage: vcfdist <query.vcf> <truth.vcf> <ref.fasta> [options]

Required:
  <STRING>	query.vcf	phased VCF file containing variant calls to evaluate 
  <STRING>	truth.vcf	phased VCF file containing ground truth variant calls 
  <STRING>	ref.fasta	FASTA file containing draft reference sequence 

Options:

  Inputs/Outputs:
  -b, --bed <STRING>
      BED file containing regions to evaluate
  -p, --prefix <STRING> [./]
      prefix for output files (directory needs a trailing slash)
  -v, --verbosity <INTEGER> [1]
      printing verbosity (0: succinct, 1: default, 2:verbose)

  Variant Filtering:
  -s, --smallest-variant <INTEGER> [1]
      minimum variant size, smaller variants ignored (SNPs are size 1)
  -l, --largest-variant <INTEGER> [5000]
      maximum variant size, larger variants ignored
  --min-qual <INTEGER> [0]
      minimum variant quality, lower qualities ignored
  --max-qual <INTEGER> [60]
      maximum variant quality, higher qualities kept but thresholded

  ReAlignment:
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

  Clustering:
  --simple-cluster
      instead of biWFA-based clustering, use gap-based clustering 

  Utilization:
  --max-threads <INTEGER> [64]
      maximum threads to use for clustering and precision/recall alignment
  --max-ram <FLOAT> [64.000GB]
      (approximate) maximum RAM to use for precision/recall alignment

  Miscellaneous:
  -h, --help
      show this help message
  -a, --advanced
      show advanced options, not recommended for most users
  -c, --citation
      please cite vcfdist if used in your analyses: thanks!
  -v, --version
      print vcfdist version (v2.1.0)
```
