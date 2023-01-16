
### Usage
```
Usage: vcfdist <query.vcf> <truth.vcf> <ref.fasta> [options]

Required:
    <FILENAME>query.vcf     phased VCF file containing variant calls to evaluate 
    <FILENAME>truth.vcf     phased VCF file containing ground truth variant calls 
    <FILENAME>ref.fasta     FASTA file containing draft reference sequence 

Options:
    -b, --bed <FILENAME>
        BED file containing regions to evaluate

    -p, --prefix <FILENAME_PREFIX> [./]
        output filepath prefix (directories should contain trailing slashes)

    -q, --min-qual <VALUE> [0]
        minimum variant quality to be considered (lower qualities ignored)

    -m, --max-qual <VALUE> [60]
        maximum variant quality (higher qualities kept, thresholded)

    -i, --max-cluster-iterations<VALUE> [4]
        maximum iterations for growing clusters

    -g, --supercluster-gap <VALUE> [50]
        minimum base gap between independent superclusters

    -x, --mismatch-penalty <VALUE> [3]
        integer Smith-Waterman substitution penalty

    -o, --gap-open-penalty <VALUE> [2]
        integer Smith-Waterman gap opening penalty

    -e, --gap-extend-penalty <VALUE> [1]
        integer Smith-Waterman gap extension penalty

    -v, --verbosity <VALUE> [0]
        printing verbosity (0: default, 1: verbose, 2:debug, 3:verbose debug)

    -r, --realign-only
        realign truth/query VCFs with Smith-Waterman params, then exit

    -h, --help
        show this help message

    -a, --advanced
        show advanced options

    --version
        show version number (vcfdist v0.0.1)
```
