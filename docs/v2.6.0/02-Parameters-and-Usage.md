# Usage
vcfdist &lt;query.vcf&gt; &lt;truth.vcf&gt; &lt;ref.fasta&gt; [optional arguments]

> [!NOTE]
> The most exact clustering algorithm (biWFA) is enabled by default. It is not recommended if your dataset includes large variants.
> If evaluating structural variants, we recommend using `--cluster size 100 --max-supercluster-size 10000`.

# Required Arguments
### query.vcf
*type: string, default: none*

Phased VCF file containing variant calls to evaluate.

### truth.vcf
*type: string, default: none*

Phased VCF file containing ground truth variant calls.

### ref.fasta	
*type: string, default: none*

FASTA file containing the draft reference sequence.

# Optional Arguments

## Inputs/Outputs:
### -b, --bed
*type: string, default: none*

BED file containing regions to evaluate. Variants located on the border of a BED region are currently excluded from the evaluation (details [here](https://github.com/TimD1/vcfdist/wiki/03-Variant-Filtering#bed-region)).

### -v, --verbosity
*type: integer, default: 1*

Printing verbosity (0: succinct, 1: default, 2: verbose).
- Succinct: Only warnings, errors, and the precision-recall summary are logged to console.
- Default: High-level info on parsed variants, superclustering, phasing, output results, and timing is additionally logged.
- Verbose: For debugging; warnings are printed each time they occur with helpful data included.

### -p, --prefix
*type: string, default: ./*

Prefix for output files (directories need a trailing slash).
For example `-p results/` will store `results/summary.vcf`, `-p test_` will store `test_summary.vcf`.

### -n, --no-output-files
*type: flag*

Skip writing output files, only print summary to console. 

## Variant Filtering/Selection:
### -f, --filter
*type: comma-separated string, default: all variants pass filtering stage*

Select just variants passing these FILTERs (OR operation).

### -l, --largest-variant
*type: integer, default: 5000*

Maximum variant size to be evaluated, larger variants are ignored.

### -sv, --sv-threshold
*type: integer, default: 50*

Variants of this size or larger are considered SVs, not INDELs. This is useful because precision-recall summary statistics are reported separately for SNPs, INDELs, and SVs.

### -mn, --min-qual
*type: integer, default: 0*

Minimum variant quality, lower quality variants are ignored.

### -mx, --max-qual
*type: integer, default: 60*

Maximum variant quality, higher quality variants are kept but their Q-score is thresholded to this value.

## Clustering:
### -c, --cluster
*type: string [and integer], default: biwfa*

Select clustering method, one of: `biwfa`, `size N`, and `gap N`
#### biwfa
Clusters are generated using bi-directional wave-front alignment, essentially an efficient algorithm for finding possible alternate alignments (and therefore if nearby variants are independent). See the papers on [BiWFA](https://academic.oup.com/bioinformatics/article/39/2/btad074/7030690) and [WFA](https://academic.oup.com/bioinformatics/article/37/4/456/5904262) for more details. This is the currently recommended (and default) vcfdist clustering algorithm because it is the most accurate; it will always find dependencies if they exist. However, when evaluating large structural variants (above 1kbp) it tends to create large clusters, which results in large memory usage and slower evaluations. For evaluating large variants, `--cluster size 100` may be preferable.

#### gap N
Gap-based clustering is the simplest and fastest clustering method: group together all variants less than N bases apart. It is also the least accurate, and will miss variant dependencies if N is too small. Conversely, as N nears the reciprocal of the background rate of genomic variation between humans (one SNP every 1000 bases), clusters will grow to be very large. We recommend 50 < N < 200, and to limit evaluations to small variants when using this option.

#### size N
This is a heuristic that compromises in terms of efficiency and accuracy, basically extending the gap N heuristic to work with larger variants. Once a variant is larger than size N, the required gap to consider it independent of an adjacent variant is the size of the variant, instead of N. 

### -s, --max-supercluster-size
*type: integer, default: 10000*

The maximum size of a supercluster (group of variants that are evaluated together). A larger limit will result in higher RAM usage and a longer runtime, but may lead to a more accurate evaluation for a few edge cases.

### -i, --max-iterations
*type: integer, default: 4*

Maximum number of iterations for expanding/merging clusters, only applicable if `--cluster biwfa` is selected (which is the default).

## Re-Alignment:
### -rq, --realign-query
*type: flag*

Realign query variants using Smith-Waterman parameters -x -o -e

### -rt, --realign-truth
*type: flag*

Realign truth variants using Smith-Waterman parameters -x -o -e

### -ro, --realign-only
*type: flag*

Standardize truth and query variant representations, then exit.

### -x, --mismatch-penalty
*type: integer, default: 3*

Smith-Waterman mismatch (substitution) penalty.

### -o, --gap-open-penalty
*type: integer, default: 2*

Smith-Waterman gap opening penalty.

### -e, --gap-extend-penalty
*type: integer, default: 1*

Smith-Waterman gap extension penalty.

## Precision and Recall:
### -ct, --credit-threshold 
*type: float, default: 0.70*

Minimum partial credit (calculated as a fractional reduction in edit distance over if the variant is omitted) to consider a query variant a true positive. 

## Phasing:
### -pt, --phasing-threshold 
*type: float, default: 0.60*

Minimum fractional reduction in edit distance over other phasing in order to consider this supercluster phased. Phased superclusters are then used to calculate switch and flip errors.

## Distance:
### -d, --distance
*type: flag*

Flag to include alignment distance calculations, which are skipped by default.

### -ex, --eval-mismatch-penalty
*type: integer, default: 3*

Mismatch penalty (`--distance` evaluation only).

### -eo, --eval-gap-open-penalty
*type: integer, default: 2*

Gap opening penalty (`--distance` evaluation only).

### -ee, --eval-gap-extend-penalty
*type: integer, default: 1*

Gap extension penalty (`--distance` evaluation only).

## Utilization:
### -t, --max-threads
*type: integer, default: 64*

Maximum threads to use for clustering and precision/recall alignment.

### -r, --max-ram
*type: float, default: 64.00*

Approximate maximum RAM (measured in GB) to use for precision/recall alignment. Evaluation of superclusters requiring RAM usage above this threshold will still occur, but with a warning.

## Miscellaneous:
### -h, --help
*type: flag*

Prints a help message listing all required and optional command-line parameters.

### -a, --advanced
*type: flag*

Prints a help message listing all command-line parameters, including advanced options that are not recommended for most users.

### -ci, --citation
*type: flag*
Prints the BibTeX and MLA formatted citations for vcfdist.

### -v, --version
Prints the current version of vcfdist.