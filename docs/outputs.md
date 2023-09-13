## vcfdist Outputs
Please see the detailed field descriptions at the bottom of this page.

### Precision/Recall Summary Metrics
#### precision-recall-summary.tsv
High-level precision/recall overview of SNP/INDEL performance. For each category (SNP/INDEL), there is one line for performance at the chosen minimum quality score, and one line for the quality score threshold that results in the best performance.

| VAR_TYPE | MIN_QUAL | TRUTH_TP | QUERY_TP | TRUTH_FN | QUERY_FP | PREC | RECALL | F1_SCORE | F1_QSCORE |
|-|-|-|-|-|-|-|-|-|-|

#### precision-recall.tsv
For each category (SNP/INDEL), there is one line reporting the precison/recall performance at each possible quality score.

| VAR_TYPE | MIN_QUAL | PREC | RECALL | F1_SCORE | TRUTH_TOTAL | TRUTH_TP | TRUTH_PP | TRUTH_FN | QUERY_TOTAL | QUERY_TP | QUERY_PP | QUERY_FP |
|-|-|-|-|-|-|-|-|-|-|-|-|-|

### Alignment Distance Summary Metrics
#### distance-summary.tsv
High-level alignment distance overview of SNP/INDEL performance. For each category (SNP/INDEL/ALL), there is one line for performance at the minimum, best, and maximum quality score threshold.

| VAR_TYPE | MIN_QUAL | EDIT_DIST | DISTINCT_EDITS | ED_QSCORE | DE_QSCORE | ALN_QSCORE |
|-|-|-|-|-|-|-|

#### distance.tsv
There is one line reporting the alignment distance performance at each possible quality score.

| MIN_QUAL | SUB_DE | INS_DE | DEL_DE | SUB_ED | INS_ED | DEL_ED | DISTINCT_EDITS | EDIT_DIST | ALN_SCORE | ALN_QSCORE |
|-|-|-|-|-|-|-|-|-|-|-|

### Detailed Results
#### edits.tsv
This file reports for each edit (where called query sequence differs from truth sequence) the contig, pos, hap, len, supercluster, and the quality range for which these edits occur.

| CONTIG | START | HAP | TYPE | SIZE | SUPERCLUSTER | MIN_QUAL | MAX_QUAL
|-|-|-|-|-|-|-|-|

#### phase-blocks.tsv
Reports the size, location, and composition of each phase block.

| CONTIG | START | STOP | SIZE | SUPERCLUSTERS |
|-|-|-|-|-|

#### superclusters.tsv
Reports the size, location, and composition of each supercluster.

| CONTIG | START | STOP | SIZE | QUERY1_VARS | QUERY2_VARS | TRUTH1_VARS | TRUTH2_VARS | ORIG_ED | SWAP_ED | PHASE | PHASE_BLOCK |
|-|-|-|-|-|-|-|-|-|-|-|-|

#### query.tsv, truth.tsv
Reports detailed information regarding each variant.

| CONTIG | POS | HAP | REF | ALT | QUAL | TYPE | ERR_TYPE | CREDIT | CLUSTER | SUPERCLUSTER | LOCATION |
|-|-|-|-|-|-|-|-|-|-|-|-|

### Output VCFs
#### orig-query.vcf, orig-truth.vcf
Original query and truth VCFs, as parsed by vcfdist.
#### query.vcf, truth.vcf
Output query and truth VCFs, standardized by vcfdist (at point C).

### Field Descriptions
| NAME | TYPE | DESCRIPTION |
|------|------|--------------------------------------------------------------|
| TYPE | string | Specific variant type (REF/SNP/INS/DEL/CPX). |
| VAR_TYPE | string | General variant type (ALL/SNP/INDEL). |
| ERR_TYPE | string | Classification error type (TP/FP/FN,PP). |
| MIN_QUAL | integer | Variants with at least this quality are retained. |
| (QUERY/TRUTH)_TOTAL | integer | Total count of query/truth variants. |
| (QUERY/TRUTH)_TP | integer | Total count of query/truth true positive variants. |
| (QUERY/TRUTH)_PP | integer | Total count of query/truth partial positive variants. |
| TRUTH_FN | integer | Total count of truth false negative variants. |
| QUERY_FP | integer | Total count of query false positive variants. |
| PREC | float | Variant calling precision [0,1]. |
| RECALL | float | Variant calling recall [0,1]. |
| F1_SCORE | float | Variant calling F1-score [0,1]. |
| (F1/ED/DE/ALN)_QSCORE | float | Quality score representation of various metrics [0,100] |
| EDIT_DIST | integer | Edit distance between query and truth after applying eligible variants (Calculated at design point B). |
| (SUB/INS/DEL)_ED | integer | Total edit distance for edits of a particular type. |
| DISTINCT_EDITS | integer | Distinct edits between query and truth after applying eligible variants (Calculated at design point B). |
| (SUB/INS/DEL)_DE | integer | Total distinct edits for edits of a particular type. |
| ALN_SCORE | integer | Total alignment score of query to truth after applying eligible variants (Calculated at design point B). |
| CONTIG | string | Contig which region is located on. |
| START | integer | 0-based index of leftmost position of current region. |
| STOP | integer | 0-based index of rightmost position of current region. |
| HAP | integer | Haplotype of current region (0/1). |
| SIZE | integer | Size of current region. |
| QUAL | float | Quality of current variant. |
| (SUPER)CLUSTER(S) | integer | Either 0-based index of or total (super)clusters in this region. |
| (QUERY/TRUTH)(1/2)_VARS | integer | Total variants on a particular haplotype within this region. |
| (ORIG/SWAP)_ED | integer | Total edit distance (minimum) of supercluster for both possible phasings. |
| PHASE | char | Character representing phasing. (=/X/?) for same, swap, unknown |
| PHASE_BLOCK | integer | 0-based index of current phase block. |
| POS | integer | 0-based index of position within contig. |
| REF/ALT | string | String of reference/alternate sequence at this position. |
| CREDIT | float | Fraction of partial credit this variant received. |
| LOCATION | string | Location of variant relative to BED regions (INSIDE/OUTSIDE/BORDER/OFF CTG)|
