
# Output `.vcf` Files
### summary.vcf 
Output summary VCF containing detailed information regarding each query and truth input variant. Please [note](https://github.com/TimD1/vcfdist/wiki/11-Implementation-Notes) that homozygous variants are split into two heterozygous variants prior to vcfdist evaluation.
| FORMAT Tag | Name | Type | Description |
| ---------- | ---- | ---- | ----------- |
| GT | GenoType | string | Genotype of the current variant. |
| BD | Benchmarking Decision | categorical | Benchmarking decision classification of the current variant (TP/FP/FN). |
| BC | Benchmarking Credit | float | Benchmarking credit allocated to the current variant in the range [0.0, 1.0]. |
| RD | Reference edit Distance | integer | Edit distance between the reference and truth sequences for this sync group. |
| QD | Query edit Distance | integer | Edit distance between the query and truth sequences for this sync group. |
| BK | Benchmarking Kategory | categorical | GA4GH category is assigned to 'gm' if credit == 1 (genotype match), 'lm' if credit > 0 (local match), or else '.' (no match). |
| QQ | Quality | integer | Variant PHRED-scaled quality score. |
| SC | SuperCluster | integer | 0-based index of supercluster in which this variant occurs, on this contig. |
| SG | Sync Group | integer | 0-based index of the current variant's sync group in this supercluster on this haplotype. Variants in the same sync group are considered complex, evaluated together, and assigned a single credit value. |
| PS | Phase Set | integer | Phase set identifier for the current variant. |
| PB | Phase Block | integer | 0-based index of current phase block in this contig. |
| BS | Block State | integer | Phasing state of the current phase block, which determines the truth-to-query haplotype mappings: 0 = T1Q1:T2Q2, 1 = T1Q2:T2Q1. |
| FE | Flip Error | integer | Whether the current supercluster's phasing is flipped (1) or not (0). |

For full false positive (FP) variants where credit is 0.0, RD and QD are missing (.) because the query variant is never aligned and there is no corresponding truth variant.

### orig-query.vcf, orig-truth.vcf
Original query and truth VCFs, as parsed by vcfdist.

### query.vcf, truth.vcf
Output query and truth VCFs, (optionally standardized by vcfdist using --realign).

# Output `.tsv` Files
`.tsv`, or tab-separated value files, are a simple text-based format used to store vcfdist's output.

## Variants
### query.tsv, truth.tsv
Reports detailed information regarding each variant.
| Column | Type | Description |
| ------ | ---- | ----------- |
| CONTIG | string | Name of the contig on which this variant is located. |
| POS | integer | 0-based index of genomic position where this variant occurs on this contig. Insertions occur before this base. |
| HAP | integer | Variant occurs on this haplotype (0 or 1 if diploid). |
| REF | string | Reference allele of this variant. |
| ALT | string | Alternate allele of this variant. |
| QUAL | float | Quality score of this variant. |
| TYPE | categorical | Type of variant (REF/SNP/INS/DEL/CPX). |
| ERRTYPE | categorical | Classification error type (TP/FP/FN). |
| CREDIT | float | Fraction of partial credit (in terms of correctness) that this variant received. --credit-threshold determines whether it is a true or false positive. |
| CLUSTER | integer | 0-based index of the cluster on the current haplotype of this contig containing this variant. |
| SUPERCLUSTER | integer | 0-based index of the supercluster on this contig containing this variant. |
| SYNC_GROUP | integer | 0-based index of credit grouping for variants on this haplotype and within this supercluster. |
| REF_DIST | integer | Edit distance between the reference and truth sequences for the current sync group. |
| QUERY_DIST | integer | Edit distance between the query and truth sequences for the current sync group. |
| LOCATION | categorical | Location of variant relative to BED regions (INSIDE/OUTSIDE/BORDER/OFF_CTG) |

For full false positive (FP) variants where credit is 0.0, REF_DIST and QUERY_DIST set to 0 because the query variant is never aligned and there is no corresponding truth variant.

## Clusters
### superclusters.tsv
Reports the size, location, and composition of each supercluster.
| Column | Type | Description |
| ------ | ---- | ----------- |
| CONTIG | string | Name of the contig on which this supercluster is located. |
| SUPERCLUSTER | integer | 0-based index of supercluster on this contig. |
| START | integer | 0-based leftmost genomic position of this supercluster on this contig, inclusive. |
| STOP | integer | 0-based rightmost genomic position of this supercluster on this contig, exclusive |
| SIZE | integer | Size of this supercluster, STOP-START. |
| QUERY1_VARS | integer | Number of variants from query haplotype 0 within this supercluster. |
| QUERY2_VARS | integer | Number of variants from query haplotype 1 within this supercluster. |
| TRUTH1_VARS | integer | Number of variants from truth haplotype 0 within this supercluster. |
| TRUTH2_VARS | integer | Number of variants from truth haplotype 1 within this supercluster. |
| ORIG_ED | integer | The minimum edit distance of this supercluster, assuming the original reported phasing is correct. |
| SWAP_ED | integer | The minimum edit distance of this supercluster, assuming the phasing of all variants is flipped. |
| PHASE_STATE | integer | Current phase block phasing state at this supercluster (0 = T1Q1:T2Q2, 1 = T1Q2: T2Q1). |
| SC_PHASE | categorical | Phase state of this supercluster, based on minimum edit distances. |
| PHASE_SET | integer | Phase set of this supercluster (`PS` tag input, genomic position of leftmost variant in phase set). |
| PHASE_BLOCK | integer | 0-based index of phase block within this contig. |
| FLIP_ERROR | integer | Whether this supercluster is considered a flip error (1) or not (0). |

## Precision and Recall
### precision-recall-summary.tsv
High-level precision/recall overview of SNP/INDEL/SV performance. For each category, there is one line for performance at the chosen minimum quality score, and one line for the quality score threshold that results in the best performance.
| Column | Type | Description |
| ------ | ---- | ----------- |
| VAR_TYPE | categorical | Variant type category. SNP for substitutions, INDEL for variants smaller than --sv-threshold, and SV for larger variants.  |
| MIN_QUAL | integer | Variants above this minimum quality score are used to calculate the performance metrics on this line. |
| TRUTH_TP | integer | Total measured true positive truth variants in this variant category. |
| QUERY_TP | integer | Total measured true positive query variants in this variant category. |
| TRUTH_FN | integer | Total measured false negative truth variants in this variant category. |
| QUERY_FP | integer | Total measured false positive query variants in this variant category. |
| PREC | float | Measured variant calling precision in the range [0.0, 1.0]: QUERY_TP / (QUERY_FP + QUERY_TP). |
| RECALL | float | Measured variant calling recall in the range [0.0, 1.0]: TRUTH_TP / (TRUTH_TP + TRUTH_FN). |
| F1_SCORE | float | Measured variant calling F1 score in the range [0.0, 1.0]: (2*PREC*RECALL) / (PREC + RECALL). |
| F1_QSCORE | float | PHRED-scaled measure of F1 score: -10*log10(1-F1_SCORE). |

### precision-recall.tsv
For each category (SNP/INDEL/SV), there is one line reporting the precison/recall performance at each possible quality score.
| Column | Type | Description |
| ------ | ---- | ----------- |
| VAR_TYPE | categorical | Variant type category. SNP for substitutions, INDEL for variants smaller than --sv-threshold, and SV for larger variants.  |
| MIN_QUAL | integer | Variants above this minimum quality score are used to calculate the performance metrics on this line. |
| PREC | float | Measured variant calling precision in the range [0.0, 1.0]: QUERY_TP / (QUERY_FP + QUERY_TP). |
| RECALL | float | Measured variant calling recall in the range [0.0, 1.0]: TRUTH_TP / (TRUTH_TP + TRUTH_FN). |
| F1_SCORE | float | Measured variant calling F1 score in the range [0.0, 1.0]: (2*PREC*RECALL) / (PREC + RECALL). |
| F1_QSCORE | float | PHRED-scaled measure of F1 score: -10*log10(1-F1_SCORE). |
| TRUTH_TOTAL | integer | Total truth variants in this variant category. |
| TRUTH_TP | integer | Total measured true positive truth variants in this variant category. |
| TRUTH_FN | integer | Total measured false negative truth variants in this variant category. |
| QUERY_TOTAL | integer | Total query variants in this variant category. |
| QUERY_TP | integer | Total measured true positive query variants in this variant category. |
| QUERY_FP | integer | Total measured false positive query variants in this variant category. |

## Phasing Analysis
### phase-blocks.tsv
Reports the size, location, and composition of each phase block.
| Column | Type | Description |
| ------ | ---- | ----------- |
| CONTIG | string | Name of the contig which contains the current phase block. |
| PHASE_BLOCK | integer | 0-based index of this phase block within the contig. |
| START | integer | 0-based index of the leftmost genomic position within current phase block, inclusive. |
| STOP | integer | 0-based index of the rightmost genomic position within current phase block, exclusive. |
| SIZE | integer | Size of the current phase block. |
| SUPERCLUSTERS | integer | Total superclusters comprising this phase block. |
| FLIP_ERRORS | integer | Total flip errors that occur within this phase block. |
| SWITCH_ERRORS | integer | Total switch errors that occur within this phase block. |

### phasing-summary.tsv
High-level summary of phasing performance.
| Column | Type | Description |
| ------ | ---- | ----------- |
| PHASE_BLOCKS | integer | Total phase blocks, across all evaluated contigs. |
| SWITCH_ERRORS | integer | Total switch errors, across all evaluated contigs. |
| FLIP_ERRORS | integer | Total flip errors, across all evaluated contigs. |
| NG_50 | integer | NG50 (break regions on new phase block), across all evaluated contigs. |
| SWITCH_NGC50 | integer | Switch NGC50 (break regions on new phase block or switch error), across all evaluated contigs. |
| SWITCHFLIP_NGC50 | integer | Switchflip NGC50 (break regions on new phase block, switch error, or flip error), across all evaluated contigs. |

### switchflips.tsv
Detailed breakdown of all switch and flip errors. 
| Column | Type | Description |
| ------ | ---- | ----------- |
| CONTIG | string | Contig on which phasing error occurred. |
| START | integer | 0-based index of leftmost genomic position where phasing error could have occurred. |
| STOP | integer | 0-based index of rightmost genomic position where phasing error could have occurred. |
| SWITCH_TYPE | categorical | Type of phasing error. |
| SUPERCLUSTER | integer | 0-based index of supercluster on this contig where the phasing error occurred. |
| PHASE_BLOCK | integer | 0-based index of phase block on this contig which the phasing error occurred within. |

## Alignment Distance
### distance-summary.tsv
High-level alignment distance overview of SNP/INDEL performance. For each category (SNP/INDEL/ALL), there is one line for performance at the minimum, best, and maximum quality score threshold.
| Column | Type | Description |
| ------ | ---- | ----------- |
| VAR_TYPE | categorical | Variant category considered on this line (ALL/SNP/INS/DEL/INDEL). |
| MIN_QUAL | integer | Minimum variant quality considered. |
| EDIT_DISTANCE | integer | Total edit distance from truth sequence after query variants above MIN_QUAL are applied. |
| DISTINCT_EDITS | integer | Total number of edits from truth sequence after query variants above MIN_QUAL are applied. |
| ED_QSCORE | float | PHRED-scaled measure of query sequence quality in terms of edit distance: -10*log10(EDIT_DISTANCE / ORIG_EDIT_DISTANCE), where ORIG_EDIT_DISTANCE is the query edit distance from the truth sequence when no query variants are applied. |
| DE_QSCORE | float | PHRED-scaled measure of query sequence quality in terms of distinct edits: -10*log10(DISTINCT_EDITS / TRUTH_TOTAL), where TRUTH_TOTAL is total number of truth variants. |
| ALN_QSCORE | float | PHRED-scaled measure of query sequence quality in terms of alignment score: -10*log10(ALN_SCORE / ORIG_ALN_SCORE), where ORIG_ALN_SCORE is the query alignment score (to the truth sequence) when no query variants are applied. |

### distance.tsv
There is one line reporting the alignment distance performance at each possible quality score.
| Column | Type | Description |
| ------ | ---- | ----------- |
| MIN_QUAL | integer | Minimum variant quality considered. |
| SUB_DE | integer | Total number of substitution edits during alignment. |
| INS_DE | integer | Total number of insertion edits during alignment. |
| DEL_DE | integer | Total number of deletion edits during alignment. |
| SUB_ED | integer | Total number of bases substituted during alignment. |
| INS_ED | integer | Total number of bases inserted during alignment. |
| DEL_ED | integer | Total number of bases deleted during alignment. |
| DISTINCT_EDITS | integer | Total number of edits from truth sequence after query variants above MIN_QUAL are applied. |
| EDIT_DIST | integer | Total edit distance from truth sequence after query variants above MIN_QUAL are applied. |
| ALN_SCORE | integer | Smith-Waterman alignment score of query sequence to truth sequence, after applying query variants above MIN_QUAL.
| ALN_QSCORE | float | PHRED-scaled measure of query sequence quality in terms of alignment score: -10*log10(ALN_SCORE / ORIG_ALN_SCORE), where ORIG_ALN_SCORE is the query alignment score (to the truth sequence) when no query variants are applied. |

### edits.tsv
This file reports for each edit (where called query sequence differs from truth sequence) the contig, pos, hap, len, supercluster, and the quality range for which these edits occur.
| Column | Type | Description |
| ------ | ---- | ----------- |
| CONTIG | string | Name of the contig on which this edit occurs. |
| START | integer | 0-based genomic position on this contig where this edit is located. |
| HAP | integer | Haplotype on which this edit occurs. | 
| TYPE | categorical | Edit category type (SNP/INS/DEL). |
| SIZE | integer | Size of the edit, in terms of reference or query bases affected. |
| SUPERCLUSTER | integer | 0-based index of the supercluster containing the edit, on this contig. |
| MIN_QUAL | integer | The edit is present at and above this quality score. |
| MAX_QUAL | integer | The edit is present below this quality score. |

# Output `.txt` Files
### parameters.txt
This file stores all internal and command-line parameters used by vcfdist (for reproducibility and debugging purposes). Each line stores a single parameter in the format `parameter = value`.
