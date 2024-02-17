### Homozygous variant counting
Prior to evaluation, all homozygous variants with genotype `1|1` are split into two heterozygous variants `0|1` and `1|0` (with identical reference and alternate alleles). We do this because variants on separate haplotypes may be realigned to a new representation, and to deal with cases where only one of the two calls is correct (or part of a complex variant). Thus, a fully-correct homozygous variant call will be counted as two true positives, whereas most other evaluation tools will count one true positive.

### Limitations
The current version of vcfdist is not designed to support:
 - somatic variants
 - unphased variants
 - overlapping variants (on the same haplotype)
 - polyploid contigs
