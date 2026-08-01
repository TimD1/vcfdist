# Synthetic integration-test fixtures

These files back the single-pass precision/recall scenarios in
`tests/integration/test-integration.yml`. Each scenario uses a small synthetic reference so it
isolates one FN-dropping / credit-assignment behavior. (This directory also holds the large
chr20 fixtures for the separate genome-scale test; those are unrelated to the scenarios below.)

The whole `data/` directory is gitignored; these inputs are force-added so the suite runs from a
clean checkout. Generated per-run outputs (`*_out_*`, `*_ct07_*`, `*_ct1_*`, …) stay ignored.

## Reference

- `synthetic.fasta` — one 400 bp contig, `sc1` (with `synthetic.fasta.fai`).
- `synthetic.bed` — `sc1  0  400`, the whole contig.

Every variant `REF` allele below matches `sc1` at its 1-based position, and every variant is
homozygous (`1/1`). Summary counts are therefore per-haplotype (doubled) and are listed as
`TRUTH_TP / QUERY_TP / TRUTH_FN / QUERY_FP` from `*precision-recall-summary.tsv`.

## Scenarios

Each scenario has `<name>_truth.vcf` and `<name>_query.vcf`.

### swallowed_snps — a large FN SV must not swallow nearby SNPs
- truth: SNPs 140 T>C, 146 G>A, 256 G>A, and a 100 bp deletion at 150.
- query: the three SNPs only.
- default `-ct`: the SNPs are recovered as TP and the deletion is FN, without the SNPs being
  folded into the deletion's sync group. SNP 6/6/0/0, SV 0/0/2/0.

### partial_match_del — a partial match across the credit threshold
- truth: 10 bp deletion at 150 (`AACAAGACGTC` > `A`).
- query: 8 bp deletion at 150 (`AACAAGACG` > `A`); query_dist 2, ref_dist 10, credit 0.8.
- `-ct 0.7`: clears the bar, TP@0.8. INDEL 2/2/0/0.
- `-ct 1.0`: only an exact match stays on-path, so truth is FN and query is FP. INDEL 0/0/2/2.

### colocated_spurious_snp — independent FN and FP, no cross-crediting
- truth: SNP 200 A>G (missed by query).
- query: an unrelated SNP 200 A>C at the same position.
- default `-ct`: truth is FN and query is FP independently (the r=0 case). SNP 0/0/2/2.

### entangled_overlap_del — entangled-overlap residual (exploratory)
- truth: 10 bp deletion at 150 (`AACAAGACGTC` > `A`).
- query: 10 bp deletion at 152 (`CAAGACGTCCT` > `C`); different start, 8 bp of overlap.
- `-ct 0.7`: both calls land in one sync group and both score TP@0.8 despite not sharing a
  start position — a nonzero-r cross-crediting residual. INDEL 2/2/0/0. This is characterized,
  not asserted to be correct behavior.

### adjacent_fn_no_overcredit — consecutive-bypass excision regression
- truth: SNP 300 A>C (reproduced), plus adjacent SNPs 301 G>T and 302 T>A (both missed).
- query: SNP 300 A>C only.
- default `-ct`: all three share one sync group, so both FN bypass spans must be excised from
  the TP's reference edit distance. The fix yields `RD=1` for the TP; leaving the inner bypass
  un-excised inflates it to `RD=2`. The call stays TP (credit 1.0) either way, so `RD` is the
  field the test pins. SNP 2/2/4/0.
