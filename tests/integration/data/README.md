# Synthetic integration-test fixtures

These files back the single-pass precision/recall scenarios in
`tests/integration/test-integration.yml`. Each scenario uses a small synthetic reference so it
isolates one FN-dropping / credit-assignment behavior. (This directory also holds the large
chr20 fixtures for the separate genome-scale test; those are unrelated to the scenarios below.)

The whole `data/` directory is gitignored; these inputs are force-added so the suite runs from a
clean checkout. This directory holds inputs only: each scenario writes its outputs into the
temporary working directory `pytest-workflow` runs it in, which is also where the `files:` checks
look for them, so a test run leaves this directory untouched.

## Reference

- `synthetic.fasta` — a 400 bp contig `sc1` and a 100 bp contig `sc2`, which is a copy of the
  first 100 bp of `sc1` (with `synthetic.fasta.fai`, written by `samtools faidx`).
- `synthetic.bed` — `sc1  0  400`, the whole of `sc1`. Used by every single-contig scenario, which
  therefore never sees `sc2`.
- `synthetic_2ctg.bed` — `sc1  0  400` and `sc2  0  100`, both contigs in full.

Every variant `REF` allele below matches the reference at its 1-based position, and every variant
is homozygous (`1/1`, or `1|1` in `record_shapes` and `preserved_fields`) except in the
`contig_start_snp` and `contig_end_snp` scenarios, which are phased heterozygous (`1|0`), and the
one het-alt call that `record_shapes` and `preserved_fields` each add. Summary counts are therefore
per-haplotype (doubled, except in those heterozygous scenarios) and are listed as
`TRUTH_TP / QUERY_TP / TRUTH_FN / QUERY_FP` from
`*precision-recall-summary.tsv`.

## Scenarios

Each single-contig scenario has `<name>_truth.vcf` and `<name>_query.vcf`.

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

### contig_start_snp — a variant at reference position 0 (#177)
- truth and query are identical: SNPs 1 A>G and 50 T>C, both phased `1|0`.
- default `-ct`: both are TP on both sides. `POS` 1 is 0-based position 0, where no one-base left
  flank node can precede the variant in the alignment graph; the aligner's fixed origin used to
  land on the variant node itself, so the call scored FP and its truth counterpart was left
  unlabeled (`Unknown error type`). The `POS` 50 SNP is the mid-contig control. SNP 2/2/0/0
  (not doubled — these are heterozygous fixtures).

### contig_end_snp — a variant on a contig's final base (#189)
- truth and query are identical: SNPs 350 T>C and 400 T>C, both phased `1|0`.
- default `-ct`: both are TP on both sides. `POS` 400 is 0-based position 399, the final base of the
  400 bp `sc1`, where the reference window would otherwise reach past the contig end and leave the
  graph's trailing node holding no bases while its coordinate span claimed two. The `POS` 350 SNP is
  the mid-contig control, whose window fits inside the contig. SNP 2/2/0/0 (not doubled).

### record_shapes — one summary-VCF record per variant

- truth and query are identical, so every call is a gm TP and only the record shape is under test.
- both: hom SNP 200 A>G, hom CPX 210 `CAAGA`>`TT`, hom deletion 220 `CAACT`>`C`, and a het-alt
  SNP 250 A>C,G called `1|2`.
- default `-ct`: each homozygous call is one record whose per-haplotype fields carry two values,
  and the CPX is split into an INS and a DEL at parse time, each of which is one such record. The
  het-alt is the exception: parsing splits it into two entries with different ALTs, which nothing
  rejoins, so it stays two co-located records, each carrying one value for its ALT allele and `.`
  for its reference allele. SNP 4/4/0/0, INDEL 6/6/0/0.

### preserved_fields — source ID/QUAL/FILTER/INFO/FORMAT on the summary VCF

- query: hom SNP 200 A>G (`rs200`, QUAL 31, `PASS`), hom CPX 210 `CAAGA`>`TT` (`rs210`, QUAL 32,
  `LowConf`), and a het-alt SNP 250 A>C,G called `1|2` (`rs250`, QUAL 33, `PASS`). Its header
  declares a `Number=1` (`INFO/DP`, `FORMAT/SDP`), a fixed `Number=2` (`INFO/SB`, `FORMAT/SAC`),
  and a `Flag` (`INFO/SOMATIC`) field, all of which survive onto the output, plus `Number=A/R/G`
  fields (`INFO/AF`, `FORMAT/AD`, `FORMAT/PL`) that do not.
- truth: the same three calls, plus SNP 256 G>A (`tv256`) that the query misses, and an
  `INFO/TRUTHSET` field the query never declares.
- default `-ct`: the three shared calls are TP and owned by the query, so they carry the query's
  columns and the truth sample writes `.` for each appended FORMAT key; the false negative at 256
  is owned by the truth and the query sample pads instead. The `LowConf` `FILTER` is preserved
  verbatim on an evaluated record rather than rewritten to `PASS`.

### one_sided_contig — a contig called by only one callset (#166, #174)

This scenario uses `synthetic_2ctg.bed` and one pair of VCFs rather than a `_truth`/`_query` pair,
because the two directions are the same inputs with the query and truth arguments swapped:

- `one_sided_contig_both.vcf` — SNP `sc1` 200 A>G and SNP `sc2` 50 T>C.
- `one_sided_contig_sc1only.vcf` — SNP `sc1` 200 A>G only.

Both files declare `sc1` and `sc2` in their headers, so the contig is known to both callsets while
only one of them calls a variant on it. `sc1` is matched in both directions and stays TP.

- query-only (`both` as query): the `sc2` call has no truth counterpart. SNP 2/2/0/2.
- truth-only (`both` as truth): the `sc2` call is missed by the query. SNP 2/2/2/0.

Both directions used to segfault while superclustering, so these pin that the run completes and
that the one-sided contig's calls are classified and counted. Both also pin `*summary.vcf`: the
truth-only direction is where `sc2`'s false negatives used to be dropped, because that output
skipped any contig the query does not call on, and the query-only direction pins the mirror. With
no query call on `sc2` there is no phase block or phase to report there, so its truth record
reports `PB=0` and `BS=.`.
