# vcfdist v3.0.0 — Design Document

|  |  |
| :-- | :-- |
| **Version:** | 0 (Draft) |
| **Authored By:** | Tim Dunn |
| **Status:** | Draft — living document; expands as deliverables land |
| **Companion SOW:** | `D0_vcfdist-v3-SOW.md` |
| **Child design docs:** | [`D2_vcfdist-v3-unit-tests.md`](./D2_vcfdist-v3-unit-tests.md) (§5.2), [`D3_retain-info-format-fields.md`](./D3_retain-info-format-fields.md) (§5.3), [`D4_ga4gh-compatibility.md`](./D4_ga4gh-compatibility.md) (§5.4), [`D5_multi-bed-stratification.md`](./D5_multi-bed-stratification.md) (§5.5), [`D6_vcfdist-v3-unphased-eval.md`](./D6_vcfdist-v3-unphased-eval.md) (§5.6), [`D7_vcfdist-v3-benchmarking-design.md`](./D7_vcfdist-v3-benchmarking-design.md) (§5.7), [`D8_vcfdist-v3-release.md`](./D8_vcfdist-v3-release.md) (§5.8), [`D9_docs-and-tutorial.md`](./D9_docs-and-tutorial.md) (§5.9), [`D10_vcfdist-v3-multiqc-module.md`](./D10_vcfdist-v3-multiqc-module.md) (§5.10), [`D11_pipeline-ecosystem-integration.md`](./D11_pipeline-ecosystem-integration.md) (§5.11) |

---

## 1. Purpose & Scope

This is the design document called for by **Deliverable #1** of the v3.0.0 SOW. It serves two audiences:

1. **Reference for the v3 architecture** — what changed from v2 and why, so the codebase is legible to future contributors and to the reviewer.
2. **Implementation sketch for every SOW deliverable** — a per-deliverable plan (§5) that turns the SOW's *what* into a *how*, surfacing the concrete files, interfaces, and open decisions before code is written.

Scope boundaries follow the SOW verbatim: v3.0.0 targets germline small-variant + SV evaluation with phasing analysis, stratification, unphased support, and ecosystem integration. Faster alignment, larger/symbolic SVs, and somatic evaluation are **Potential Future Work** (SOW #12–#14); large SVs in general (>1 kb) beyond the #13 target and general runtime optimization beyond avoiding explosion are explicitly deferred to v3.1.

---

## 2. Background: what v2 does

vcfdist is a distance-based germline variant-calling evaluation tool that, unlike `hap.py`/`vcfeval`, **requires locally-phased** truth and query VCFs and, in return, reports **flip and switch phasing errors** and produces more accurate SNP/INDEL/SV precision-recall curves — particularly around complex variants — by aligning candidate haplotypes rather than matching variant representations.

The v2 pipeline, per `src/main.cpp` runs roughly:

```
  → parse VCFs
  → normalize VCFs (v2 only)
  → cluster (per-haplotype WFA clustering)
  → recluster (cluster across haplotypes)
  → supercluster (merge overlapping truth/query clusters)
  → precision/recall (align candidate haplotypes, assign credit)
  → edit distance (v2 only)
  → phasing (phase blocks, switch/flip errors)
  → write outputs
```

Core source modules: `bed` (region I/O), `cluster` (clustering), `dist` (precision/recall + credit), `edit` (edit-distance metrics), `phase` (phasing + `summary.vcf`), `variant` (VCF model), `print`, `globals`, `fasta.h`, `defs.h`.

### v2 output files

| File | Contents |
| :-- | :-- |
| `precision-recall-summary.tsv` | Headline P/R/F1 per variant type at the chosen quality threshold |
| `precision-recall.tsv` | Full P/R curve: metrics across all quality thresholds |
| `phasing-summary.tsv` | Switch/flip error counts and rates |
| `switchflips.tsv`, `phase-blocks.tsv`, `superclusters.tsv` | Per-event / per-block / per-supercluster detail |
| `distance-summary.tsv`, `distance.tsv`, `edits.tsv` | Edit-distance metrics (**removed in v3**) |
| `query.tsv`, `truth.tsv` | Per-variant classification tables |
| `summary.vcf` | Two-sample `TRUTH`/`QUERY` GA4GH-style annotated VCF (`BD`/`BK`/`BC`/…) |
| `query.vcf`, `truth.vcf`, `orig-{query,truth}.vcf` | Normalized / original input round-trips |

---

## 3. v3 Architecture

> This section describes the v3 engine **as it exists today**, recorded as settled context. The
> capabilities not yet built are summarised per deliverable in §5.

### 3.1 Pipeline overview

v3 keeps the cluster → supercluster → precision/recall → phasing → write shape. The edit-distance stage has been removed; the alignment core moved to a graph/WFA formulation; and phasing decisions moved from per-supercluster to per-variant.

### 3.2 Major v2 → v3 changes

- **Per-variant phasing (was per-supercluster).** v2 assigned a single phase decision to each supercluster; v3 decides phase per variant, giving finer-grained flip/switch attribution. *Known constraint:* phase groups still require adjacent variants (§4).
- **Removal of the edit-distance code path.** The `edit` module and its `distance*.tsv`/`edits.tsv` outputs are gone. Edit distance was a secondary metric that overlapped with the credit-based precision/recall and added maintenance surface. Callers of those outputs move to the P/R metrics, which the v3.0.0 release notes cover (§5.8).
- **Removal of the VCF-normalization code path.** v3 now assumes normalized input (handled upstream) rather than carrying vcfdist's own normalization; reduces representation-handling divergence.
- **Graph-based / WFA alignment.** The candidate-haplotype alignment moved to a graph/WFA formulation — the foundation the future-work SV scaling (SOW #12–#13) uses.
- **Genotype and allele-count error summaries.** New per-variant summaries beyond TP/FP/FN were added (the `GE`, genotype-error, and allele-count fields already visible in `summary.vcf`).
- **Bounded superclusters.** Supercluster size is now bounded to prevent the runtime/RAM explosion seen on dense/large regions (clusters >20 kb were already problematic).

### 3.3 Assumptions carried into v3

- Sequence-resolved REF/ALT (no symbolic/BND alleles until SOW #13).
- Diploid genotypes with local phasing available (until SOW #6 relaxes the phasing requirement for the query callset).
- Single confident/evaluation BED (until SOW #5 adds multi-BED stratification).

---

## 4. Known Limitations

Limitations of the v3 engine as it stands today (§3). The *Addressed by* column names the deliverable
that closes each one; the change itself is summarised in §5. Entries marked future work or
documented-only are not slated for v3.0.0.

| Limitation | Addressed by |
| :-- | :-- |
| Requires locally-phased query input | #6 Unphased evaluation ([`D6`](./D6_vcfdist-v3-unphased-eval.md)) |
| Requires locally-phased truth input | Documented; not fixed in v3 ([`D6`](./D6_vcfdist-v3-unphased-eval.md) §9.1) |
| Single BED only | #5 Multi-BED stratification ([`D5`](./D5_multi-bed-stratification.md)) |
| Drops annotations / unevaluated variants | #3 Retain fields + mark `N` ([`D3`](./D3_retain-info-format-fields.md)) |
| Overlapping variants dropped | #6 for the query callset ([`D6`](./D6_vcfdist-v3-unphased-eval.md) §4.1.3); truth side documented only (D6 §9.3) |
| Phase groups require adjacent variants | Documented; not fixed in v3 |
| 10 kb practical SV limit | #12 Faster alignment (future work) |
| Sequence-resolved alleles only (no `<DEL>`/BND) | #13 Large & symbolic SV (future work) |
| Diploid/germline only | #14 Somatic (future work) |

---

## 5. Per-Deliverable Changes

One subheading per SOW deliverable, in SOW order, summarising the **changes, additions, removals, and
fixes** each one makes. Each child design doc is the detailed specification for its deliverable.

### 5.1 Deliverable #1 — Design Document

This document: §2–§4 as the architecture reference, §5 as the per-deliverable summary. Each deliverable
updates its own subsection to describe what it implemented.

### 5.2 Deliverable #2 — Comprehensive Unit Tests

Detailed design: [`D2_vcfdist-v3-unit-tests.md`](./D2_vcfdist-v3-unit-tests.md), which enumerates the
suite test-by-test — name, input, expected result, and the branch, boundary, or bug each case exists to
catch — for every function in every `src/*.cpp`.

- **Additions:** one test file `tests/unit/src/test_<name>.cpp` per `src/<name>.cpp` (`bed`, `variant`,
  `cluster`, `dist`, `phase`, `globals`, `print`, `timer`); shared in-memory builders in
  `test_helpers.{h,cpp}` plus a `write_tmp_vcf` emitter; a `GlobalsGuard` RAII fixture; two tiny on-disk
  fixtures (`tiny.fasta`, `tiny.bed`); and a `bcftools` VCF-conformance integration test over
  `summary.vcf`.
- **Changes:** GoogleTest for unit tests launched by `pytest`, with `pytest-workflow` retaining the
  integration tests; the `Makefile` `OBJS` list grows so every test object links into the existing
  `test_vcfdist` binary. `.github/workflows/test.yml` gains `bcftools`.
- **Refactor:** the P/R/F1 arithmetic, currently inlined twice in `print.cpp`, is extracted into a pure
  `compute_pr_f1(query_tp, query_fp, truth_tp, truth_fn)` helper in `print.{h,cpp}` (#94) — the single
  source of truth for the per-stratum rows of #5 and the MultiQC module of #10.
- **`main.cpp`** gets no unit test; its string tables move to `strings.cpp` and are tested there.

### 5.3 Deliverable #3 — Retain Original INFO/FORMAT Fields

Detailed design: [`D3_retain-info-format-fields.md`](./D3_retain-info-format-fields.md). This deliverable
establishes the `summary.vcf` record layout that #4 annotates, #5 stratifies, and #6 extends.

- **Changes:** `summary.vcf` emits **one record per variant** rather than one per haplotype copy — a
  homozygous variant becomes a single `GT=1|1` record, while het-alt (`1|2`) records stay split as two
  co-located biallelic records. Per-haplotype fields (`BC`, `RD`, `QD`, `SG`, and `BD`/`BK` until #4
  collapses them) become `Number=.` lists in emitted-`GT` order. `ID`, `QUAL`, `FILTER`, `INFO`, and
  sample `FORMAT` come from the record owning each output line instead of hardcoded placeholders.
  `Number=A`/`R`/`G` fields are subset to the allele each line represents, with their header cardinality
  rewritten to match.
- **Additions:** per-variant `rec_idxs`, `alt_idxs`, and `ploidies` vectors on `ctgVariants`; a sideline
  container holding retained-but-unevaluated variants outside `ctgVariants`; `BD=N` for those variants;
  six `VCFDIST_*` `FILTER` tags naming why each was not evaluated.
- **Removals:** the per-contig ploidy value, its chrX warning suppression, and the non-spec `ploidy=`
  attribute on `##contig` lines — each record's `GT` now carries its own ploidy.
- **Out of scope:** overlapping variants and unphased heterozygous genotypes, both handled by #6.

### 5.4 Deliverable #4 — GA4GH Benchmarking VCF Compatibility

Detailed design: [`D4_ga4gh-compatibility.md`](./D4_ga4gh-compatibility.md). Brings `summary.vcf` into
conformance with the GA4GH benchmarking (intermediate) VCF contract, so a consumer of that format can
read vcfdist's per-variant decisions. vcfdist's own metrics remain its native TSVs.

- **Changes:** **counting becomes per site** — one unit per variant per callset, matching `hap.py` and
  `vcfeval`, rather than one per haplotype copy. This changes every count and every derived rate. `BD`
  and `BK` collapse from #3's `Number=.` lists to a single value per site, assigned by an explicit
  match-tier ladder. `FORMAT/SC` becomes `INFO/BS`, the benchmarking superlocus ID; `FORMAT/BS`
  (phase-block state) becomes `FORMAT/PBS`.
- **Additions:** `am` (allele match) and `pm` (phased match) in the `BK` vocabulary, the latter
  implementing **GA4GH Comparison Method #4** — the first implementation of that method, since neither
  incumbent reports phasing; `--stringency {lm,am,gm,pm}` (default `gm`) and `--require-phase` /
  `--no-require-phase`; and `FP_GT`, `FP_AL`, `QUERY_UNK`, `FRAC_NA`, TiTv and het:hom columns in the
  precision-recall summaries.
- **Removals:** the `AC_ERR_2_TO_1` truth-counter correction in `print.cpp`, which existed only because
  per-haplotype counting mis-attributed allele-count errors.
- **Not emitted:** `BVT`, `BLT`, `BI`, `Regions`, and `VTC`, all of which a GA4GH consumer derives itself
  and overwrites.
- **Caveat:** `BD` is a hard integer label, so vcfdist's fractional `BC` credit cannot survive a GA4GH
  export. Exported counts reflect vcfdist's thresholded decisions, not its partial-credit result.
- **Out of scope:** `som.py`, which makes its own TP/FP/FN decisions rather than consuming a pre-decided
  VCF; and stratification, which is supplied externally rather than embedded in the VCF.

### 5.5 Deliverable #5 — Multi-BED Stratification Support

Detailed design: [`D5_multi-bed-stratification.md`](./D5_multi-bed-stratification.md). N named
stratification region sets that **label** variants without removing any; `-b` remains the sole filter.

- **Additions:** `-st`/`--stratification <manifest.tsv>`, a hap.py-compatible two-column manifest whose
  relative paths resolve against its own directory; a flat `strata_bits` bitset on `ctgVariants`; a
  `STRATUM` leading column on `precision-recall.tsv`, `precision-recall-summary.tsv`, and
  `genotype-errors.tsv`, `*` first; a set-valued `STRATA` column on `query.tsv`/`truth.tsv`; a
  `bedData::merge()` sorting and merging third-party region sets at load; a pure `classify()` extracted
  from `contains()` so the monotonic cursor sweep and `contains()` share one decision tree; and a `WARN`
  naming any stratum with zero contig overlap.
- **Changes:** `-s` becomes `-sc`/`--max-supercluster-size`, with **no `-s` alias retained**, so the
  short form can mean stratification. `-b` and stratum BEDs accept gzip and bgzip via htslib
  `hts_open`/`hts_getline`. `contains()` no longer consults the `g.bed_exists` global. Best-F1 selection
  becomes per stratum. `phasing-summary.tsv` splits into `phasing-variants-summary.tsv`
  (variant-attributed, stratified) and `phasing-blocks-summary.tsv` (genome-level contiguity), and their
  rate columns become bare-float fractions with a terminating newline.
- **Fixes:** a coordinate-unsorted input VCF is now a hard `ERROR` rather than being silently degraded
  record-by-record by the overlap filter.
- **Out of scope:** stratifying `summary.vcf`, `switchflips.tsv`, or `phase-blocks.tsv`; per-stratum
  NG50/NGC50; and region-set algebra over strata.

### 5.6 Deliverable #6 — Unphased Variant Evaluation

Detailed design: [`D6_vcfdist-v3-unphased-eval.md`](./D6_vcfdist-v3-unphased-eval.md). Removes the hard
local-phasing requirement for the **query** VCF. No haplotype-pairing search is needed: the v3 alignment
core is already phase-blind on the query side — `Graph` builds the query as a DAG containing every
unevaluated query variant with no haplotype filter, and query phase is read nowhere in the alignment,
credit, or P/R tally.

- **Changes:** the haplotype merge moves **ahead of** clustering, so clustering runs per callset on the
  merged variant list and no longer depends on which haplotype an unphased heterozygote was assigned to;
  `wf_swg_cluster` operates on a single merged `ctgVariants` and the thread fan-out drops from
  `HAPS × contigs` to `contigs` per callset. Unphased heterozygotes are forced to `PHASE_NONE`, which the
  phasing DP already treats as zero-cost in both orientations, so they are transparent rather than
  fabricating switch and flip observations. The switch/flip denominator becomes `ASSESSED_PAIRS` —
  adjacent phaseable-heterozygote pairs per phase block — rather than every merged query variant.
  Unphased genotypes are written with `/` rather than `|`.
- **Additions:** an `is_phased` field on `ctgVariants`, distinct from `phase_sets`; `ASSESSED_PAIRS`,
  `PHASED_HET_VARIANTS`, `TOTAL_HET_VARIANTS`, and `PHASED_HET_FRACTION` in
  `phasing-variants-summary.tsv`; `PHASED_REGION_FRACTION` in `phasing-blocks-summary.tsv`.
- **Removals:** `-i`/`--max-iterations` — cluster iteration is fixed at 1 and becomes a correctness
  precondition; the cluster-merging half of `load_and_merge_callset_vars_across_haps`, and with it the
  cross-haplotype reach heuristic; the overlapping-variant filter for the `QUERY` callset, whose
  overlaps the query graph already represents as mutually exclusive alternative paths.
- **Fixes:** the stale contig index in `fix_phase_set_tags`'s no-phase-set early exit; the
  genome-spanning phase block an unphased query reports as near-perfect phasing; and an out-of-bounds
  `phase_sets` read in the phasing DP's backward pass.
- **Out of scope:** unphased **truth** VCFs and overlapping **truth** variants — the truth side of the
  alignment graph is a linear per-haplotype chain by construction — and a graph-based cluster-reach
  redesign.

### 5.7 Deliverable #7 — Benchmarking v3 vs. v2 / vcfeval / hap.py

Detailed design: [`D7_vcfdist-v3-benchmarking-design.md`](./D7_vcfdist-v3-benchmarking-design.md).

- **Additions:** a reproducible pixi + Snakemake harness running vcfdist v3, vcfdist v2.6.4, `vcfeval`,
  and `hap.py` (`xcmp`) over the T2T-Q100 HG002 truth set with the HPRC, PAV, and GIAB-TR query
  call-sets; precision-recall curves with all tools on one plot; runtime and peak-RAM figures;
  per-variant concordance matrices cross-tabulating each tool's TP/FP/FN decision against vcfdist's; and
  switch/flip validation against `whatshap compare` as an independent phasing oracle.
- **Changes:** all published figures are generated on one on-demand AWS `x86-64` instance, so `hap.py`
  runs natively and every tool's runtime is measured on the same ISA. `--credit-threshold` is pinned
  identically for v2 and v3, since their defaults differ.
- **Removals:** Truvari and its WFA/MAFFT/POA refine variants.
- **Reporting:** the counting-convention change from #4 and the supercluster-boundary shift from #6 are
  each measured and reported separately, so neither is presented as an accuracy result.

### 5.8 Deliverable #8 — v3.0.0 Release & Launch Announcement

Detailed design: [`D8_vcfdist-v3-release.md`](./D8_vcfdist-v3-release.md).

- **Additions:** `CHANGELOG.md` at the repository root, with v3.0.0 as the first entry; a
  `3.0.0-rc1` pre-release with a two-week feedback window and a stated output-format freeze; a
  `cppcheck`-based dead-code sweep; and a check that `src/globals.h`'s `VERSION` and `README.md` cannot
  diverge.
- **Changes:** the release procedure moves to `docs/v3.0.0/13-Release-Instructions.md`, with the Docker
  image built against a pinned tag rather than whatever `master` HEAD was, and its `ubuntu:20.04` base
  and HTSlib 1.17 bumped. Default parameter values are chosen from the #7 benchmarking evidence, each
  recorded with its basis or explicitly noted as conventional and unmeasured.
- **Announcement:** a Fulcrum blog post drawing on the #7 figures, concise GitHub release notes, and
  LinkedIn/BlueSky posts. The post states the per-site counting change before showing any figure, since
  v3's numbers are not comparable to v2's or to the prior publications' without a convention-matched run.

### 5.9 Deliverable #9 — Documentation & sandbox.bio Tutorial

Detailed design: [`D9_docs-and-tutorial.md`](./D9_docs-and-tutorial.md).

- **Removals:** the `04-VCF-Normalization` and `08-Alignment-Distance` wiki pages, whose code paths no
  longer exist.
- **Additions:** `Unphased-Evaluation` and `GA4GH-Compatibility` wiki pages; an interactive sandbox.bio
  tutorial written by the platform maintainers from our outline; and a small self-contained demo dataset
  — a synthetic contig plus phased truth/query VCFs and a BED, built by a fixed-seed script — that
  exercises exact TPs, a genotype error, a partial-credit complex variant, an FP, an FN, and a flip plus
  a switch.
- **Changes:** `02-Parameters-and-Usage` regenerated from v3 `--help`; `09-Outputs` rewritten against
  the v3 file set; `01-Overview`, `05-Variant-Clustering`, `06-Precision-and-Recall`, and
  `07-Phasing-Analysis` updated for the v3 engine; `README.md`'s version string, demo command, and
  expected output regenerated. `demo/pr_plot.py` filters to `STRATUM == "*"`.
- **WebAssembly:** in-browser execution needs a biowasm build of vcfdist, which the sandbox.bio
  maintainers own.

### 5.10 Deliverable #10 — MultiQC Module

Detailed design: [`D10_vcfdist-v3-multiqc-module.md`](./D10_vcfdist-v3-multiqc-module.md). A native
module contributed upstream to MultiQC, with a standalone pip-installable plugin as the fallback if
upstream review stalls. Both build the same module class; the plugin adds packaging and hooks.

- **Additions:** `multiqc/modules/vcfdist/`, content-anchored search patterns, and a `module_order`
  entry. The report carries headline SNP/INDEL F1 and recall plus switch+flip counts in General
  Statistics; a precision-recall bargraph and curve with a variant-type switcher; a switch/flip
  bargraph, which is the metric no competing tool reports; a merged phasing table; per-variant-type
  stratified F1 heatmaps; and a weakest-strata table. Parsed values are written to `multiqc_data/`,
  including the full sample × stratum × variant-type matrix.
- **Reads:** `precision-recall-summary.tsv`, `precision-recall.tsv`,
  `phasing-variants-summary.tsv`, and `phasing-blocks-summary.tsv`, routing rows by `STRATUM` so
  genome-wide and per-stratum figures never mix.
- **Out of scope:** the per-variant detail files, per-stratum precision-recall curves, and per-stratum
  phasing errors.

### 5.11 Deliverable #11 — nf-core & snakemake-wrappers Integration

Detailed design:
[`D11_pipeline-ecosystem-integration.md`](./D11_pipeline-ecosystem-integration.md). Three upstream
contributions; no vcfdist source changes.

- **Additions:** an `nf-core/modules` `VCFDIST` module — `main.nf`, `meta.yml`, `environment.yml`, and
  `nf-test` cases — with one named output per file vcfdist writes; vcfdist as a first-class `method` in
  `nf-core/variantbenchmarking`, including a `get_vcfdist_results()` reader in `bin/merge_reports.py`
  and a datavzrd template; and a `bio/vcfdist` snakemake-wrappers entry with `wrapper.py`,
  `environment.yaml`, `meta.yaml`, a lockfile, and example rules that double as the rendered
  documentation.
- **Contract:** the three positional arguments come first; the reference FASTA must be uncompressed;
  VCF and FASTA indexes are not read and are not declared as inputs; truth input must be locally phased;
  `--max-threads` and `--max-ram` are derived from the executor's allocation; and
  `-n`/`--no-output-files` is unsupported, since it clears the flag the file writers check.
- **Out of scope in the first round:** the pipeline's per-decision `tagged_variants` channel, phasing
  metrics, stratification wiring, and the structural-variant benchmarking path.

### 5.12 Deliverable #12 (future work) — Faster Alignment Algorithm

Speed up the existing sequence-resolved WFA / graph-based alignment so denser and larger clusters are
tractable without runtime or RAM explosion — via truth-graph and query-reference node merging during
WFA, plus mixed exact/heuristic clustering. Raises the practical SV size limit toward 100 kb within the
current sequence-resolved approach, without changing what kinds of variants can be represented.

### 5.13 Deliverable #13 (future work) — Large & Symbolic (BND) SV Support

Handle large SVs and a distinct representation path for symbolic and breakend alleles (`<DEL>`,
`<DUP>`, `<INV>`, BND mate pairs) that carry no resolved REF/ALT sequence, via a new algorithm — these
cannot be evaluated by sequence alignment.

### 5.14 Deliverable #14 (future work) — Somatic Variant Support

Evaluate somatic calls, where the diploid-genotype and local-phasing assumptions underpinning germline
evaluation do not hold: variant allele fractions, sub-clonal calls, and tumor/normal comparison rather
than haplotype-resolved genotypes. Validated against SEQC2. Requires its own design pass.

---

## 6. Open Decisions

**Sample identity in outputs (#3/#10).** How each vcfdist run declares its sample name for MultiQC
merging. The MultiQC module derives it by stripping the known suffix from `-p/--prefix`, which works but
is a convention rather than a declaration; an explicit CLI flag or embedded field would be firmer.

---

## 7. References

- Companion SOW: `D0_vcfdist-v3-SOW.md`
- Deliverable #2 unit-test design: [`D2_vcfdist-v3-unit-tests.md`](./D2_vcfdist-v3-unit-tests.md) (expands §5.2)
- Deliverable #3 field-retention / record-shape design: [`D3_retain-info-format-fields.md`](./D3_retain-info-format-fields.md) (expands §5.3)
- Deliverable #4 GA4GH-compatibility design: [`D4_ga4gh-compatibility.md`](./D4_ga4gh-compatibility.md) (expands §5.4)
- Deliverable #5 stratification design: [`D5_multi-bed-stratification.md`](./D5_multi-bed-stratification.md) (expands §5.5)
- Deliverable #6 unphased-evaluation design: [`D6_vcfdist-v3-unphased-eval.md`](./D6_vcfdist-v3-unphased-eval.md) (expands §5.6)
- Krusche et al. 2019, *Best practices for benchmarking germline small-variant calls in human genomes*, Nat. Biotechnol. (GA4GH benchmarking standard)
- Dunn & Narayanasamy 2023, *vcfdist: Accurately benchmarking phased small variant calls in human genomes*, Nat. Commun.
- MultiQC documentation — writing modules & plugins (pin to the targeted release)
