# vcfdist v3.0.0 — MultiQC Module (Design)

|  |  |
| :-- | :-- |
| **Version:** | 1 (Draft) |
| **Authored By:** | Tim Dunn |
| **Status:** | Draft — plan of record for SOW Deliverable #10 |
| **Companion docs:** | `D0_vcfdist-v3-SOW.md`, [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md), [`D6_vcfdist-v3-unphased-eval.md`](./D6_vcfdist-v3-unphased-eval.md), [`D5_multi-bed-stratification.md`](./D5_multi-bed-stratification.md), [`D4_ga4gh-compatibility.md`](./D4_ga4gh-compatibility.md), [`D7_vcfdist-v3-benchmarking-design.md`](./D7_vcfdist-v3-benchmarking-design.md) |

---

## Purpose & Scope

This document is the implementation plan for **Deliverable #10** of the v3.0.0 SOW.
It describes a **native MultiQC module** that ingests vcfdist's per-run TSV outputs and renders
them in an aggregated MultiQC report. The goal is to let users compare many vcfdist runs (e.g. across
samples, callers, or parameter sweeps) in a single interactive HTML report, in the same way the existing
`happy` (hap.py) module aggregates hap.py `summary.csv` files.

The module is developed as a native module inside a fork of **MultiQC v1.35** (branch
`feature/vcfdist-module`), modeled directly on the in-tree `happy` module so it is PR-ready for
upstreaming to `MultiQC/MultiQC`.

**In scope:** the four vcfdist aggregate summary TSVs, genome-wide headline metrics, and a per-stratum
accuracy breakdown.

**Out of scope (YAGNI):** the per-variant detail files (`switchflips.tsv`, `phase-blocks.tsv`,
`query.tsv`, `truth.tsv`); `genotype-errors.tsv`, which is a fifth *aggregate* summary rather than a
detail file and is simply not parsed yet; per-stratum precision-recall *curves*; per-stratum *phasing*
errors; and parsing `parameters.tsv` for a software-version string. *Deferred Work*, below, records why
each is deferred and what it would cost.

### Schema changes from D4, D5 and D6

The module reads only vcfdist's aggregate TSVs, and three deliverables shape them:

- [**D4**](./D4_ga4gh-compatibility.md) adds `FP_GT`, `FP_AL`, `QUERY_UNK`, `FRAC_NA`, and TiTv /
  het-hom ratio columns to the precision-recall summaries, and — the change that matters most for a
  report aggregating several runs — **switches counting from per-haplotype to per-site**, so a
  homozygous variant contributes 1 to TP/FP/FN instead of 2. Every count and rate moves.
- [**D5**](./D5_multi-bed-stratification.md), in its output-schema section, prepends a `STRATUM` key
  column to `precision-recall.tsv`, `precision-recall-summary.tsv`, and `genotype-errors.tsv`, and
  **splits** `phasing-summary.tsv` into `phasing-variants-summary.tsv` (per-variant, stratified) and
  `phasing-blocks-summary.tsv` (genome-level contiguity).
- [**D6**](./D6_vcfdist-v3-unphased-eval.md) **redefines the switch/flip error-rate denominator** from
  every merged query variant to `ASSESSED_PAIRS`, **removing `VARIANTS`** in the process, and adds
  `PHASED_HET_VARIANTS`, `TOTAL_HET_VARIANTS`, and `PHASED_HET_FRACTION` to the per-variant phasing file
  plus `PHASED_REGION_FRACTION` to the genome-level one.

The module is written against the **v3 schema**, while remaining readable on v2 output (see *Version
tolerance*, below), so one report can aggregate runs from either.

## Inputs Consumed

vcfdist writes files as `<out_prefix><name>.tsv`, where `<out_prefix>` is the user-supplied `-p/--prefix`
(e.g. `out/56a8e1.`). The **sample name** is the shared prefix with the known suffix stripped, so all
files from one run collapse to a single sample.

| File | Key columns | Used for |
| :-- | :-- | :-- |
| `*precision-recall-summary.tsv` | `STRATUM`, `VAR_TYPE`, `THRESHOLD` (`NONE`/`BEST`), `MIN_QUAL`, `PREC`, `RECALL`, `F1_SCORE`, `F1_QSCORE`, `TRUTH_TP`, `QUERY_TP`, `TRUTH_FN`, `QUERY_FP` | General stats, PR summary bargraph, stratified heatmaps, weakest-strata table |
| `*precision-recall.tsv` | `STRATUM`, `VAR_TYPE`, `MIN_QUAL`, `PREC`, `RECALL`, `F1_SCORE` (swept over quality) | PR-curve line plot (`*` only) |
| `*phasing-variants-summary.tsv` | `STRATUM`, `SWITCH_ERRORS`, `FLIP_ERRORS`, `SWITCH_ERROR_RATE`, `FLIP_ERROR_RATE`, `ASSESSED_PAIRS`, `PHASED_HET_VARIANTS`, `TOTAL_HET_VARIANTS`, `PHASED_HET_FRACTION` (`VARIANTS` on post-D5, pre-D6 output only) | Phasing table + general stats |
| `*phasing-blocks-summary.tsv` | `PHASE_BLOCKS`, `NG_50`, `SWITCH_NGC50`, `SWITCHFLIP_NGC50`, `PHASED_REGION_FRACTION` | Phasing table (genome-level; no `STRATUM` by construction) |

The two phasing files are merged into a **single row per sample**, reconstituting the pre-D5 phasing
table for display. This is the one place the module deliberately undoes a vcfdist-side split: D5 split
the file because per-stratum contiguity metrics are not meaningful, which is a *stratification*
concern, not a display concern — at the `*` level the columns belong in one table.

`PHASED_REGION_FRACTION` sits in the genome-level file rather than the per-variant one (D6 §5.3), so
the merged row is where it appears next to `PHASED_HET_FRACTION`. The two answer different questions and
disagree informatively where phasing is sparse, which is a reason to show them adjacently — with the
denominators named in the column descriptions, since NG50 normalizes by contig length while
`PHASED_REGION_FRACTION` uses evaluated BED bases.

### `STRATUM` handling

Every stratified file is read in one pass that routes rows by `STRATUM`:

- `STRATUM == "*"` → the genome-wide store, feeding general stats, the PR summary bargraph, the PR
  curve, and the phasing table.
- anything else → the stratified store, feeding the stratified heatmaps and the weakest-strata table.
- **column absent** → the row is treated as `*` (pre-D5 output).

**This filter is not optional.** `csv.DictReader` tolerates the added column silently, so without it
the module would aggregate all 182 strata into one sample's numbers and report wrong precision and
recall with no error. D5's downstream-consumers section names this module as the consumer that fails
this way, and it is the dangerous kind of failure: a plausible wrong number rather than a traceback.
At GIAB's 181 strata,
`precision-recall-summary.tsv` grows 8 → 1,456 rows and `precision-recall.tsv` 244 → ~44,408, so the
mis-aggregation is large, not marginal.

### Version tolerance

All parsing is by **column header name** (`csv.DictReader`), so column-set differences across vcfdist
versions are tolerated — missing columns are simply omitted from the tables. Concretely, the module
reads all of: pre-D6/D5 output (a single `phasing-summary.tsv`, no `STRATUM`, no rate columns at all
before commit `a1a2852`), post-D5-only output, and post-D6/D5 output. Search patterns for the legacy
`*phasing-summary.tsv` are retained alongside the two new ones, and the three globs cannot collide since
neither new filename ends in `phasing-summary.tsv`. `VARIANTS` is read when present — post-D5, pre-D6
output has it as the rate denominator — and ignored once `ASSESSED_PAIRS` is there.

**Two version differences are *not* silently tolerable**, and both are cases where a column keeps its
name while its meaning changes — exactly what header-name parsing cannot catch:

1. **D6's rate denominator.** `SWITCH_ERROR_RATE`/`FLIP_ERROR_RATE` change meaning without a rename, so a
   report aggregating pre- and post-D6 runs would place non-comparable rates in one column. The presence
   of `ASSESSED_PAIRS` is the discriminator: if some samples carry it and others do not, the module emits
   a `log.warning` naming the affected samples and leaves both rate columns blank for the pre-D6 ones.
   Counts (`SWITCH_ERRORS`, `FLIP_ERRORS`) are unaffected and always shown.
2. **D4's counting convention.** Every count and rate in `precision-recall*.tsv` changes when counting
   goes from per-haplotype to per-site, again with no column rename. Aggregating runs from both sides of
   that change puts non-comparable numbers in one table and one heatmap. The discriminator is the
   presence of any D4-added column (`QUERY_UNK` is the most reliable, since `FP_GT`/`FP_AL` can be
   legitimately zero): mixed presence triggers a `log.warning` naming the affected samples. Unlike the
   rate case the module does **not** blank the columns — precision and recall are still individually
   valid, they are simply computed on different denominators — so the warning is the whole mitigation,
   and it must say which samples used which convention.

Both checks share one implementation: compare the presence of a discriminator column across samples,
warn with names, and either blank or annotate. Neither can be inferred from the data itself.

## Report Content

The report is ordered **genome-wide first, per-stratum second**: every section before the stratified
ones shows only the `*` row, so a run without `-st` produces exactly today's report and a stratified
run adds sections rather than changing existing ones.

### General Statistics

Headline numbers per sample from the `BEST`-threshold `*` PR summary rows plus phasing: SNP F1, SNP
recall, INDEL F1, INDEL recall, and total switch+flip errors. D6's `PHASED_HET_FRACTION` is present but
hidden by default.

Switch+flip errors are **shown, not hidden**: they are the metric no competing tool reports, and the
SOW names them as vcfdist's differentiator. A metric hidden behind the column-configuration modal is
not a headline metric.

### Precision-Recall summary (`bargraph.plot`)

Grouped Precision/Recall/F1 bars per sample at the `BEST` threshold, with a variant-type switcher
(SNP/INDEL/SV/ALL) as separate datasets (`data_labels`).

A tabbed *table* was the original intent, but MultiQC tables do not support the dataset switcher (the
violin/table engine keeps a single dataset), so a grouped bargraph is used to get per-variant-type
tabs. Exact TP/FN/FP counts remain in the exported `multiqc_vcfdist_pr_summary` data file. The
upstream PR should call this trade-off out explicitly, since it is the first thing a reviewer will
question.

### Precision-Recall curve (`linegraph.plot`)

Recall (x) vs precision (y), one line per sample, from the quality-swept `precision-recall.tsv`
`*` rows, with the same variant-type switcher.

### Phasing errors (`bargraph.plot`)

Switch and flip error **counts** per sample, stacked, from the `*` row of
`phasing-variants-summary.tsv` — one bar per sample, two categories.

This is the plot the SOW asks for by name, and it earns headline placement for a reason the other
sections do not: precision and recall are directly comparable to `hap.py`/`vcfeval` output, whereas
switch and flip errors have no counterpart in either tool. A reader scanning the report should not have
to reach the final table to find the metric that distinguishes vcfdist.

Counts, not rates, are plotted. After D6 the rates share a denominator (`ASSESSED_PAIRS`) that varies
per sample, so bar heights would not be comparable across samples in the way a stacked count plot is;
the rates remain in the phasing table below, next to the denominator that produced them.

### Phasing summary (`table.plot`)

One row per sample from the merged `*` phasing rows: switch/flip errors, error rates, D6's
`PHASED_HET_FRACTION` and `PHASED_REGION_FRACTION`, `NG_50`, and switch/switchflip NGC50.
D6's denominator columns (`PHASED_HET_VARIANTS`, `TOTAL_HET_VARIANTS`, `ASSESSED_PAIRS`) and
`PHASE_BLOCKS` are present but hidden by default.

`PHASED_HET_FRACTION` is shown, not hidden: after D6 a switch-error rate is uninterpretable without
knowing what fraction of heterozygotes were phaseable at all. Rate columns need a **units-aware** read rather than a plain `%`-strip: pre-D5 output emits them as
e.g. `0.020930%`, a percent *string*, while [`D5`](./D5_multi-bed-stratification.md) §7.4's split writers
emit bare-float fractions (`0.00020930`) under the same column names. The `%` sigil is therefore the
discriminator — strip it and divide by 100 when present, take the value as-is otherwise. Unlike the two
cases below, the data itself resolves this one, so it needs no cross-sample warning.

### Stratified F1 (`heatmap.plot`, one per variant type)

**Emitted only when the input carries strata beyond `*`** (i.e. vcfdist was run with `-st`).

- Rows = samples, columns = strata, cell = F1 at that stratum's `BEST` threshold.
- One heatmap per variant type present in the data, each its own sub-section
  (`vcfdist-strat-f1-snp`, `…-indel`, `…-sv`, `…-all`).
- `pconfig`: `xcats_samples: False` (stratum names are **not** sample names — without this, MultiQC
  applies sample-name cleaning and `--ignore-samples` to strata), `ycats_samples: True`,
  `cluster_rows: False`, `cluster_cols: False` (deterministic column order — reports must be
  reproducible, and manifest order is meaningful to the user), `min: 0`, `max: 1`, `zlab: "F1"`,
  `square: False`, `display_values: False`.
- Columns follow **manifest order**, which is the order D5 emits after the leading `*` row.
- A stratum with no evaluated variants is `None`, not `0` — it renders as a gap. D5's
  zero-overlap-warning section expects empty strata as a normal outcome, and a `0` would read as total
  failure.

A heatmap is the only plot type in MultiQC 1.35 that scales visually to 181 columns. The two
alternatives fail concretely: a stratum-keyed dataset switcher (`data_labels`) would render 182
buttons, and a full sample × stratum table exceeds `config.max_table_rows` (default 500,
`config_defaults.yaml:110`), above which `table.plot()` silently degrades to a violin plot
(`plots/violin.py:936`). The cost of the heatmap is that `HeatmapConfig` has **no** `data_labels`
field (`plots/heatmap.py:48-71`), so one plot is one metric × one variant type — hence F1 only, and
one section per variant type. Exact per-stratum values live in the weakest-strata table and the
exported stratified matrix.

### Weakest strata (`table.plot`)

The exact-numbers companion to the heatmaps: one table, defaulting to the **20 worst
(stratum, variant-type) cells**, worst F1 first.

- **Row key** is the composite `"<sample> · <stratum> · <var_type>"` with the sample also carried as
  its own column. This follows the in-tree `happy` module, which keys rows
  `f"{s_name}_{row['Type']}_{row['Filter']}"` with a `sample_id` field (`happy.py:94`) — stratified
  composite row keys are upstream-precedented.
- **Columns:** sample, `STRATUM`, `VAR_TYPE`, `F1_SCORE`, `PREC`, `RECALL`, `MIN_QUAL`, and
  `TRUTH_TP`/`TRUTH_FN`/`QUERY_FP` (hidden by default). MultiQC tables are user-sortable, so
  worst-F1-first is the default order only.
- **Selection:** rank (stratum, variant-type) cells by their *minimum* F1 across samples, then emit
  **all** samples' rows for each selected cell, so samples stay comparable within a row group. Cells
  with `TRUTH_TP + TRUTH_FN == 0` are excluded from ranking — otherwise the empty strata D5 expects
  would occupy every slot with `F1 = 0` and crowd out real weakness. The excluded count is logged.
- **Row budget:** `20 × n_samples` rows. If that would exceed `config.max_table_rows`, the cell count
  is reduced to fit and a `log.warning` names the truncation — a silent cap would read as "these are
  the only weak strata". Configurable via `vcfdist_max_strata_cells`.

## Exported Data Files

Written via `write_data_file` so values land in `multiqc_data/` for reuse:

| File | Contents |
| :-- | :-- |
| `multiqc_vcfdist_pr_summary` | `*` rows per sample × variant type, all `BEST`-threshold columns |
| `multiqc_vcfdist_phasing` | merged `*` phasing row per sample |
| `multiqc_vcfdist_stratified` | *new* — the **full** sample × stratum × variant-type matrix, including strata omitted from the heatmaps and the weakest-strata table (JSON, mirroring `happy`'s nested exports) |

The third file is the escape hatch that makes the weakest-strata top-N honest: nothing the report
elides is unavailable.

## Module Structure & Registration

Modeled on `multiqc/modules/happy/`:

- `multiqc/modules/vcfdist/__init__.py` — `from .vcfdist import MultiqcModule`
- `multiqc/modules/vcfdist/vcfdist.py` — `MultiqcModule(BaseMultiqcModule)`; discovers files via
  `find_log_files`, parses each TSV, calls `add_data_source` / `is_ignore_sample` /
  `add_software_version(None)`, writes parsed data with `write_data_file`, and raises
  `ModuleNoSamplesFound` when nothing matches.
- `multiqc/search_patterns.yaml` — content-anchored sub-patterns: `vcfdist/prsummary`,
  `vcfdist/prcurve`, `vcfdist/phasing_errors`, `vcfdist/phasing_metrics`, and `vcfdist/phasing`
  (legacy pre-D5 single file).
- `multiqc/config_defaults.yaml` — add `- vcfdist` to `module_order` next to `happy`.
- `pyproject.toml` — entry point `vcfdist = "multiqc.modules.vcfdist:MultiqcModule"`.

Pin to MultiQC **1.35**: the module/plot API moved across versions, so third-party plugin examples
cannot be trusted verbatim, as [`D1`](./D1_vcfdist-v3-design-doc.md) notes in its MultiQC section.

## Design Decisions

| # | Decision | Rationale |
| :-: | :-- | :-- |
| 1 | Genome-wide `*` sections first; stratified sections appended | A run without `-st` yields today's report unchanged; stratification adds, never rewrites |
| 2 | Filter `STRATUM == "*"` for all headline metrics; absent column ⇒ `*` | Without it, 182 strata silently aggregate into wrong P/R — the failure D5 documents for its other downstream consumers, applied to this module |
| 3 | Merge the two split phasing files into one display row | D5's split serves per-stratum correctness; at `*` the original eight columns are one table |
| 4 | Switch/flip errors get a headline bargraph and a visible general-stats column | The SOW names them as vcfdist's differentiator; unlike P/R they have no `hap.py`/`vcfeval` counterpart, so they must not be buried in the last table |
| 5 | Plot error **counts**, not rates | After D6 the rate denominator (`ASSESSED_PAIRS`) varies per sample, so bar heights would not be comparable; rates stay in the table beside their denominator |
| 6 | Do not parse `switchflips.tsv` or `phase-blocks.tsv` | The aggregates already carry every displayed number, and phase blocks are *not* derivable from `switchflips.tsv` — see *Deferred Work* |
| 7 | Heatmap for the per-stratum overview | Only plot type that scales to 181 columns; a switcher would need 182 buttons, and a table degrades to a violin past 500 rows |
| 8 | F1 only in the heatmap, one heatmap per variant type | `HeatmapConfig` has no `data_labels` (`plots/heatmap.py:48-71`) — one plot is one metric × one type |
| 9 | Worst-20 cells table alongside the heatmap | Heatmap shows shape but no values; the table gives auditable numbers where they matter |
| 10 | Exclude zero-truth-variant cells from worst-N ranking; log the count | Empty strata are expected by D5 and would otherwise fill every slot at `F1 = 0` |
| 11 | Warn on any truncation or mixed rate-denominator input | A silent cap reads as completeness; a silent denominator change reads as a regression |
| 12 | Full stratified matrix always exported to `multiqc_data/` | Makes the top-N display cap lossless |
| 13 | Composite `sample · stratum · var_type` table row keys | Precedented in-tree by `happy.py:94`; a `STRATUM` key column cannot exist in a sample-keyed table |
| 14 | Per-stratum PR curves and per-stratum phasing deferred | No usable multi-series idiom at 181 strata — see *Deferred Work* |

## Deferred Work

- **Per-stratum precision-recall curves.** D5's decision record puts the per-stratum QUAL sweep in
  scope specifically to enable these, and `precision-recall.tsv` carries the data. They are
  nonetheless deferred: the only switcher available to `linegraph.plot` is `data_labels`, which at
  181 strata is the same unusable 182-button control rejected for the heatmaps, and one line per
  sample × stratum in a single plot is unreadable. A viable follow-on is a curve plot restricted to
  the strata named in the weakest-strata table. **D5's expectation is deliberately unmet, not
  overlooked.**
- **Per-stratum phasing errors.** `phasing-variants-summary.tsv` is stratified, but the display would
  need a second heatmap family for a metric with no natural [0,1] scale. Deferred.
- **`genotype-errors.tsv`.** Also gains `STRATUM` in D5; not parsed at all today.
- **`switchflips.tsv` and `phase-blocks.tsv`.** Not parsed: the headline switch/flip bargraph and the
  phasing table are both fully served by the aggregate summaries, so neither per-event file is needed
  for anything in this design. Two properties of `switchflips.tsv`, measured on
  `out/56a8e1.switchflips.tsv` (chr20), are recorded here so a future follow-on does not rediscover
  them the hard way:

  - **Phase blocks cannot be derived from it.** Its columns are `CONTIG`, `START`, `STOP`,
    `SWITCH_TYPE`, `VARIANT`, `PHASE_BLOCK`, and only blocks that *contain an error* appear at all —
    60 distinct `PHASE_BLOCK` ids against 676 blocks in `phase-blocks.tsv`, i.e. 9%. `max(PHASE_BLOCK)`
    is a lower bound on the count, not the count, and the file carries no block coordinates, so block
    sizes, `NG_50`, and the NGC50 metrics are underivable in principle. Block count is already a
    column of `phasing-blocks-summary.tsv`, so nothing needs deriving.
  - **One flip is two rows.** `SWITCH_TYPE` takes `FLIP_BEG`, `FLIP_END`, and `SWITCH_ERR`; the chr20
    file has 46/46/14, whose 106 rows correspond to the summary's 46 flips and 14 switches. Counting
    rows would report 106 errors instead of 60.

- **Other per-variant detail files** (`query.tsv`, `truth.tsv`), which additionally gain D5's
  set-valued `STRATA` column.
- **Software version from `parameters.tsv`.** D5 adds `stratification_tsv` and `n_strata` rows there,
  so it would also replace the `ASSESSED_PAIRS` heuristic under *Version tolerance* with an exact
  version comparison. The module calls `add_software_version(None)` until then.

## Verification

Dev-install the fork (`pip install -e .`) and run `multiqc` against:

1. **`out/`** — the checked-in pre-D5 outputs (prefixes `56a8e1.`, `df5b32.`), confirming two samples
   parse, all sections render, and the legacy `phasing-summary.tsv` path still works.
2. **A stratified fixture** — a D5 run over the `strat.tsv` manifest from D5's integration-test
   section, confirming the `*` sections are numerically identical to the same run without `-st`, the
   heatmaps appear, the empty-but-on-contig stratum renders as a gap rather than `0`, and the
   weakest-strata table's numbers are consistent with `multiqc_vcfdist_stratified`.

Both are added as a small CI fixture set.

**Do not verify against `demo/results/`.** That directory holds two generations of output — `dev.*`
files alongside older unprefixed `phasing-summary.tsv` / `precision-recall-summary.tsv` and a
`superclusters.tsv` that v3 no longer writes. The prefix-stripping rule reads the unprefixed set as a
degenerate sample name, and the stale files no longer match current column sets. Either regenerate
that directory as part of [`D9`](./D9_docs-and-tutorial.md) or leave it out of the fixture set.

## Relationship to D4

MultiQC ships a core `happy` module that consumes a GA4GH quantifier's `*.summary.csv`. That is not a path
to vcfdist's metrics: [`D4`](./D4_ga4gh-compatibility.md) §1.2 keeps the native TSVs as what vcfdist
reports, so there is no vcfdist-produced `summary.csv` for the `happy` module to read. A *user* with an
existing GA4GH pipeline can do that round-trip themselves, but it is not a path this project provides.

Two consequences:

- **This module is the only path by which vcfdist's own metrics reach a MultiQC report**, which is what
  makes the core-module route worth the upstream review cycle and the standalone-plugin fallback worth
  keeping genuinely ready.
- **The `happy` module could not substitute for this one anyway**, and the reasons belong in the upstream PR
  description: it cannot show phasing, which is vcfdist's differentiator and the reason this deliverable
  exists; and a GA4GH consumer counts hard integer `BD` labels, so vcfdist's fractional `BC` partial credit
  cannot survive the export.

Where D4 *does* touch this module is the schema: it adds the `FP_GT`/`FP_AL`/`QUERY_UNK`/`FRAC_NA` and ratio
columns the tables surface, and it changes the counting convention under fixed column names — see *Version
tolerance*, which is where that lands as real work.
