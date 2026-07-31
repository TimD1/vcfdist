# D4: GA4GH Benchmarking VCF Compatibility

- **Status:** Design — complete, pending review. Decisions recorded in §11.
- **Issue:** [#49](https://github.com/TimD1/vcfdist/issues/49) (D4: GA4GH Benchmarking VCF Compatibility)
- **Branch:** `49_td_D4-ga4gh-compatibility`
- **Area:** `src/phase.cpp`, `src/variant.{cpp,h}`, `src/print.cpp`, `src/globals.{cpp,h}`,
  `src/defs.h`
- **Owns:** the GA4GH annotations on the `summary.vcf` record layout, the per-site counting convention,
  and the match-tier decision model.

## 1. Background

### 1.1 What the GA4GH intermediate VCF is for

The GA4GH benchmarking framework ([Krusche et al. 2019](https://doi.org/10.1038/s41587-019-0054-x))
separates a *comparison engine* from a *quantifier*. The engine decides, per variant, whether truth
and query agree; the quantifier turns those decisions into stratified counts, precision/recall, and
ROC curves. The interchange between them is the **intermediate VCF (VCF-I)**: a two-sample
`TRUTH`/`QUERY` VCF carrying per-sample decision annotations.

vcfdist already emits something close to this shape in `summary.vcf` and describes it as
GA4GH-compatible. It is not yet consumable by a GA4GH quantifier; §3–§5 enumerate why.

### 1.2 Who consumes it

VCF-I is a **one-hop** format: the only tools that read `BD`/`BK` are GA4GH quantifiers. Two producers
exist today — hap.py's built-in `xcmp`, and `rtg vcfeval` in its GA4GH build. Conforming makes vcfdist
a third producer of the same interchange format, following an established pattern.

**Conformance buys exactly one thing: a user who already has a GA4GH pipeline can point it at
`summary.vcf`.** vcfdist's own metrics are its native TSVs, and every downstream path this project
builds reads those directly — D10's MultiQC module, D11's nf-core module and snakemake wrapper, and
D7's benchmarking harness. A GA4GH export re-denominates the counts and discards fractional `BC`
credit (§4), so it is not the reporting path, and **numbers derived from one are never quoted as
vcfdist results** — not in the release notes, not in the blog post, not in the benchmarking report.

### 1.3 What "compatible" means here

The design targets what a GA4GH quantifier *reads*, established from the reference implementations'
source rather than the prose spec alone. Since nothing else consumes VCF-I, **the standard is
agreement with those implementations, not spec purity** — where the spec and both implementations
disagree, the implementations win, because they are what parses our output (§3.4).

Conformance is asserted against vcfdist's own output (§9), which cannot catch a behavior neither the
spec nor the source made apparent; §9's conformance-by-inspection item is what covers that.

Primary sources: `ga4gh/benchmarking-tools` `doc/ref-impl/{intermediate,inputs,outputs}.md`
and `doc/standards/GA4GHBenchmarkingPerformanceMetricsDefinitions.md`; hap.py
`src/python/Haplo/quantify.py`, `src/c++/main/quantify.cpp`,
`src/c++/lib/quantify/{GA4GHQuantify,BlockQuantify}.cpp`,
`src/c++/lib/variant/VariantStatistics.cpp`, `src/c++/lib/tools/Roc.cpp`; rtg-tools
`src/main/java/com/rtg/vcf/eval/Ga4ghEvalSynchronizer.java`; `samtools/hts-specs` `VCFv4.4.tex`.

## 2. Goals / Non-Goals

**Goals**

- Emit a `summary.vcf` that satisfies the GA4GH benchmarking-VCF contract, so any GA4GH quantifier
  parses it without error.
- Align vcfdist's own counting unit with hap.py and vcfeval — one unit per variant per callset, not
  one per haplotype copy — because divergent denominators are an incompatibility, not a capability
  (§4).
- Collapse the per-haplotype `BD`/`BK` lists on each record into one decision per site, via an
  explicit match-tier ladder (§6).
- Extend the match-kind vocabulary to cover phasing, implementing GA4GH Comparison Method #4 (§6).
- Report natively the metrics a GA4GH quantifier would derive from our own output — the
  `FP.gt`/`FP.al` split and the unassessed-call fraction — so vcfdist's TSVs and the export never
  disagree about what is measurable (§7).

**Non-Goals**

- **The `summary.vcf` record shape.** Un-splitting homozygous records, keeping het-alt records split,
  per-variant ploidy, and the `Number=.` per-haplotype fields are all
  [`D3`](./D3_retain-info-format-fields.md)'s (D3 §4.2–§4.5, summarised in §5's table). This
  deliverable annotates that shape; it does not change it.
- `som.py`. It makes its own TP/FP/FN decisions rather than consuming a pre-decided VCF.
- Changing the alignment or credit engine. This deliverable changes how per-haplotype credit is
  *aggregated into a per-site decision*, not how credit is computed.
- Stratification, inside the VCF or in the TSVs. Regions are supplied to a GA4GH consumer externally
  (§3.3); [`D5`](./D5_multi-bed-stratification.md) gives the §7.1 summary columns their stratum axis.

## 3. The `quantify` contract

### 3.1 Tags vcfdist must emit

| Tag | Kind | Number | Type | Value | Read at |
| --- | ---- | ------ | ---- | ----- | ------- |
| `GT` | FORMAT | 1 | String | the sample's reported genotype | `GA4GHQuantify.cpp:173`, `VariantStatistics.cpp:489` |
| `BD` | FORMAT | 1 | String | `TP` / `FP` / `FN` / `N` | `GA4GHQuantify.cpp:103`, `BlockQuantify.cpp:498-499` |
| `BK` | FORMAT | 1 | String | `.` / `lm` / `am` / `gm`, plus vcfdist's `pm` (§6.4) | `GA4GHQuantify.cpp:104` |
| `QQ` | FORMAT | 1 | Float | variant quality, for ROC | `BlockQuantify.cpp:326,351,358,505,523` |
| `BS` | INFO | . | Integer | benchmarking superlocus ID | `BlockQuantify.cpp:444`, `quantify.cpp:595` |

That is the entire required set. Both reference producers emit exactly these
(`Ga4ghEvalSynchronizer.java:123-129`), independently confirming the list is complete.

`INFO/BS` is the item neither the SOW nor D1 anticipated, and it is not optional in practice.
`quantify` groups records by `BS` to propagate query `QQ` onto truth records
(`BlockQuantify.cpp:321-383`: each truth record's ROC quality derives from query TPs sharing its
superlocus) and to set the `TS_boundary` / `TS_contained` region flags. It also splits its parallel
work blocks on `BS` boundaries (`quantify.cpp:595`). Absent `BS`, every record becomes its own
superlocus: quantify still runs, but truth-side ROC degrades and boundary flagging is lost.
vcfdist's supercluster index is exactly the intended value.

Boundaries are detected by *change*, not uniqueness — `vbs != current_bs || vbs < 0 || vchr !=
current_chr` (`BlockQuantify.cpp:456-457`) — so IDs need not be globally unique, but all records of
one superlocus must be adjacent. A per-contig monotonic index satisfies this, since the chromosome
is compared too.

**`Number=.` for `BS`, not `Number=1`.** `intermediate.md` specifies `Number=1`, but both
implementations declare `.` — `Ga4ghEvalSynchronizer.java:123` uses `VcfNumber.DOT` and
`GA4GHQuantify.cpp:59` appends `Number=.`. Only the first integer is ever read.

### 3.2 Tags vcfdist must *not* emit

These are computed by `quantify` and written back over whatever the input contained:

| Tag | Kind | Overwritten at |
| --- | ---- | -------------- |
| `BVT` | FORMAT | `GA4GHQuantify.cpp:208`, derived from REF/ALT + `GT` |
| `BLT` | FORMAT | `GA4GHQuantify.cpp:209`, derived from `GT` |
| `BI` | FORMAT | `GA4GHQuantify.cpp:207`, from `VariantStatistics::extraCountsToBI` |
| `Regions` | INFO | `QuantifyRegions`; header appended at `quantify.cpp:382` |
| `VTC` | INFO | `GA4GHQuantify.cpp:217`, only under `--output-vtc` |

A GA4GH quantifier **derives** variant type from REF/ALT rather than reading `BVT`/`BLT`/`BI` as
input. `vcfeval` emits none of them either.

There is an active reason not to emit them: `GA4GHQuantify::updateHeader` appends the `BVT` and
`BLT` header lines **unconditionally**, with no `bcf_hdr_get_hrec` guard (unlike `BD`, `BK`, `BI`,
`QQ`, and `BS`), so emitting our own would produce duplicate declarations.

### 3.3 Site-level requirements, and what a consumer must do

The first three items are obligations on **our output**. The last three are things the *user* does
when feeding a quantifier; they are recorded so D9 can document the export path, and none of them
is a step this deliverable performs (§1.2).

- **Sample columns must be named `TRUTH` and `QUERY`, in that order.** Indices 0 and 1 are
  hardcoded (`BlockQuantify.cpp:498-501`), and `rocEvaluate` returns early unless there are exactly
  two samples (`:493-497`). Already satisfied; keep stable.
- **The absent sample gets `GT=.`** A missing genotype yields `BVT=NOCALL`, which `rocEvaluate`
  skips (`:503`, `:521`), so an FP query record contributes nothing to the truth side. This is what
  `print_var_empty` already does, and what `vcfeval` does (`Ga4ghEvalSynchronizer.java:214-215`).
- **Non-`PASS` `FILTER` on an evaluated record is demoted by the consumer.** `addROCValue` turns
  filter-failed TPs into FNs and filter-failed FPs/UNKs into Ns (`BlockQuantify.cpp:143-159`).
  vcfdist's own `VCFDIST_*` tags land only on `BD=N` records, which are non-assessed anyway — but
  [`D3`](./D3_retain-info-format-fields.md) §4.6 preserves the *input's* `FILTER` verbatim, so a user
  who accepts a non-`PASS` filter via `--filter` will see those calls demoted by a quantifier.
  Deliberately not prevented: rewriting the caller's `FILTER` to `PASS` would misreport the input,
  Documented in D9 as a caveat of the export path.
- **bgzip + tabix are the user's step.** `quantify.cpp:317` sets `reader->require_index = 1` and
  errors with "Failed to open or file not indexed", so a plain `.vcf` cannot be quantified. vcfdist
  continues to write uncompressed `summary.vcf`; compressing and indexing it is one `bgzip` +
  `tabix` away and is not made a vcfdist output.
- **Invocation (user-side):** `--type ga4gh` (defaults to `xcmp`) and `--roc QQ` (defaults to `QUAL`).
- **Confident regions (user-side):** pass vcfdist's `--bed` as the quantifier's confident-region file. It is named `CONF` and
  relabels anything outside it as `UNK` (`GA4GHQuantify.cpp:106-108`). Stratification BEDs go via
  `--stratification` or `--stratification-region NAME:file.bed`, never embedded in the VCF.

### 3.4 Where the spec and the implementations disagree

Recorded so future readers do not "fix" our output toward the prose:

| Point | Spec text | Implementations | We follow |
| ----- | --------- | --------------- | --------- |
| `INFO/BS` cardinality | `Number=1` | `Number=.` (both) | `Number=.` |
| `BVT`/`BLT`/`BI` | listed as VCF-I annotations | computed by quantify, input ignored | do not emit |
| Confident regions | applied *after* comparison (`inputs.md`) | quantify applies them | vcfdist's own — filtered at parse time |

### 3.5 `BK` does not affect totals

`rocEvaluate` branches on `BD` and `BVT` only. `BK` reaches the counts by two narrower routes: it
is part of quantify's count-map key (`GA4GHQuantify.cpp:123`), and `makeObservationFlags`
(`Roc.cpp:49-65`) maps it to observation flags where `am` → `FP.gt` and `lm` → `FP.al`
(`Roc.hh:126-137`).

Omitting `am` is harmless *for TP/FP/FN totals*, but it leaves hap.py's `FP.gt` / `FP.al` breakdown
wrong. vcfdist currently never emits `am`, so `FP.gt` is always zero and
genotype errors are misattributed to `FP.al`. `ac_errtype` supplies the missing signal, and
`vcfeval` confirms the intended pairing: `BD=FP` with `BK=am` for a genotype mismatch
(`Ga4ghEvalSynchronizer.java:245-246`).

This section is also what licenses the `pm` extension of §6.4: since no total depends on `BK`, adding a
value cannot move a count. Only the two narrower routes are affected, and §6.4 works through both.

## 4. Counting convention

**vcfdist counts per haplotype copy.** `write_precision_recall` tallies inside
`for (int hi = 0; hi < HAPS; hi++) { if (var_on_hap(vi, hi)) ... }` (`print.cpp:288-300`, truth loop
at `:318-331`), so a homozygous variant contributes **2** to TP/FP/FN. The per-haplotype split in
`summary.vcf` faithfully mirrors this: two lines, two counted alleles.

**hap.py and vcfeval count per site.** `VariantStatistics::add` treats a homalt genotype as a single
allele (`VariantStatistics.cpp:491-500`) and increments its bucket by exactly one per record per
sample (`:557`). A `1|2` record is likewise one location. GA4GH Comparison Method #2 Note #1 is
explicit: *"A homozygous variant shall be considered a single ALT allele."* vcfeval is also
per-record, reconciling many-to-many matches with an `INFO/CALL_WEIGHT` float rather than by
changing record shape.

So this is a **genuine methodological difference in denominators**, not a formatting artifact: under
vcfdist's convention a homozygous variant carries twice the weight of a heterozygous one relative to
hap.py.

**Decision: adopt the per-site convention.** vcfdist counts one unit per variant per callset.

1. The SOW's goal is to be "at least as capable on every axis" as the incumbents; divergent
   denominators are not a capability, they are an incompatibility.
2. It makes D7's concordance analysis interpretable. Under divergent conventions every homozygous
   site looks like a disagreement, which would swamp the cross-tool comparison and mask the
   complex-variant differences that are the actual finding. (This reason stands independently of the
   GA4GH export, since D7 compares per-variant decisions directly, not through a quantifier.)
3. It removes a fragile special case. `print.cpp:305-315` currently detects `AC_ERR_2_TO_1` on a
   *non-carried* haplotype and hand-corrects the truth counters (`truth_tp--`, `truth_fn++`)
   precisely because per-haplotype counting mis-attributes allele-count errors. A site-level
   decision from the match tier (§6) handles this directly, and the special case is deleted.

   Stratification ([`D5`](./D5_multi-bed-stratification.md)) adds a stratum axis to the same counting
   loops, and stratifies one uniform per-site count rather than this correction.

**Consequences, all of which must be release-noted:**

- **Every published vcfdist number changes.** Prior-publication figures used per-haplotype counting.
  D7 must regenerate baselines and report the delta; comparisons against v2 need convention-matched
  runs.
- **Genotype errors become one FP/FN instead of a TP plus an FP.** Under Method #3 a genotype error
  is FN on truth and FP on query at that site, which is also what `FP.gt` reports. This is alignment
  with hap.py's default, not a new opinion.
- **Multiallelic sites remain a residual discrepancy** of +1 per hetalt site, measured and reported
  rather than assumed negligible (§5.3, measured in D7).
- **Fractional credit cannot survive any GA4GH export.** `BD` is a hard label, and vcfeval's
  `CALL_WEIGHT` — the natural home for a fractional weight — has no reader anywhere in hap.py
  (`quantify.cpp`, `GA4GHQuantify.cpp`, `BlockQuantify.cpp`, `VariantStatistics.cpp`,
  `Roc.cpp` all ignore it); it is vcfeval-internal. So a quantifier fed `summary.vcf` reports
  vcfdist's *thresholded* decisions. vcfdist's native TSVs remain the only place its partial-credit
  scoring is visible, which is the concrete reason §1.2 keeps them as the reporting path.

## 5. Phase 1 — conformant output

The record shape is [`D3`](./D3_retain-info-format-fields.md)'s. What D3 §4.2–§4.5 establishes, and
what this deliverable builds on it:

| Established in D3 | Consequence here |
| ----------------- | ---------------- |
| One record per variant; homozygotes un-split | A record is a site, so a site-level `BD`/`BK` is expressible at all (§6) |
| Het-alt (`1\|2`) records stay split | A permanent count residual at het-alt sites vs. tools that count one location per record; measured by D7 (§5.3) |
| Per-variant ploidy; `GT` from `orig_gts` + `ploidies` | `quantify` derives `BVT`/`BLT` and counts alleles from `GT` (`VariantStatistics.cpp:486-543`), so a haploid call written `1\|0` would be counted `CT_HET` instead of `CT_HEMI`. D3 fixes that; nothing further is needed here |
| Per-haplotype fields as `Number=.` lists, including `BD`/`BK` | §6 replaces the `BD`/`BK` lists with a single value each, declared `Number=1`. `BC`, `RD`, `QD`, and `SG` stay `Number=.` |

The remaining Phase-1 work is tag naming and the header (§5.1–§5.2), plus the residual that
D3's het-alt decision creates (§5.3).

### 5.1 Tag renames

Two tags swap names, because `BS` is exactly the name GA4GH gives to vcfdist's supercluster and
vcfdist currently spends it on something else:

- **`SC` → `BS`.** `FORMAT/SC` is the supercluster index (`phase.cpp:47`), which *is* hap.py's
  benchmarking superlocus. It moves to `INFO/BS` (§3.1) — a rename in effect, but also a kind change
  from FORMAT to INFO, since that is where `quantify` reads it (`BlockQuantify.cpp:444`). No
  `FORMAT/SC` remains: one value per site needs no per-sample column.
- **`BS` → `PBS`.** `FORMAT/BS` currently means "Block Phase" — the phase block's keep/swap state
  (`phase.cpp:51`) — and must vacate the name. It becomes `FORMAT/PBS`, "Phase Block State",
  adjacent to the existing `PB` ("Phase Block"). The two could legally coexist as different header
  types, but `BS` meaning both a superlocus ID and a phase state in one file is a permanent footgun.

Rename now, while the record shape is changing anyway. Both are breaking output changes and belong in
the release notes (#8) and D9's docs; `PBS` is the only tag in `summary.vcf` whose *meaning* is
unchanged from v2 under a new name, so it is the one most likely to be silently mis-parsed by an
existing consumer.

#### 5.1.1 `SC` is the only field that changes column

`SC` moving to INFO invites the same question of every other per-supercluster and phasing field —
`PB`, `PBS`, `VP`, `FE`, `SG`. **Decided: they all stay in FORMAT.** Only `SC` moves. Recorded here
because the tempting argument for moving them is wrong in an instructive way.

**Being `.` on TRUTH is evidence of the opposite.** `print_var_sample` writes `PBS`, `VP`, and `FE`
as `query ? <value> : "."` (`variant.cpp:532-534`), and it is tempting to read a permanently empty
TRUTH column as a sign the value is really a property of the locus. It is the reverse: INFO has no
per-sample slot at all, so "a value for one callset and none for the other" is exactly the
distinction FORMAT exists to draw. A field being unpopulated for a sample says nothing about whether
it is per-sample; it says that sample has no value for it.

**These are phasing properties of the variant, not of the locus.** A variant's phase is which
haplotype its alleles sit on — a property of the call, adjacent to `GT`, and VCF places it in FORMAT
accordingly: the `|` separator lives inside `GT`, and `PS` is a **reserved FORMAT key**
(`VCFv4.4.tex`). That is the direct precedent, because a phase set already *is* a locus-spanning
grouping that the spec nonetheless makes per-sample. `PB` is vcfdist's output analogue of the input
`PS` (phase blocks are cut from query phase sets, `phase.cpp:250-259`), so if `PS` is FORMAT then
`PB` is too — and keeping `PS` in FORMAT while moving `PB` to INFO would have been internally
inconsistent. `PBS` is that block's state, `VP` the variant's phase within it, and `FE` whether the
two disagree; all four sit on the same axis.

| Tag | Meaning | Column |
| --- | ------- | ------ |
| `SC` → `BS` | supercluster index | **INFO** — the sole move (§5.1) |
| `PB` | phase block index in contig | FORMAT — output analogue of `PS` |
| `BS` → `PBS` | phase block keep/swap state | FORMAT — state of that sample's block |
| `VP` | variant phase | FORMAT — phase of that sample's call |
| `FE` | flip error | FORMAT — that sample's phase vs. its block |
| `SG` | sync group index | FORMAT — per-haplotype, see below |
| `BC`, `RD`, `QD` | credit, edit distances | FORMAT — genuinely differ |
| `PS` | input phase set | FORMAT — reserved key |
| `GE` | genotype error | FORMAT — see below |
| `GT`, `BD`, `BK`, `QQ` | GA4GH decision set | FORMAT (§3.1) |

**Why `SC` is different.** It is not a property of either callset's variant but of a *grouping across
both*: a truth variant and a query variant belong to one supercluster for the same reason, and the
index means the same thing from either side. That is a locus property in the way a phase is not. GA4GH
independently requires it in INFO (`BlockQuantify.cpp:444`), so both arguments point the same way —
but the semantic one would hold even if the spec were silent.

**`SG` also stays, despite sharing an ID space across callsets.** `sync_group` is assigned to truth
and query from the same counter under the same `truth_hap` (`dist.cpp:550`, `:564`), so a matched pair
genuinely carries the same group ID — the one field here with a real claim to being per-locus. But it
is **per-haplotype**, and `phase.cpp:389` swaps query's `HAP1`/`HAP2` entries when the computed
genotype is swapped, so query's list is in query-haplotype order while truth's stays in
truth-haplotype order. Under D3 §4.5 a per-haplotype field's value order follows the emitted `GT`, which
is per-sample; one INFO list cannot carry both orderings.

**`GE` stays, and §6.5 is why it is not even close.** It reads `.` on TRUTH today only because
`ac_errtype` is never set for truth (`variant.cpp:137` defaults it to `AC_UNKNOWN`, rendering `.` via
`ac_strs`, `globals.cpp:25`) — the gap §6.5 commits to closing. Afterwards a truth record's value is
the inverse of the query record's (`-` where the query reads `+`).

**One residual inconsistency.** `print_var_empty` populates `PB` (and `SC`) even for a
sample with no call, while `PBS`, `VP`, and `FE` are `.` on TRUTH (`variant.cpp:492`). Since `PB` is
per-sample, a truth-side value is arguably wrong for the same reason `PBS`'s is — though it is also
defensible, as a truth variant does fall within the block's interval. Not changed here: it is
cosmetic, and this commit should stay a pure move-and-rename. Noted so the populated TRUTH
column is not later mistaken for the per-locus tell it is not.

### 5.2 Header

`##fileformat=VCFv4.2` is retained. Following D3 §4.5's choice of `Number=.` over `Number=P`, nothing
here uses a post-4.2 feature — `Number=.`, `Number=1`, and an Integer `INFO` field are all 4.2-era —
so there is no htslib floor for consumers either, including the ancient htslib bundled with hap.py.

Emission was never constrained by the linked htslib in any case: vcfdist writes VCFs with `fprintf`
(`phase.cpp`, `variant.cpp:221`), not htslib's writer.

### 5.3 The het-alt count residual

[`D3`](./D3_retain-info-format-fields.md) §4.3 settles that a `1|2` record stays two co-located
biallelic records, for reasons that are about normalization and credit rather than about GA4GH
conformance. The consequence *here* is a permanent count residual at het-alt sites relative to
hap.py and vcfeval, which count a `1|2` record as one location: vcfdist reports +1 per het-alt site.

That residual is measured and reported by D7 rather than assumed negligible. It is the one place
where §4's counting alignment is deliberately incomplete, and stating its size is what keeps the
alignment claim honest.

### 5.4 What changes in `summary.vcf`

Relative to the shape [`D3`](./D3_retain-info-format-fields.md) leaves behind — not relative to v2,
which D3 §4 already documents:

| | After D3 | After this deliverable |
| - | -------- | ---------------------- |
| `BD` / `BK` | `Number=.`, one value per haplotype | `Number=1`, one decision per site (§6) |
| `BK` vocabulary | `gm` / `lm` / `.` | adds `am` and `pm` (§3.5, §6.4) |
| `INFO` | `.` | `BS=<supercluster>` — the only INFO field (§5.1.1) |
| `FORMAT/SC` | supercluster index | renamed and moved → `INFO/BS` (§5.1) |
| `FORMAT/BS` | Block Phase (keep/swap) | renamed → `FORMAT/PBS`, same meaning and column (§5.1) |
| `BC`, `RD`, `QD`, `SG` | `Number=.` | unchanged — genuinely per-haplotype |
| `fileformat` | VCFv4.2 | VCFv4.2, unchanged (§5.2) |

## 6. Phase 2 — the match-tier ladder

### 6.1 Rationale: GA4GH Comparison Methods #1–#4

`GA4GHBenchmarkingPerformanceMetricsDefinitions.md` defines four comparison methods of increasing
stringency, and the extended match vocabulary maps onto them exactly:

| Tier | GA4GH Comparison Method |
| ---- | ----------------------- |
| `lm` | #1 Loose Regional Comparison |
| `am` | #2 Allele Match Required |
| `gm` | #3 Genotype Match Required |
| `pm` | #4 Genotype Match **and Local Phasing** Required |

Method #4 is fully specified — *"sites for which phasing of the variants in the site is incorrect"*
count as FN, and phasing errors count as FP — and is implemented by no tool, because neither
incumbent reports phasing. Method #5 is marked *"To be developed."*

`--stringency pm` therefore makes vcfdist the first implementation of Method #4, and places
vcfdist's phasing differentiator *inside* the standard's framework rather than beside it. That is the
framing for the release announcement (#8).

### 6.2 Tier definitions

Each tier adds one criterion. Reaching tier *k* means every criterion up to *k* is satisfied, so the
tiers are a conjunction and `pm ⟹ gm ⟹ am ⟹ lm`.

| Tier | Criterion | Default accepted | Meaning |
| ---- | --------- | ---------------- | ------- |
| `lm` | `MaxAlleleCredit` | `CREDIT_NONZERO`, `CREDIT_PASS` | some partial credit somewhere |
| `am` | `MaxAlleleCredit` | `CREDIT_PASS` | ≥1 allele fully matched |
| `gm` | `AlleleCount` | `ALLELE_COUNT_EQUAL` | matched count == claimed count |
| `pm` | `PhasingMatch` | `PHASE_CORRECT`, `PHASE_HOMOZYGOUS` | phase correct, or not applicable |

**Why `MaxAlleleCredit` appears twice and `MinAlleleCredit` not at all.** `Min ≤ Max` always, so
`Min ≥ x ⟹ Max ≥ x`: Min is the *stronger* predicate. Assigning Min to the looser tier breaks
monotonicity either way — with equal accepted sets, `lm ⟹ am` unconditionally and the `am` rung
never rejects anything; with a stricter set on `am`, a query `1|1` against a truth `0|1`
(Max = `PASS`, Min = `ZERO`) satisfies `am` while failing `lm`, ranking a genotype error *below* a
weak partial match. Two thresholds on `Max` are monotone instead, and they discriminate on
heterozygotes too, where a single carried haplotype makes `Min == Max`.

Min is then redundant, because `AlleleCount` already encodes it: `calc_gts` is the count that
*matched* and `orig_gts` the count *claimed*, so `AC_ERR_2_TO_2` already means "both copies matched"
and `AC_ERR_1_TO_2` means "one matched, two claimed".

### 6.3 Criterion enums

```
MaxAlleleCredit : CREDIT_ZERO | CREDIT_NONZERO | CREDIT_PASS   (PASS boundary = --credit-threshold)
AlleleCount     : ALLELE_COUNT_LOSS | ALLELE_COUNT_EQUAL | ALLELE_COUNT_GAIN
PhasingMatch    : PHASE_CORRECT | PHASE_INCORRECT | PHASE_UNPHASED | PHASE_HOMOZYGOUS
```

**These are accepted *sets*, not thresholds.** Only `MaxAlleleCredit` is ordered. `AlleleCount` is
categorical with the *middle* value correct — `EQUAL` is right, `LOSS` and `GAIN` are errors in
opposite directions. `PhasingMatch` is correct/incorrect plus two not-applicable values. A minimum
comparison is meaningless for the latter two, and would silently accept `GAIN` whenever `EQUAL` was
requested.

`AlleleCount` derives from the existing 8-value `ac_errtype` (`defs.h:105-112`), which needs an
explicit mapping:

| `ac_errtype` | `AlleleCount` |
| ------------ | ------------- |
| `AC_ERR_1_TO_1`, `AC_ERR_2_TO_2` | `EQUAL` |
| `AC_ERR_1_TO_2` | `GAIN` |
| `AC_ERR_2_TO_1` | `LOSS` |
| `AC_ERR_0_TO_1`, `AC_ERR_0_TO_2`, `AC_ERR_1_TO_0`, `AC_ERR_2_TO_0`, `AC_UNKNOWN` | *no tier* — nothing matched, so `am` was not reached |

`PhasingMatch` derives from existing per-variant state: `PHASE_UNPHASED` from
`phases[vi] == PHASE_NONE`, `PHASE_HOMOZYGOUS` from `orig_gts[vi] ∈ {GT_ALT1, GT_ALT1_ALT1}`, and
`PHASE_CORRECT`/`PHASE_INCORRECT` from `phases[vi]` against `pb_phases[vi]` — the comparison the
writer already performs as `flip_error` (`phase.cpp:75-88`).

`PHASE_HOMOZYGOUS` also covers **haploid** variants, which have no phase to be right or wrong about;
its description should read "not applicable (homozygous or haploid)".

### 6.4 Decision and reporting

`--stringency {lm,am,gm,pm}`, **default `gm`**, sets the tier at which a call is correct: `BD=TP` if
the variant reached the stringency tier, else `FP` on query and `FN` on truth. `BD=N` for D3's
retained-but-unevaluated variants.

**The tier is always computed; `--stringency` only moves the TP/FP boundary.** It does not gate which
criteria are evaluated — a correctly-phased variant reaches `pm` under `--stringency gm` exactly as
it does under `pm`; what differs is whether failing `pm` costs it its `TP`. At the default `gm`, the
`pm` tier is therefore reported without affecting precision/recall unless the user asks for it.

**`BK` carries the true tier, `pm` included.** There is one home for the tier, the spec's own field,
extended by one value. At the default `--stringency gm`:

| Variant | `BK` | `BD` (truth / query) |
| ------- | ---- | -------------------- |
| genotype match, phase correct | `pm` | TP / TP |
| genotype match, phase wrong (flip) | `gm` | TP / TP |
| allele match, genotype wrong | `am` | FN / FP |
| partial credit below threshold | `lm` | FN / FP |
| no credit | `.` | FN / FP |

Under `--stringency pm` every `BK` value above is unchanged; only row 2's `BD` becomes FN / FP. That
is the whole of what the flag does, and the whole of the §6.4.1 irregularity.

**`pm` is a vocabulary extension, and the cost is confined to `extended.csv` bucketing.**
`makeObservationFlags` (`Roc.cpp:49-65`) tests only `lm`/`am`/`gm`, so `pm` falls through to flags
`0`; and `BK` is part of quantify's count-map key (`GA4GHQuantify.cpp:123`), so a novel value opens
new buckets. Neither is an error — nothing validates `BK` against a vocabulary — but note the scope:
because the tier is computed at every stringency, **`BK=pm` lands on the majority of TP records at the
default `gm`**, not only under `--stringency pm`. Two consequences follow, and both are acceptable:

- **Totals are untouched.** `rocEvaluate` branches on `BD` and `BVT` only (§3.5), so every
  precision/recall figure and every `TOTAL`/`TP`/`FP`/`FN` count is exactly what it would be with
  `pm` collapsed to `gm`. This is the property that makes the extension safe rather than merely
  tolerable.
- **`extended.csv` gains `pm` rows**, and the flags-`0` fall-through means a `pm` record contributes
  to neither `FP.gt` nor `FP.al` — which costs nothing, since a `pm` record is a *TP*, and those two
  columns sub-classify FPs. The FP sub-classification of §6.4.1 is unaffected: a phase-error FP
  reaches tier `gm`, not `pm`, so it still carries `BK=gm`.

**No `MT` field.** The tier lives in `BK` alone, not in a second vcfdist-native field. `BK` is the
field the format already has for the match kind, extending it is what a fourth comparison method
requires, and hap.py tolerates the value; a parallel field would duplicate the tier and make `BK` a
lossy view of vcfdist's own decision. The trade is that `BK` is no longer drawn solely from the
published vocabulary — see the documentation requirement in §6.4.1.

Nothing further is needed to explain *why* a variant stopped where it did: the genotype criterion is
already visible in `GE` (`ac_strs`, `globals.cpp:25`) and the phase criterion in `VP`. No
per-criterion tags are added.

#### 6.4.1 How each stringency lands in `quantify`

The tiers were designed against `intermediate.md`'s published `BK → BD` tables and match them
exactly where those tables exist:

| `--stringency` | GA4GH method | Published table | Emitted `(BD, BK)` pairs | `quantify` agrees |
| -------------- | ------------ | --------------- | ------------------------ | ----------------- |
| `lm` | #1 Loose Regional | yes | `lm`/`am`/`gm`/`pm` → TP/TP; `.` → FN/FP | ✓ |
| `am` | #2 Allele Match | yes | `am`/`gm`/`pm` → TP/TP; `.`/`lm` → FN/FP | ✓ |
| `gm` | #3 Genotype Match | yes | `gm`/`pm` → TP/TP; `am`/`lm`/`.` → FN/FP | ✓ **(default)** |
| `pm` | #4 Genotype + Local Phasing | none | `pm` → TP/TP; `gm` → FN/FP | ✓ totals, see below |

`pm` appears at every stringency, since the tier is always computed (§6.4). It is not in any published
table, so **a consumer applying the Method #3 `BK → BD` table to vcfdist output must treat `pm` as it
treats `gm`** — a match. That mapping belongs in D9's docs and the release notes (#8), and is the one
obligation the extension creates for downstream readers.

**Under `pm` the counts are still right; only the FP sub-classification is coarser.** A phase-error
record carries `BD=FP` with `BK=gm`, and `BD` is what all counting runs off (§3.5), so `QUERY.FP` and
every precision/recall figure are correct — the entire point of Method #4. The one consequence is in
the `FP.gt`/`FP.al` breakdown: `makeObservationFlags` sets `OBS_FLAG_GM`, but `addObs`
(`Roc.hh:123-138`) consults only `OBS_FLAG_AM` (→ `FP.gt`) and `OBS_FLAG_LM` (→ `FP.al`), so a
`gm`-flagged FP falls into neither sub-column.

That is a pre-existing category, not something vcfdist introduces: a plain FP carries `BK=.`, which
also sets neither flag, so `FP.gt + FP.al ≤ QUERY.FP` already holds in hap.py's own output. Phase
errors join the residue — arguably correctly, since a phase error is neither a genotype nor an allele
mismatch — but they are then invisible in the sub-columns, recoverable only from `BK`/`VP` in
`summary.vcf` or vcfdist's own phasing outputs. There is no better option: GA4GH has no `BK` value for
"genotype matched, phase wrong", `am` would misreport it as a genotype error and inflate `FP.gt`, and
`.` would claim no match at all. Note the asymmetry this creates with `pm`: the *correct* case gets a
new `BK` value, while the *error* case has to reuse `gm`. That is the right way round — a new value
for a passing record cannot be mistaken for a different verdict, whereas a new value for a failing one
would need every consumer to know it means "not a match" before its FPs counted correctly.

**`--stringency pm` emits a `WARN`** noting that phase-error FPs are counted in `QUERY.FP` but not
sub-classified in `FP.gt`/`FP.al`, that Method #4 has no published `BK` mapping so a consumer applying
the Method #3 table will read a phase error's `BK=gm` as a match, and that `--no-require-phase` is
available if unphased input is being penalised unintentionally (§6.6). The output stays fully
parseable, so nothing downstream would surface any of this on its own.

### 6.5 Known gaps in this phase

**`ac_errtype` is query-only, so the `gm` tier has no truth-side input.** `set_allele_errtype` reads
`calc_gts` (`variant.cpp:171-195`), and `calc_gts` is only ever written for query — every assignment
in `phase.cpp:383-425` is on `qvars`, and `variant.h:83` documents it as "only set for query". The
truth sample needs the mirrored comparison. This is a concrete work item, not a detail: without it
the truth-side `BD` cannot distinguish a genotype error from a miss.

**`pm` sees flip errors but not switch errors.** `PHASE_CORRECT` compares a variant's phase against
its phase-block state, which is the flip test. After a switch error, downstream variants are
consistent with the *new* block state and so all read `PHASE_CORRECT` while being globally
mis-phased. `--stringency pm` therefore penalises flips and is blind to switches. If Method #4 is
meant to capture both, `pm` needs a block-level input; scoped here as a documented limitation with
its own issue.

### 6.6 CLI

**Two flags are exposed.** The four accepted-value sets of §6.2 remain the internal model, so the
tier logic stays uniform and every row of §6.2/§6.3 is directly testable — but only `pm`'s is wired
to the CLI, as a negatable boolean matching the `--no-output-files` convention already in the
codebase:

```
--stringency {lm,am,gm,pm}     default: gm    tier at which a call counts as correct
--require-phase                default: on    PHASE_UNPHASED does not satisfy the pm tier
--no-require-phase                            PHASE_UNPHASED satisfies the pm tier
```

**Why the credit and genotype sets are not exposed.** Every alternative value is either redundant
with a flag that already exists or breaks the ladder:

| Tier | Alternative set | Effect |
| ---- | --------------- | ------ |
| `lm` | add `CREDIT_ZERO` | tier becomes vacuous; `--stringency lm` gives precision = recall = 1 |
| `lm` / `am` | swap their sets | the two rungs collapse into one, or invert (§6.2) |
| `gm` | `{LOSS, EQUAL, GAIN}` | forgives genotype errors — already `--stringency am` |
| `gm` | `{EQUAL, GAIN}` | asymmetric; produces figures nobody else can reproduce |

For the credit tiers the knob users actually want is the threshold, and that already exists as
`--credit-threshold`. Suppressing the rest also removes the need for start-up validation: exposing
both credit sets would let a user invert the ladder and break every assumption downstream of §6.2,
whereas `--require-phase` cannot, since `pm` is the top tier and its criterion is orthogonal to those
below — `pm ⟹ gm ⟹ am ⟹ lm` holds for either setting.

**What `--require-phase` changes, precisely.** It only ever moves a genotype-matched *unphased*
variant between the `gm` and `pm` tiers:

| | `--require-phase` (default) | `--no-require-phase` |
| - | -------------------------- | -------------------- |
| `BK` for a genotype-matched unphased variant | `gm` | `pm` |
| `BD` at `--stringency lm`/`am`/`gm` | unchanged | unchanged |
| `BD` at `--stringency pm` | FN / FP | TP / TP |

So it affects `TP` **only at `--stringency pm`**; at every other stringency it is reporting-only,
visible in `BK`. That is what makes it safe on by default — the strict reading is recorded without
ever costing a call its `TP` unless Method #4 was explicitly requested. It matters most once D6 lands
and `PHASE_UNPHASED` goes from rare to common.

## 7. Native summary metrics

The match tier (§6) computes signals vcfdist's own TSVs do not currently report — the same signals a
GA4GH quantifier would derive from `summary.vcf`. Reporting them natively is what keeps vcfdist's
metrics complete without a quantifier in the loop. Three are cheap and land here; three larger ones
are deferred to §7.2.

**These columns are stratified by [`D5`](./D5_multi-bed-stratification.md).** `FP_GT`, `FP_AL`, and
`QUERY_UNK` are counts and take the stratum axis exactly as `TRUTH_TP` and `QUERY_FP` do; `FRAC_NA`
and the two ratios are derived per stratum from those counts. This deliverable adds them as
genome-wide columns; D5 §7.3 widens them.

### 7.1 Added to `precision-recall-summary.tsv`

| Column | Definition | Source |
| ------ | ---------- | ------ |
| `FP_GT` | query FPs where an allele matched but the genotype did not | tier `am` (§6.2) |
| `FP_AL` | query FPs with only a partial match, below the credit threshold | tier `lm` (§6.2) |
| `QUERY_UNK` | query variants retained but not assessed | `BD=N` (D3) |
| `FRAC_NA` | `QUERY_UNK / QUERY_TOTAL` | derived |
| `TRUTH_TITV`, `QUERY_TITV` | transition/transversion ratio, per callset | SNP REF/ALT |
| `TRUTH_HET_HOM`, `QUERY_HET_HOM` | het:hom ratio, per callset | `orig_gts` + `ploidies` (D3 §4.4) |

`FP_GT` and `FP_AL` are the two columns hap.py exposes as `FP.gt` and `FP.al`, and the mapping is
already fixed by §3.5: `makeObservationFlags` reads `am` → `FP.gt` and `lm` → `FP.al`. Emitting them
natively is what makes a GA4GH export checkable on these columns rather than only on totals. **They
do not sum to `QUERY_FP`** — an FP with `BK=.` matched nothing and belongs to neither, and under
`--stringency pm` a phase-error FP carries `BK=gm` and also lands in neither (§6.4.1). So
`FP_GT + FP_AL ≤ QUERY_FP`, exactly as in hap.py's own output; the documentation must say so, since a
reader will otherwise treat the shortfall as a bug.

`QUERY_UNK` and `FRAC_NA` are the same quantity as a count and a fraction, and both are needed.
`QUERY_UNK` makes the row reconcile — `QUERY_TP + QUERY_FP + QUERY_UNK = QUERY_TOTAL` — while
`FRAC_NA` is the trust signal, because **precision excludes unassessed calls from its denominator**.
A run at `FRAC_NA = 0.4` reports precision over 60% of the query calls and is not comparable to one at
`0.01`; the fraction is what catches a mis-specified `--bed`. **This depends on D3.** vcfdist filters
on `--bed` at parse time (§3.4), so absent D3's retention there is nothing to count and `QUERY_UNK`
would be structurally zero — misreporting "nothing was skipped" rather than "skipping is not tracked".

The two ratios are pure QC and affect no other figure. `HET_HOM` earns its place by construction: it
is the direct symptom of the record-splitting this deliverable removes, since a per-haplotype split
reports every homozygote as two heterozygotes. It is only meaningful once §4's per-site counting and
D3 §4.4's per-variant ploidy are in place.

`FP_GT`, `FP_AL`, and `QUERY_UNK` are counts and also go in `precision-recall.tsv`, whose rows are a
`MIN_QUAL` sweep. The ratios stay in the summary only — they are composition checks on a callset, and
recomputing them per quality threshold invites reading noise as signal.

### 7.2 Deferred, with the enabling work noted

Recorded as future options rather than omissions, since each has a natural home:

| Option | What it adds | Belongs to |
| ------ | ------------ | ----------------------- |
| Stratified long-format summary | hap.py's `extended.csv` axes: indel size bins (`I1_5`/`I6_15`/`I16_PLUS`, `D…`), genotype (het/homalt/hetalt), per-subset rows, and region sizes for per-Mb rates | The subset axis is **D5**'s; vcfdist has only SNP/INDEL/SV at `--sv-threshold` |
| ROC on an alternate score field | `--score-field`, matching vcfeval's `--vcf-score-field` and hap.py's `--roc`; callers often rank better on `GQ` or a model score than `QUAL` | **D3** retains the INFO/FORMAT fields that make it reachable |
| Single JSON summary | one machine-readable file (≈ hap.py's `metrics.json.gz`) instead of N TSVs | **D10** — MultiQC modules much prefer one file |

The stratified table is the largest of the three and the reason it is not attempted here: adding those
axes turns `precision-recall-summary.tsv` from a wide table into a long one keyed
`(VAR_TYPE, SUBTYPE, SUBSET, GENOTYPE, MIN_QUAL)`. That is a bigger break for existing parsers than
D4's other output changes combined, and it should not ride along with them.

**Separated `tp`/`fp`/`fn` VCFs are declined outright**, not deferred. They are vcfeval's most-used
convenience output, but `BD` already carries the decision, so they are one `bcftools` filter away.
The right response is documentation: **D9's docs and the README should show the one-liners** —

```bash
bcftools view -i 'FORMAT/BD=="FP"' summary.vcf   # query false positives
bcftools view -i 'FORMAT/BD=="FN"' summary.vcf   # truth false negatives
bcftools view -i 'FORMAT/BD=="TP"' summary.vcf   # true positives, both samples
```

— together with the two facts that make them correct: `FP` is only ever assigned to `QUERY` and `FN`
only to `TRUTH` (`variant.cpp:514-522`), so neither needs a sample subscript, while `TP` is assigned to
both and takes `-s QUERY` / `-s TRUTH` to pick a side. Worth adding to the README's usage section
because the absence of these files is the most likely reason a vcfeval user assumes vcfdist cannot
produce them.

## 8. Data flow

```
parse_variants ──▶ cluster ──▶ prec_recall_aln ──▶ phase
                                     │                │
                          credit[hi][vi]      phases[vi], pb_phases[vi]
                                     │                │
                                     └──── match tier (§6.2) ────┐
                                                                 │
                    ┌── per-site BD/BK ◀── --stringency ◀────────┘
                    │
  write_summary_vcf ┤── INFO/BS   ◀── superclusters[vi]   (the only INFO field, §5.1.1)
                    ├── Number=.  ◀── BC, RD, QD, SG  (FORMAT, GT-allele order)
                    └── GT        ◀── orig_gts[vi] + ploidies[vi]  (not calc_gts)

  write_precision_recall ──▶ one count per variant per callset (§4)
                         └──▶ FP_GT/FP_AL ◀ match tier; QUERY_UNK/FRAC_NA ◀ BD=N; ratios (§7.1)

  summary.vcf ──▶ vcfdist's own consumers: D10 MultiQC module, D11 nf-core/snakemake, D7 harness
              └──▶ [user-side, not run here] bgzip + tabix ──▶ any GA4GH quantifier
```

## 9. Testing plan

Conformance is asserted against **our own output**, since no quantifier is run (§1.3). Each item below
is checkable with `bcftools` and the existing integration harness.

Record shape, per-variant ploidy, and per-haplotype value counts are
[`D3`](./D3_retain-info-format-fields.md) §9's tests and are not repeated here; the cases below assume
they pass.

- **Structure:** `bcftools view` the output clean, including a haploid contig; no version floor
  (§5.2).
- **`BD`/`BK` are single-valued:** exactly one value each on every record — the assertion that this
  deliverable actually collapsed D3's `Number=.` lists — and their `##FORMAT` lines declare
  `Number=1`, while `BC`/`RD`/`QD`/`SG` still declare `Number=.`.
- **`INFO/BS`:** monotonically non-decreasing within a contig, and all records of one supercluster
  adjacent, since a quantifier splits work blocks on `BS` boundaries.
- **Tag renames (§5.1):** neither `SC` nor a FORMAT-column `BS` appears in the header or the records;
  `INFO/BS` carries the value `FORMAT/SC` used to, and `FORMAT/PBS` the value `FORMAT/BS` used to,
  asserted on the same fixture so the swap cannot half-land.
- **Columns (§5.1.1):** `BS` is the *only* `INFO` field, and the FORMAT key set is `SC`-free with
  `PBS` in place of `BS` — `PB`, `PBS`, `VP`, and `FE` all still in FORMAT, and no `MT` (§6.4). Pinned
  so a later "tidy" does not migrate the phasing fields to `INFO`, where a per-sample value cannot be
  expressed at all.
- **`SG` stays per-sample (§5.1.1):** on a merged record whose computed genotype is swapped, TRUTH and
  QUERY carry the same `SG` values in *opposite* order — asserted, since that ordering is the reason
  a shared ID space still cannot become one `INFO` list.
- **Match tiers:** one record per row of the §6.2 and §6.3 tables, including each `ac_errtype`
  mapping and each `PhasingMatch` value, asserted against `BK`. One case per `--stringency` level
  confirming `BD` changes only where the tier boundary moves, and that `BK` does not change with
  stringency at all. A haploid variant takes `PHASE_HOMOZYGOUS` and is not scored as a phasing error
  under `--stringency pm`.
- **`BK=pm` at the default stringency:** a correctly-phased genotype-matched variant reports
  `BK=pm` with `BD=TP` under `--stringency gm` — the case that makes `pm` common rather than rare, so
  it is asserted at the default and not only under `--stringency pm`.
- **`--require-phase`:** the §6.6 table asserted directly — a genotype-matched unphased variant
  reports `BK=gm` by default and `BK=pm` under `--no-require-phase`, with `BD` identical at
  `lm`/`am`/`gm` and differing only at `pm`.
- **`--stringency pm`:** emits the §6.4.1 `WARN`; a flip-error record is counted in `QUERY_FP` and in
  neither `FP_GT` nor `FP_AL`. Also assert the baseline that makes this benign — a plain `BK=.` FP is
  likewise in neither sub-column, so `FP_GT + FP_AL ≤ QUERY_FP` at every stringency.
- **Sample invariants:** columns named `TRUTH`, `QUERY`, in that order; absent sample `GT=.`. A record
  carrying a preserved non-`PASS` input `FILTER` is emitted as-is, with the demotion caveat documented
  rather than tested (§3.3).
- **Counting (§4):** a fixture with known hom/het/het-alt composition, asserting the per-site totals in
  `precision-recall-summary.tsv` directly. A homozygous variant contributes **1**, not 2 — the
  regression test for the convention change — and the het-alt residual is asserted at its known value
  rather than left implicit (§5.3).
- **`BD=N` inertness:** D3's retained variants change no `TP`/`FP`/`FN` count in any TSV — but they
  *do* move `QUERY_UNK` and `FRAC_NA` (§7.1), which is the one intended effect and must be asserted
  rather than assumed absent.
- **New summary columns (§7.1):** a fixture with a genotype error and a below-threshold partial match
  asserts `FP_GT` and `FP_AL` land in the right column and that `FP_GT + FP_AL ≤ QUERY_FP` holds —
  including the two cases that make it a strict inequality, a `BK=.` FP and, under `--stringency pm`, a
  phase-error FP. `QUERY_TP + QUERY_FP + QUERY_UNK == QUERY_TOTAL` on every row, and `FRAC_NA` matches
  the ratio of those columns. `HET_HOM` on a known hom/het fixture pins §4's counting change: it would
  have read 2× under per-haplotype counting.
- **Conformance-by-inspection, recorded not automated:** the tag set of §3.1/§3.2 is checked against
  the consuming source at implementation time and the source revision is recorded in the PR. This is
  the compensating control for not running a quantifier (§1.3) — a re-read on the pinned revisions,
  not a round-trip.

## 10. Out of scope

| Item | Why | Where |
| ---- | --- | ----- |
| The `summary.vcf` record shape (un-splitting, ploidy, `Number=.` fields) | Owned by D3, which rewrites the writer | §5, D3 §4.2–§4.5 |
| Multiallelic record merging | D3 §4.3's decision; the het-alt count residual is reported rather than removed | §5.3 |
| Polyploid input (ploidy > 2) | `ERROR` retained | D3 §4.4 |
| Switch errors in the `pm` tier | Needs a block-level input, not per-variant phase | §6.5, own issue |
| `som.py` compatibility | Somatic; makes its own decisions | §2 |
| Stratification in the VCF | Supplied externally to quantify; D5's concern | §3.3 |
| Stratified long-format summary (indel size bins, genotype, per-subset rows) | Changes the TSV from wide to long — a bigger parser break than D4's other output changes combined; subset axis is D5's | §7.2 |
| ROC on a score field other than `QUAL` | Needs D3's retained INFO/FORMAT fields to reach the score | §7.2 |
| Single JSON summary file | D10's MultiQC module is the consumer that wants it | §7.2 |
| Separated `tp`/`fp`/`fn` VCFs | Declined, not deferred — `BD` already carries the decision, so a `bcftools` filter suffices; documented in the README and D9 instead | §7.2 |

## 11. Decision record

| Decision | Section | Rationale |
| -------- | ------- | --------- |
| Target what `quantify` reads; agreement with it, not spec purity, is the standard | §1.3, §3.4, §6.4.1 | Nothing but `quantify` consumes VCF-I, and the implementations diverge from the spec in three places |
| Record shape delegated to D3 | §5 | D3 rewrites the writer for field retention anyway; retaining the per-haplotype split only for D4 to remove it would mean writing the field-duplication logic twice |
| `BD`/`BK` collapse from D3's `Number=.` lists to `Number=1` here | §5.4, §6 | The aggregation policy *is* the match tier; separating it from the mechanical un-split keeps each deliverable reviewable on its own |
| Delete the `AC_ERR_2_TO_1` correction before stratification lands | §4 | D5 adds a stratum axis to the same loops; deleting first means the per-stratum version of a doomed special case is never written |
| The §7.1 columns are added genome-wide here and stratified by D5 | §7.1 | They are counts, so they take the stratum axis exactly as the existing count columns do; D5 §7.3 owns that widening |
| Emit exactly `GT`, `BD`, `BK`, `QQ`, `INFO/BS` | §3.1 | Verified complete against both reference producers |
| `INFO/BS` from the supercluster index, as `Number=.` | §3.1 | Not optional in practice — truth-side ROC and boundary flagging depend on it; both implementations declare `.` despite the spec's `Number=1` |
| Do not emit `BVT`/`BLT`/`BI`/`Regions`/`VTC` | §3.2 | Computed by quantify and overwritten; `BVT`/`BLT` headers are appended unguarded, so emitting ours duplicates them |
| Adopt per-site counting | §4 | Divergent denominators would swamp D7's concordance analysis and make the round-trip a re-denomination; also deletes the `AC_ERR_2_TO_1` special case |
| Per-haplotype tags as `Number=.`, not `Number=P` | D3 §4.5, §5 | `P` is semantically correct but needs htslib ≥ 1.23 to read; these fields are informational and D11's adoption path pins older htslib. Byte-identical output; only validation is given up, which §9 supplies |
| `##fileformat` stays `VCFv4.2` | §5.2 | The bump existed only for `Number=P`; nothing else here is post-4.2 |
| Swap the two names: `SC` → `INFO/BS`, and `BS` (Block Phase) → `PBS` | §5.1 | `BS` is GA4GH's name for the superlocus, which is exactly what `SC` holds; the incumbent `BS` must vacate the name rather than coexist as a second meaning in one file |
| `SC` is the only field that changes column; `PB`, `PBS`, `VP`, `FE` stay in FORMAT | §5.1.1 | Phase is a property of the *variant*, not the locus — VCF puts it in FORMAT (`GT`'s `\|`, and `PS` as a reserved FORMAT key, itself a locus-spanning grouping the spec still makes per-sample). `PB` is the output analogue of `PS`, so moving it while keeping `PS` would be inconsistent. A field being `.` on TRUTH is what a per-sample field looks like when that sample has no value, not evidence of a locus property — INFO has no per-sample slot to express the difference. `SC` differs because it indexes a grouping *across* both callsets, meaning the same thing from either side |
| `SG` stays in FORMAT despite sharing an ID space across callsets | §5.1.1 | It is per-haplotype, and query's list is re-ordered to query-haplotype order at `phase.cpp:389`, so TRUTH and QUERY hold the same values in different orders. One INFO list cannot carry both, and D3 §4.5 pins per-hap order to the emitted `GT` |
| `GE` stays in FORMAT | §5.1.1 | Query-only only because `ac_errtype` is unset for truth; §6.5 closes that, after which truth's value is the inverse of query's |
| Multiallelic stays split — settled, not pending D3 | §5.3, §10 | Provenance alone is insufficient: merging needs a fallback to un-normalized source alleles and a cross-entry credit rule. Residual reported instead |
| `MaxAlleleCredit` at both `lm` and `am`; no `MinAlleleCredit` | §6.2 | `Min ≥ x ⟹ Max ≥ x`, so Min at a looser tier makes the ladder non-monotone; `AlleleCount` already encodes Min |
| Criteria are accepted sets, not thresholds | §6.3 | `AlleleCount` is categorical with the middle value correct; `PhasingMatch` has two not-applicable values |
| `PHASE_HOMOZYGOUS` also covers haploid | §6.3 | Both mean "phase not applicable"; avoids a fourth not-applicable value |
| Default `--stringency gm` | §6.4 | Matches hap.py's default and GA4GH Method #3; keeps `pm` reported but not scoring |
| `BK` carries the true tier including `pm`; no `MT` field | §6.4 | `BK` is the field the format already has for the match kind, and Method #4 needs a fourth value. Nothing validates `BK`, and `rocEvaluate` counts off `BD`, so totals are identical to a `pm`-collapsed export; the cost is confined to extra `extended.csv` buckets. A second field would have duplicated the tier and made `BK` a lossy view of vcfdist's own decision |
| A phase-error FP keeps `BK=gm` rather than getting its own value | §6.4.1 | A new value on a *passing* record cannot be misread as a different verdict; a new value on a *failing* one would need every consumer to know it means "not a match" before its FPs counted correctly |
| No per-criterion tags added | §6.4 | `GE` already exposes the genotype criterion and `VP` the phase criterion |
| Add `FP_GT`/`FP_AL`, `QUERY_UNK`/`FRAC_NA`, and TiTv/het-hom ratios natively | §7.1 | Otherwise hap.py fed by vcfdist's own export reports metrics vcfdist itself does not. The signals already exist internally: the match tier, D3's `BD=N`, and the genotypes |
| `FP_GT + FP_AL ≤ QUERY_FP`, documented rather than forced to sum | §7.1 | A `BK=.` FP matched nothing and a phase-error FP carries `BK=gm`, so neither belongs to a sub-column — the same strict inequality hap.py's own output has |
| Both `QUERY_UNK` and `FRAC_NA`, not one | §7.1 | The count makes the row reconcile; the fraction is the trust signal, since precision excludes unassessed calls from its denominator |
| Stratified long-format summary, alternate ROC score field, and JSON summary deferred | §7.2 | Each has a natural owner — D5, D3, D10 — and the long-format change alone is a bigger parser break than the rest of D4 |
| Separated `tp`/`fp`/`fn` VCFs declined, documented as `bcftools` one-liners | §7.2 | `BD` already carries the decision; `FP` is QUERY-only and `FN` TRUTH-only, so the filters need no sample subscript |
| Only `--stringency` and `--require-phase` exposed; the four sets stay internal | §6.6 | Every alternative value is redundant with `--credit-threshold`/`--stringency` or breaks monotonicity; suppressing them removes the need for start-up ladder validation, while keeping each §6.2/§6.3 row testable |
| `--require-phase` is a negatable boolean, default on | §6.6 | Matches `--no-output-files` already in the codebase, keeps the §6.2 default set, and affects `TP` only at `--stringency pm` |
| Each change lands as its own commit | — | Four changes here can each move P/R; a squash makes them indistinguishable when D7's numbers shift |
