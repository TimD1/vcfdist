# vcfdist v3.0.0 — Unphased Variant Evaluation Design

|  |  |
| :-- | :-- |
| **Version:** | 0 (Draft) |
| **Authored By:** | Tim Dunn |
| **Status:** | Draft — design of record for Deliverable #6 |
| **Companion docs:** | `D0_vcfdist-v3-SOW.md`, [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md), [`D2_vcfdist-v3-unit-tests.md`](./D2_vcfdist-v3-unit-tests.md) |
| **Companion designs:** | [`D3`](./D3_retain-info-format-fields.md) (record shape, `BD=N`), [`D4`](./D4_ga4gh-compatibility.md) (per-site counting, match tiers), [`D5`](./D5_multi-bed-stratification.md) (stratum axis, phasing-output split). This deliverable's metric changes land in D5's split, stratified phasing files. |
| **Implements:** | Deliverable #6 (SOW), [`D1`](./D1_vcfdist-v3-design-doc.md) §5.6 |
| **Issue:** | [#46](https://github.com/TimD1/vcfdist/issues/46) (D6: Unphased Variant Evaluation) |

---

## 1. Purpose & Scope

This document is the design for unphased variant evaluation, called for by **Deliverable #6** and sketched in [§5.6 of the design doc](./D1_vcfdist-v3-design-doc.md#56-deliverable-6--unphased-variant-evaluation). The SOW names this the highest research-risk item in the project; the analysis in §2 below shows why that risk is smaller than it appears — the v3 alignment core is already phase-blind — and where the real work actually sits, which is in parsing, clustering, and metric definitions.

**In scope.** Removing the hard local-phasing requirement for the **query** VCF: unphased and partially-phased query calls are evaluated for precision/recall exactly as phased calls are, and phasing metrics are reported over the subset of query variants that carry usable phase, with the phased fraction of the callset reported explicitly.

**Out of scope.** Unphased **truth** VCFs (§9.1), the graph-based clustering redesign (§9.2), and overlapping **truth** variants (§9.3). All three are specified far enough to show that this deliverable does not foreclose them, and §9.1 and §9.3 share a root cause: the truth side of the alignment graph is a linear per-haplotype chain.

---

## 2. Where phase is load-bearing today

The premise this design rests on is that v3's precision/recall is already independent of query phasing. That is very nearly true, and the exceptions matter, so both are established here with evidence.

### 2.1 The alignment core is phase-blind on the query side

`Graph::Graph` builds the query side as a DAG containing *every* unevaluated query variant in the supercluster, with no haplotype filter whatsoever (`dist.cpp:912-950`). Only the truth side is a linear per-haplotype chain, gated on `tvars->var_on_hap(var_idx, truth_hap)` (`dist.cpp:1003`). `evaluate_variants` then runs each truth haplotype independently, and `errtypes` are indexed per haplotype, so both truth haplotypes see the complete query variant set. A query variant's `calc_gt` is *derived* from which truth haplotype it matched (`dist.cpp:555`), never read from the input genotype.

The precision/recall tally reads `orig_gt` only to map results back to the input representation:

```cpp
int calc_hi = hi ^ qvars->calcgt_is_swapped(vi);      // print.cpp:290
if (qvars->var_on_hap(vi, hi)) { ... }                // print.cpp:292
```

For a heterozygous variant exactly one `hi` satisfies `var_on_hap`, and the XOR undoes the orientation, so the totals are invariant to whether the input said `0|1` or `1|0`. **Query phase orientation is used nowhere in the alignment, credit, or tally math.**

*Code shape note:* [`D4`](./D4_ga4gh-compatibility.md) §4 replaces this per-haplotype tally with one unit per variant per callset, so by the time this deliverable is implemented the loop above no longer exists in this form. The conclusion is unaffected and in fact strengthened — a per-site decision cannot depend on which haplotype an allele was written on, so there is no orientation left for query phase to influence. The snippet is retained because it is the evidence for the claim in the code as it stands today.

### 2.2 The two exceptions

1. **Clustering is per-haplotype.** `main.cpp:69` and `:86` spawn `HAPS × contigs` `wf_swg_cluster` threads over `variants[hap][ctg]`; reaches are computed per haplotype and merged afterwards (`cluster.cpp:238-317`). Cluster reaches determine supercluster boundaries, which determine *which variants are co-evaluated*, which affects credit assignment. P/R is invariant to phase given fixed superclusters — but superclusters are not fixed under a change in per-haplotype variant layout. §4.2 removes this dependence.

2. **The phasing DP feeds back into genotype assignment.** `fix_allele_counts` breaks credit *ties* using `pb_phases`, the output of the phasing DP (`phase.cpp:406-413`, `:421-428`). That sets `calc_gt`, which drives `calcgt_is_swapped`, which selects which haplotype's `errtype` is counted. On exact credit ties, therefore, the phasing result influences P/R. This is rare but real, and is documented rather than removed — the tie-break must pick *something*, and current phasing is a defensible choice.

### 2.3 What actually blocks unphased input

`parse_variants` discards unphased heterozygous calls before any of the above runs:

```cpp
// skip unphased heterozygous variants (1/1 is allowed, 0/1 is not)
if (ngt == 2 && !same && !bcf_gt_is_phased(gt[HAP2])) {   // variant.cpp:1031
    ...
    unphased_gt_total += 1;
    continue;
}
```

This applies to both callsets and also discards unphased `1/2` compound heterozygotes. Today an unphased query therefore evaluates to approximately zero heterozygous variants. Removing this skip is the substance of the deliverable.

---

## 3. Design decisions

| Decision | Resolution | Rationale |
| :-- | :-- | :-- |
| Callset scope | Query only; truth must remain locally phased | Truth graph is a linear per-hap chain by construction (`dist.cpp:1003`); GIAB/HPRC truth sets are phased. See §9.1. |
| Clustering | Per-callset on the merged variant list; merge across haps *before* clustering | Removes the arbitrary haplotype assignment entirely rather than choosing one. §4.2. |
| Cluster iterations | Fixed at 1; `-i` closed off | Becomes a correctness precondition, not a tuning knob, once clustering runs on a merged list. §4.2.3. |
| Phased-status tracking | New `is_phased` field on `ctgVariants`, surfaced per-variant in `summary.vcf` | `phase_sets` cannot serve as the proxy — it is overwritten by backfill. §4.1.2. |
| Overlapping variants | Filter removed for `QUERY`, retained for `TRUTH` | The query graph already represents overlaps as mutually exclusive alternative paths; the truth chain is linear and silently corrupts. §4.1.3. |
| Phasing-rate denominator | Adjacent phaseable-heterozygote pairs, per phase block | Follows from switches being inter-variant regions. §5.1. |
| Summary output | Explicit denominator columns alongside counts; both variant- and region-based phased fractions | Lets any consumer recompute any convention. §5.3. |

---

## 4. Implementation

### 4.1 Parsing (`variant.cpp`)

#### 4.1.1 Admit unphased heterozygous query calls

Remove the skip at `variant.cpp:1030-1038` (the guard itself is at `:1031`) for the query callset. For the truth callset the skip is retained (§9.1), but its warning is escalated: `unphased_gt_total` is currently reported only through the aggregate warning at `:1208-1210`, with per-variant detail behind `verbosity > 1`. For truth, an unphased heterozygote now indicates the input is unusable for phasing analysis, so the summary warning must state that phasing metrics will be degraded, not merely that variants were skipped.

No change is needed to haplotype placement. The existing loop iterates `hap` over `[0, |ngt|)` and selects the allele via `bcf_gt_allele(gt[hap])` (`variant.cpp:1011-1028`), so allele *k* lands on haplotype *k*. An unphased `1/2` compound heterozygote therefore places its two alternate alleles on opposite haplotypes exactly as a phased `1|2` would, and needs no special handling.

#### 4.1.2 Add per-variant phased status

Add `std::vector<bool> is_phased` to `ctgVariants`, populated by both `add_var` overloads (the copy overload must propagate it, alongside the other per-variant fields at `variant.cpp:53-58`). Set it from `bcf_gt_is_phased` at parse time, with homozygous and haploid calls recorded as phased — they carry no phase ambiguity — so that `is_phased == false` means specifically "heterozygous and unphased".

This field must be distinct from `phase_sets`. `fix_phase_set_tags` overwrites zero-valued phase sets with the running phase set (`phase.cpp:328-331`), which destroys any signal carried in `phase_sets` by the time `phase()` runs; that overwrite is deliberate and must be preserved (§4.3.2). The `// TODO: only set phase sets for 1|1 variants (for when we add unphased eval)` at `phase.cpp:329` is resolved by this separation rather than by conditionalizing the backfill.

`load_and_merge_callset_vars_across_haps` must carry `is_phased` through the merge. A merged homozygous record (both haplotypes agree, `cluster.cpp:134-150`) is phased by definition; a merged heterozygous record inherits the flag from the contributing haplotype.

#### 4.1.3 Remove the overlap filter for the query callset

`variant.cpp:1131-1139` skips any variant overlapping the previous one *on the same haplotype*, for both callsets. Its purpose is to guarantee that downstream consumers receive a well-formed linear haplotype. §4.2 retires that requirement for clustering, and what remains of it differs by callset.

**The query side already supports overlapping variants and the filter is purely conservative.** In `Graph`'s query construction (`dist.cpp:912-950`), `ref_pos` advances only via the reference-filling loop, and only as far as `var_pos`; adding a variant node does not advance it, because the variant's end is pushed onto the `qnode_ends` priority queue instead (`:948-949`). Since variants are position-sorted, `ref_pos <= var_pos` holds invariantly and the `assert(ref_pos == var_pos)` at `:941` cannot fire.

Overlapping query variants therefore become parallel branches in the DAG. For a deletion at position 100 of reference length 10 followed by a substitution at 105: processing the substitution, `qnode_ends.top() == 110` is not `< 105`, so a reference node `[100, 105)` is emitted and `ref_pos` becomes 105. The deletion node `[100, 110)` and the reference node `[100, 105)` both begin at 100, and since connectivity is coordinate-based (`qbegs[n1] == qends[n2]`, `:983-987`), **no path can traverse both** — 110 does not equal 105. That is the correct semantics: overlapping calls are mutually exclusive alternatives and the aligner selects whichever matches the truth.

**The truth side cannot support them, and fails silently.** `Graph`'s truth construction (`dist.cpp:999-1024`) builds a linear chain: it emits a reference node only when `var_pos > ref_pos` (`:1006`), and it *does* advance `ref_pos` past each variant (`:1022`). With the same two variants on one truth haplotype, the deletion sets `ref_pos = 110`; the substitution's `105 > 110` test fails so no reference node is emitted, yet its alternate allele is still appended to `this->truth`; the trailing segment then resumes at `ref_pos = 106` (`:1028`), re-emitting reference bases the deletion already removed. The constructed truth haplotype does not exist, and there is no assertion, error, or warning. The overlap filter is currently the only thing preventing this.

**Accordingly, this deliverable removes the filter for the `QUERY` callset entirely — phased and unphased alike — and retains it for `TRUTH`.** Conditioning on callset rather than on `is_phased` is both safer and simpler to document, and it recovers overlapping phased query calls that are silently dropped today. Truth-side support requires making the truth side a DAG, which is out of scope (§9.3).

This matters concretely for unphased input even though it is not specific to it: because callers overwhelmingly emit `0/1` rather than `1/0`, newly-admitted unphased heterozygotes concentrate on HAP2, so a heterozygous substitution falling inside a heterozygous deletion — which phased input would have placed in trans on opposite haplotypes — would otherwise collide and be dropped with only a verbosity-2 warning.

Two consequences to carry through:

- The `simple_gt` downgrade at `variant.cpp:1141-1144`, which converts a homozygous call to heterozygous when its counterpart on the other haplotype was skipped by the filter, becomes unreachable for query and must be retained only for truth.
- Admitting overlaps anywhere depends on clusters remaining singletons, since `generate_str` aborts the program on overlapping variants within its range (`dist.cpp:120-130`). This is a second, independent reason for the iteration lock in §4.2.3.

### 4.2 Clustering (`cluster.cpp`, `main.cpp`)

#### 4.2.1 Move the haplotype merge ahead of clustering

Today the pipeline is: cluster per `(hap, contig)` → merge haplotypes and merge their clusters → supercluster. It becomes: merge haplotypes → cluster per `(callset, contig)` → supercluster.

`wf_swg_cluster` changes signature from `(variantData*, ctg_idx, hap, ...)` to operate on a single merged `ctgVariants`, and the thread fan-out at `main.cpp:69` and `:86` drops from `HAPS × contigs` to `contigs` per callset.

The cluster-merging half of `load_and_merge_callset_vars_across_haps` (`cluster.cpp:238-317`) is deleted outright, since clusters no longer exist at the time the merge runs. That removes the cross-haplotype reach heuristic, including the acknowledged `// TODO: why no std::min here?` at `cluster.cpp:257` — a line that assigns `next_left_reach` from one haplotype without taking a minimum, and which appears to be a genuine defect. The variant-merging half (`:123-236`), which sets `orig_gt` and detects homozygotes, is retained and now runs against unclustered per-haplotype lists.

`sort_superclusters` estimates per-supercluster RAM from per-haplotype lengths via `var_on_hap` on `orig_gt` (`cluster.cpp:592-599`). This remains correct — it is estimating the cost of the *alignment* stage, which is still per truth haplotype — and is unchanged.

#### 4.2.2 What this fixes

- **Unphased variants need no haplotype to be clustered.** Because clustering now operates on the merged list, the arbitrary placement of an unphased heterozygote no longer perturbs cluster reaches, supercluster boundaries, or which variants are co-evaluated.
- **Overlapping variants become representable in the clustered list**, which is what permits §4.1.3.
- **The cross-haplotype reach heuristic disappears**, along with its suspected defect.

#### 4.2.3 Iterations remain fixed at 1

`max_cluster_itrs` defaults to `1` (`globals.h:44`) and `-i` is commented out of the usage text (`globals.cpp:494`). Tracing the loop at `cluster.cpp:666-670`: pass 1 runs with `prev_clusters[i] = i`, so **every cluster is a single variant**; the merge passes then run; `iter = 2 > 1` breaks. Every reach computed in production today is therefore computed for one variant applied to the reference.

Under the restructure this stops being a default and becomes a precondition, for two independent reasons:

1. **`generate_str` aborts the program on overlapping variants.** With a deletion at position 100 of reference length 10 followed by a variant at position 105, the deletion advances `ref_pos` to 110; the next iteration finds `ref_pos != poss[var_idx]`, computes `ref_end = min(end_pos, 105) = 105`, hits `ref_end < ref_pos`, and calls `ERROR(...)` → `std::exit(1)` (`dist.cpp:120-130`). A singleton range cannot self-overlap, so one pass is safe; multi-variant clusters drawn from a list that now admits overlaps are not.
2. **Merged-list clusters are not haplotypes.** After merging, a multi-variant cluster can contain two heterozygous variants in *trans*. `generate_str` applies all variants in its range to one sequence, synthesizing a chimera present on neither haplotype, and any reach derived from it is meaningless. Singletons are immune.

Iteration would additionally amplify the "assumes all variants are true positives" property documented at `cluster.cpp:633`: a cluster containing a false positive synthesizes an incorrect haplotype, inflating its alignment score, lengthening its reach, and merging more aggressively on the next pass. That property is harmless per-variant and dangerous per-cluster.

Accordingly `-i` is **removed** from `parse_args`, and `max_cluster_itrs` becomes an internal constant rather than a configurable field. Removal is preferred over retaining a flag that errors, since the flag is already undocumented (`globals.cpp:494`) and no released version exposed it. The iterate-to-fixpoint scaffolding (`prev_active`, the `while` loop) is left in place at no cost, so the graph-based reach work (§9.2) can re-enable it once it is sound.

#### 4.2.4 Accepted limitation

With one pass, a merged cluster's reach is approximated as the minimum and maximum of its members' reaches (`cluster.cpp:862-878`, `:897-908`). A *set* of variants can jointly be re-represented over a wider window than any member alone — compensating indels being the canonical case — so this can under-estimate reach and therefore under-merge. This is a genuine accuracy gap and is recorded here as the concrete motivation for §9.2, where alternative paths and variant skipping are native and joint reach falls out correctly rather than being approximated by a linear haplotype whose validity degrades as it merges.

### 4.3 Phasing (`phase.cpp`)

#### 4.3.1 Exclude unphased variants from the DP

`phase()` classifies each variant as `PHASE_ORIG`, `PHASE_SWAP`, or `PHASE_NONE` purely by comparing `orig_gt` against `calc_gt` (`phase.cpp:533-541`). An unphased heterozygote carries an arbitrary orientation and would be classified as a real ORIG or SWAP observation, fabricating switch and flip errors.

The fix is to force `PHASE_NONE` wherever `is_phased == false`. No change to the DP itself is required: `PHASE_NONE` already costs 0 in *both* orientations (`phase.cpp:558-561`), so such variants neither anchor a phasing nor penalize a swap, and they are already excluded from flip detection (`phase.cpp:608`). They become correctly transparent.

#### 4.3.2 Preserve the phase-set backfill

It is tempting to mark unphased variants by giving them a distinct phase set. This must not be done. `cost_swap` is zero whenever `phase_sets[i] != phase_sets[i+1]` (`phase.cpp:571-574`), and a swap is only recorded as a switch error when the phase set matches across the boundary (`phase.cpp:601-604`). A distinct phase set on an unphased variant therefore creates two cost-free switch boundaries around it, silently erasing genuine switch errors and shattering phase blocks.

Unphased variants must retain the surrounding backfilled phase set while being independently marked unphaseable — hence the separate `is_phased` field (§4.1.2).

#### 4.3.3 Defects fixed in scope

Three existing defects are latent today and become the common path under this deliverable, so all three are **fixed here**. Two are specific to unphased input; the third is not, and is additionally **filed on its own** (§6) because it affects every input:

- **Stale contig index in the no-phase-set early exit.** `fix_phase_set_tags` pushes `this->lengths[ctg_idx]` and then `continue`s (`phase.cpp:302-305`), skipping the `ctg_idx++` at `:337`. Every contig lacking phase sets therefore reports the length of whichever contig the index is stuck on. A fully unphased query makes this fire on every contig.
- **An unphased VCF reports a genome-spanning phase block.** The `phaseblockData` constructor builds phase blocks from raw `phase_sets` (`phase.cpp:250-262`) *before* `fix_phase_set_tags` runs at `:265`. With no `PS` in the header every variant gets `phase_set = 0` (`variant.cpp:984-990`), yielding exactly one phase block per contig; `fix_phase_set_tags` then pushes the full contig length as that block's size. The result is an NG50 approximately equal to contig length — an unphased callset reporting perfect phasing. Phase-block construction must treat "no phase information" distinctly from "one phase set".
- **Out-of-bounds read on the phase-set boundary check.** The backward pass begins at `i = qvars->n` and reads `qvars->phase_sets[i]` (`phase.cpp:601`), but `phase_sets` has size `n`. This fires whenever `ptrs[phase][n] == PHASE_PTR_SWAP`, which is common. The forward pass guards correctly (`if (i < qvars->n-1 && ...)`, `:572`); the backward pass does not. Unlike the other two this is not specific to unphased input and is filed separately (§6), but the fix lands here because this deliverable rewrites the surrounding code.

---

## 5. Metric definitions

### 5.1 Switch and flip error rates

The DP records phase changes as *inter-variant regions*: `mat`/`ptrs` are sized `n+1`, the forward pass writes `mat[·][j]` while consuming variant `j-1` (`phase.cpp:576-585`), `cost_swap` keys off the boundary `phase_sets[i] != phase_sets[i+1]` (`:571-574`), and the backward pass guards on that same boundary before recording a switch (`:601-604`). The header comment at `:567-568` states this directly.

It follows that the opportunity count is boundaries, not variants — and that boundaries adjacent to `PHASE_NONE` variants are not distinguishable opportunities, since a swap anywhere within a run of zero-cost variants is the same event. The denominator is therefore:

> **`ASSESSED_PAIRS`** — summed over phase blocks, `max(0, k - 1)` where *k* is the number of phaseable heterozygous query variants in that block (`phases[i] != PHASE_NONE`).

```
SWITCH_ERROR_RATE = SWITCH_ERRORS / ASSESSED_PAIRS
FLIP_ERROR_RATE   = FLIP_ERRORS   / ASSESSED_PAIRS
```

Both are **fractions**, which is also what the split writers emit:
[`D5`](./D5_multi-bed-stratification.md) §7.4 replaces the current `%.6f%%` percent-string with a bare
float and a terminating newline, so this deliverable inherits the corrected format rather than
changing it again.

Blocks containing a single phaseable heterozygote contribute zero, correctly reflecting that they present no switch opportunity. This quantity coincides with WhatsHap's `all_assessed_pairs`, so v3 rates are directly quotable alongside `whatshap compare` output in the Deliverable #7 benchmarking — a convenience, not the derivation.

The `phases[i] != PHASE_NONE` population excludes homozygous variants, unphased-input heterozygotes, and heterozygotes whose `calc_gt` returned `GT_REF_REF`. That is precisely the set the DP is able to penalize, so it is the correct denominator population regardless of the pairs-versus-sites question.

**This replaces the current denominator**, `qvars->n` — every merged query variant on the contig, homozygotes included (`phase.cpp:639`, consumed at `:808-809`). No comparability is lost: `SWITCH_ERROR_RATE` and `FLIP_ERROR_RATE` were added by commit `a1a2852` (2024-11-10) and are **not an ancestor of any release tag**. The latest release, v2.6.4, emits `PHASE_BLOCKS  SWITCH_ERRORS  FLIP_ERRORS  NG_50  SWITCH_NGC50  SWITCHFLIP_NGC50` and no rates, matching the documented v2.5.3 schema (`docs/v2.5.3/09-Outputs.md:124-133`). All published vcfdist figures use counts and NG50/NGC50, which this change does not touch.

### 5.2 Phased fraction

Two figures are reported, because they answer different questions and disagree informatively where phasing is sparse:

> **`PHASED_HET_FRACTION`** — phaseable heterozygous query variants ÷ evaluated heterozygous query variants.
>
> **`PHASED_REGION_FRACTION`** — evaluated BED bases covered by a phase block spanning at least two phaseable heterozygotes ÷ total evaluated BED bases.

A callset can be high on the first and low on the second (phasing concentrated in variant-dense regions), or the reverse (long blocks anchored by few heterozygotes).

Note a pre-existing inconsistency this exposes: `calculate_ng50` normalizes by total *contig* length (`phase.cpp:829-832`), not evaluated BED bases. `PHASED_REGION_FRACTION` deliberately uses the BED denominator, since it describes the evaluated region. The two denominators must be documented explicitly rather than left for the reader to infer; harmonizing NG50 onto BED bases is out of scope here.

### 5.3 Output changes

By the time this deliverable lands, [`D5`](./D5_multi-bed-stratification.md) §7.1 has split
`phasing-summary.tsv` into a variant-attributed file and a genome-level one. The new columns land in
whichever file matches what they measure:

**`phasing-variants-summary.tsv`** — stratified, `STRATUM` first column:

| Column | Type | Description |
| :-- | :-- | :-- |
| `ASSESSED_PAIRS` | integer | **Replaces `VARIANTS`** as the denominator of the switch and flip error rates (§5.1). |
| `PHASED_HET_VARIANTS` | integer | Phaseable heterozygous query variants (`phases[i] != PHASE_NONE`). |
| `TOTAL_HET_VARIANTS` | integer | Evaluated heterozygous query variants. |
| `PHASED_HET_FRACTION` | float | §5.2 — derived per stratum from the two counts above. |

**`phasing-blocks-summary.tsv`** — genome-level, no `STRATUM`:

| Column | Type | Description |
| :-- | :-- | :-- |
| `PHASED_REGION_FRACTION` | float | §5.2. Genome-level because it measures phase-block *span* coverage; it could in principle be restricted to a stratum's bases, but that needs per-stratum base accounting rather than the variant attribution everything else in the stratified file uses. Left genome-level, with the option noted. |

`VARIANTS` is removed from the stratified file when `ASSESSED_PAIRS` replaces it: it was exposed by
D5 only as the denominator it then divided by (`sum(qvars->n)`, homozygotes included), and keeping a
denominator nothing divides by would invite exactly the confusion the rate change is meant to fix.
D5 §7.1 records the hand-off from its side.

**Stratifying a pair count is an attribution, not a partition.** `ASSESSED_PAIRS` counts boundaries
*between* consecutive phaseable heterozygotes, and a boundary can span two variants in different
strata. Each pair is attributed to the strata of the **second** variant — the same variant
`switchflips.tsv` reports in its `VARIANT` column and the same one D5 §7.4 attributes switch errors
to — so a stratum's rate is always errors-attributed-to-it over pairs-attributed-to-it, and numerator
and denominator can never disagree about which stratum an event belongs to. What this does *not* give
is a decomposition: pair counts over complementary strata do not sum to the genome-wide total, for
the same boundary-straddling reason D5 §3 gives for variants. It is approximate in exactly that
sense, and the wiki must say so rather than presenting per-stratum rates as a partition of the
genome-wide rate.

Emitting the denominators rather than only the rates lets any consumer recompute any convention,
keeps the D10 MultiQC module robust to a later convention change, and documents the denominator in
the data rather than only in the wiki.

**No new `summary.vcf` field for phased status.** A per-variant phased-status FORMAT field would let a
user audit which calls were excluded from the phasing metrics, but it is unnecessary: the information
is already on the record twice over.

- **`GT` carries it directly.** Unphased calls are written with `/` and phased calls with `|` (below),
  so `0/1` *is* the statement "heterozygous, unphased, excluded from the phasing metrics".
- **`VP` carries the resolved phase** for variants that have one, and
  [`D4`](./D4_ga4gh-compatibility.md) §6.3 derives its `PHASE_UNPHASED` criterion from
  `phases[vi] == PHASE_NONE` — the same predicate — so a `BK` of `gm` rather than `pm` on a
  genotype-matched call is the second visible signal.

A third field would duplicate both and add a schema element for no new information.

**Genotype separators.** `print_variant` emits `|`-separated genotypes unconditionally
(`variant.cpp:562`, `:569`). Round-tripping unphased input as phased output is misleading; unphased
calls must be written with `/`. This applies to `summary.vcf`, which after
[#96](https://github.com/TimD1/vcfdist/issues/96) is the only VCF vcfdist emits. Because this is what
makes `GT` load-bearing for the point above, it is a requirement of this deliverable rather than a
cosmetic fix.

---

## 6. Defects filed separately

Two items surfaced during this design that are not caused by, and not confined to, unphased evaluation. Both are *filed* separately; the first is nonetheless *fixed* as part of this deliverable, since it sits in code this work rewrites:

1. **`phase.cpp:601` out-of-bounds read** (described in §4.3.3). Fixed as part of this work because it sits in rewritten code, but filed independently since it affects all inputs.
2. **The graph-based clustering redesign** (§9.2), which is a deliverable in its own right.

Per repository convention, neither issue references this document or any planning artifact.

---

## 7. Testing

Deliverable #2 owns the unit-test suite; this section names the cases this deliverable adds to that enumeration, to be folded into [`D2`](./D2_vcfdist-v3-unit-tests.md).

**`test_variant.cpp`** — `pv-unphased-het-kept` (replaces the existing `pv-unphased-het-skipped`), `pv-unphased-hom-phased-flag`, `pv-unphased-compound-het-trans`, `pv-phased-het-flag-set`, `pv-is-phased-copy-roundtrip`, `pv-query-overlap-admitted`, `pv-query-phased-overlap-admitted`, `pv-truth-overlap-still-skipped`, `pv-truth-simple-gt-downgrade-retained`, `pv-unphased-gt-slash-separator`.

**`test_dist.cpp`** — these refine [`D2`](./D2_vcfdist-v3-unit-tests.md) §5.4's single
`graph-overlapping-query-vars` case, which should be replaced by them rather than kept alongside:
`graph-overlapping-query-del-sub` (deletion `[100,110)` and substitution `[105,106)` both appear as nodes, a reference node `[100,105)` is emitted, and no path traverses both), `graph-overlapping-query-assert-holds`, `graph-overlapping-query-same-pos`, `graph-ins-then-sub-same-pos-cis` (guards that a zero-`rlen` insertion is not pushed to `qnode_ends` and remains on the cis path).

**`test_cluster.cpp`** — `merge-preserves-is-phased`, `cluster-on-merged-list-single-callset`, `cluster-overlapping-vars-no-crash`, `cluster-unphased-het-no-hap-dependence` (same superclusters whether a het is presented as `0|1`, `1|0`, or `0/1`), and removal of the cases covering the deleted cross-haplotype reach merge.

**`test_phase.cpp`** — `phase-unphased-forced-none`, `phase-unphased-transparent-to-dp` (a run of unphased variants between two ORIG variants yields zero switches and zero flips), `phase-unphased-keeps-backfilled-ps` (guards §4.3.2), `phase-assessed-pairs-single-het-block-zero`, `phase-assessed-pairs-multi-block`, `phase-rate-denominator-excludes-hom`, `fps-ctg-idx-not-stale` (pins §4.3.3), `pbd-all-zero-ps-not-one-block` (pins §4.3.3), `phase-backward-pass-no-oob` (pins §6, item 1).

**Integration.** A fully-unphased query fixture and a partially-phased query fixture, each asserting that P/R summary counts match the phased equivalent of the same callset, and that phasing metrics degrade to zero-with-explicit-denominators rather than silently reporting perfect phasing.

---

## 8. Risks and expected deltas

| Risk | Assessment |
| :-- | :-- |
| **Supercluster boundaries shift** | Per-callset variant density differs from per-haplotype density, so clustering produces different superclusters even for phased input, moving P/R numbers. Magnitude is **unmeasured**; Deliverable #7 quantifies it against its baseline and reports the delta. |
| **Runtime and RAM** | Clustering thread fan-out halves, but each thread handles roughly twice the variants; per-cluster work is unchanged since reaches remain per-variant. Net effect expected to be small, unverified. |
| **Credit-tie coupling to phasing** | §2.2(2) remains. P/R is not strictly independent of phasing on exact credit ties. Documented, not removed. |
| **Overlap filter removal (query)** | Admitting overlapping query variants exercises `Graph` paths that real input could not previously reach, since the filter guaranteed they never occurred. The `qnode_ends` priority-queue logic (`dist.cpp:908`, `:922-928`, `:956-959`) is untested for the overlap case and needs direct coverage. Also increases query node counts in dense regions, with a corresponding wavefront cost; unmeasured. Affects phased callsets too, so it contributes to the D7 delta independently of unphased support. |
| **Metric change is user-visible** | The rate denominator changes and new columns appear. Low risk given no rate column has shipped (§5.1), but it belongs in the v3.0.0 release notes (Deliverable #8) regardless. |

---

## 9. Out of scope

### 9.1 Unphased truth VCFs

The truth side of the alignment graph is a linear per-haplotype chain, selected by `var_on_hap` (`dist.cpp:1003`). Supporting an unphased truth requires evaluating both haplotype pairings per supercluster and taking the best-scoring assignment, which is a new branch in truth-graph construction and interacts with the retry loop in `evaluate_variants`. It is a substantially larger change with genuine research risk, and it is not needed for the target use case: GIAB, HPRC, and PAV truth sets are locally phased.

Unphased truth heterozygotes therefore continue to be skipped at parse, with an escalated warning (§4.1.1).

### 9.2 Graph-based cluster reach

The clustering restructure in §4.2 fixes three of the four known problems with the current algorithm: it handles unphased variants, permits overlapping variants, and removes the cross-haplotype merge heuristic. It does not fix the fourth — that reaches are computed by applying variants to the reference without allowing them to be skipped (`cluster.cpp:633`) — nor the joint-reach under-merge described in §4.2.4.

The sound fix is to compute reach on the variant graph rather than on a synthesized linear haplotype. Most of the machinery already exists: `Graph` plus `calc_prec_recall_aln` is already a graph-to-linear wavefront alignment that represents overlapping variants as alternative paths and permits skipping any variant, because the reference path is always present. What is missing is the *reach* primitive — `wf_swg_max_reach` (`dist.cpp:1135`) keys wavefront offsets by diagonal, and a graph version must key them by *(node, diagonal)* and traverse the graph in reverse for the left reach. `calc_prec_recall_aln`'s packed `idx4` provides a template for the per-node bookkeeping.

This work is closer in character to SOW #12 (faster alignment, already Potential Future Work) than to this deliverable, and it would move every supercluster boundary immediately before the v3.0.0 release. It is specified here so that §4.2 does not foreclose it: the iterate-to-fixpoint scaffolding is retained, and the merged-list clustering is the correct substrate for it.

### 9.3 Overlapping truth variants

§4.1.3 removes the overlap filter for the query callset and retains it for truth, because `Graph`'s truth side is a linear chain that silently constructs a nonexistent haplotype when fed overlapping variants (`dist.cpp:999-1024`).

Supporting overlapping truth variants means making the truth side a DAG, as the query side already is. The linearity is load-bearing beyond the construction loop: `get_truth_pos` (`dist.cpp:857-863`) derives a position by summing `tseqs` lengths in order, and `calc_prec_recall_aln` indexes `this->truth` as a flat string, so both would need to become node-aware. This overlaps §9.2's graph work and is also a prerequisite for SOW #13, where symbolic and breakend alleles cannot be laid out on a linear chain either.

Until then, overlapping truth calls continue to be dropped at parse. Because truth sets are curated and generally normalized, this affects far less data than the query-side filter did — but it should be counted and reported rather than passed over at verbosity 2. With [`D3`](./D3_retain-info-format-fields.md) landed, the natural home for that report is the same drop-counter machinery D3 §5 already standardizes.

---

## 10. References

- Companion SOW: `D0_vcfdist-v3-SOW.md`
- Parent design doc: [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md) §5.6
- Unit-test design: [`D2_vcfdist-v3-unit-tests.md`](./D2_vcfdist-v3-unit-tests.md)
- WhatsHap `compare` metric definitions — https://whatshap.readthedocs.io/en/latest/guide.html
- Martin et al. 2016, *WhatsHap: fast and accurate read-based phasing*, bioRxiv
- Dunn & Narayanasamy 2023, *vcfdist: Accurately benchmarking phased small variant calls in human genomes*, Nat. Commun.
