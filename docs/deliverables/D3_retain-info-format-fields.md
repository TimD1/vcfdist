# D3: `summary.vcf` — Record Shape, Original INFO/FORMAT Fields, and Unevaluated Variants

- **Status:** Design — complete, pending review. Decisions recorded in §10.
- **Issue:** [#48](https://github.com/TimD1/vcfdist/issues/48) (D3: Retain Original INFO/FORMAT Fields)
- **Branch:** `48_td_D3-retain-info-format-fields`
- **Area:** `src/variant.{cpp,h}`, `src/phase.cpp`, `src/cluster.cpp`, `src/bed.cpp`,
  `src/globals.{h,cpp}`
- **Owns:** the `summary.vcf` **record shape** (§4.2) as well as field retention — the record layout
  that [`D4`](./D4_ga4gh-compatibility.md) annotates, [`D5`](./D5_multi-bed-stratification.md)
  stratifies, and [`D6`](./D6_vcfdist-v3-unphased-eval.md) extends.

## 1. Background

`summary.vcf` is not a copy of the input VCFs with annotations added. It is *synthesized*
from vcfdist's internal `ctgVariants` arrays by `phaseblockData::write_summary_vcf`
(`phase.cpp:20`), which `fprintf`s a merged two-sample (TRUTH, QUERY) record set. The
site-level columns are hardcoded placeholders in `ctgVariants::print_var_info`
(`variant.cpp:461`, `473`): `ID=.`, `QUAL=.`, `FILTER=PASS`, `INFO=.`. Only CHROM/POS/REF/ALT
carry data, and those are reconstructed from *normalized* alleles — prefix/suffix-trimmed and
re-anchored on the preceding reference base (`variant.cpp:472-475`).

Separately, `parse_variants` (`variant.cpp:673`) silently discards variants at ten distinct
`continue` sites. A discarded variant never enters `ctgVariants`, so nothing survives but an
aggregate counter reported as an `INFO`/`WARN` line at `variant.cpp:1173-1239`. The output is
therefore a filtered subset of the input with no per-record accounting of what was removed
or why.

Two structural facts constrain the design:

1. **One input record does not correspond to one output line.** A record becomes 0, 1, or 2
   `ctgVariants` entries. Homozygous variants are split per haplotype at parse time, re-joined
   during superclustering (`cluster.cpp:134-159`), then re-split into per-haplotype lines by
   `write_summary_vcf`. `TYPE_CPX` variants are decomposed into an INS + DEL pair
   (`variant.cpp:1150-1153`). Multiallelic `1|2` records yield two entries with different ALTs.
2. **The hap merge is inferential, not provenance-based.** `cluster.cpp:134-139` matches
   haplotype copies by `(pos, ref, alt)` equality, so it over-merges two distinct het records
   sharing a locus and alleles, under-merges when one copy was dropped by the overlap filter
   (`variant.cpp:1141-1144` rewrites `simple_gt` to het in that case), and cannot merge a `1|2`
   record at all.

## 2. Goals / Non-Goals

**Goals**

- **Emit one record per variant, not one per haplotype copy.** A homozygous variant becomes a
  single `GT=1|1` record instead of two `1|0` / `0|1` records; het-alt (`1|2`) variants stay split
  (§4.2, §4.3).
- Track **ploidy per variant** rather than per contig, so haploid and mixed-ploidy input renders a
  correct `GT` (§4.4).
- Preserve each source record's `ID`, `QUAL`, `FILTER`, `INFO`, and sample `FORMAT` fields on
  the corresponding `summary.vcf` line(s).
- Retain variants vcfdist cannot evaluate, marked `BD=N` per GA4GH convention, each carrying a
  new `FILTER` tag naming the reason, and excluded from every analysis stage.
- Correctly subset allele-indexed (`Number=A`/`R`/`G`) fields to the allele each output line
  represents.

**Non-Goals**

- **Merging multiallelic records.** A `1|2` record remains two co-located biallelic records
  (§4.3). Only the *homozygous* split is undone.
- **Aggregating per-haplotype `BD`/`BK` into a single per-site decision.** This deliverable makes
  every per-haplotype field multi-valued, `BD`/`BK` included, which needs no aggregation policy at
  all; [`D4`](./D4_ga4gh-compatibility.md) §6 then collapses `BD`/`BK` to one value per site via
  its match-tier ladder (§4.5). Splitting the two is what keeps this deliverable a mechanical
  record-shape change.
- **The counting convention in `precision-recall*.tsv`.** vcfdist still counts one unit per
  haplotype copy after this deliverable; [`D4`](./D4_ga4gh-compatibility.md) §4 changes that.
- The `query.tsv` / `truth.tsv` tables. Scoped to `summary.vcf`. (`query.vcf` and `truth.vcf` are
  removed outright by [#96](https://github.com/TimD1/vcfdist/issues/96), so `summary.vcf` is the
  only VCF vcfdist emits.)
- Changing precision/recall, credit, or phasing semantics. Retained `BD=N` variants are inert.
- Overlapping variants and unphased heterozygous genotypes — both explicitly deferred, see §3.

## 3. Deferred Work (Out of Scope)

Two drop reasons are **not** addressed by this deliverable. For both, D3 leaves the current
behavior untouched: the variant continues to be discarded at parse time with only its existing
counter and `WARN`.

- **Overlapping variants** (`variant.cpp:1131-1139`). The intended end state is to keep *and
  evaluate* these, not to mark them `BD=N`, and that is not VCF I/O plumbing — so it is **deferred
  to D6**, not to a standalone issue. [`D6`](./D6_vcfdist-v3-unphased-eval.md) §4.1.3 found the
  query side needs less than this deliverable assumed: `Graph`'s query construction already
  represents overlapping calls as mutually exclusive alternative paths, so D6 drops the filter for
  the `QUERY` callset outright. The truth side genuinely does need graph-construction work — the
  truth chain is linear and silently synthesizes a haplotype that does not exist when fed overlaps
  — and stays out of scope for v3 (D6 §9.3).
- **Unphased heterozygous genotypes** (`variant.cpp:1030-1038`). Deferred to **D6**
  (issue #46, unphased variant evaluation), which evaluates these rather than marking them
  not-assessed.

**Known limitation this creates.** Issue #48 frames the output as "a complete accounting of the
input rather than a silently-filtered subset." With both reasons still discarding variants, that
accounting is *not* yet complete after D3. **D6 closes both** — unphased heterozygotes become
evaluated (D6 §4.1.1) and the query-side overlap filter is removed (D6 §4.1.3) — so the accounting
is complete for the query callset once D6 lands. Overlapping *truth* variants remain dropped in v3
(D6 §9.3). This should be stated in the user-facing documentation so the guarantee is not
overclaimed.

## 4. Design Overview

Four independent mechanisms:

1. **Record shape.** One output record per variant rather than one per haplotype copy, with
   per-haplotype fields becoming multi-valued (§4.2–§4.5).
2. **Provenance.** Each internal variant gains a back-reference to its source record and the
   original ALT ordinal it came from, so the writer can recover the record's fields and subset
   allele-indexed ones correctly.
3. **Sideline container.** Variants retained but not evaluated are stored *outside*
   `ctgVariants`, so "excluded from all analysis" holds by construction rather than by
   auditing every loop. The writer merges them back into position order at output time.
4. **Owning-callset rule.** Site-level columns come from the record owning the output line —
   which is what the code already does implicitly for POS/REF/ALT.

### 4.1 The owning-callset rule

`write_summary_vcf` has four emitting branches, each already choosing one callset to supply
the site-level columns:

| Branch | Line | Owner |
| ------ | ---- | ----- |
| Query matches truth (same pos, ref, alt) | `phase.cpp:160` | QUERY |
| Same position, different alleles | `phase.cpp:182` | QUERY |
| Query only | `phase.cpp:195` | QUERY |
| Truth only (pure FN) | `phase.cpp:208` | TRUTH |

Preserved fields follow the same owner. Query wins on matched lines; the truth record supplies
pure-FN lines, where no query record exists to source from. **Accepted loss:** on a matched
line the truth record's `INFO`/`FORMAT` is not preserved. No prefixing or namespacing is
introduced, and no field is demoted from `INFO` to `FORMAT`.

`BD=N` records are never matched (they are excluded from evaluation), so their lines always
have exactly one populated callset and the new `FILTER` tags never collide.

### 4.2 One record per variant — the homozygous split is undone

`write_summary_vcf` currently emits one record per haplotype per variant, via four
`for (int qhi = 0; qhi < HAPS; qhi++)` loops (`phase.cpp:151`, `:177`, `:192`, `:205`). Each loop
becomes a single emission per variant.

By the time the writer runs, a homozygous variant is **already a single entry**:
`load_and_merge_callset_vars_across_haps` collapsed its two parse-time copies by exact
`(pos, ref, alt)` equality and recorded `GT_ALT1_ALT1` (`cluster.cpp:133-159`). So this is not a
merge — it is *declining to re-split*. That is what makes the change mechanical, and it is why this
deliverable can own it without any of D4's decision machinery.

The genotype comes from `gt_strs[orig_gts[vi]]` for query and `gt_strs[orig_gts[tvi]]` for truth.
Use `orig_gts`, **not** `calc_gts`: `calc_gts` is vcfdist's *inferred* genotype, and reporting it in
the `GT` column would present vcfdist's inference as the caller's claim. That rendering is correct
for diploid variants only, so §4.4 is a prerequisite, not an independent improvement.

**Why the writer owns this.** The per-haplotype split is a property of the *writer*, and this
deliverable rewrites the writer end-to-end (§6.4) to thread preserved fields through it. Un-splitting
there means the field-duplication logic is written once, for a single line, rather than per split line.

### 4.3 Het-alt records stay split

A `1|2` record is parsed into two entries with different ALTs, `simple_gt` `1|0` and `0|1`
(`variant.cpp:1011-1022`), and the cross-haplotype merge cannot rejoin them because their ALTs
differ (`cluster.cpp:133-139`). They remain **two co-located biallelic records** on output.

Provenance does not change this, for two reasons:

- **The two alleles need not share a position after normalization.** vcfdist trims and re-anchors per
  allele (`pos += lm` for insertions, `variant.cpp:1050-1060`), so `REF=T ALT=G,TT` can yield a SNP
  and an INS at different normalized positions. Emitting one record would mean falling back to the
  source record's original `POS`/`REF`/`ALT` — available via §6.2 — but that record would then carry
  un-normalized alleles while every other record in the file carries normalized ones.
- **Credit would have to be aggregated across two independently evaluated entries**, which may sit in
  different superclusters. Unlike the homozygous case, where a single entry already exists and only
  the writer splits it, this needs a new cross-entry rule.

The consequence is a permanent count residual at het-alt sites relative to tools that count a `1|2`
record as one location; D7 measures its size.

### 4.4 Per-variant ploidy

`GT_ALT1` (`defs.h:79`, rendered `"1"`) is assigned only to the **record-level** `orig_gt` used for
`GT_counts` bookkeeping (`variant.cpp:939`, `:943`). The value *stored per variant* is `simple_gt`,
which takes only three values — `hap ? GT_REF_ALT1 : GT_ALT1_REF`, overwritten to `GT_ALT1_ALT1`
when `same` (`variant.cpp:1017-1019`) — and that is what reaches `add_var` (`:1150-1157`). For a
haploid record the hap loop runs once with `same == false`, so a haploid call is stored as
`GT_ALT1_REF`, i.e. `"1|0"`.

The only thing correcting this today is the writer's **per-contig** ploidy:
`ploidy == 1 ? "1" : (thi ? "0|1" : "1|0")` at `phase.cpp:163`, and again at `:170`, `:185`, `:198`,
`:210`. `variantData::ploidy` holds one value per contig, inferred from the first record's GT length
(`variant.cpp:931-932`); later disagreement only warns, and the warning is suppressed outright when
the contig name ends in `X` (`variant.cpp:923`). So on a chrX carrying both PAR (diploid) and
non-PAR (haploid) calls, every output genotype is forced into whichever shape the first record had —
and chrY, which gets no exemption, warns on every record after the first.

Design:

- **Add `std::vector<uint8_t> ploidies;` to `ctgVariants`**, sized `n` alongside the other
  parsed-data vectors (`variant.h:66-79`), set from `std::abs(ngt)` at parse time.
- **Thread it through `add_var`** — the explicit overload (`variant.cpp:95`), the copying overload
  (`:45`), and all four `merged_vars->add_var` call sites in
  `load_and_merge_callset_vars_across_haps` (`cluster.cpp:140`, `173`, `196`, `217`). This is the
  same plumbing §6.1 adds for `rec_idxs`/`alt_idxs`, so all three land in one pass.
- **The writer consults `ploidies[vi]`** to render `GT`: one allele for ploidy 1, the
  `gt_strs[orig_gts[vi]]` rendering for ploidy 2.

**`simple_gt` is deliberately left alone.** Extending it to take `GT_ALT1` for haploid records would
look tidier, but `var_on_hap` returns true for *both* haplotypes on `GT_ALT1` (`variant.cpp:365`,
`:367`), so a haploid variant would begin counting twice; and `set_allele_errtype` has no `GT_ALT1`
branch (`variant.cpp:171-195`), so haploid variants would fall through to `AC_UNKNOWN`. A separate
vector is additive and requires no audit of the clustering, alignment, or phasing consumers of the
genotype enum.

**What the per-contig value is reduced to:**

- **Delete the chrX special case** at `variant.cpp:923`. Mixed ploidy within a contig is now
  legitimate, so there is nothing to suppress; chrY stops warning too.
- **Keep the polyploid `ERROR`** at `variant.cpp:977`. Ploidy > 2 remains unsupported.
- **Replace the per-contig equality check** in `bed.cpp:329-338`, which warns when truth and query
  disagree on a contig's ploidy. The check has real value — it catches a haploid truth paired with a
  diploid query — but a single value per contig can no longer express the input. Track the *set* of
  ploidies observed per contig per callset and warn when the two sets differ.
- **Drop the `ploidy=` attribute from `##contig`** lines (`phase.cpp:36-37`). It is not a VCF-spec
  contig attribute, and under mixed ploidy a single value is wrong; each record's `GT` now carries
  its own. `superclusterData::ploidy` (`cluster.cpp:352`, `:363`) and `phaseblockData::ploidy`
  (`phase.cpp:239`, read at `:62`) exist only to carry the per-contig value to the writer and are
  removed with it.

### 4.5 Per-haplotype fields become multi-valued

A merged record must carry one value per haplotype for every field vcfdist computes per haplotype,
and that count varies with the variant's ploidy (§4.4). Four declarations could express this:

| Declaration | Values | Ploidy-adaptive | Consumer requirement |
| ----------- | ------ | --------------- | -------------------- |
| `Number=2` | fixed 2 | no — wrong on haploid variants | universal |
| `Number=G` | 3 for a biallelic diploid | yes, but enumerates *genotypes*, not haplotypes | universal |
| `Number=P` | one per `GT` allele | yes, exactly | htslib ≥ 1.23 |
| `Number=.` | unbounded | yes, by not constraining | universal |

**Decided: `Number=.`** `Number=P` — *"The field has one value for each allele value defined in
`GT`"* (`VCFv4.4.tex:197`, FORMAT-only) — is semantically correct, but it was introduced in VCF 4.4
and `BCF_VL_P` first appears in **htslib 1.23** (absent in 1.19 through 1.22.1), so consumers on
older `bcftools` or `pysam` would report cardinality errors. These fields are informational, and the
pipeline-integration path (D11) runs through tools that routinely pin older htslib, so compatibility
wins. `Number=.` produces byte-identical records; only declared cardinality, and hence validation, is
given up — which is why §9 asserts the value count directly. `Number=P` is the right end state once
htslib ≥ 1.23 is widespread: a header-only change.

| Tag | Source | Number after this deliverable |
| --- | ------ | ----------------------------- |
| `BC` | `credit[hi][vi]` | `.` |
| `RD` | `ref_ed[hi][vi]` | `.` |
| `QD` | `query_ed[hi][vi]` | `.` |
| `SG` | `sync_group[hi][vi]` | `.` |
| `BD`, `BK` | per haplotype | `.` — **one value per haplotype, not yet aggregated** |
| `QQ`, `PS`, `PB`, `BS`, `VP`, `FE`, `GE`, `SC` | already per-variant | `1` |

**`BD`/`BK` are lists here, and only here.** Making them multi-valued is the mechanical consequence
of un-splitting; it needs no policy, because no information is combined. D4 §6 then replaces the
list with a single per-site value derived from its match-tier ladder, which is where the "which
haplotype's decision wins" question actually belongs.

**Value order is the order of alleles in the emitted `GT`**, which may be swapped relative to
internal haplotype indices (`calcgt_is_swapped`). With `Number=.` nothing in the format conveys
either the count or the order, so both must be stated in the header descriptions.

### 4.6 `FILTER` on evaluated records

Input `FILTER` values are preserved verbatim (§6.5), including on evaluated records. One consequence
to document rather than design around: a GA4GH consumer treats a non-`PASS` `FILTER` on an evaluated
record as a filtered call and demotes it — filtered TPs become FNs, filtered FPs become Ns. So a user
who accepts a non-`PASS` filter via `--filter` and then exports `summary.vcf` into a GA4GH pipeline
will see those calls demoted.

Preserving the caller's `FILTER` is the stated goal of this deliverable, and silently rewriting it to
`PASS` would misreport the input. The interaction belongs in D9's documentation of the GA4GH export
path.

## 5. Per-Drop-Site Disposition

This table is the authoritative decision record. "Drop" means the variant does not appear in
`summary.vcf` at all; "Keep (`BD=N`)" means it appears, is annotated, and is inert.

| # | Reason | Line | Scope | Disposition |
| - | ------ | ---- | ----- | ----------- |
| 1 | FILTER not in `--filters` | `869-873` | record | **Keep (`BD=N`)** + `VCFDIST_FAILED_FILTER` |
| 2 | `QUAL < --min-qual` | `875-880` | record | **Keep (`BD=N`)** + `VCFDIST_LOW_QUAL` |
| 3 | Variant larger than `--largest-variant` | `1121-1128` | allele | **Keep (`BD=N`)** + `VCFDIST_TOO_LARGE` |
| 4a | Fully outside all BED regions (`BED_OUTSIDE`) | `1105-1111` | allele | **Keep (`BD=N`)** + `VCFDIST_BED_OUTSIDE` |
| 4b | Straddles a BED region edge (`BED_BORDER`) | `1105-1111` | allele | **Keep (`BD=N`)** + `VCFDIST_BED_BORDER` |
| 4c | Contig absent from BED (`BED_OFFCTG`) | `1105-1111` | allele | **Keep (`BD=N`)** + `VCFDIST_BED_OFF_CTG` |
| 5 | Spanning deletion (`alt == "*"`) | `1041-1045` | allele | **Drop** — no variation. Counter + `WARN` total |
| 6 | ALT identical to REF (`A→A`) | `1073-1078` | allele | **Drop** — no variation. Counter + `WARN` total |
| 7 | Both alleles missing (`.\|.`) | `1020-1025` | record | **Drop**. Counter + `WARN` total |
| 8 | One allele missing (`1\|.`, `.\|1`) | `1020-1025` | allele | **Unchanged** — the known allele is already kept and evaluated. Counter + `WARN` total retained; see §5.2 |
| 9 | Overlapping variants | `1131-1139` | allele | **Out of scope** (§3) — unchanged, still dropped |
| 10 | Unphased heterozygous GT | `1030-1038` | record | **Out of scope**, deferred to D6 (§3) — unchanged, still dropped |

Cases 5, 6, and 7 keep today's behavior of discarding the variant, but must gain (or retain) a
counter and a `WARN` reporting the total, consistent with `variant.cpp:1173-1239`.

Case 8 is **not** a behavior change. The `alt_idx < 0` `continue` sits *inside* the `for (hap)`
loop (`variant.cpp:1011`), so it drops only
the haplotype carrying the `.`; the known allele on the other haplotype still reaches `add_var`. A
`1|.` record therefore already yields one variant on HAP1 with `simple_gt == GT_ALT1_REF` today —
exactly what [`D2`](./D2_vcfdist-v3-unit-tests.md) §5.2 pins as `pv-half-missing-hap1-kept`. What
distinguishes case 8 from case 7 is a property of the current code rather than a change D3 makes:
`.|.` loses both haplotypes and so drops the record, while `1|.` loses one. What is *not* settled
about half-calls is recorded in §5.2.

Cases 4a-4c share a single `continue` site today; splitting them into three tags requires only
branching on the `loc` value already computed at `variant.cpp:1105`, which is also already
stored per variant in `ctgVariants::locs`.

### 5.1 Which drops can be partial

A record is *partially* dropped when one iteration of the `for (hap)` loop (`variant.cpp:1011`)
hits a `continue` while another reaches `add_var`. Which tests can diverge depends on the GT:

- **Multiallelic (`1|2`)** — the haps carry different ALTs, so cases 3, 4, 5, 6, and 9 can all
  split the record. Case 4 splits only narrowly: `bedData::contains` receives the original
  `rec->pos`/`reflen`, which are allele-independent, so its result varies with `type` alone,
  and `type` matters only at `bed.cpp:145` (an INS exactly at a region end becomes
  `BED_BORDER`).
- **One allele missing (`1|.`)** — case 8.
- **Homozygous `1|1` / het `0|1`** — both iterations resolve the same `alt_idx`, so every test
  is a pure function of `(record, alt_idx)` and fires identically on both haplotypes, *except*
  the overlap check (case 9), the only test reading per-haplotype state (`prev_end[hap]`,
  `prev_type[hap]`). `variant.cpp:1141-1144` already compensates by demoting `simple_gt` from
  `GT_ALT1_ALT1` to het.
- **Never partial** — cases 1 and 2 precede the loop. Case 10 is counted per allele but its
  condition is allele-independent, so it affects every non-REF allele of the record.

`prev_end`/`prev_type` are updated only after a successful `add_var` (`variant.cpp:1160-1161`),
so a dropped allele never cascades into dropping the next variant.

### 5.2 Open question — half-call fidelity

Two things about half-missing genotypes are **not** settled by the §5 table, and neither is resolved
in this revision.

1. **`GT` rendering.** Because `simple_gt` becomes `GT_ALT1_REF`, a `1|.` record round-trips through
   `summary.vcf` as `1|0` — an explicit reference call on a haplotype the caller made no call on.
   §4.2 rejects `calc_gts` for the `GT` column precisely because it "would present vcfdist's
   inference as the caller's claim," and rendering a no-call as a ref call is the same category of
   misreport. Fixing it needs either a fourth `simple_gt` value or the half-call carried alongside
   `ploidies`, both of which touch §4.4's design rather than being additive to it.
2. **Half-call reporting.** [`D2`](./D2_vcfdist-v3-unit-tests.md) §5.2 pins three defects on this
   path that are still live in `dev`: the `GT_counts` histogram labels `1|.` as `.|.`
   (`variant.cpp:948` sets `GT_MISSING` when *either* allele is missing, so a half-call is
   indistinguishable from a no-call); the summary `WARN` says the variants were "skipped"
   (`variant.cpp:1204`) when a half-call was in fact kept; and `unknown_allele_total` increments per
   haplotype rather than per record (`variant.cpp:1024`), so a single `.|.` reports 2.

Both sit in code this deliverable already rewrites, which argues for absorbing them here; both are
also independent of field retention, which argues for a separate issue. **Unresolved** — recorded so
the call is made deliberately rather than by omission.

## 6. Detailed Changes

### 6.1 Provenance fields (`variant.{h,cpp}`, `cluster.cpp`)

Add to `ctgVariants` (`variant.h:66-79`), sized `n` alongside the existing parsed data:

- `std::vector<int> rec_idxs;` — ordinal of the source record within its input VCF.
- `std::vector<int> alt_idxs;` — original ALT ordinal (`alt_idx`, `variant.cpp:1019`), 1-based.
- `std::vector<uint8_t> ploidies;` — the variant's own ploidy (§4.4).

All three land in one pass over the same call sites, since the plumbing is identical. Later
deliverables add their own per-variant vectors the same way — one at a time, when each needs it
(D5 §6.3's stratum bits, D6 §4.1.2's `is_phased`) — rather than via a shared refactor.

All three are set by the explicit `add_var` overload (`variant.cpp:95`), copied by
`add_var(other_vars, idx)` (`variant.cpp:45`), and threaded through all four `merged_vars->add_var`
call sites in `load_and_merge_callset_vars_across_haps` (`cluster.cpp:140`, `173`, `196`, `217`).

The CPX INS and DEL halves (`variant.cpp:1150-1153`) share one `alt_idx`, so both receive
identical field subsets, which is correct — they derive from the same original allele. (Before
§4.2, the two haplotype lines of a `1|1` did too; there is now only one such line.)

### 6.2 Source record retention

**Decided: retain in memory.** `parse_variants` retains, per contig and callset, the fields
needed to reconstruct output columns for record ordinal `rec_idx`: `ID`, `QUAL`, the original
`FILTER` list, and the raw `INFO` and sample `FORMAT` data. These are held for the lifetime of
the run and read back by the writer.

Two alternatives were rejected:

- **Storing reconstructed text per *variant*** — would duplicate strings across the two entries of
  a CPX record and the two co-located entries of a het-alt record. The `rec_idx` indirection stores
  one copy per record instead.
- **Re-streaming the input VCFs at write time** — zero retention, but `write_summary_vcf`
  interleaves two callsets, so it would require two readers synchronized against the internal
  arrays, and would reread and re-parse both inputs. Rejected as added complexity for memory
  vcfdist already budgets for (`--max-ram`, default 64 GB).

The cost scales with input record count, not genome size; it should still be measured on a
WGS-scale VCF (~5M records) during implementation to confirm it is a small fraction of the
alignment working set, but it is not a design risk.

### 6.3 Sideline container for `BD=N` variants

Unevaluated variants (cases 1-4c) never enter `ctgVariants`. They are appended to a separate
per-contig, per-callset container keyed by **(record ordinal, haplotype)** — not by record,
because cases 3 and 4 are per-allele and can leave one haplotype evaluated while the other is
`BD=N`. Each entry stores the position needed to interleave it into output order and the
reason code that selects its `FILTER` tag.

The alternative — an `ignored` flag on `ctgVariants` — is rejected: it would require auditing
every loop over `poss`/`n` in `cluster.cpp`, `dist.cpp`, `phase.cpp`, and `print.cpp`, where one
missed check silently corrupts supercluster indices or precision/recall counts.

### 6.4 Writer changes (`phase.cpp`, `variant.cpp`)

- `write_summary_vcf` (`phase.cpp:89`) drops the four `for (qhi)` loops (§4.2), emitting one record
  per variant, and extends its merge walk over query and truth variants to also draw from the two
  sideline containers, so `BD=N` records interleave in position order.
- `print_var_info` (`variant.cpp:456`) takes the owning record's fields and emits real `ID`,
  `QUAL`, `FILTER`, and `INFO` in place of the `.`/`.`/`PASS`/`.` literals, and appends the
  record's original FORMAT keys to the fixed key list.
- `print_var_sample` (`variant.cpp:508`) renders `GT` from `orig_gts[vi]` and `ploidies[vi]`
  (§4.2, §4.4), emits the per-haplotype fields as `.`-cardinality lists in emitted-`GT` order
  (§4.5), and appends the source sample's original FORMAT values, subset per §7.
- The header block (`phase.cpp:31-55`) gains: propagated `##INFO`/`##FORMAT` lines for every
  preserved field, `##FILTER` lines for each of the six new reason tags, a `BD` description
  updated to include `N`, `Number=.` declarations for the per-haplotype fields with their order
  documented, and `##contig` lines stripped of `ploidy=` (§4.4).

### 6.5 FILTER tag semantics

Per the VCF spec, `FILTER` lists the filters a record *failed*; `PASS` means all passed. The
new reason tag is appended to the owning record's original `FILTER` list, except that a lone
`PASS` is replaced rather than appended to.

Six tags, named after the `defs.h` constants they correspond to so they stay traceable to
code, and `VCFDIST_`-prefixed to avoid collision with input FILTER IDs. **Decided:** uppercase
`VCFDIST_` prefix, matching the convention that VCF FILTER IDs are uppercase and making
vcfdist-added tags unmistakable against a caller's own filters.

These strings live in a table in `globals.cpp` alongside `region_strs` and `error_strs`. Note
`region_strs` (`globals.cpp:27`) uses space-padded display values — `"OFF CTG"` — so it cannot
be reused verbatim as a FILTER ID source; the new table is separate.

| Tag | Trigger | Constant |
| --- | ------- | -------- |
| `VCFDIST_FAILED_FILTER` | FILTER not in `--filters` | — |
| `VCFDIST_LOW_QUAL` | `QUAL < --min-qual` | — |
| `VCFDIST_TOO_LARGE` | exceeds `--largest-variant` | — |
| `VCFDIST_BED_OUTSIDE` | fully outside all BED regions | `BED_OUTSIDE` |
| `VCFDIST_BED_BORDER` | straddles a BED region edge | `BED_BORDER` |
| `VCFDIST_BED_OFF_CTG` | contig absent from BED file | `BED_OFFCTG` |

## 7. Allele-Indexed Field Subsetting

Because output lines carry normalized, split alleles while `Number=A`/`R`/`G` fields are
positionally tied to the *original* ALT list, those fields must be subset to the allele the
line represents. Length class comes from
`bcf_hdr_id2length(hdr, BCF_HL_INFO | BCF_HL_FMT, id)`:

| Class | htslib constant | Action for original ordinal `alt_idx` |
| ----- | --------------- | ------------------------------------- |
| `Number=A` | `BCF_VL_A` | take element `alt_idx - 1` |
| `Number=R` | `BCF_VL_R` | take elements `0` and `alt_idx` |
| `Number=G` | `BCF_VL_G` | keep the three diploid entries for alleles `{0, alt_idx}`, via index `k(k+1)/2 + j` |
| `Number=1`, fixed, `.`, Flag | `BCF_VL_FIXED`, `BCF_VL_VAR` | pass through verbatim |

Subsetting changes a field's cardinality, so propagated `##INFO`/`##FORMAT` header lines must
be rewritten: an `A`/`R`/`G` field becomes `Number=1`/`2`/`3` on output, since every output
record is biallelic — un-splitting homozygotes (§4.2) does not change that, and het-alt records
stay split precisely so that it holds (§4.3). Not rewriting them produces a file that fails
`bcftools` validation.

## 8. Data Flow

```
parse_variants
  ├─ evaluated allele        ──▶ ctgVariants (+ rec_idx, alt_idx) ──▶ cluster ──▶ P/R ──▶ phase
  ├─ retained, unevaluated   ──▶ sideline container (rec_idx, hap, pos, reason)
  ├─ dropped (cases 5,6,7)   ──▶ counter only ──▶ WARN total
  ├─ deferred (cases 9,10)   ──▶ counter only ──▶ WARN total   [unchanged, §3]
  └─ all records             ──▶ retained record fields, indexed by rec_idx

write_summary_vcf   (one record per variant, §4.2)
  └─ position-ordered merge of { query vars, truth vars, query sideline, truth sideline }
       ├─ site columns  ◀── owning record's fields (§4.1)
       ├─ GT            ◀── orig_gts[vi] + ploidies[vi]  (not calc_gts) (§4.2, §4.4)
       ├─ per-hap tags  ◀── BC/RD/QD/SG/BD/BK as Number=. lists, GT-allele order (§4.5)
       ├─ A/R/G fields  ◀── subset by alt_idx (§7)
       └─ BD            ◀── credit-derived decision per haplotype, or N for sideline entries
```

## 9. Testing Plan

- **Record shape (§4.2, §4.3):** a hom SNP, a hom INDEL, and a CPX each yield **one** record with
  `GT=1|1` and two values per per-haplotype field; a `1|2` record still yields **two** co-located
  records — asserted explicitly, so the asymmetry is pinned rather than incidental.
- **Per-variant ploidy (§4.4):** a wholly haploid contig emits `GT=1` and one value per
  per-haplotype field; a chrX carrying a PAR diploid call *and* a non-PAR haploid call emits both
  correctly, asserted with the haploid record first *and* with the diploid record first, since the
  current bug is "whichever came first wins"; a chrY input produces no ploidy-mismatch warning; a
  haploid truth against a diploid query still warns, via the observed-ploidy-set comparison
  replacing `bed.cpp:329-338`; `##contig` lines carry no `ploidy=`; polyploid input still `ERROR`s.
- **Per-haplotype value counts (§4.5):** because `Number=.` is unvalidated, the number *and order*
  of values must be asserted directly — two for a diploid variant, one for a haploid one, in
  emitted-`GT` order on a record whose computed genotype is swapped.
- **Field preservation:** a query VCF with `Number=1`, `A`, `R`, `G`, and Flag `INFO`/`FORMAT`
  fields; confirm values survive onto the merged record of a `1|1` and both halves of a
  CPX, and that `A`/`R`/`G` fields are subset to the right allele on a `1|2` record.
- **Header validity:** `bcftools view` the output and confirm no cardinality complaints,
  particularly for the rewritten `A`/`R`/`G` declarations.
- **Per-drop-site coverage:** one targeted record per row of the §5 table, asserting presence
  or absence in `summary.vcf`, the `BD` value, and the `FILTER` contents. Include the partial
  cases from §5.1: `1|2` with one oversized allele, `1|.`, and a `1|1` half-dropped by overlap.
- **BED tag discrimination:** three records — one wholly outside a region, one straddling a
  region edge, one on a contig absent from the BED — must receive three distinct tags.
- **Inertness:** confirm `BD=N` records change no `precision-recall*.tsv` or phasing output
  relative to the same run on a pre-filtered input — i.e. retaining them is provably output-only.
- **Counting unchanged:** `precision-recall-summary.tsv` is byte-identical before and after this
  deliverable on a fixture with hom, het, and het-alt variants. Un-splitting the *output records*
  must not change the *counts*; the counting convention is D4's change, not this one.
- **Original FILTER preservation:** a record failing an input filter *and* dropped for a
  vcfdist reason must show both tags.
- **Deferred reasons unchanged:** an overlapping variant and an unphased het record must still
  be absent from `summary.vcf`, with their counters and `WARN`s intact.

## 10. Decision Record

**One open question remains** — half-call fidelity (§5.2). Every other decision is recorded below,
with the reasoning where it is not obvious from the change itself:

| Decision | Section | Rationale |
| -------- | ------- | --------- |
| One record per variant; the homozygous split is undone here, not in D4 | §4.2 | The hom entry is already single — only the writer splits it — and this deliverable rewrites the writer anyway, so retaining the split would mean writing the field-duplication logic twice |
| Het-alt (`1\|2`) records stay split | §4.3 | Merging needs a fallback to un-normalized source alleles and a cross-entry credit rule; the count residual is reported instead |
| Per-haplotype fields, `BD`/`BK` included, become `Number=.` lists — no aggregation policy here | §4.5 | Making them multi-valued is the mechanical consequence of un-splitting and combines no information; collapsing `BD`/`BK` to one per-site value is D4's match-tier decision |
| `Number=.` rather than `Number=P` | §4.5 | `P` is semantically correct but needs htslib ≥ 1.23 to read; output is byte-identical, and §9 asserts the value count directly |
| Ploidy tracked per variant, via a new `ploidies` vector | §4.4 | `orig_gts` cannot express it — `GT_ALT1` never reaches the stored genotype — so §4.2 is wrong for haploid input without it |
| `simple_gt` left unchanged rather than extended to `GT_ALT1` | §4.4 | `var_on_hap` is true for both haplotypes on `GT_ALT1`, and `set_allele_errtype` has no branch for it |
| Per-contig ploidy retired; chrX special case deleted; `ploidy=` dropped from `##contig` | §4.4 | Mixed ploidy within a contig is legitimate, so there is nothing to suppress; a single per-contig value cannot express the input and is not a spec attribute |
| Input `FILTER` preserved verbatim on evaluated records, GA4GH demotion documented not prevented | §4.6 | Preserving the caller's `FILTER` is this deliverable's goal; rewriting it to `PASS` would misreport the input |
| Site columns come from the line's owning callset | §4.1 | Already what the code does for POS/REF/ALT; pure-FN lines have no query record to source from |
| Truth `INFO`/`FORMAT` not preserved on matched lines | §4.1 | Accepted loss; query annotations are the ones users need |
| Unevaluated variants live in a sideline container, not a flag on `ctgVariants` | §6.3 | Makes "excluded from all analysis" true by construction rather than by auditing every loop |
| Provenance via `rec_idxs` + `alt_idxs` rather than reconstructing groupings | §6.1 | Inference on `(pos, ref, alt)` cannot recover `1\|2` records and mis-groups co-located variants |
| Subset `Number=A`/`R`/`G` fields and rewrite their header cardinality | §7 | Passing them through verbatim yields values that are silently wrong for the emitted allele |
| Three distinct BED FILTER tags rather than one | §5, §6.5 | `BORDER` is a materially different condition from "not in any region"; `loc` is already computed and stored |
| `VCFDIST_` uppercase prefix for FILTER tags | §6.5 | Matches VCF convention for FILTER IDs; unmistakable against a caller's own filters |
| Retain source record fields in memory | §6.2 | Re-streaming would need two readers synchronized against the internal arrays; memory is within what `--max-ram` already budgets |
| Overlapping variants and unphased het GT deferred to D6, not to standalone issues | §3 | Neither is VCF I/O plumbing. D6 §4.1.1 evaluates unphased hets; D6 §4.1.3 drops the overlap filter for QUERY, since the query graph already represents overlaps as alternative paths. Only overlapping *truth* variants need graph work, and stay out of scope (D6 §9.3) |
| Accounting is incomplete until D6 lands | §3 | Recorded so the "complete accounting" claim in #48 is not overstated in user-facing docs; complete for the query callset after D6, with truth-side overlaps still dropped |
| Case 8 (`1\|.`) is unchanged | §5, §5.2 | The `alt_idx < 0` `continue` is per-haplotype, so the known allele already reaches `add_var`; D2 §5.2 pins this behavior |
| Half-call `GT` fidelity and the three half-call reporting defects left **open** | §5.2 | Both touch §4.4's ploidy design and are independent of field retention, so neither is silently folded in nor silently dropped |
