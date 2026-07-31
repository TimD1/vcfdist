# D5: Multi-BED Stratification Support

- **Status:** Design — complete, pending review. Decisions recorded in §11.
- **Issue:** [#47](https://github.com/TimD1/vcfdist/issues/47) (D5: Multi-BED Stratification Support)
- **Branch:** `47_td_D5-multi-bed-stratification`
- **Area:** `src/bed.{cpp,h}`, `src/globals.{cpp,h}`, `src/variant.{cpp,h}`, `src/cluster.cpp`,
  `src/print.cpp`, `src/phase.{cpp,h}` (no downstream consumers — see §10)

## 1. Background

vcfdist accepts exactly one BED file today, via `-b/--bed`, and that BED is a **filter**, not a
label. Three things follow from that:

1. **It removes variants.** `parse_variants` calls `g.bed.contains()` at `variant.cpp:1105` and
   `continue`s on `BED_OUTSIDE`, `BED_BORDER`, and `BED_OFFCTG` (`variant.cpp:1107-1111`). Only
   `BED_INSIDE` variants ever enter `ctgVariants`, so `locs[]` (`variant.h:71`) is a constant
   column in practice and `query.tsv`'s `LOCATION` field always reads `INSIDE`.
2. **It selects contigs.** `intersect_contigs` (`bed.cpp:187`) prunes query, truth, and reference
   contigs to those named in `-b`, and errors if a BED contig is absent from the reference.
3. **It is validated strictly.** `bedData::check()` (`bed.cpp:66`) *errors* on unsorted or
   overlapping intervals and warns on adjacent ones.

Stratification is the opposite operation: N named region sets that **label** variants without
removing any, without touching contig selection, and supplied by third parties whose sort/merge
state we do not control. GIAB ships 181 strata in
`GRCh38/v3.1-GRCh38-all-stratifications.tsv`, every one a `.bed.gz` at a path relative to the
manifest.

Today the only stratified reporting path is external: `docs/v2.5.3/10-Variant-Stratification.md`
tells users to bgzip `summary.vcf` and stratify it with a GA4GH quantifier. That works but requires a
hap.py install and a Python 2 virtualenv, and it means vcfdist's own
`precision-recall-summary.tsv` can never answer a stratified question.

Two structural facts constrain the design:

1. **A variant belongs to many strata at once.** GIAB's strata overlap heavily — one variant is
   simultaneously in `alldifficultregions`, `AllTandemRepeats`, `lowmappabilityall`, a GC-content
   bin, and an ancestry tract. Membership is therefore set-valued, which determines the output
   shape of every per-record file (§7.2).
2. **Evaluation results live in a *merged* container, not the parsed one.** `write_precision_recall`
   reads `ctg_scs->callset_vars[QUERY]` (`print.cpp:282`), which `cluster.cpp:140,173,196,217`
   builds by copying variants through `ctgVariants::add_var(other, idx)` (`variant.cpp:45`). Any
   new per-variant field must be threaded through that copy or it is silently lost between parse
   time and reporting time.

## 2. Goals / Non-Goals

**Goals**

- Accept N named stratification BEDs via a hap.py-compatible manifest TSV, including `.bed.gz`.
- Report precision/recall/F1 per stratum per variant type, at both the `NONE` and per-stratum
  `BEST` quality thresholds, in the GA4GH/GIAB reporting idiom (`*` = all regions).
- Report the full per-QUAL sweep and the allele-count error matrix per stratum.
- Record each variant's stratum set in `query.tsv` / `truth.tsv`.
- Stratify the switch/flip **error counts**, and split the genome-wide phasing **contiguity**
  metrics into their own file so the two are not conflated (§7.1).
- Warn loudly when a stratum matches nothing, so 181 silently-empty rows cannot read as results.
- Turn "the input VCF is coordinate-sorted" from an unchecked assumption into an enforced
  precondition (§6.4), since the membership sweep depends on it and vcfdist currently degrades
  silently without it.

**Non-Goals**

- **Changing which variants are evaluated.** Strata never filter. `-b` remains the sole filter,
  and passing `-st` must not change a single number in the `*` rows.
- **Stratifying `summary.vcf`.** A GA4GH consumer takes strata externally
  ([`D4`](./D4_ga4gh-compatibility.md) §3.3), and embedding a stratum field would duplicate the
  manifest and bloat an already-18 MB file. (`query.vcf` and `truth.vcf` no longer exist —
  [#96](https://github.com/TimD1/vcfdist/issues/96) removes them before rc1, so `summary.vcf` is the
  only VCF vcfdist emits.)
- **Stratifying `switchflips.tsv` or `phase-blocks.tsv`.** Both are interval/event records rather
  than variant records; a phase block spans many strata. Deferred, not blocked (§11, decision 16).
- **Per-stratum NG50 / NGC50 / phase-block counts.** These are contiguity metrics over genomic
  span; computed over phase blocks intersected with a fragmented stratum they are not meaningful.
  Explicitly genome-wide only (§7.1).
- **Region-set algebra.** No union/intersection/complement operators over strata. Users compose
  BEDs with `bedtools` before handing them to vcfdist.

## 3. Membership Semantics

A variant is a member of a stratum if its interval **overlaps** any region in that stratum, i.e.
`bedData::contains()` returns `BED_INSIDE` or `BED_BORDER`. `BED_OUTSIDE` and `BED_OFFCTG` are
non-members.

The interval queried is the same one the `-b` filter uses — the **original, pre-normalization**
representation `(rec->pos, rec->pos + reflen, type)` from `variant.cpp:1105` — so that a variant's
`-b` disposition and its stratum memberships are always derived from one interval.

**Divergence from hap.py, deliberately accepted.** hap.py assigns a variant by its *trimmed*
effective reference position, so a variant whose anchor falls outside a region is not a member
even if the variant body overlaps it.[^happy] Under vcfdist's rule a deletion straddling a tandem
repeat's edge *is* counted in that repeat stratum, which is usually the case the stratum exists to
explain. Two consequences must be documented rather than hidden:

- Per-stratum counts will not match hap.py's exactly on boundary-straddling variants. D7's
  concordance work should treat this as a known, characterized difference, not a bug.
- **Strata do not partition.** Two complementary BEDs covering the `-b` span sum to *at least* the
  `*` row, not exactly — anything straddling their shared boundary is counted in both. The
  integration test asserts `>=`, with equality only when nothing straddles (§9.2).

vcfdist's existing insertion rule is retained unchanged: `contains()` returns `BED_BORDER` for an
insertion at exactly a region's last base (`bed.cpp:144-146`), which mirrors hap.py's requirement
that a region contain both the padding base and the following base. Under the overlap rule such an
insertion is a member.

[^happy]: <https://github.com/Illumina/hap.py/blob/master/doc/happy.md> — "Variants are assigned to
stratification regions based on their effective reference position… an insertion at position 110
would be captured if the region spans positions 109-111."

## 4. Command-Line Interface

### 4.1 Flag rename (breaking)

`-s` is currently `--max-supercluster-size` (`globals.cpp:339`, usage at `globals.cpp:492`). It is
renamed so that the stratification flag can take the mnemonic short form:

| Before | After | Long form |
| ------ | ----- | --------- |
| `-s` | `-sc` | `--max-supercluster-size` |
| — | `-st` | `--stratification` |

`-sc`/`-st` match the existing multi-character convention (`-sv`, `-mq`, `-ct`, `-md`, `-ci`).
**No `-s` alias is retained** — a silent alias would let an old `-s 15000` invocation be
reinterpreted as a stratification path.

No in-repo caller passes vcfdist's `-s`. Three prose references to `--max-supercluster-size` in
`src/cluster.cpp:608,1001,1033` and the validation message at `globals.cpp:427` still read correctly,
since only the short form changes. The rename is user-facing and belongs in D8's release notes.

### 4.2 `-st, --stratification <STRING>`

Takes one manifest TSV, hap.py-compatible:

```
ancestry_AFR	ancestry/GRCh38_ancestry_AFR.bed.gz
refseq_cds	FunctionalRegions/GRCh38_refseq_cds.bed.gz
```

- Two tab-separated columns: stratum name, then BED path. Extra columns ignored.
- **Relative paths resolve against the manifest's own directory**, not the working directory, so
  GIAB's shipped TSVs work unmodified. Absolute paths are used as-is.
- Blank lines and lines beginning with `#` are skipped.
- Load order is manifest order; that is also output row order after `*`.

Validation, all `ERROR`:

| Condition | Message intent |
| --------- | -------------- |
| Manifest unreadable | Name the manifest path |
| A line has fewer than 2 fields | Name the manifest and line number |
| Duplicate stratum name | Name the collision; ambiguous output rows otherwise |
| Stratum named `*` | Reserved for the all-regions row |
| Stratum BED unreadable | Name **both** the manifest-relative and the resolved absolute path |

And one `WARN`, per §8.3: a stratum whose contigs do not intersect the reference at all.
A manifest that parses to zero strata is a `WARN`, not an error, and behaves as if `-st` were
absent.

New `Globals` members:

```cpp
std::string strat_tsv_fn;                 ///< Stratification manifest TSV filename (empty if none)
std::vector<std::string> strat_names;     ///< Stratum names, in manifest order
std::vector<bedData> strata;              ///< Parsed stratum regions, parallel to strat_names
int nstrata = 0;                          ///< Number of strata (0 if --stratification unused)
```

## 5. `bed.{cpp,h}` Changes

### 5.1 Compressed input (also fixes `-b`)

`bedData::bedData` reads with `std::ifstream` (`bed.cpp:27`), so a gzipped BED yields binary
garbage, `std::stoi` throws, and `parse_args` reports the misleading
`"Invalid BED filename provided"` (`globals.cpp:154`). Since every GIAB stratum is `.bed.gz`,
compressed input is a prerequisite, not a nicety.

Replace the `ifstream` with htslib's `hts_open` + `hts_getline` (`kstring_t`), which handles
plain, gzip, and bgzip transparently. htslib is already a link dependency — `bcf_open` is used
for VCFs at `globals.cpp:115`. This fixes `-b regions.bed.gz` as a side effect.

### 5.2 Lenient load for strata

`check()` errors on unsorted and overlapping intervals. That strictness is right for `-b`, where a
malformed evaluation region silently changes the denominator of every metric, and wrong for
third-party strata we neither author nor control. It is also load-bearing: `contains()`'s two
binary searches (`bed.cpp:129-136`) assume sorted, non-overlapping intervals, so an unmerged
stratum would not merely be untidy — it would return wrong answers.

Add a second constructor parameter and a merge step:

```cpp
bedData(const std::string & bed_fn, bool merge_overlaps = false);
void merge();  ///< Sorts intervals by start and merges overlapping/adjacent ones per contig
```

`-b` keeps `merge_overlaps = false` and the existing `check()`. Strata pass `true`, which sorts
then merges, recomputes `size`, and at verbosity ≥ 2 reports how many intervals were coalesced
per stratum. Merging is safe here precisely because strata are used only for membership tests —
`size` is the only aggregate derived from them, and merging makes it *more* correct by removing
double-counted bases.

### 5.3 Decouple `contains()` from the `-b` global

`bedData::contains()` opens with `if (!g.bed_exists) return BED_INSIDE;` (`bed.cpp:115`) — a
method on one `bedData` consulting a global flag about a *different* `bedData`. For strata this is
outright wrong: membership would depend on whether `-b` was supplied.

Remove that line and short-circuit at the sole `-b` call site instead (`variant.cpp:1105`):

```cpp
uint8_t loc = g.bed_exists
        ? g.bed.contains(ctg, rec->pos, rec->pos + reflen, type)
        : BED_INSIDE;
```

This is a small latent-coupling fix that D5 requires and that also makes `contains()` unit-testable
without global setup (§9.1).

## 6. Membership Computation and Storage

### 6.1 Where

In `parse_variants`, immediately after the `-b` gate passes (`variant.cpp:1105-1111`) and before the
size and overlap filters. Membership is computed once per source record-allele and reused for both
halves of a `TYPE_CPX` split (`variant.cpp:1150-1153`) and for both haplotype copies of a
homozygous variant, since all of them derive from one interval.

### 6.2 Storage

A flat bitset on `ctgVariants`:

```cpp
std::vector<uint64_t> strata_bits;  ///< Flat per-variant stratum membership, strata_words per variant
int strata_words = 0;               ///< ceil(g.nstrata / 64.0); 0 when --stratification is unused
```

`strata_words = (g.nstrata + 63) / 64` — 3 words for GIAB's 181 strata, so 24 B per variant, or
roughly 120 MB at 5M variant-alleles. Flat rather than `vector<vector<bool>>` to avoid a heap
allocation per variant. When `-st` is absent, `strata_words == 0` and the vector stays empty, so
the feature costs nothing when unused.

Accessors on `ctgVariants` keep the word arithmetic in one place:

```cpp
bool in_stratum(int vi, int si) const;
void set_stratum(int vi, int si);
```

### 6.3 Threading through the hap merge — the critical path

Per §1 fact 2, `write_precision_recall` reads the **merged** container. Both
`ctgVariants::add_var` overloads must carry the bits:

- `add_var(pos, rlen, ...)` (`variant.cpp:95`) appends `strata_words` words. Rather than adding 3
  words to an already 23-parameter signature, it appends zeros and the caller sets bits via
  `set_stratum()` immediately after — the parse-time path already knows the membership set.
- `add_var(other, idx)` (`variant.cpp:45`) must **copy** `other`'s words for `idx`. This is the
  overload the four merge sites use (`cluster.cpp:140,173,196,217`). Missing this is the one
  failure mode that produces plausible-looking all-zero strata with no error anywhere, so §9.2
  asserts it end-to-end rather than trusting review.

Because the hap merge matches copies by `(pos, ref, alt)` (`cluster.cpp:134-139`), both copies of a
homozygous variant carry identical bits and the merge is a copy of an identical set — no OR
semantics or conflict policy needed.

### 6.4 Lookup strategy — a monotonic cursor, with ordering established up front

The naive approach is `nstrata` calls to `contains()` per variant-allele: 181 strata × two binary
searches × ~5M alleles ≈ 2×10⁹ comparisons, against a reading stage that currently takes 8.5 s on
chr20. A sweep-style cursor per stratum makes this O(1) amortized instead — but only if both sides
of the sweep are ordered. Rather than implement both and switch, settle both requirements now.

**Requirement 1 — stratum intervals sorted and non-overlapping. Satisfied by construction.**
Not an assumption: §5.2's `merge()` sorts and merges every stratum at load, precisely because
`contains()`'s binary searches already require it. vcfdist controls this side of the sweep outright.

**Requirement 2 — the variant query stream non-decreasing within a contig. Enforced, see below.**
Two properties already hold, neither sufficient on its own:

- **Contigs never interleave.** `variant.cpp:837-850` errors with `"Unsorted %s VCF '%s', contig
  '%s' already parsed"` when a record returns to a contig already left. So a contig's records are
  contiguous, and cursors are only ever needed for the *current* contig.
- **Accepted records are non-decreasing per haplotype.** The overlap filter
  (`variant.cpp:1131-1139`) skips a record unless `pos >= prev_end[hap]`, and
  `prev_end[hap] == prev_pos + prev_rlen >= prev_pos`, so any *accepted* record has
  `pos >= prev_pos`.

That second property is an accident of a filter written for another purpose, it is per-haplotype
while the cursor is per-contig, and it holds only for accepted records — the membership query at
`variant.cpp:1105` runs *before* the overlap filter at 1131. Relying on it would make stratum
counts silently depend on an unrelated filter's semantics.

**So D5 adds the missing precondition check.** In the record loop, track the previous record
position within the current contig and `ERROR` when `rec->pos` moves backwards, reusing the
existing `"Unsorted %s VCF"` phrasing and reset alongside `prev_end`/`prev_type` at
`variant.cpp:848-849`. One integer comparison per record.

This is a genuine robustness improvement independent of stratification: today an unsorted VCF is
not rejected, it is silently degraded — records are dropped one at a time by the overlap filter
with only a `WARN`, quietly changing every denominator. It is also a **behavior change** (runs that
previously produced misleading numbers now fail loudly) and belongs in D8's release notes.

With both requirements settled, there is one implementation and no fallback path:

- Cursor state is `nstrata` integers, not `nstrata × ncontigs`, since contigs are contiguous. On a
  contig change, re-point each stratum at that contig's interval array and reset its index to 0.
- Strata with no intervals on the current contig are skipped wholesale for that contig — membership
  is trivially false.
- **Advance is position-only; `type` is used solely for classification.** The cursor advances while
  `stops[cursor] <= variant_start`; the insertion rule (`bed.cpp:144-146`) is then applied against
  the interval the cursor landed on. So the repeated equal-position queries produced by the
  haplotype loop (`variant.cpp:1011`) are free, and the `type`-dependent branch never affects
  cursor state.

**The classification logic must not be duplicated.** A cursor that reimplements
`INSIDE`/`BORDER`/`OUTSIDE` alongside `contains()` gives two implementations of one rule, free to
drift — and the `contains()` copy is the one D2 §5.1 tests. `contains()` currently interleaves two
concerns in `bed.cpp:128-156`: locating the variant (`start_idx` via `upper_bound` on `starts`,
`stop_idx` via `lower_bound` on `stops`) and classifying from how those indices relate. Split them:

- `classify(contig, start, stop, type, start_idx, stop_idx)` — the decision tree at
  `bed.cpp:138-156` (insertion rule, between-regions case, multi-region span) **plus the
  before-all / after-all `BED_OUTSIDE` early returns at `bed.cpp:123-126`**, which are load-bearing:
  a variant entirely left of the first region yields `start_idx == -1` and would otherwise be
  misclassified `BED_BORDER` by line 139. **Pure, no search.** Only the `BED_OFFCTG` contig-presence
  check (`bed.cpp:121`) stays with the caller, since the sweep already skips contigs where a
  stratum has no intervals.
- `contains()` keeps its two binary searches and calls `classify()`, so its behavior and its D2
  tests are unchanged.
- The stratum sweep derives the same two indices from its cursor and calls the same `classify()`.

That leaves one decision tree with two ways of locating the variant, which is also what makes
§9.1's brute-force oracle worth writing: it can compare *locations*, knowing classification is
shared.

**One entry point, confirmed.** `parse_variants` is the only site that computes membership.
`variantData::add_variants` (`variant.cpp:594`, declared `variant.h:120`), the CIGAR-derived path,
**has no callers anywhere in `src/`**. The `cluster.cpp` merge
sites copy already-computed bits (§6.3) and perform no lookups.

## 7. Output Schema Changes

Aggregate tables are widened **in place** with `STRATUM` as the leading column, `*` first,
matching hap.py's `Subset` convention. Without `-st`, each file contains only its `*` rows — the
numbers are identical to today's with one added constant column. The console summary table
(`print.cpp:413-457`) prints `*` only and is otherwise unchanged, so interactive output stays
readable at 181 strata.

### 7.1 Files

| File | Change | Rows without `-st` → with 181 strata |
| ---- | ------ | ------------------------------------ |
| `precision-recall.tsv` | `STRATUM` first column | 244 → ~44,408 |
| `precision-recall-summary.tsv` | `STRATUM` first column | 8 → 1,456 |
| `genotype-errors.tsv` | `STRATUM` first column | 4 → 728 |
| `phasing-summary.tsv` | **Split into the two below** | — |
| `phasing-variants-summary.tsv` | *new* — variant-attributed, `STRATUM` first column | 1 → 182 |
| `phasing-blocks-summary.tsv` | *new* — genome-level, no `STRATUM` | 1 → 1 |
| `query.tsv`, `truth.tsv` | `STRATA` set column appended | unchanged |
| `parameters.tsv` | `stratification_tsv`, `n_strata` rows | +2 |
| `switchflips.tsv`, `phase-blocks.tsv` | unchanged (§2) | — |
| `summary.vcf` | unchanged (§2) | — |
| `runtime.tsv` | unchanged | — |

**The `precision-recall*.tsv` columns being stratified include D4's.**
[`D4`](./D4_ga4gh-compatibility.md) §7.1 adds `FP_GT`, `FP_AL`, `QUERY_UNK`,
`FRAC_NA`, `TRUTH_TITV`/`QUERY_TITV`, and `TRUTH_HET_HOM`/`QUERY_HET_HOM`. The three counts take the
stratum axis exactly as `TRUTH_TP` and `QUERY_FP` do; `FRAC_NA` and the four ratios are derived per
stratum from per-stratum counts. Nothing here is special-cased — they are added to the same counter
arrays as §7.3 describes — but they must not be left genome-wide by accident, which is what §9.2's
"every stratum name appears in every widened file" assertion catches.

**Counting is per site.** D4 §4 replaces per-haplotype counting with one unit per variant per
callset, so a variant contributes to its strata once, not once per carried haplotype.
The stratum loop is therefore over a variant's stratum set, with no haplotype dimension.

**The `phasing-summary.tsv` split.** The split is along one axis: **can this quantity be restricted
to a subset of variants at all?** Everything variant-attributed goes in the stratifiable file;
everything that measures contiguity over genomic span stays genome-level.

| New file | `STRATUM` | Columns after this deliverable | Why |
| -------- | --------- | ------------------------------ | --- |
| `phasing-variants-summary.tsv` | yes, first column | `VARIANTS`, `SWITCH_ERRORS`, `FLIP_ERRORS`, `SWITCH_ERROR_RATE`, `FLIP_ERROR_RATE` | Each error is attributable to a single query variant, hence to that variant's strata |
| `phasing-blocks-summary.tsv` | no | `PHASE_BLOCKS`, `NG_50`, `SWITCH_NGC50`, `SWITCHFLIP_NGC50` | Contiguity over genomic span; an NG50 computed over phase blocks intersected with a fragmented stratum is not meaningful |

The names say what distinguishes them — per-variant quantities versus per-block ones — rather than
"errors" versus "metrics", because [`D6`](./D6_vcfdist-v3-unphased-eval.md) adds non-error columns
(phased-variant counts and fractions) to the first file, at which point an `-errors-` name would be
wrong.

`VARIANTS` is newly exposed. Today it is only a divisor (`phase.cpp:808-809`) computed as
`sum(qvars->n)` over contigs (`phase.cpp:639`); per-stratum rates are unauditable without the
per-stratum denominator, so it becomes a column.

**D6 then replaces that denominator.** D6 §5.1 substitutes `ASSESSED_PAIRS` — adjacent
phaseable-heterozygote pairs per phase block — for `sum(qvars->n)`, and adds `PHASED_HET_VARIANTS`,
`TOTAL_HET_VARIANTS`, and `PHASED_HET_FRACTION` alongside it. The phaseable subset does not exist
until D6 marks unphased variants, so this deliverable exposes the denominator it actually divides by,
under its true name.

**Stratifying a *pair* count is approximate, and the caveat belongs in the docs.** `ASSESSED_PAIRS`
counts boundaries *between* consecutive phaseable heterozygotes, and a boundary can span two variants
in different strata. D6 attributes each pair to the strata of the second variant, matching the
attribution `switchflips.tsv` and the per-stratum switch counts already use (§7.4), so numerator and
denominator stay consistent — a stratum's rate is always errors-attributed-to-it over
pairs-attributed-to-it. What it is *not* is a partition: pair counts across complementary strata do
not sum to the genome-wide total the way variant counts do, for the same boundary-straddling reason
§3 gives for variants. Document it as a stratified rate over an attributed denominator, not as a
decomposition of the genome-wide rate.

Nothing **in this repo** parses `phasing-summary.tsv` — the only in-repo references are the docs
pages (`docs/v2.*/09-Outputs.md:124`). One consumer outside the repo does: the MultiQC module of
[`D10`](./D10_vcfdist-v3-multiqc-module.md), whose search pattern is anchored on the filename, so the
split renames the file out from under it (see *Downstream Consumers*, below). That fix is owned by D10
and costs no work here, but the split is not consumer-free.

### 7.2 `STRATA` on `query.tsv` / `truth.tsv`

A **set-valued** column, not a key column: comma-separated stratum names in manifest order, `.`
when empty. Appended last, after `LOCATION` (`print.cpp:535-536`, `576-577`).

A `STRATUM` *key* column here would require emitting one row per (variant, stratum) pair. At
5.8 MB for chr20's `query.tsv`, whole-genome is roughly 30× that (~170 MB) before multiplying by
the mean strata-per-variant — order of a gigabyte per file. The set column adds one field and no
rows.

This also makes `LOCATION` genuinely informative for the first time: today only `BED_INSIDE`
variants survive parsing, so that column is a constant. [`D3`](./D3_retain-info-format-fields.md)
has already changed that by retaining `BD=N` variants, so `STRATA` reads `.` for
retained-but-unevaluated variants, which are excluded from every count.

### 7.3 Aggregation changes in `print.cpp`

`write_precision_recall` (`print.cpp:267`) gains a leading stratum axis on both counter arrays
(`print.cpp:271-276`):

```cpp
// ax0: stratum (0 = "*"), ax1: SNP/INDEL/SV/ALL, ax2: TP/FP/FN, ax3: QUAL
std::vector<std::vector<std::vector<std::vector<float>>>> query_counts(g.nstrata+1, ...);
```

At 182 strata that is 182 × 4 × 3 × 61 floats ≈ 2.1 MB per array — negligible.

Four things must change together, and each is a place where a stratified version can silently
diverge from the `*` version:

1. **The accumulation loops** (`print.cpp:286-338`) gain an inner loop over the variant's set
   strata, always also writing index 0.
2. **The truth FN carry above the call quality** (`print.cpp:333-336`) must likewise be
   per-stratum.

3. **`max_f1_score` / `max_f1_qual`** (`print.cpp:354-355`) become `[nstrata+1][VARTYPES]`, so
   `BEST` is chosen per stratum. A stratum's best-F1 threshold genuinely differs from the
   genome-wide one; reusing the global threshold would report a non-optimal row as `BEST`.

The "no variants pass all filters" warnings (`print.cpp:429-435`) must be **restricted to
stratum 0**. Unrestricted they would fire up to 182 × 3 times, and an empty stratum is expected,
not anomalous.

**The `compute_pr_f1` helper.** [`D2`](./D2_vcfdist-v3-unit-tests.md) §5.7 extracts the P/R/F1
arithmetic — duplicated at `print.cpp:371-373` and `438-440` — into a pure
`compute_pr_f1(query_tp, query_fp, truth_tp, truth_fn)` under issue #94. The stratum axis adds one
loop around that one helper rather than multiplying divergent copies.

### 7.4 Per-stratum switch/flip attribution

`ctg_pbs->switches[]` and `flips[]` hold **query variant indices** (`phase.cpp:602,609`;
confirmed by their use as indices at `print.cpp:516,520`). An error is attributed to the strata of
that variant, and `switchflips.tsv` already reports the same index in its `VARIANT` column
(`phase.cpp:734,745,766`), so the two files agree by construction.

One convention to document: a switch error is localized to the *interval* between two consecutive
phased variants (`beg`/`end` in `switchflips.tsv`), and that interval can straddle a stratum
boundary. Attributing to the variant at `next_vi` is a choice, made for consistency with the
existing `VARIANT` column rather than because the interval has a single correct stratum.

Per-stratum `VARIANTS` is the count of query variants in that stratum, mirroring the genome-wide
`sum(qvars->n)` at `phase.cpp:639`. `write_phasing_summary` (`phase.cpp:795`) is split into
`write_phasing_variants_summary` and `write_phasing_blocks_summary`, declared in `phase.h` — named
for the files they write and for the per-variant/per-block axis of §7.1, **not** `…_errors_…` /
`…_metrics_…`, for the same reason the filenames avoid that split: D6 adds non-error columns to the
first file.

**The split also fixes the rate format.** `write_phasing_summary` currently emits
`SWITCH_ERROR_RATE`/`FLIP_ERROR_RATE` with `%.6f%%` (`phase.cpp:810`) — a *string* with a trailing
percent sign in a numeric column — and writes its single data row **without a terminating newline**.
Both are worked around by every consumer today ([`D10`](./D10_vcfdist-v3-multiqc-module.md) carries a
`%`-stripping helper; [`D11`](./D11_pipeline-ecosystem-integration.md) §8.2 records them as parse
hazards). The split rewrites this writer anyway, so it emits **bare floats as true fractions** —
`0.00020930`, not `0.020930%` — and terminates the row. Note what this does *not* change: the column
names. A consumer comparing a pre- and post-split run therefore reads the same column with values
100× apart, which is the "name kept, meaning changed" hazard D10 warns about generally. The `%` sigil
in legacy output is self-describing, so D10 can discriminate on the value rather than a column name;
that is the mitigation, and it belongs in D8's release notes as a breaking output change.

## 8. Interaction With Existing Behavior

### 8.1 Strata must not affect contig selection

`intersect_contigs` (`bed.cpp:187`) reads only `g.bed`, so this holds by construction today. The
plan adds no stratum reference there, and §9.2 asserts it: a run whose strata name contigs absent
from `-b` must produce byte-identical `*` rows.

### 8.2 Strata must not affect which variants are evaluated

Membership is computed after the `-b` gate and never consulted by a `continue`. The strongest
statement of this invariant is the regression test in §9.2: `*` rows identical with and without
`-st`.

### 8.3 Zero-overlap warning

The realistic failure mode is an assembly/naming mismatch — a `chr20`-prefixed stratum against
`20`-named reference contigs, or GRCh37 strata against a GRCh38 run — yielding 181 all-zero rows
that read as genuine results.

At load, compare each stratum's contig set against the reference FASTA's contigs and `WARN` by
name for any stratum with **zero** intersection. If *every* stratum has zero overlap, escalate the
message to name the likely cause (contig naming or assembly mismatch), since that pattern is
diagnostic rather than incidental.

A stratum that overlaps contigs but contains no evaluated variants is *not* warned — that is a
legitimate result, and it emits a full zero row rather than being omitted, so a consumer joining on
stratum name never sees a missing key.

## 9. Testing Plan

### 9.1 Unit — `tests/unit/src/test_bed.cpp` (gtest)

**This file is not new: D2 §5.1 already specifies it**, calling `bedData::contains` "the
highest-value function in this file." D5 *extends* that enumeration rather than replacing it, and
two of D2's planned tests collide with D5's changes — see §9.3.

**Manifest parsing:** relative-path resolution against the manifest directory; absolute paths;
comment and blank lines; extra columns ignored; duplicate name → ERROR; name `*` → ERROR; missing
BED → ERROR naming both paths; short line → ERROR with line number; zero-stratum manifest → WARN.

**Compressed input:** the same intervals parse identically from plain, gzip, and bgzip encodings
of one BED.

**`merge()`:** already-merged input is unchanged; unsorted; overlapping; adjacent; nested;
identical duplicates; single interval; empty contig. Assert both the interval list and the
recomputed `size`.

**`contains()` — membership derivation.** D2 §5.1 already covers every `BED_*` return value
(`contains-middle-inside`, `contains-partial-overlap-border`, `contains-between-regions-outside`,
`contains-ins-at-region-end-border`, `contains-spans-multiple-border`, `contains-unknown-contig-offctg`,
and the left-of-all / right-of-all cases). D5 adds only what is new: for each of those inputs,
assert the **derived membership boolean** (§3), so the `INSIDE`-or-`BORDER` rule is pinned
independently of the raw return code. Plus two inputs D2 does not cover: a variant abutting a
region start, and a zero-length query interval.

**Cursor correctness (§6.4):** membership from the cursor sweep matches a brute-force linear scan
over a randomized-but-sorted interval set, across contig transitions, on a contig where a stratum
has no intervals, and with repeated queries at an identical position (the haplotype-loop case).
The brute-force scan is a test oracle only — it is not a shipped code path.

**VCF sort precondition (§6.4)** — these belong in `tests/unit/src/test_variant.cpp` rather than this
file, since the check lives in `parse_variants`: a within-contig backwards position `ERROR`s; the
existing contig-re-entry `ERROR` still fires; a correctly sorted VCF is unaffected.

### 9.2 Integration — `tests/integration/test-integration.yml` (pytest-workflow)

Fixtures under `tests/integration/data/`, all scoped to the existing `chr20.bed` span: a
`strat.tsv` manifest plus (a) a stratum equal to `chr20.bed`, (b) two complementary strata
partitioning that span, (c) an empty-but-on-contig stratum, (d) a `.bed.gz` stratum, (e) a
`20`-named (no `chr` prefix) stratum to trip the §8.3 warning.

The invariant assertions are the substance of this deliverable:

| Assertion | What it protects |
| --------- | ---------------- |
| `*` rows with `-st` are numerically identical to today's output without it | §8.1, §8.2 — the schema change and strata are both inert on existing numbers |
| Stratum (a) reproduces the `*` row exactly | End-to-end membership plumbing, incl. the §6.3 merge copy |
| Strata (b) counts sum to **≥** the `*` row, with equality iff no variant straddles their shared boundary | §3 — the overlap rule, asserted in its true form rather than as a false partition |
| Stratum (c) emits a full zero row with `PREC`/`RECALL` = 1 | §8.3 — matches the existing `query_tot == 0 ? 1` convention (`print.cpp:371-372`); no missing keys |
| Stratum (e) triggers the zero-overlap `WARN` naming that stratum | §8.3 — the guardrail itself |
| Every stratum name in the manifest appears in every widened file | No silent row dropping |
| `STRATA` column in `query.tsv` is consistent with the per-stratum counts | §7.2 vs §7.3 agree |
| `phasing-variants-summary.tsv` `*` row matches the old `phasing-summary.tsv` variant-attributed columns; `phasing-blocks-summary.tsv` matches its contiguity columns | §7.1 — the split loses nothing |

CLI cases: `-st` without `-b`; `-st` with a nonexistent manifest; `-st` with a zero-stratum
manifest; `-sc` accepted; **bare `-s` rejected** (pinning §4.1's no-alias decision).

One extra fixture: a deliberately position-unsorted VCF, asserting the §6.4 precondition `ERROR`
rather than today's silent record-by-record dropping. Worth an integration case and not just a
unit test, because the current failure mode is a *successful* run with wrong denominators.

Verification is `make` in `src/` and `tests/unit/build/`, then `pytest -vv` from `tests/`
(`tests/conftest.py` rebuilds both automatically).

### 9.3 Coordination with D2

Three D2-planned tests interact with D5's changes, and are updated here:

| D2 test | D2 loc | Conflict | Resolution |
| ------- | ------ | -------- | ---------- |
| `contains-no-bed-inside` — asserts `g.bed_exists=false` → `BED_INSIDE` | §5.1 | **Direct contradiction.** D5 §5.3 deletes exactly this early return | D5 replaces it with an equivalent test at the `variant.cpp` call site, plus a `contains()` test asserting the value is now computed from the intervals |
| `bedData::bedData` — classified FIXTURE, "reads file, ERROR path" | §5.1 | D5 §5.1 changes the reader from `ifstream` to htslib | Extend with the plain/gzip/bgzip equivalence cases; the ERROR path is unchanged |
| `add_var` (copy overload) — `copy-roundtrip` "(every field)" | §5.2 | "Every field" silently stops being every field once `strata_bits` exists | Extend `copy-roundtrip` to assert the copied stratum words. This is the unit-level guard for §6.3's critical path |

`bedData::check` (D2 §5.1, PURE) is unaffected: D5 leaves `check()` itself untouched and adds
`merge()` alongside it.

## 10. Downstream Consumers

Widening in place breaks parsers. Three consumers matter, none of them fixed here:

| Consumer | Failure mode | Action |
| -------- | ------------ | ------ |
| MultiQC `vcfdist` module | **Silent**, twice over. `csv.DictReader` tolerates the new `STRATUM` column and then aggregates all 182 strata into one sample, so the report shows wrong precision and recall; and the `*phasing-summary.tsv` search pattern matches neither new file, so the phasing section vanishes without an error while the module otherwise runs normally | Owned by [`D10`](./D10_vcfdist-v3-multiqc-module.md), not by D5: route rows by `STRATUM` treating an absent column as `*`, and register search patterns for both new phasing files |
| `docs/v2.*/09-Outputs.md` | Documents the old schemas | D5 updates the current-version page; [`D9`](./D9_docs-and-tutorial.md) owns the full v3 docs refresh |
| `demo/pr_plot.py` | Aggregates all 182 strata into one curve | Owned by [`D9`](./D9_docs-and-tutorial.md) §2.7 — one `STRATUM == "*"` filter, alongside the demo-output regeneration it already does |

The MultiQC failure is the dangerous one because it is silent in both directions — a wrong number and
a missing section, neither of which raises. It is tracked in D10 rather than here, but the split and
the widening are what create it, so it is recorded on this side too.

## 11. Decision Record

| # | Decision | Rationale |
| - | -------- | --------- |
| 1 | Manifest TSV only; no repeated `--strat NAME=file` flag | GIAB ships the manifests; one code path to validate. Adding flags later is trivial |
| 2 | Overlap (`INSIDE` or `BORDER`) = membership | A variant straddling a repeat edge belongs to that repeat stratum. Divergence from hap.py documented (§3) |
| 3 | Widen aggregate tables in place, `*` row first | One schema per metric, one parse target for D10; matches hap.py's `Subset`. Accepts the downstream breakage in §10 |
| 4 | Per-stratum QUAL sweep in scope | Counters are already per-QUAL; output plumbing only. Enables per-stratum P/R curves in D10 |
| 5 | `STRATA` set column on `query.tsv`/`truth.tsv` | Only way to explain a surprising per-stratum number without re-deriving membership externally |
| 6 | Per-stratum `genotype-errors.tsv` | Mechanically identical to the P/R aggregation |
| 7 | Split `phasing-summary.tsv` into `phasing-variants-summary.tsv` (stratified) + `phasing-blocks-summary.tsv` (genome-level) | The split axis is "can this be restricted to a subset of variants at all". Named for per-variant vs. per-block rather than errors vs. metrics, because D6 adds non-error columns to the first file |
| 7a | `VARIANTS` here, `ASSESSED_PAIRS` in D6 | This deliverable cannot compute the phaseable subset; it exposes the denominator it actually divides by, and D6 swaps it. Both land before rc1, so no user sees the intermediate |
| 7c | The split writers emit bare-float **fractions** and terminate the data row | The current `%.6f%%` (`phase.cpp:810`) puts a string in a numeric column and omits the trailing newline; every consumer works around both. Column names are unchanged, so pre/post values differ 100× — the legacy `%` sigil is the discriminator D10 keys on, and the change is release-noted |
| 7b | A stratified pair-count denominator is attributed, not partitioned | A pair spans two variants that may be in different strata; attributing to the second variant matches `switchflips.tsv`, keeping numerator and denominator consistent, but per-stratum pair counts do not sum to the genome-wide total. Documented, not hidden |
| 8 | `-s` → `-sc`, new `-st`; no `-s` alias | An alias would reinterpret an old `-s 15000` as a path. No in-repo caller affected |
| 9 | gzip via htslib `hts_open`/`hts_getline` | Required for GIAB `.bed.gz`; fixes `-b *.bed.gz`, which fails today with a misleading error |
| 10 | Lenient sort+merge for strata, strict `check()` for `-b` | Third-party BEDs aren't guaranteed sorted; `contains()`'s binary searches *require* merged input |
| 11 | Remove `g.bed_exists` from `contains()` | A method must not consult a global about a different object; membership can't depend on `-b` |
| 12 | Reuse [`D2`](./D2_vcfdist-v3-unit-tests.md) §5.7's `compute_pr_f1` helper (#94) | The stratum axis adds one loop around one helper rather than multiplying duplicated arithmetic |
| 12a | D4's new summary columns are stratified here | They are counts on the same counter arrays; leaving them genome-wide would make a stratified row internally inconsistent |
| 13 | Complementary strata assert `sum >= *`, not `==` | Direct consequence of decision 2; the honest invariant |
| 14 | `WARN` on zero contig overlap | 181 silently-empty strata is the realistic failure mode |
| 15 | No stratification of `summary.vcf` | A GA4GH consumer takes strata externally (D4 §3.3); `query.vcf`/`truth.vcf` no longer exist (#96) |
| 16 | `switchflips.tsv`/`phase-blocks.tsv` unstratified | Interval/event records, not variant records. Deferred |
| 17 | Monotonic cursor only — no binary-search fallback | Both sweep inputs are *established*, not assumed: strata are merged at load (§5.2) and the VCF sort order becomes an enforced precondition (decision 18). Two paths would mean two behaviors to keep in agreement |
| 17a | Extract a pure `classify()` from `contains()`; cursor and `contains()` share it | Two lookup strategies must not mean two copies of the `INSIDE`/`BORDER` decision tree. Keeps `contains()`'s behavior and its D2 §5.1 tests unchanged |
| 18 | `ERROR` on a within-contig backwards VCF position | Makes requirement 2 a precondition instead of an accident of the overlap filter. Independently fixes silent degradation on unsorted input. Behavior change → D8 release notes |

## 12. Implementation Order

Each step is independently reviewable. Step 1 is a pure refactor with no behavior change; steps 2–3
are behavior changes, kept separate so each is reviewable on its own.

1. **`bed.{cpp,h}` I/O** — htslib reader (§5.1), `merge()` + lenient constructor (§5.2), remove the
   `g.bed_exists` guard and short-circuit at the call site (§5.3). Unit tests from §9.1 land here,
   as does the `contains-no-bed-inside` reconciliation in §9.3.
2. **VCF sort precondition** — `ERROR` on a within-contig backwards position (§6.4, decision 18).
   Deliberately ahead of the cursor that relies on it, so the precondition is enforced before
   anything assumes it, and so its release-note-worthy behavior change is reviewable on its own.
3. **CLI rename** — `-s` → `-sc` (§4.1), plus the integration case pinning bare `-s` as rejected.
4. **Manifest loading** — `-st`, `Globals` members, validation, zero-overlap warning (§4.2, §8.3).
5. **Membership** — cursor sweep, bitset storage, and the `add_var` copy path (§6). The §9.2
   "stratum (a) reproduces `*`" assertion is what confirms this step, plus the `copy-roundtrip`
   extension in §9.3.
6. **P/R aggregation and output** — stratum axis, per-stratum `BEST`, widened
   `precision-recall{,-summary}.tsv` and `genotype-errors.tsv` (§7.3).
7. **Phasing split** — `phasing-variants-summary.tsv` / `phasing-blocks-summary.tsv` (§7.1, §7.4).
8. **`STRATA` column** — `query.tsv` / `truth.tsv` (§7.2); `parameters.tsv` provenance rows.
9. **Measurement** — one run against GIAB's full 181-stratum manifest, recording peak RSS, added
    wall-clock in the reading stage, and output file sizes. This measures **memory only**; the
    lookup strategy is settled by §6.4 and is not contingent on the result. The `.bed.gz` files are
    distributed via the GIAB FTP rather than committed to the git repo, so no interval counts are
    available to estimate from, and §6.2's ~120 MB figure covers only the per-variant bitset, not
    the resident interval arrays. If resident size proves unacceptable, the fallback is a second
    pass with one stratum resident at a time — which requires re-reading the VCFs, so measure
    before designing it.
10. **Docs** — `print_usage()` (§4), `docs/v2.5.3/09-Outputs.md` schemas, and
    `docs/v2.5.3/10-Variant-Stratification.md` rewritten around native usage, with the external
    GA4GH-quantifier path ([`D4`](./D4_ga4gh-compatibility.md) §3.3) retained as the interop option a
    user can drive themselves.
