# vcfdist v3.0.0 — Unit-Test Design

|  |  |
| :-- | :-- |
| **Version:** | 0 (Draft) |
| **Authored By:** | Tim Dunn |
| **Status:** | Draft — enumeration of the target unit-test suite |
| **Companion docs:** | `D0_vcfdist-v3-SOW.md`, [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md), [`D8_vcfdist-v3-release.md`](./D8_vcfdist-v3-release.md) |
| **Implements:** | Deliverable #2 (SOW), §5.2 of the design doc |

---

## 1. Purpose & Scope

This document is the design for the comprehensive unit-test suite called for by **Deliverable #2** and sketched in [§5.2 of the design doc](./D1_vcfdist-v3-design-doc.md#52-deliverable-2--comprehensive-unit-tests). It **enumerates** the tests that should exist — it does not write them. For every function in every `src/<name>.cpp`, it lists one minimal test per distinct input case, branch, boundary, or likely bug, so the implementation work is reduced to writing already-specified tests.

**Convention.** For each `src/<name>.cpp` there is one test file `tests/unit/src/test_<name>.cpp`. Each test file `#include`s the corresponding header and exercises that translation unit's functions. This mirrors the existing `tests/unit/src/test_dist.cpp`.

**Guiding principles** (from the task brief):

- **One minimal test per distinct case.** Each test targets a single branch, boundary value, error path, or plausible bug — not a broad scenario.
- **Smallest possible fixtures.** Reference sequences of 4–12 bp; 1–3 variants; single contig unless the test is specifically about multiple contigs. Prefer in-memory construction over on-disk files wherever the API allows.
- **Assert behavior, not implementation.** Each test pins a function's current observable behavior rather than asserting on internal implementation details.

**Out of scope.** End-to-end pipeline behavior, threading correctness, and file-format round-trips remain the province of the existing `pytest`/`pytest-workflow` integration tests (`tests/integration/`). This suite targets *functions*, not the assembled pipeline.

**This suite pins current behavior, and later deliverables amend it.** The four output-affecting
deliverables — #3 record shape and field retention, #4 GA4GH conformance, #5 stratification, #6 unphased
evaluation — each update the cases they change as part of their own work. The amendments known in advance:

| Named case(s) here | Amended by | How |
| :-- | :-- | :-- |
| `pv-unphased-het-skipped` (§5.2) | D6 §7 | Replaced by `pv-unphased-het-kept`; unphased query heterozygotes are evaluated rather than dropped |
| `pa-itrs-ok`, `pa-itrs-too-small-errors` (§5.6) | D6 §4.2.3 | Deleted — `-i` is removed from `parse_args` |
| `pa-supercluster-ok` (§5.6) | D5 §4.1 | Retargeted to `-sc`, plus a new case pinning bare `-s` as **rejected** |
| `contains-no-bed-inside` (§5.1) | D5 §5.3, §9.3 | Replaced — the `g.bed_exists` early return is deleted, so the equivalent assertion moves to the `variant.cpp` call site |
| `copy-roundtrip` "(every field)" (§5.2) | D3 §6.1, D5 §9.3, D6 §7 | Extended each time a per-variant vector is added (`rec_idxs`/`alt_idxs`/`ploidies`, then `strata_bits`, then `is_phased`) — "every field" silently stops being every field otherwise |
| `wsv-header-fields`, `pvi-format-key-string` (§5.2, §5.5) | D3 §4.5, D4 §5.1 | The FORMAT key list changes: `Number=.` per-haplotype fields, `SC` → `INFO/BS`, `BS` → `PBS` |
| `wsv-*` per-haplotype record cases (§5.5) | D3 §4.2 | One record per variant instead of two, so the expectations change rather than the cases disappearing |
| `wps-header-and-row`, `wps-zero-variants-div` (§5.5) | D5 §7.1, D6 §5.3 | `write_phasing_summary` is split in two, then its denominator changes to `ASSESSED_PAIRS` |
| `pa-*` default-value assertions (§5.6) | D8 §4.3 | Any default G12 moves must be updated where a test asserts it |
| `-mr` cases (§5.6) | D8 §2.2 | Coupled to PR #92; if single-pass merges, the flag may go the way of `-i` |
| `compute_pr_f1` extraction + its five cases (§5.7) | — **owned here** | Issue [#94](https://github.com/TimD1/vcfdist/issues/94). Extracted here; [`D5`](./D5_multi-bed-stratification.md) §7.3 reuses the helper for its per-stratum rows |

Each of those documents names its own amendments in its testing section, so the obligation is recorded on
both sides rather than only here.

**The defect annotations below are verified against `dev`.** The defects surfaced while reading the source
were filed as #59–#74. Where a fix has landed, the case asserts the **corrected** behavior —
`set_allele_errtype` (§5.2), `add_variants`'s CIGAR `default:` (§5.2), `split_large_supercluster` (§5.3),
`idx4::operator<` (§5.4), `fix_allele_counts`'s truth loop bound and `write_phasing_summary`'s division
guard (§5.5), `--max-ram`'s zero guard (§5.6), `write_params`'s `fopen` check, and `get_ptr_repr`'s cell
widths (§5.7). The rows marked **bug** are live. Cited line numbers drift as fixes land; re-read each
function before writing its test.

**One exception.** §6 specifies a single *integration* test — a `bcftools` conformance check on
`summary.vcf`. It lives here because it is test work, and it asserts a property no unit test can: that the
assembled pipeline emits a valid VCF.

---

## 2. Framework, Build & CI

### 2.1 Unit Test Framework — GoogleTest

The repository uses **pytest** for launching tests, **GoogleTest** for unit tests, and **pytest-workflow** for integration tests.

### 2.2 File & build layout

```
tests/unit/
  src/
    test_bed.cpp        test_cluster.cpp   test_dist.cpp     (exists)
    test_globals.cpp    test_phase.cpp     test_print.cpp
    test_timer.cpp      test_variant.cpp
    test_helpers.h      test_helpers.cpp   (new — shared fixtures/builders, §4)
  build/
    Makefile            (extend OBJS to add the new test_*.o + test_helpers.o)
```

The current `Makefile` compiles a single `test_dist.o`. It must be extended so `OBJS` includes every `test_<name>.o` and `test_helpers.o`, and a compile rule is added per test file (mirroring the existing `test_dist.o` rule). All test objects link into the one `test_vcfdist` binary run by CI.

### 2.3 CI

`tests/unit/test-unit-tests.yml` already runs `./unit/build/test_vcfdist` under `pytest-workflow`; no change is needed beyond the binary growing more `TEST()` cases. The GitHub Actions `test.yml` invokes the same path. New test files are picked up automatically once added to `OBJS`.

---

## 3. Testability Taxonomy & Prioritization

Every function is classified to set expectations for how much test value it offers and how much scaffolding it needs:

| Class | Meaning | Test value |
| :-- | :-- | :-- |
| **PURE** | No I/O; callable with plain values or hand-built structs. | **Highest** — test exhaustively. |
| **FIXTURE** | Needs a constructed `ctgVariants` / `superclusterData` / `fastaData` / small VCF/BED file, but is deterministic. | High — test the distinct cases with shared builders (§4). |
| **HEAVY** | Threading, full-pipeline orchestration, or RAM scheduling. | Low as a *unit* — one smoke/guard test at most; covered by integration. |
| **FORMATTING** | Produces only console/file text. | Low — assert schema/headers and the one or two pieces of embedded logic; otherwise skip. |

**Priority order for implementation:** PURE first (they need almost no scaffolding and catch the most arithmetic/boundary bugs), then FIXTURE (once the shared builders in §4 exist), then a thin layer of FORMATTING/HEAVY guard tests.

The genuinely high-value PURE targets, worth calling out up front:

- **Alignment primitives** (`dist.cpp`): `wf_ed`, `wf_swg_align`, `wf_swg_max_reach`, `generate_str`, `calc_ng50`.
- **Genotype/allele-count logic** (`variant.cpp`): `get_vartype`, `set_allele_errtype`, `var_on_hap`, `set_var_calcgt_on_hap`, `calcgt_is_swapped`, CIGAR parsing in `add_variants`.
- **Clustering index/range math** (`cluster.cpp`): `get_min/max_ref_pos`, `get_supercluster_range`, `get_supercluster_split_location`, `split_cluster`, `get_next_variant_info`.
- **BED interval classification** (`bed.cpp`): `contains`, `check`.
- **Scalar helpers**: `qscore`, `get_ptr_repr` (`print.cpp`); `parent_path` (`globals.cpp`); the phasing DP `phase()` and `calculate_ng50` (`phase.cpp`).

---

## 4. Shared Test Infrastructure

Placed in `tests/unit/src/test_helpers.{h,cpp}` and reused across files.

### 4.1 In-memory builders (no files on disk)

Most of `ctgVariants`/`ctgSuperclusters` state is public, so tests build structures directly rather than parsing VCFs:

- **`make_fasta(ctg, seq)` → `shared_ptr<fastaData>`** — populate the public `fasta[ctg]`/`lengths[ctg]` maps on a heap `fastaData` (bypassing the `FILE*` ctor). Sequences are already uppercase in real use; include a lowercase variant for the mixed-case tests.
- **`make_ctgVariants(ctg, {var…})` → `shared_ptr<ctgVariants>`** — construct via `ctgVariants(ctg)` then `add_var(...)` per variant; a variant descriptor carries `{pos, rlen, type, ref, alt, gt, qual, phase_set, supercluster}`. Also exposes setters for the per-hap lanes (`errtypes[HAP][i]`, `sync_group`, `callq`, `ref_ed`, `query_ed`, `credit`) and cluster metadata (`clusters`, `left_reaches`, `right_reaches`, `nc`).
- **`make_ctgSuperclusters(qvars, tvars)` → `shared_ptr<ctgSuperclusters>`** — sets `callset_vars[QUERY]`/`[TRUTH]`.
- **`make_superclusterData(...)`** — assembles `contigs`, `lengths`, `ploidy`, `superclusters[ctg]`; used by `supercluster()`, `sort_superclusters()`, and the phasing tests.
- **`make_graph(scs, ref, ctg, truth_hap)`** — thin wrapper over the `Graph` ctor; a subsequent `calc_prec_recall_aln` call yields a valid `ptrs` map for `calc_prec_recall`/`get_ptr_repr` tests.
- **`alloc_reach_offs(qlen, tlen, x, o, e)`** — allocates the `MATS*(max(x,o+e)+1)*(qlen+tlen-1)` int buffer initialized to `-2` that `wf_swg_max_reach` requires from its caller.

### 4.2 Minimal on-disk fixtures (only where htslib/file parsing is unavoidable)

Two stable, suite-wide inputs are checked in under `tests/unit/data/` (kept tiny):

- **`tiny.fasta`** — one contig `chr1` (~12 bp, e.g. `ACGTACGTACGT`), a second contig `chrX` for the ploidy-exception path, and coordinates chosen so a variant can sit at position 0 (contig-start off-by-one).
- **`tiny.bed`** — a central interval on `chr1` (e.g. `chr1  2  8`) so `contains` returns each of INSIDE/OUTSIDE/BORDER at known positions, plus a contig absent from the BED for OFFCTG. Separate malformed BEDs (flipped, zero-length, unsorted, overlapping, shared-boundary) for `check`.

The many single-sample VCFs that `parse_variants` needs are **generated programmatically at test time, not checked in**. `parse_variants` must read a real file because it calls htslib `bcf_open`, so a helper emits one:

- **`write_tmp_vcf(records, opts) → std::string path`** — writes a minimal valid single-sample VCF (configurable `##contig` / `##FILTER` / `##FORMAT` header lines and `TRUTH`/`QUERY` sample name via `opts`; one text line per record) to a unique path under the temp dir (§4.3) and returns it. Each `parse_variants` case (§5.2) builds only the header lines and records it needs inline and passes the returned path straight to `parse_variants`. This keeps the one genotype, filter, or malformed header a case actually targets visible next to its assertion, instead of scattering ~40 near-identical `.vcf` files across `tests/unit/data/`. Written files live in the temp dir and are removed by its teardown.

### 4.3 Cross-cutting harness

- **Death tests.** `ERROR(...)` calls `std::exit(1)`; the `argc<4` help/version/citation paths in `parse_args` call `std::exit(0)`. Negative cases use GoogleTest `EXPECT_EXIT` / `EXPECT_DEATH` asserting the exit code and a message substring.
- **stderr / stdout capture.** `WARN`/`INFO` write to stderr; `print_*`/`print_version`/`print_usage`/`print_citation` write to stdout. Use `testing::internal::CaptureStdout/Stderr` (already available with gtest) and assert on substrings/anchors, never full-text equality for large blocks.
- **Globals reset.** Nearly every function reads the global `g` (`g.max_qual`, `g.sv_threshold`, `g.credit_threshold`, `g.max_size`, `g.min_qual`, `g.max_supercluster_size`, `g.reach_min_gap`, `g.thread_nsteps`, `g.ram_steps`, …). A `GlobalsGuard` RAII fixture saves and restores `g` around each test and sets `g.verbosity = 0` to silence logging. Tests that depend on a threshold set it explicitly.
- **Temp output dir.** File-writers point `g.out_prefix` at the session scratch dir; tests read the file back and assert header + row counts/columns.

---

## 5. Per-File Test Enumeration

Legend for the "Targets" column: the branch, boundary, or bug the case exists to catch. Line numbers reference the files as read for this document.

### 5.1 `test_bed.cpp` (from `bed.cpp`)

**`bedData::bedData(const std::string&)` — FIXTURE (reads file, ERROR path).**

| Test | Input | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `open-fail-errors` | nonexistent path | `ERROR` exit | file-not-found guard |
| `parses-three-cols` | `chr1  2  8` | one region, `size==6` | happy parse |
| `ignores-extra-cols` | `chr1  2  8  name  0  +` | region parsed, extras dropped | 3-column-only contract |

**`bedData::add` — PURE.**

| Test | Input | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `add-new-contig` | add on unseen `chr1` | `contigs=={chr1}`, one region | new-contig branch |
| `add-appends-same-contig` | two regions same contig | `regions[chr1].n==2` | append branch |
| `add-accumulates-size` | regions len 6 + len 4 | `size==10` | size accumulation |

**`bedData::check` — PURE (ERROR/WARN).**

| Test | Input | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `check-flipped-errors` | region `8..2` | `ERROR` | `stop<start` |
| `check-zero-length-errors` | region `5..5` | `ERROR` | `stop==start` |
| `check-unsorted-errors` | `10..20` then `1..5` | `ERROR` | out-of-order |
| `check-overlap-errors` | `2..8` then `5..10` | `ERROR` | overlap |
| `check-shared-boundary-warns` | `2..5` then `5..8` | `WARN` (not fatal) | mergeable-adjacent |
| `check-valid-passes` | `2..5` then `7..9` | no error/warn | happy path |

**`bedData::contains` — PURE (the highest-value function in this file).**

| Test | Input | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `contains-no-bed-inside` | `g.bed_exists=false` | `BED_INSIDE` | early-return short-circuit |
| `contains-flipped-errors` | `stop<start` | `ERROR` | interval validation |
| `contains-unknown-contig-offctg` | contig not in BED | `BED_OFFCTG` | missing-contig |
| `contains-before-all-outside` | var ends ≤ first start | `BED_OUTSIDE` | left-of-all |
| `contains-after-all-outside` | var starts ≥ last stop | `BED_OUTSIDE` | right-of-all |
| `contains-middle-inside` | var strictly inside one region | `BED_INSIDE` | central case |
| `contains-ins-at-region-end-border` | `TYPE_INS` at `stop-1` | `BED_BORDER` | insertion-at-boundary special case |
| `contains-between-regions-outside` | var in gap between two regions | `BED_OUTSIDE` | inter-region gap |
| `contains-partial-overlap-border` | var straddling one boundary | `BED_BORDER` | `start_idx<0` / `stop_idx` off-end |
| `contains-spans-multiple-border` | var covering ≥2 regions | `BED_BORDER` | multi-region span |

**`bedData::operator std::string` — FORMATTING.** One test: non-empty bed renders each contig + `start-stop` lines.

**`intersect_contigs` — FIXTURE/HEAVY.**

| Test | Input | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `ic-bed-drops-extraneous` | query has contig absent from BED | contig erased from query | BED-driven pruning |
| `ic-bed-missing-in-fasta-errors` | BED contig absent from FASTA | `ERROR` | fasta-required guard |
| `ic-bed-adds-empty-contig` | BED contig absent from a VCF | empty `ctgVariants` added, ploidy 0 | empty-contig injection |
| `ic-nobed-query-only-contig-warns` | no BED, contig only in query | `WARN`, added empty to truth | query-only path |
| `ic-nobed-truth-only-contig-warns` | no BED, contig only in truth | `WARN`, added empty to query | truth-only path |
| `ic-nobed-fasta-missing-errors` | no BED, truth contig absent from FASTA | `ERROR` | fasta-required-without-BED |
| `ic-ploidy-mismatch-warns` | same contig ploidy 1 vs 2 | `WARN` | ploidy-consistency check |

### 5.2 `test_variant.cpp` (from `variant.cpp`)

**`ctgVariants` ctor — PURE:** `ctor-sets-ctg`, `ctor-zeroes-n`, `ctor-allocates-two-phase-lanes` (asserts all six per-hap lane vectors sized 2), `ctor-empty-ctg`.

**`add_var` (copy overload) — PURE:** `copy-roundtrip` (every field), `copy-preserves-hap1-hap2-distinct` (guards HAP1/HAP2 argument transposition), `copy-second-of-two`, `copy-appends-not-overwrites`.

**`add_var` (full overload) — PURE:** `add-all-fields`, `add-qual-capped` (`var_qual>max_qual`→`max_qual`), `add-qual-below-cap`, `add-qual-negative` (no lower clamp — pins behavior), `add-pushes-phase-defaults` (`PHASE_NONE/PHASE_NONE/AC_UNKNOWN`), `add-header-defaults` (`supercluster==-1`, `calc_gt==GT_REF_REF`, `errtypes==ERRTYPE_UN`), `add-lane-lengths-track-n`, `add-ins-rlen-zero`, `add-del-empty-alt`.

**`get_vartype` — PURE (SNP/INDEL/SV classification):**

| Test | Input | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `vt-sub-snp` | `TYPE_SUB` | `VARTYPE_SNP` | SUB branch |
| `vt-small-ins-indel` | INS, alt len `sv_threshold-1` | `VARTYPE_INDEL` | INS/INDEL |
| `vt-boundary-ins-sv` | INS, alt len `== sv_threshold` | `VARTYPE_SV` | strict-`<` off-by-one |
| `vt-small-del-indel` | DEL, ref len `sv_threshold-1` | `VARTYPE_INDEL` | DEL/INDEL |
| `vt-boundary-del-sv` | DEL, ref len `== sv_threshold` | `VARTYPE_SV` | DEL boundary |
| `vt-ins-uses-alt` | INS tiny alt, large ref | `VARTYPE_INDEL` | INS reads `alts` size |
| `vt-del-uses-ref` | DEL large ref, tiny alt | `VARTYPE_SV` | DEL reads `refs` size |
| `vt-cpx-falls-to-sv` | `TYPE_CPX` | `VARTYPE_SV` | else fall-through |

**`set_allele_errtype` — PURE (all AC_ERR_* transitions):** one test per code — `ac-0-to-1`, `ac-0-to-2`, `ac-1-to-0`, `ac-1-to-1`, `ac-1-to-2`, `ac-2-to-0`, `ac-2-to-1`, `ac-2-to-2`, plus `ac-refref-refref-unknown` (returns **and stores** `AC_UNKNOWN` — since #59, `variant.cpp:195` writes `ac_errtype[vi]` on every path) and `ac-haploid-unknown`.

**`calcgt_is_swapped` — PURE (ERROR path):** `swap-equal-false`, `swap-orig-hom-false`, `swap-orig-ref-false`, `swap-calc-ref-false`, `swap-het-opposite-01-10-true`, `swap-het-opposite-10-01-true`, `swap-orig01-calc-hom-credit-hap1-true`, `swap-orig01-calc-hom-credit-hap2-false`, `swap-orig10-calc-hom-credit-hap2-true`, `swap-credit-tie-false` (strict `>`), `swap-unexpected-errors` (death).

**`var_on_hap` — PURE (ERROR path):** `onhap-haploid-alt-both`, `onhap-10`, `onhap-01`, `onhap-hom-both`, `onhap-refref-none`, `onhap-haploid-ref-none`, `onhap-calc-flag-selects`, `onhap-hap-gt-1-errors` (death).

**`set_var_calcgt_on_hap` — PURE (state machine, many ERROR paths):** the 18 transitions enumerated — every valid set/unset from REF_REF, REF_ALT1, ALT1_REF, ALT1_ALT1 on each hap, each invalid transition's `ERROR` (death), the `GT_MISSING`/haploid default `ERROR`, and `setgt-errors-suppressed-with-ignore` (all invalid transitions become no-ops when `ignore_errors=true`).

**`add_variants` (CIGAR parsing) — no tests; slated for deletion.** `variantData::add_variants`
(`variant.cpp:594`, declared `variant.h:120`) **has no callers anywhere in `src/`** — it is residue of the
removed VCF-normalization/realign path, as [`D5`](./D5_multi-bed-stratification.md) §6.4 independently
establishes when ruling it out as a stratum-membership entry point. No cases are written for it: the
function ships to no one, and [`D8`](./D8_vcfdist-v3-release.md) §4.1 records it as a known target for the
G9 dead-code sweep. If it is ever revived, its CIGAR stride and run-merge logic is worth enumerating then.

**`variantData` ctor — PURE:** `vd-default-callset-query`, `vd-two-hap-maps`, `vd-empty-members`.

**`write_vcf` / `print_variant` / `print_var_info` / `print_var_empty` / `print_var_sample` — FORMATTING/FIXTURE.** Assert the merge logic and coordinate conventions rather than full text:

- `write_vcf`: `wv-homozygous-merge` (single `1|1`), `wv-het-two-records` (`1|0`+`0|1`), `wv-haploid-gt` (`"1"`), `wv-hap1-only` (`1|0`), `wv-hap2-only` (`0|1`), `wv-merge-order` (min-pos first), `wv-empty` (header only).
- `print_var_info`: `pvi-sub-pos-plus-1`, `pvi-ins-flanking-base` (index `pos-1`), `pvi-del-flanking-base`, `pvi-format-key-string` (regression-lock the `GT:BD:BC:…` key list), `pvi-unknown-type-errors` (death), `pvi-ins-at-pos0` (off-by-one at contig start).
- `print_variant`: `pv-sub`, `pv-ins-flanking` (index `pos`, contrast with `print_var_info`'s `pos-1`), `pv-del-flanking`, `pv-missing-contig-errors` (death), `pv-unknown-type-errors` (death).
- `print_var_sample`: `pvs-credit-1-tp-gm`, `pvs-credit-0-query-fp`, `pvs-credit-0-truth-fn`, `pvs-credit-at-threshold-tp-lm` (`>=` boundary), `pvs-partial-query-fp-lm`, `pvs-partial-truth-fn-lm`, `pvs-ref-ed-zero-dots` (**both RD and QD gated on `ref_ed==0`**, pin behavior), `pvs-nonzero-eds`, `pvs-query-switch-flip`, `pvs-truth-phase-dots`.
- `print_var_empty`: `pve-query-newline`, `pve-truth-tab`, `pve-embeds-sc-pb`.

**`parse_variants` — FIXTURE (needs tiny VCFs; many ERROR/WARN paths).** Grouped by branch; each row is one minimal VCF emitted inline by `write_tmp_vcf` (§4.2):

*Header/validation:* `pv-invalid-callset-errors`, `pv-contig-missing-idx-errors`, `pv-multisample-errors`, `pv-selected-filter-absent-warns`, `pv-unsorted-contig-errors`, `pv-seqnames-null-errors` (goto-cleanup path).

*Filter & quality:* `pv-filter-fail-skipped`, `pv-no-filters-pass`, `pv-below-minqual-skipped`, `pv-qual-nan-zero`, `pv-gq-int`, `pv-gq-float-fallback`, `pv-gq-missing-zero`.

*Genotype:* `pv-no-gt-header-warns-monoploid`, `pv-haploid-alt`, `pv-haploid-ref-skipped`, `pv-0-1`, `pv-1-0`, `pv-1-1-hom`, `pv-1-2-compound`, `pv-2-1`, `pv-0-2-other`, `pv-polyploid-errors`, `pv-ploidy-mismatch-warns`, `pv-ploidy-mismatch-chrX-silent`. Missing (`.`) alleles are enumerated separately below.

*Phase set:* `pv-ps-not-in-header-warns`, `pv-ps-missing-on-het-warns`, `pv-ps-present`.

*Allele filtering:* `pv-unphased-het-skipped`, `pv-unphased-hom-allowed`, `pv-spanning-deletion-skipped`.

*Missing (`.`) alleles:* two independent code paths touch a missing allele — the genotype-histogram tally (`variant.cpp:943`, `:948`) and the per-haplotype parse loop (`:1020`) — and they disagree, so the cases are enumerated individually. VCFs come from `write_tmp_vcf` (§4.2); the histogram and summary-warning rows use `testing::internal::CaptureStderr` (§4.3) and assert on substrings only, since these counters are function-local and print-only (`ntypes` likewise, consumed at `:1247`). The four rows marked **bug** follow this document's document-don't-enforce convention for known defects — cf. `pvs-ref-ed-zero-dots` (§5.2) and `gpr-query-node-wide` (§5.7) — pinning current output so CI stays green until each is fixed under its own issue.

| Test | Input GT | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `pv-half-missing-hap1-kept` | `1\|.` | one variant on HAP1, `orig_gt==GT_ALT1_REF` | `.` acts as REF on the other hap; record kept (`:1014`, `:1156`) |
| `pv-half-missing-hap2-kept` | `.\|1` | one variant on HAP2, `orig_gt==GT_REF_ALT1` | mirror case; a phased `.` carries the phase bit, so it clears the unphased-het guard (`:1031`) |
| `pv-half-missing-unphased-skipped` | `./1` | no variants, counted in `unphased_gt_total` | unphased-het guard fires before the present allele can be kept |
| `pv-both-missing-skipped` | `.\|.` | no variants | full no-call dropped at `:1020` |
| `pv-haploid-missing-skipped` | `.` (`ngt==1`) | no variants | haploid no-call dropped |
| `pv-half-missing-histogram-label` | `1\|.` | histogram row is `.\|.` | **bug:** `:948` sets `GT_MISSING` if *either* allele is missing, so a half-call is indistinguishable from a no-call |
| `pv-half-missing-warn-says-skipped` | `1\|.` | warning claims `skipped` though the variant was kept | **bug:** message wording at `:1204` |
| `pv-both-missing-warn-double-counts` | one `.\|.` record | warning reports `2`, not `1` | **bug:** `unknown_allele_total` counts haplotypes, not records (`:1024`) |
| `pv-haploid-missing-histogram-alt` | `.` (`ngt==1`) | histogram row is `1` (haploid ALT), not `.\|.` | **bug:** `bcf_gt_allele()` returns a truthy `-1` at `:943`, labeling a haploid no-call `GT_ALT1` |
| `pv-missing-vs-star-ref-tally` | `1\|.` vs `1\|*` | `*` bumps `ntypes[hap][TYPE_REF]`, `.` does not | inconsistent per-hap type accounting between the `*` branch (`:1042`) and the `.` branch (`:1020`) |

*Classification (prefix/suffix trim + CPX split):* `pv-snp`, `pv-refcall-skipped`, `pv-insertion-trim`, `pv-deletion-trim`, `pv-mnp-suffix-sub`, `pv-cpx-equal-len-split`, `pv-cpx-ins-like-split`, `pv-cpx-del-like-split`, `pv-lowercase-uppercased` (mixed-case).

*Region & size & overlap:* `pv-outside-skipped`, `pv-border-skipped`, `pv-offctg-skipped`, `pv-inside-kept`, `pv-too-large-skipped` (`>g.max_size`), `pv-overlap-skipped`, `pv-two-ins-same-pos-skipped`, `pv-homref-downgrade-on-skip`, `pv-het-imbalance-warns`.

### 5.3 `test_cluster.cpp` (from `cluster.cpp`)

**`ctgSuperclusters` ctor — PURE:** `ctor-size` (`callset_vars.size()==CALLSETS`), `ctor-nulls`.

**`get_min_ref_pos` — PURE:** `min-both-empty` (`INT_MAX-1`), `min-query-only`, `min-truth-only`, `min-both-query-smaller`, `min-both-truth-smaller`, `min-uses-start-index-only` (does *not* scan — pins the sorted assumption), `min-tie`.

**`get_max_ref_pos` — PURE:** `max-both-empty` (`INT_MAX`, asymmetric — no `+1`), `max-query-single` (`pos+rlen+1`), `max-scans-full-range`, `max-cross-callset`, `max-rlen-zero-ins`, `max-truth-only`.

**`get_supercluster_range` — PURE:** `range-both-empty-degenerate` (`{INT_MAX,-1}`), `range-query-single` (both off-by-one ends), `range-multi-cluster`, `range-cross-callset-min-max`, `range-start-idx-oob-errors` (death, `>` not `>=`), `range-end-idx-oob-errors` (death), `range-one-callset-empty`.

**`get_next_variant_info` — PURE:** `next-both-exhausted` (`callset==-1`), `next-query-only`, `next-truth-only`, `next-query-earlier`, `next-truth-earlier`, `next-tie-prefers-query` (strict `<`), `next-end-pos-of-winner`.

**`get_supercluster_split_location` — PURE (complex arithmetic — high value):** `split-fewer-than-two-empty`, `split-no-gaps-empty` (dense region scores 0 → returns `{}`), `split-single-gap`, `split-prefers-central-gap`, `split-edge-gap-low-score`, `split-log2-div-zero-guard` (`size_reduction_factor==1.0`), `split-index-after-curr`, `split-cross-callset-gap`.

**`split_cluster` — mutating, plain-vector:** `sc-at-existing-boundary-noop`, `sc-midcluster-inserts` (`nc++`, one new entry per lane), `sc-reach-reassignment` (which cluster keeps which reach), `sc-later-breakpoints-incremented`, `sc-both-callsets-independent`, `sc-clust-idx-minus-one-safety`.

**`split_large_supercluster` — mutating, plain-vector:** `sls-already-small-noop`, `sls-one-split`, `sls-multiple-splits` (monotonic breakpoints), `sls-no-valid-split-bails` (retains the oversized supercluster and **surfaces** it — #64 replaced the silent path with an explicit report), `sls-mutates-end-indices` (documented side effect), `sls-breakpoint-shift-after-insert`.

**`sort_superclusters` — FIXTURE:** `sort-empty-query-skipped` (keyed on QUERY only), `sort-nscs-count` (`max(q,t)+1`), `sort-small-low-bucket`, `sort-large-last-bucket-warn` (`>g.max_ram`), `sort-ctg-sc-paired`, `sort-len-lower-upper-bound`, `sort-empty-callset-nc-zero-guard`.

**`load_and_merge_callset_vars_across_haps` — FIXTURE:** `merge-empty-contig-skipped`, `merge-hom-variant`, `merge-het-hap1-first`, `merge-het-tie-prefers-ins`, `merge-het-tie-default-hap1`, `merge-hap1-only`, `merge-hap2-only`, `merge-two-clusters-separate`, `merge-two-clusters-overlap`, `merge-sentinel-appended`.

**`supercluster` — FIXTURE:** `scl-empty-skipped`, `scl-single-cluster`, `scl-two-far-two-scs`, `scl-query-truth-overlap-one-sc`, `scl-chained-overlap`, `scl-oversized-split-invoked`, `scl-var-count-stat-guard` (empty-callset `nc==0` guard), `scl-brks-advance-remainder`.

**`superclusterData` ctor — HEAVY:** `sd-contig-union-dedup`, `sd-samples-filenames-order`, `sd-query-only-contig`, `sd-truth-only-contig`, `sd-end-to-end-smoke`.

**`wf_swg_cluster` — HEAVY/FIXTURE:** `wsc-no-vars-returns`, `wsc-single-variant`, `wsc-two-adjacent-merge`, `wsc-two-far-separate`, `wsc-max-itrs-break`, `wsc-reach-min-gap-boundary`, `wsc-contig-edge-clamp`.

### 5.4 `test_dist.cpp` (from `dist.cpp` — extends the existing file)

**Existing coverage:** `TestNG50.TestNG50Calc` covers `calc_ng50` with four cases (exact-50%, two-block tie, single over-threshold, never-reaches). **Gaps to add** below.

**`contains` (set & map templates) — PURE:** `contains-set-present/absent/empty`, `contains-set-idx4` (real key type + hash), `contains-map-present/absent/empty`, `contains-map-idx4`.

**`calc_ng50` — PURE (add the missing cases):** `ng50-empty`, `ng50-total-zero`, `ng50-empty-and-zero`, `ng50-single-below`, `ng50-exact-even-ties`, `ng50-unsorted-input`.

**`generate_str` — PURE (haplotype synthesis):** `gs-no-variants` (ref slice), `gs-snp`, `gs-ins`, `gs-del`, `gs-cpx`, `gs-minqual-filters-sub` (filtered var keeps ref base), `gs-minqual-filters-del` (filtered DEL keeps ref bases), `gs-prefix-skip` (skip `< beg_pos`), `gs-boundary-end` (exclusive end), `gs-missing-contig-errors` (death).

**`wf_ed` — PURE (unit-cost edit distance, top primitive):** `ed-identical` (0), `ed-single-sub`, `ed-single-ins`, `ed-single-del`, `ed-empty-query` (=len), `ed-empty-truth` (=len), `ed-both-empty` (0), `ed-two-subs`, `ed-all-different`, `ed-symmetry` (`wf_ed(a,b)==wf_ed(b,a)`).

**`wf_swg_align` — PURE (gap-affine score; keep both strings ≥1):** `swg-identical` (0), `swg-single-sub` (`=x`), `swg-gap-len1` (`=o+e`), `swg-gap-len2-affine` (`=o+2e`, not `2(o+e)`), `swg-two-separate-gaps` (`=2(o+e)`), `swg-sub-vs-gap-choice`, `swg-penalty-sensitivity` (varies with `x`).

**`wf_swg_max_reach` — PURE (needs `alloc_reach_offs`):** `reach-full-match`, `reach-partial-budget` (truncated by `max_score`), `reach-score-allows-sub`, `reach-reverse-flag`, `reach-gap-extension`, `reach-main-diag-block` (the `off+1<main_diag_off` guard).

**`idx4` operators & `std::hash<idx4>` — PURE:** `idx4-eq-neq`, `idx4-hash-usable`, `idx4-operator-less-strict-weak-ordering` (assert irreflexivity, asymmetry and transitivity over the `qni`/`tni`/`qi`/`ti` lexicographic order at `dist.h:100-105` — #66 made this a valid strict weak ordering, so the case now **enforces** rather than documents).

**`Graph` ctor / `get_truth_pos` — FIXTURE:** `graph-ref-only`, `graph-single-query-snp`, `graph-skip-evaluated`, `graph-query-pointers`, `graph-truth-hap-filter`, `graph-truth-linear-chain`, `graph-overlapping-query-vars`, `graph-ref-string`; `truthpos-first-node`, `truthpos-multi-node` (per-node `-1` for `'_'` prefix), `truthpos-accumulate`.

**`calc_prec_recall_aln` — FIXTURE:** `pra-identical-ref` (0), `pra-single-sub` (1), `pra-single-indel` (1), `pra-matching-variant-zero`, `pra-endpoint-reached`, `pra-empty-queue-errors` (death).

**`calc_prec_recall` — FIXTURE:** `cpr-default-all-fp`, `cpr-tp-perfect`, `cpr-credit-boundary` (`>=g.credit_threshold`), `cpr-fn-truth-only`, `cpr-ref-dist-zero-guard` (div-by-zero), `cpr-callq-min`, `cpr-sync-group-increment`.

**`evaluate_variants` — FIXTURE:** `ev-max-retries-zero` (default single pass), `ev-simple-tp`, `ev-fn-indel-retry` (`max_retries≥1`), `ev-exclude-size-sort`, `ev-sub-never-excluded`.

**HEAVY (guard tests only):** `precision_recall_wrapper` → `prw-empty-range-early-return`; `precision_recall_threads_wrapper` and `Graph::print` → no unit tests (integration).

### 5.5 `test_phase.cpp` (from `phase.cpp`)

**`phase()` — FIXTURE (the switch/flip DP — highest value in this file):**

| Test | Input pattern (one phase set unless noted) | Expect | Targets |
| :-- | :-- | :-- | :-- |
| `phase-classify-orig` | orig `1\|0`, calc `1\|0` | `phases[i]==PHASE_ORIG` | same-orientation |
| `phase-classify-swap` | orig `1\|0`, calc `0\|1` | `PHASE_SWAP` | opposite |
| `phase-classify-none` | `1\|1` | `PHASE_NONE` | unphaseable |
| `phase-perfect-block` | all ORIG | 0 switch, 0 flip | baseline |
| `phase-all-swap-block` | all SWAP | 0 switch, 0 flip, backtrace starts SWAP | uniform-swap ≠ N switches |
| `phase-single-flip` | `ORIG,SWAP,ORIG` | 1 flip at middle, 0 switch | isolated flip |
| `phase-single-switch` | `ORIG,ORIG,SWAP,SWAP` | 1 switch, 0 flip | mid-block switch |
| `phase-switch-plus-flip` | mixed | 1 switch + 1 flip | independent accounting |
| `phase-boundary-free` | `ORIG,ORIG‖SWAP,SWAP` across PS boundary | swap costs 0, not a switch | boundary exemption |
| `phase-multiple-blocks` | two consistent PS | resolved independently | multi-block DP |
| `phase-empty-contig` | `n==0` | no crash, 0/0 | empty guard |
| `phase-none-not-flipped` | NONE between ORIGs | not counted as flip | homozygous exclusion |
| `phase-flips-switches-ascending` | ≥2 each | ascending index order after reverse | output ordering |
| `phase-invalid-phase-errors` | inject bad `phases[i]` | `ERROR` (death) | DP-cost guard |
| `phase-pb-phases-assignment` | any | `pb_phases[i-1]` matches trace | off-by-one |

**`fix_phase_set_tags()` — FIXTURE (must run before `phase`/`fix_allele_counts`):** `fps-no-phase-sets-contig`, `fps-backfill-leading-zeros`, `fps-propagate-unphased-middle`, `fps-new-ps-span-reset`, `fps-same-ps-extends-end`, `fps-final-span-pushed`, `fps-both-callsets`, `fps-ng50-reported`.

**`fix_allele_counts()` — FIXTURE:** `fac-unknown-errors` (death), `fac-1-to-1-tallied`, `fac-1-to-2-tallied`, `fac-0-to-2-tallied`, `fac-force-1-1-keeps-gt`, `fac-force-1-1-swaps-hap-data`, `fac-2-to-1-hap1-better`, `fac-2-to-1-hap2-better`, `fac-2-to-1-tie-orig`, `fac-2-to-1-tie-swap`, `fac-0-to-1-hap1-better`, `fac-0-to-1-hap2-better`, `fac-0-to-1-tie-orig`, `fac-0-to-1-tie-swap`, `fac-truth-fn-2-to-0`, `fac-truth-fn-1-to-0-hap2`, `fac-truth-fn-1-to-0-hap1`, `fac-truth-loop-bound` (truth loop iterates `vi < tvars->n`, `phase.cpp:435` — pins the corrected bound).

**`calculate_ng50(break_on_switch, break_on_flip)` — FIXTURE:** `ng-no-breaks`, `ng-switch-break`, `ng-switchflip-break`, `ng-flip-ignored-when-off`, `ng-switch-ignored-when-off`, `ng-empty-returns-zero`, `ng-threshold-selection`, `ng-never-reaches-half-zero`, `ng-final-block-appended`, `ng-switch-and-flip-same-sc`, `ng-next-vi-regression-errors` (death).

**`phaseblockData` ctor — FIXTURE/HEAVY:** `pbd-copies-metadata`, `pbd-boundaries-single-ps`, `pbd-boundaries-multi-ps`, `pbd-boundary-at-index-zero`, `pbd-empty-query-contig`, `pbd-pipeline-order`.

**FORMATTING (schema/logic only):**
- `write_summary_vcf`: `wsv-header-fields` (lock the 15 FORMAT IDs + `#CHROM…TRUTH\tQUERY`), `wsv-empty-query-skipped`, `wsv-flip-error-swap-orig`, `wsv-flip-error-orig-swap`, `wsv-phase-none-no-flip`, `wsv-matched-truth-query`, `wsv-positional-tie-diff-vars`, `wsv-truth-only-hi-mapping`, `wsv-haploid-gt`, `wsv-indel-pos-decrement`.
- `write_switchflips`: `wsf-header`, `wsf-empty-skipped`, `wsf-flip-two-rows` (FLIP_BEG/FLIP_END), `wsf-switch-err-row`, `wsf-clean-switch-no-row`, `wsf-switch-and-flip-priority`, `wsf-next-vi-regression-errors` (death).
- `write_phasing_summary`: `wps-header-and-row`, `wps-zero-variants-guarded` (`variants == 0` yields a rate of `0`, `phase.cpp:808-809`; the unguarded division was fixed under #79).
- `write_genotype_error_summary`: `wges-header-columns`, `wges-row-per-vartype`.

### 5.6 `test_globals.cpp` (from `globals.cpp`)

**`parent_path` — PURE:** `pp-nested` (`a/b/c`→`a/b/`), `pp-single-dir`, `pp-bare-filename` (`""`), `pp-trailing-slash`, `pp-root` (`/`→`/`), `pp-absolute-file` (`/file`→`/`), `pp-empty` (`""`), `pp-dotslash`, `pp-dotdot`.

**`create_directory` — FIXTURE (scratch dir; death on failure):** `cd-single-dir` (final component **not** created), `cd-nested-trailing-slash`, `cd-already-exists-ok` (EEXIST tolerated), `cd-absolute-skips-root`, `cd-mkdir-failure-errors` (death), `cd-no-slash-noop`.

**`init_timers` — FIXTURE (mutates `g.timers`):** `it-populates` (6 named timers), `it-empty-input`, `it-appends-not-clears`, `it-writes-global-not-this` (writes `g.timers`, not `this->timers`).

**`parse_args` — FIXTURE (argv arrays; opens files; many death paths).** Grouped:

*Short-circuit (`argc<4`, all `exit(0)`):* `pa-argc1-usage`, `pa-help-short`, `pa-help-long`, `pa-version-short` (`-v` means version here, verbosity in main loop), `pa-version-long`, `pa-citation`, `pa-unknown-short-argc`.

*Mandatory positional:* `pa-optional-before-mandatory-warns`, `pa-query-open-fail-errors`, `pa-truth-open-fail-errors`, `pa-ref-open-fail-errors`, `pa-all-mandatory-ok`.

*Verbosity pre-pass:* `pa-verbosity-valid`, `pa-verbosity-missing-errors`, `pa-verbosity-non-numeric-errors`, `pa-verbosity-oob-low-errors`, `pa-verbosity-oob-high-errors`.

*Per-flag (happy + missing-value + bad-numeric + bound violation):* `-b` (`pa-bed-ok`, `pa-bed-missing-errors`, `pa-bed-bad-file-errors`), `-p` (`pa-prefix-relative-dotslash`, `pa-prefix-absolute-kept`, `pa-prefix-dotslash-kept`, `pa-prefix-dotdot-kept`, `pa-prefix-missing-errors`), `-f` (`pa-filter-single`, `pa-filter-comma-list`, `pa-filter-trailing-comma`, `pa-filter-missing-errors`), `-l` (`pa-largest-ok`, `pa-largest-missing-errors`, `pa-largest-non-numeric-errors`), `-sv` (`pa-sv-ok`, `pa-sv-too-small-errors`, `pa-sv-missing/non-numeric-errors`), `-q` (`pa-minq-ok`, `pa-minq-negative-errors`, `pa-minq-missing/non-numeric-errors`), `-mq` (`pa-maxq-ok`, `pa-maxq-missing/non-numeric-errors`), `-n` (`pa-no-output`), `-h`/`--version`/`-ci` in main loop (`pa-help-main-loop`, `pa-version-main-loop`, `pa-citation-main-loop` — none exit), `-x`/`-o`/`-e` (`pa-*-ok`, `pa-*-negative-errors`), `-i` (`pa-itrs-ok`, `pa-itrs-too-small-errors`), `-md`/`-mr` (`pa-*-ok`, `pa-*-missing/non-numeric-errors`), `-s` (`pa-supercluster-ok`, `pa-supercluster-too-small-errors`), `-t` (`pa-threads-ok`, `pa-threads-too-small-errors`), `-ct` (`pa-credit-ok`, `pa-credit-zero-errors`, `pa-credit-above-one-errors`, `pa-credit-upper-inclusive`), `-r` (`pa-ram-ok`, `pa-ram-trailing-units-warns`, `pa-ram-zero-errors`, `pa-ram-negative-errors` — `globals.cpp:399` rejects `max_ram <= 0`, so zero and negative now share one guard), `pa-unknown-option-errors`, `pa-verbosity-skipped-in-main-loop`.

*Cross-field validation:* `pa-maxq-lt-minq-errors`, `pa-maxsize-lt-sv-warns`, `pa-supercluster-lt-maxsize-plus2-errors`, `pa-maxsize-1-warns-snps-only`.

*Thread/RAM scheduler:* `pa-thread-steps-default-64` (`{64,32,…,1}`, 7 steps), `pa-thread-steps-single`, `pa-thread-steps-non-power-of-two` (integer halving).

**FORMATTING (stdout capture):** `print_version` → `pv-format`; `print_usage` → `pu-required-section`, `pu-lists-documented-flags`, `pu-omits-commented-flags` (`-i/-x/-o/-e/-mr` hidden), `pu-interpolates-defaults`; `print_citation` → `pc-both-formats`.

### 5.7 `test_print.cpp` (from `print.cpp`)

**`qscore` — PURE (highest value in this file):** `qs-p-1` (0), `qs-p-0p1` (10), `qs-p-0p01` (20), `qs-p-0p001` (30), `qs-clamp-high` (`1e-11`→100), `qs-p-0` (`log10(0)`→+inf→clamped 100 — pin the clamp value), `qs-negative-input` (NaN domain — document), `qs-p-gt-1` (→0 via lower clamp), `qs-rounding` (float, not integer Phred). **Note the code clamps to 100.0** (not 60) — the test locks that.

**`get_ptr_repr` — PURE (9 branches):** `gpr-not-found` (`"  ."`), `gpr-up` (`"  |"`), `gpr-left` (`"  _"`), `gpr-diag` (`"  \\"`), `gpr-invalid-same-matrix` (`" ?1"`), `gpr-query-node` (`"^NN"` zero-pad), `gpr-query-node-wide` (a node id ≥ 100 widens to 4 chars — `std::setw(2)` pads but never truncates, so **this** is the genuine grid-misalignment case), `gpr-truth-node` (`"<NN"`), `gpr-invalid-other` (`" ?2"`), plus `gpr-width-consistency` (every branch returns exactly 3 characters for node ids < 100, the invalid markers included, so the grid stays aligned).

**Color wrappers `GREEN/RED/BLUE/YELLOW/PURPLE` (int/char/string) — PURE (trivial).** These are **not** `isatty()`-dependent — they unconditionally emit escape codes; the `isatty` check lives only in the `COLOR_*` macros. Test **one representative per color** to lock the code constant (31/32/33/34/35) and guard copy-paste errors; the int/char/string overloads add little. `green-str-empty` covers the empty-string edge.

**`write_params` — FORMATTING:** `wp-filters-join` (`A,B,C`), `wp-filters-empty` (guards `g.filters[0]` OOB), `wp-filters-single` (no trailing comma), `wp-fopen-fail` (NULL check present at `print.cpp:132-134` — death test on an unwritable `g.out_prefix`).

**HEAVY / metrics extraction:** `write_precision_recall` computes the P/R/F1 arithmetic inline, and does so *twice* — once in the per-quality curve loop and again in the summary loop (`print.cpp:371–373` and `:438–440`), the two copies differing only cosmetically (`(precision+recall)` vs. `precision+recall > 0`). This duplicated, file-local math is refactored out into a single pure helper so it can be unit-tested directly and both call sites collapse onto one implementation:

```cpp
struct prec_recall_f1 { float precision; float recall; float f1; };
prec_recall_f1 compute_pr_f1(int query_tp, int query_fp, int truth_tp, int truth_fn);
```

The helper takes all four counts: precision keys off the query TP/FP total and recall off the truth TP/FN total. It lives in `print.{h,cpp}` beside its callers and becomes the single source of truth for P/R/F1 (reused by the #5 per-stratum rows and the #10 MultiQC module); `write_precision_recall`/`write_results` retain only formatting and I/O. Tested in `test_print.cpp` — PURE:

- `pr-precision-zero-query` — `query_tp+query_fp==0` → precision `1` (vacuous-perfect convention).
- `pr-recall-zero-truth` — `truth_tp+truth_fn==0` → recall `1`.
- `pr-f1-zero-denominator` — `precision+recall==0` (zero TP on both sides) → f1 `0`.
- `pr-f1-normal` — prec `1.0`, recall `0.5` (query 1 TP / 0 FP, truth 1 TP / 1 FN) → f1 `0.667`.
- `pr-both-empty` — no query and no truth variants → prec `1`, recall `1`, f1 `1` (pins the doubly-vacuous case).

`print_wfa_ptrs`/`print_graph_ptrs` are debug dumps — skip (their only real logic, `get_ptr_repr`, is tested directly). `b2s` is file-local and trivial — test only if `print.cpp` is compiled into the binary (it is).

### 5.8 `test_timer.cpp` (from `timer.cpp`)

**`timer` — PURE (ERROR paths):** `t-get-name` (ctor name), `t-default-name` (`"default"`), `t-start-twice-errors` (death), `t-stop-not-running-errors` (death), `t-total-while-running-errors` (death), `t-total-accumulates` (two start/stop intervals sum; assert monotonic non-negative rather than an exact wall-clock value), `t-total-seconds-scale` (ns→s division). **FORMATTING:** `timer::print` and `write_runtime` → one schema/no-crash test each.

### 5.9 `main.cpp` — no dedicated test file

`main()` is the pipeline entry point and is exercised end-to-end by the integration tests, not unit tests. Its only unit-testable content is the **string tables**. Once those move to `strings.cpp`, add `test_strings.cpp` asserting each table's size matches its constant count and its index/label mapping (e.g. `error_strs[ERRTYPE_TP]=="TP"`, `type_strs.size()==TYPES`, `gt_strs[GT_ALT1_ALT1]` correct). These guard the string-lookup arrays every other module indexes into.

---

## 6. VCF-Conformance Integration Test (bcftools)

The one integration test this document specifies (§1). It asserts that `summary.vcf` is a *valid* VCF, which no unit test can establish — the file is produced by the assembled pipeline, and the failures worth catching are format-level rather than function-level.

**Motivation.** Once [#96](https://github.com/TimD1/vcfdist/issues/96) lands, `summary.vcf` is the **only** VCF vcfdist emits (`src/main.cpp:127`), and [`D3`](./D3_retain-info-format-fields.md) rewrites how it is built. `summary.vcf` is not a copy of the input with annotations added — it is synthesized from vcfdist's internal arrays by `phaseblockData::write_summary_vcf` (`phase.cpp:20`), with site-level columns hardcoded in `print_var_info` and alleles reconstructed from *normalized*, re-anchored representations. A synthesized VCF undergoing a format rewrite is exactly where a parse-level contract test pays for itself.

**Placement.** `tests/integration/`, in the existing `pytest-workflow` harness, against the checked-in chr20 fixture (`tests/integration/data/`).

**Checks.** There is no `bcftools check` subcommand; the test composes four existing ones, each chosen for a specific failure mode:

| Command | Failure mode it catches |
| :-- | :-- |
| `bcftools view -Ou summary.vcf > /dev/null` | Malformed header, undeclared INFO/FORMAT key, wrong field arity or type — the direct risk from D3's field pass-through and GA4GH tag emission |
| bgzip + `bcftools index` | Records out of coordinate order, or otherwise not indexable |
| `bcftools norm --check-ref e -f <ref>` | REF that disagrees with the reference FASTA — the re-anchoring bug class described above |
| `bcftools query -f` over `BD`/`BK`/`BC` | Tags present as text but not properly declared/typed in the header |

**CI change required.** Unlike the unit tests (§2.3), this one does not come for free. `.github/workflows/test.yml` installs HTSlib only, so `bcftools` must be added to the workflow. `src/Dockerfile` already builds bcftools 1.17 from source, so the version to match is established.

---

## 7. Rollout Plan

Sequence the work so the suite compiles and adds value incrementally:

1. **Infrastructure** — extract `strings.cpp`, add `test_helpers.{h,cpp}` (§4, incl. the `write_tmp_vcf` emitter), extend the `Makefile` `OBJS` + per-file rules, confirm the empty new test files build and CI stays green.
2. **PURE tier** — `test_timer`, `test_bed` (`add`/`check`/`contains`), the `dist` primitives (`wf_ed`, `wf_swg_align`, `wf_swg_max_reach`, `generate_str`, `calc_ng50` gaps), `variant` genotype/type logic, `cluster` index/range math, `qscore`/`get_ptr_repr`, `parent_path`. Highest bug-per-line, minimal scaffolding.
3. **FIXTURE tier** — once the builders are proven: `Graph`/`calc_prec_recall*`, `supercluster`/`load_and_merge`/`sort_superclusters`, the `phase()` DP + `fix_*` + `calculate_ng50`, `parse_variants` (with the tiny VCF fixtures), `parse_args` (with death tests).
4. **FORMATTING/HEAVY guard tier** — schema-lock the writers, add the handful of smoke/guard tests, extract `compute_pr_f1` (§5.7) and test it.
5. **VCF-conformance integration test** (§6) — independent of the tiers above, and written once [`D3`](./D3_retain-info-format-fields.md) has settled the `summary.vcf` format it asserts against.
