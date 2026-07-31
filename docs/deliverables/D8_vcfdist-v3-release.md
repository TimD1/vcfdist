# vcfdist v3.0.0 — Release & Launch Announcement (Plan)

|  |  |
| :-- | :-- |
| **Version:** | 0 (Draft) |
| **Authored By:** | Tim Dunn |
| **Status:** | Draft — plan of record for SOW Deliverable #8 |
| **Issue:** | [#51](https://github.com/TimD1/vcfdist/issues/51) (D8: v3.0.0 Release & Launch Announcement) |
| **Companion docs:** | `D0_vcfdist-v3-SOW.md` #8, [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md) §5.8, [`D7_vcfdist-v3-benchmarking-design.md`](./D7_vcfdist-v3-benchmarking-design.md) §10 |

---

## 1. Purpose & Scope

This document is the plan for **Deliverable #8** of the v3.0.0 SOW: cutting the
v3.0.0 release and announcing it. It answers three questions — *what must be true before we tag*, *how the
tag is cut*, and *how the release is announced* — and records the decisions behind each.

**In scope.** The pre-release checklist and its exit criteria (§2); the version sequence and format freeze (§3); the
pre-release engineering work this deliverable owns (§4); `CHANGELOG.md` and release notes (§5); the release
procedure, which replaces [`docs/v2.5.3/13-Release-Instructions.md`](../v2.5.3/13-Release-Instructions.md)
(§6); distribution via GitHub/Zenodo, Docker Hub, and bioconda (§7); and the launch announcement (§8).

**Out of scope.** Deliverables #9 (documentation & sandbox.bio), #10 (MultiQC module), and #11 (nf-core &
snakemake-wrappers) all consume a *released* artifact, and each is planned in its own document. The
engineering content of #2–#7 is likewise out of scope — §2 states only the conditions this deliverable
checks, not their design.

---

## 2. Pre-Release Checklist

### 2.1 What must be true before tagging

Deliverables #2–#7 are complete by the time this one starts, so the conditions below are the ones this
deliverable owns. Each is stated as an observable condition, so readiness is checkable rather than argued.

| # | Condition | Exit criterion |
| :-- | :-- | :-- |
| C1 | **Extension deliverables merged and green** | #2–#6 merged, CI green, and the defect issues #59–#74 closed with each test asserting the corrected behavior. The `*` rows of every stratified file are numerically identical with and without `-st`, and `precision-recall-summary.tsv` reconciles: `QUERY_TP + QUERY_FP + QUERY_UNK == QUERY_TOTAL`. |
| C2 | **Benchmarking complete** | **Genome-wide** run finished on the `linux-64` AWS host ([`D7`](./D7_vcfdist-v3-benchmarking-design.md) §9) — the chr20 smoke run is not sufficient. Produces the figures §8 draws on, the concordance evidence that v3 is safe to adopt, and the parameter evidence C6 decides on. The per-site counting delta (D4 §4) and the supercluster-boundary shift (D6 §4.2) are each measured and reported separately. |
| C3 | **Output-format issue cleanup** | #94 (`compute_pr_f1`), #95 (`LOCATION` trailing space), and #96 (remove `query.vcf`/`truth.vcf`) closed — all three change output the freeze locks in. |
| C4 | **bcftools VCF conformance** | The integration test specified in [`D2`](./D2_vcfdist-v3-unit-tests.md) §6 green in CI against the frozen output format. |
| C5 | **Dead-code sweep** | Sweep run and its findings resolved — deleted or justified in writing (§4.1). |
| C6 | **Defaults chosen from benchmarking evidence** | Every default the freeze locks is either measured or explicitly justified as unmeasured (§4.3). |
| C7 | **Repository & version hygiene** | §4.2 checklist complete. |
| C8 | **PR #92 disposition** | Decided per §2.2 — merged, deferred, or merged behind a default-off flag. Not left open across the tag. |

C3–C6 close **before rc1**, since the freeze (§3) locks output, code surface, and defaults alike, and C4
must be green *on* the frozen format.

### 2.2 PR #92 (single-pass precision/recall) — path to merge

[PR #92](https://github.com/TimD1/vcfdist/pull/92) replaces the greedy retry-based FN-dropping loop with a
single-pass mechanism (truth-graph reference-allele bypass nodes, `skip_cost = ceil((1-ct)*len)`, weighted
wavefront forward pass). All nine plan tasks are implemented and per-task reviewed.

**The chr20 evidence points one way.** Measured against `dev` at default settings, single-pass is net
slightly worse: **SNP −32 TP, INDEL −9 TP** (~0.04% lower recall *and* precision), a 42:1 loss:gain ratio.
Every loss traces to one mechanism — the entangled-loci credit-splitting residual, where a truth-only
neighbor indel dilutes an adjacent clean SNP or indel, demoting `gm` to `lm`. Three fix attempts (a
coincident parallel query-reference node; restricting the un-swallow drop to truth variants the min-cost
path actually bypassed; a bypass-vs-align rounding tie-break) each failed to recover the losses, and two
made chr20 marginally worse. The predicted broad insertion-handling win is **not visible on chr20 at all**.

The open question is therefore not whether the implementation is correct — it is whether the predicted win
exists at genome scale, which the C2 run answers:

1. **Rebase `91_D12_td_single-pass-prec-recall` onto `dev`.** `dev` has moved substantially — 13 defect
   fixes plus the default changes in #93, which altered `credit_threshold`. Because `skip_cost` is defined
   in terms of `ct`, a default change is not neutral to this branch, so the chr20 comparison is re-taken
   after the rebase.
2. **Park the uncommitted experiments.** The merge candidate is the reviewed Task 1–9 implementation only;
   the failed attempts stay as artifacts (`task-10-fix.patch`) referenced from the PR. The
   `select_fn_drop_candidate()` extraction and `skip_cost_saved()` are independently unit-tested and can be
   kept on their own merits, separately from the behavior change they were written to support.
3. **Register `vcfdist-v3-singlepass` as its own tool id in the D7 harness**, so one genome-wide run yields
   the baseline and the single-pass comparison together across all three query sets, rather than paying for
   a second AWS run.
4. **Decide against a rule fixed in advance:**
   - **Merge** if genome-wide SNP+INDEL F1 is neutral-or-better on all three query sets *and* the predicted
     insertion/SV win appears.
   - **Merge behind a default-off flag** if the SV/insertion win is real but the small-variant regression
     persists — this preserves the work without changing default behavior.
   - **Defer to v3.1** if neither holds. The residual then needs core forward-pass/backtrack surgery, which
     is #12-scale work and correctly outside a v3.0.0 tag.

## 3. Version Sequence, Cadence & Freeze

**Sequence:** `3.0.0-b0` (current, in `src/globals.h:84`) → `3.0.0-rc1` → `3.0.0`.

**Why an rc.** The v2→v3 output surface changes substantially — edit-distance and normalization outputs
removed, `query.vcf`/`truth.vcf` removed, new stratification columns, rewritten `summary.vcf`, changed
defaults. Those are exactly the changes that break downstream consumers in ways internal benchmarking and CI
do not surface, and they are cheap to adjust before a stable tag and expensive after. An rc buys real-user
signal at the point where it is still actionable.

**What rc1 is:**

- A GitHub release marked **pre-release**, tagged `v3.0.0-rc1`.
- A Docker Hub tag `timd1/vcfdist:3.0.0-rc1`. **The `latest` tag is not moved.**
- **No bioconda submission.** BiocondaBot watches releases; keeping rc1 flagged as a pre-release avoids a
  recipe update for a version that is not the recommended install.
- Announced only where it reaches people who will exercise it — the GitHub release, the relevant open
  issues, and direct outreach to known users. **No blog or social post at rc**; those belong to the stable
  tag (§8).

**Feedback window: two weeks** from the rc1 tag. Short enough not to stall the release, long enough for a
user to run a real callset and report back.

**Freeze policy, effective at the rc1 tag:**

| Change class | Allowed after rc1? |
| :-- | :-- |
| Bug fix with no output-schema change | Yes |
| Documentation, release notes | Yes |
| New CLI flag, renamed flag, changed default | No |
| New/renamed/removed output file, column, or VCF tag | No |
| Changed numeric results from an engine change | No |

Any change in a "No" row after rc1 forces **rc2** and restarts the two-week window. Stating this up front is
the entire value of the freeze — without a defined cost, "one more small format tweak" is always locally
reasonable.

---

## 4. Pre-Release Engineering Work

Unlike §6–§8, which are procedure, this section is work that must be built. The bcftools conformance test
(C4) is the other addition to the SOW's #8 scope, but it is an integration test rather than release work, so
it is specified in [`D2`](./D2_vcfdist-v3-unit-tests.md) §6 and appears here only as a checklist item.

### 4.1 Dead C++ code sweep (C5)

**Motivation.** v3 removed two whole code paths (edit distance, VCF normalization), and four commits
(`df5b321`, `e029965`, `df54be1`, `d33d3a1`) already removed dead code found incidentally. That is evidence
there is more, and that finding it has so far been ad hoc. The sweep converts it into a repeatable step, run
before the freeze because deletions change the surface the freeze locks.

**Tooling — one primary, one cross-check, one supporting signal:**

- **Primary: `cppcheck --enable=unusedFunction`** over `src/`. Cross-translation-unit unused-function
  detection, apt-installable, and cheap enough to add to CI as a non-blocking job so the sweep does not
  regress.
- **Cross-check: linker garbage collection.** Build with `-ffunction-sections -fdata-sections
  -Wl,--gc-sections -Wl,--print-gc-sections`; the linker reports exactly which functions were collected.
  This is precise where cppcheck is heuristic. GNU `ld` only — macOS `ld64` uses `-dead_strip` with
  different reporting — so it runs on the Linux CI runner, not the M3 dev host.
- **Supporting only: gcov/lcov** over the D2 suite. This finds code that is *reachable but never exercised*,
  which is a test-gap signal, not dead code. It feeds D2; it is not grounds for deletion.

**One target is already known**, so it does not need rediscovering: `variantData::add_variants`
(`variant.cpp:594`, declared `variant.h:120`) has **no callers anywhere in `src/`** — residue of the
removed VCF-normalization/realign path, confirmed independently by
[`D5`](./D5_multi-bed-stratification.md) §6.4 while ruling it out as a stratum-membership entry point.
[`D2`](./D2_vcfdist-v3-unit-tests.md) §5.2 writes no tests for it on the strength of this. Delete it under C5 unless there is a reason to keep it, in which case the
annotation below applies and D2's enumeration should be restored.

**Resolution rule.** Every finding is either deleted or annotated with why it is retained (e.g. debug-only
helpers reached by a flag no test sets). "Reported and ignored" does not satisfy C5.

### 4.2 Repository & version hygiene (C7)

- `src/globals.h:84` `VERSION` and the version string in `README.md` must agree. **`README.md:108`
  currently reads `vcfdist v2.5.0`** against a source version of `3.0.0-b0` — a full major behind, and
  four minor releases behind the latest v2 tag (v2.6.4). Add a check so they cannot diverge again.
- Prune the tracked developer scripts `src/run` and `src/profile`, or move them out of the build directory.
- Confirm the `.gitignore` covers the build artifacts and stray outputs sitting in `src/` (`*.o`, `vcfdist`,
  `summary.vcf`, `query.vcf`, `truth.vcf`, `diff.txt`, `vcfdist.dSYM`) so a release build from a clean clone
  is reproducible.
- Regenerate `demo/output.txt` against the release build. The release procedure's first step is comparing
  demo output against the previous release; with v3's output changes, that diff is expected to be large and
  should be reviewed deliberately rather than accepted wholesale.
- Regenerate Doxygen (`cd src && doxygen Doxyfile`) and confirm no warnings.

### 4.3 Choose default parameter values from benchmarking evidence (C6)

**Motivation.** The freeze (§3) makes a changed default an rc2-forcing event, so v3.0.0 ships whatever
defaults are in `globals.h` at rc1 and is stuck with them for the major version. Several were last set
by judgement rather than measurement, and #93 already changed a batch of them without a published
basis. D7 produces exactly the evidence needed, on the datasets the release claims accuracy on, so the
decision costs a configuration sweep rather than a separate study.

**The parameters in scope**, with what each decision needs:

| Default | Current | What the evidence has to show |
| ------- | ------- | ----------------------------- |
| `credit_threshold` (`-ct`) | `0.98` | The P/R operating point across the sweep. This is the most consequential default in the tool — it decides what counts as a TP — and D7 already pins it explicitly for the v2-vs-v3 comparison, so the sweep is a matter of adding points rather than new machinery. |
| `max_dist` (`-md`) | `100` | Runtime/RAM against accuracy loss on the datasets where alignment is expensive. |
| `max_retries` (`-mr`) | `0` | Whether retries buy recall at acceptable cost — and, if PR #92 merges (C8), whether the flag survives at all. |
| `max_supercluster_size` (`-sc`) | `15000` | The splitting threshold's effect on runtime/RAM tail behavior versus accuracy on dense regions. |
| `sv_threshold` (`-sv`) | `50` | Convention, not measurement: 50 bp matches the community indel/SV boundary and D7's own size classes. Justify and keep rather than sweep. |
| `--stringency` | `gm` | Convention: GA4GH Comparison Method #3, matching hap.py's default (D4 §6.4). Justify and keep. |

**The harness currently runs at non-default `-md 1000 -mr 10`** (`workflows/config/config.yml`), against
defaults of `100` and `0`. That has to be resolved before the figures are published, because it is the
difference between "these are vcfdist's numbers" and "these are vcfdist's numbers at settings a user
would not get". Two acceptable outcomes: the run uses defaults, or the run sweeps and the sweep *is* the
evidence for changing them. Either way the published operating point and the shipped defaults must
agree, and D7 §2's `credit_threshold` pin sets the precedent for how to record it.

**Exit criterion (C6).** For every row above: either a measured basis recorded with the figure it came
from, or an explicit written note that the value is conventional and unmeasured. "It has always been
that" does not satisfy it; "it matches the community boundary and D7's size classes, so it is not
swept" does.

---

## 5. CHANGELOG & Release Notes

**Create `CHANGELOG.md` at the repository root.** None exists today; release history lives only in GitHub
release descriptions. [Keep a Changelog](https://keepachangelog.com) format, newest first, with v3.0.0 as the
first entry — earlier versions can be backfilled from the GitHub releases later, and that backfill does not
gate the tag.

The v3.0.0 entry's **Removed** and **Changed** sections are the migration path for v2 users, and per
[`D1`](./D1_vcfdist-v3-design-doc.md) §5.8 this is ordinary release-note content, not a separate migration
document. Each breaking change gets a one-line before/after. The list below is limited to changes verified in
the repository; D3–D6's entries are written as each lands.

| Class | Change | Source |
| :-- | :-- | :-- |
| **Removed** | Edit-distance metrics and the `distance-summary.tsv`, `distance.tsv`, `edits.tsv` outputs — migrate to the P/R metrics | D1 §3.2 |
| Removed | VCF normalization and the `-rq`/`--realign-query`, `-rt`/`--realign-truth`, `-ro`/`--realign-only` flags; v3 assumes normalized input | D1 §3.2 |
| Removed | `query.vcf` and `truth.vcf` outputs, and with them the `-p` prefix hazard where a colliding prefix overwrote the input VCFs — these were written at `src/main.cpp:53-56` before evaluation ran | #96 |
| Removed | `-i`/`--max-iterations`. Cluster iteration is fixed at 1 and is now a correctness precondition rather than a tuning knob; the flag was undocumented and no released version exposed it | D6 §4.2.3 |
| Removed | `phasing-summary.tsv`, replaced by `phasing-variants-summary.tsv` (per-variant, stratified) and `phasing-blocks-summary.tsv` (genome-level contiguity) | D5 §7.1 |
| Removed | `FORMAT/SC` in `summary.vcf` → moved to `INFO/BS`, the GA4GH benchmarking-superlocus field | D4 §5.1 |
| **Changed — CLI** | `-s` is now `-sc`/`--max-supercluster-size`; **`-s` is not retained as an alias**, because a silent alias would reinterpret an old `-s 15000` as a stratification file path. An invocation using `-s` fails loudly | D5 §4.1 |
| Changed — CLI | New `-st`/`--stratification <manifest.tsv>` for multi-BED stratification | D5 §4.2 |
| Changed — CLI | New `--stringency {lm,am,gm,pm}` (default `gm`) and `--require-phase`/`--no-require-phase` (default on) | D4 §6.6 |
| Changed — CLI | `-b`/`--bed` now accepts gzipped/bgzipped BED files | D5 §5.1 |
| Changed — CLI | Default parameter values, now chosen from the benchmarking evidence rather than by judgement — each documented with its basis | #93, §4.3 |
| **Changed — behavior** | A coordinate-unsorted input VCF is now a hard `ERROR` instead of being silently degraded record-by-record by the overlap filter. Previously such a run *succeeded* with quietly wrong denominators | D5 §6.4 |
| Changed — behavior | **Counting is per site, not per haplotype copy.** A homozygous variant contributes 1 to TP/FP/FN, not 2, matching `hap.py` and `vcfeval`. **Every count and every derived rate changes**; figures from v2 and from prior publications are not directly comparable. The het-alt residual (+1 per `1\|2` site vs. those tools) is reported in the benchmarking report | D4 §4 |
| Changed — output | Switch/flip error rates are emitted as **bare-float fractions** with a terminating newline, replacing the `%.6f%%` percent-*string* and missing final newline of `phasing-summary.tsv`. The column names are unchanged, so a value that read `0.020930%` now reads `0.00020930` — a 100× difference on the same column name. Consumers discriminate on the legacy `%` sigil | D5 §7.4 |
| Changed — behavior | Switch/flip error rates are now taken over `ASSESSED_PAIRS` — adjacent phaseable-heterozygote pairs per phase block — rather than over every merged query variant. No released version emitted a rate column, so no published figure changes | D6 §5.1 |
| Changed — behavior | Overlapping query variants are evaluated instead of dropped; overlapping *truth* variants are still dropped | D6 §4.1.3 |
| Changed — behavior | Unphased genotypes are written with `/` rather than `\|`, so unphased input no longer round-trips as phased output | D6 §5.3 |
| Changed — output | `summary.vcf` emits **one record per variant**: a homozygous variant is a single `GT=1\|1` record rather than two per-haplotype records. Het-alt (`1\|2`) records remain split | D3 §4.2 |
| Changed — output | Per-haplotype `summary.vcf` fields (`BC`, `RD`, `QD`, `SG`) are `Number=.` lists in emitted-`GT` order; `BD`/`BK` are one value per site | D3 §4.5, D4 §5.4 |
| Changed — output | `FORMAT/BS` (phase-block state) renamed to `FORMAT/PBS` — same meaning, new name, because `BS` is GA4GH's superlocus field. **The tag most likely to be silently mis-parsed by an existing consumer** | D4 §5.1 |
| Changed — output | Ploidy is tracked per variant, so mixed-ploidy contigs render correct genotypes; the non-standard `ploidy=` attribute is dropped from `##contig` lines | D3 §4.4 |
| Changed — output | `STRATUM` is the leading column of every aggregate TSV, `*` (all regions) first. Without `-st` the numbers are unchanged with one added constant column | D5 §7.1 |
| Changed — output | `LOCATION` in `query.tsv`/`truth.tsv` no longer emits a trailing space | #95 |
| **Added** | Genotype and allele-count error summaries | D1 §3.2 |
| Added | Unphased query evaluation, with the phased fraction of the callset reported explicitly | D6 |
| Added | Multi-BED stratification via a hap.py-compatible manifest | D5 |
| Added | INFO/FORMAT pass-through, and `BD=N` retention of variants vcfdist does not evaluate, each carrying a `VCFDIST_*` `FILTER` tag naming the reason | D3 |
| Added | GA4GH-conformant `summary.vcf`: `INFO/BS`, the `am` match kind, and a `pm` (phased match) extension implementing GA4GH Comparison Method #4 — the first implementation of that method. A consumer applying the Method #3 `BK`→`BD` table must treat `pm` as a match | D4 |
| Added | `FP_GT`, `FP_AL`, `QUERY_UNK`, `FRAC_NA`, and TiTv / het-hom ratio columns in the precision-recall summaries | D4 §7.1 |
| Added | `STRATA` set column in `query.tsv`/`truth.tsv`, naming each variant's strata | D5 §7.2 |

**Caveat to state once, prominently:** a GA4GH consumer fed `summary.vcf` reports vcfdist's *thresholded*
decisions — `BD` is a hard label, so fractional `BC` credit cannot survive the export. Those counts are a
compatibility export, not the vcfdist result (D4 §1.2).

Also list the 16 defect fixes (#59–#74) under **Fixed**, grouped rather than enumerated one per line.

---

## 6. Release Mechanics

[`docs/v2.5.3/13-Release-Instructions.md`](../v2.5.3/13-Release-Instructions.md) is the current procedure. It
still works, but has four gaps that matter for v3:

1. **No CHANGELOG step** — the file does not exist yet (§5).
2. **The Docker image is not reproducible.** `src/Dockerfile` runs `git clone https://github.com/TimD1/vcfdist`
   with no tag or ref, so the image contains whatever `master` HEAD was at build time. The procedure works
   only because the merge-to-master step happens first. Fix: build with a `--build-arg` pinning the tag, and
   verify `vcfdist --version` inside the built image before pushing.
3. **The base image and dependencies are stale.** `ubuntu:20.04` is past end of standard support, and HTSlib
   1.17 is several releases behind. Bumping both is low-risk and appropriate for a major release; it needs a
   build-and-test pass, so it is scheduled before rc1, not between rc1 and the tag.
4. **No version-consistency check** — nothing prevents `globals.h` and `README.md` from disagreeing, which is
   how the README came to be two majors stale (§4.2).

**Updated procedure**, to be written to `docs/v3.0.0/13-Release-Instructions.md` as part of the #9 wiki
refresh:

1. Verify every condition in §2.1, including C6 — the shipped defaults match the operating point the published
   figures were generated at (§4.3).
2. Run `demo/demo.sh` and the benchmarking harness against the previous release; review and update
   `demo/output.txt`.
3. Bump `VERSION` in `src/globals.h` **and** the version string in `README.md`; confirm they agree.
4. Finalize `CHANGELOG.md` for the version.
5. Cache the wiki: `cd docs && git clone https://github.com/TimD1/vcfdist.wiki.git v3.0.0 && rm -rf v3.0.0/.git`.
   (Depends on #9 having refreshed the wiki for v3 first.)
6. Commit (`git add src/globals.h README.md CHANGELOG.md docs/`), merge `dev` → `master`, push.
7. Tag and create the GitHub release. For rc1, mark it **pre-release**.
8. Build and push the Docker image, pinned to the tag. Move `latest` **only** for the stable release.
9. For the stable release only: verify the BiocondaBot recipe PR (§7).

---

## 7. Distribution

**GitHub release.** The tag is the primary artifact. The repository carries a
[Zenodo DOI badge](https://zenodo.org/badge/latestdoi/472945373), so a Zenodo record and DOI are minted
automatically on publish — confirm the record's metadata and authorship after the stable tag. Attach the D7
benchmarking figures to the release, or link them into the repository, so the release notes' accuracy claims
are checkable.

**Docker Hub.** `timd1/vcfdist:3.0.0-rc1` at rc; `timd1/vcfdist:3.0.0` plus `latest` at the stable tag.
Verify `vcfdist --version` in the built image before pushing.

**bioconda.** BiocondaBot opens a recipe-update PR automatically on a new (non-pre-)release. It is not
fully hands-off for a major version — verify before merging:

- The build/test command in the recipe still exists in v3. Any recipe test invoking a removed flag
  (e.g. `--realign-query`, `-i`, or a bare `-s`) fails and must be updated. Note the current recipe's
  only test is `vcfdist --version`, which passes regardless — so this is a manual review, not a
  check the bot performs (§3, P3).
- The compiler and HTSlib pins match what v3 actually requires.
- The recipe's `summary`/`about` text is not describing v2 behavior.

---

## 8. Launch Announcement

**Two surfaces, one narrative** (decision: full post on the Fulcrum blog, concise notes on GitHub linking to it).

**Fulcrum Genomics blog — the primary post.** Draws directly on the #7 figures. Outline:

1. **What vcfdist is and the problem it solves** — representation-tolerant matching, and phasing-error
   reporting that neither incumbent provides.
2. **What v3 changes** — unphased query support (removing the adoption barrier), multi-BED stratification,
   GA4GH-conformant output, annotation pass-through. Framed as *capability parity plus vcfdist's advantage*,
   which is the SOW's actual thesis: strictly better, not merely different.
3. **How it compares** — the D7 accuracy, runtime, and RAM figures against v2, `vcfeval`, and `hap.py`.
4. **Concordance evidence** — near-total agreement with `hap.py` on simple variants, with residual
   disagreements characterized. This is the "safe to adopt" section, and it is the most persuasive part of
   the post for anyone with an existing pipeline.
5. **How to get it** — bioconda, Docker, source; a link to the demo.

Three honesty constraints, all established upstream and all load-bearing for credibility:

1. **Every count changed, and the post says so before it shows a figure.** D4 §4 adopts per-site
   counting, so v3's numbers are not comparable to v2's or to the prior publications' without a
   convention-matched run. Leading with improved figures while omitting a denominator change would be
   the single most damaging thing this post could do.
2. **The GA4GH export reflects vcfdist's *thresholded* decisions** — `BD` is a hard label, so fractional
   `BC` credit cannot survive it. It is a compatibility export, not the vcfdist result (D4 §1.2).
3. **Any P/R shift from D6's clustering restructure is reported, not smoothed over**, and is attributed
   separately from (1), which is what C2's separate baselines make possible.

**GitHub release notes.** Concise: headline changes, the breaking-change list from §5, install instructions,
and a link to the blog post.

**Social — LinkedIn and BlueSky.** One post each, pointing at the blog post, leading with a single figure.
Publish after the blog post is live and the bioconda package is available, so every link in the post resolves
to something installable.

---

## 9. References

- SOW: `D0_vcfdist-v3-SOW.md` #8
- Design doc: [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md) §5.8
- Gating deliverable designs: [`D2`](./D2_vcfdist-v3-unit-tests.md), [`D6`](./D6_vcfdist-v3-unphased-eval.md),
  [`D5`](./D5_multi-bed-stratification.md), [`D3`](./D3_retain-info-format-fields.md),
  [`D7`](./D7_vcfdist-v3-benchmarking-design.md)
- Current release procedure: [`docs/v2.5.3/13-Release-Instructions.md`](../v2.5.3/13-Release-Instructions.md)
- Krusche et al. 2019, *Best practices for benchmarking germline small-variant calls in human genomes*,
  Nat. Biotechnol. — GA4GH benchmarking standard
- [Keep a Changelog](https://keepachangelog.com) — CHANGELOG format
