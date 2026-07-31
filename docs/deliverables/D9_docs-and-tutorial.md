# vcfdist v3.0.0 — Documentation & sandbox.bio Tutorial (Design)

|  |  |
| :-- | :-- |
| **Version:** | 0 (Draft) |
| **Authored By:** | Tim Dunn |
| **Status:** | Draft — plan of record for SOW Deliverable #9 |
| **Companion docs:** | `D0_vcfdist-v3-SOW.md`, [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md) (§5.9), [`D10_vcfdist-v3-multiqc-module.md`](./D10_vcfdist-v3-multiqc-module.md) |
| **Tracking issue:** | [#52](https://github.com/TimD1/vcfdist/issues/52) |

---

## 1. Purpose & Scope

Plan of record for **Deliverable #9** of the v3.0.0 SOW, which has two independent
workstreams:

- **A. v3 documentation refresh** — bring the README and the GitHub wiki in line with the v3 CLI and
  outputs, removing the v2-only features that no longer exist.
- **B. sandbox.bio tutorial** — an in-browser, no-setup tutorial written by the sandbox.bio maintainers
  from our outline, running on a small demo dataset we produce.

**In scope:** the wiki/README delta below, the tutorial outline, the demo dataset, and
`demo/pr_plot.py` (§2.7).
**Out of scope:** everything on the sandbox.bio side — the tutorial prose *and* the biowasm/WebAssembly
build (§3.2), both of which Robert owns — plus the release notes / CHANGELOG (D8) and the MultiQC report
docs (D10).

---

## 2. Workstream A — v3 documentation refresh

### 2.1 Mechanism

The **live GitHub wiki is the source of truth**; `docs/v<version>/` is a snapshot cached from it at
release time (`git clone …/vcfdist.wiki.git vX.X.X && rm -rf vX.X.X/.git`, step 4 of
[`13-Release-Instructions`](../v2.5.3/13-Release-Instructions.md)). So the work is edited against a wiki
clone and pushed there; `docs/v3.0.0/` falls out of the D8 release process, not out of this deliverable.

### 2.2 Page-by-page delta (from the `v2.5.3` snapshot)

| Page | Action |
| :-- | :-- |
| `Home` | Update TOC for deleted/added pages |
| `01-Overview` | Rewrite pipeline stages: no normalization, no edit-distance; graph/WFA alignment; per-variant phasing (D1 §3) |
| `02-Parameters-and-Usage` | Regenerate from v3 `--help` (§2.3) |
| `03-Variant-Filtering` | Review; largely unchanged |
| `04-VCF-Normalization` | **Delete** — code path removed in v3 |
| `05-Variant-Clustering` | Update: supercluster size bounding and splitting (`-sc`, renamed from `-s` — §2.3) |
| `06-Precision-and-Recall` | Update; add genotype / allele-count error summaries |
| `07-Phasing-Analysis` | Update: per-variant (not per-supercluster) phase decisions |
| `08-Alignment-Distance` | **Delete** — `edit` module removed in v3 |
| `09-Outputs` | Rewrite against the actual v3 file set (§2.4) |
| `10-Variant-Stratification` | Rewrite for native multi-BED, keeping the external GA4GH path as the interop option a user can drive themselves |
| `11-Implementation-Notes` | Review |
| `12-Installation-Help` | Review; bump version strings |
| `13-Release-Instructions` | Update to match the v3 release flow (coordinate with **D8**) |
| *new* `Unphased-Evaluation` | New page, per [`D6`](./D6_vcfdist-v3-unphased-eval.md) |
| *new* `GA4GH-Compatibility` | New page: retained INFO/FORMAT, the record shape, `BD`/`BK`/`N`, the match-tier ladder, and the three export caveats in §2.6 |

### 2.3 CLI delta to document

Comparing the flags documented in `v2.5.3/02-Parameters-and-Usage` against `Globals::print_usage`, and
including the flags the extension deliverables add:

- **Gone:** `--advanced`, `--cluster`, `--distance`, `--eval-mismatch-penalty`,
  `--eval-gap-open-penalty`, `--eval-gap-extend-penalty`, `--phasing-threshold`, `--realign-only`,
  `--realign-query`, `--realign-truth`.
- **Gone:** `-i`/`--max-iterations` ([`D6`](./D6_vcfdist-v3-unphased-eval.md) §4.2.3) — cluster iteration
  is fixed at 1 and is a correctness precondition once clustering runs on the merged variant list, not a
  tuning knob.
- **Renamed (breaking):** `-s` → `-sc`/`--max-supercluster-size`
  ([`D5`](./D5_multi-bed-stratification.md) §4.1). **No `-s` alias is retained**, deliberately: a silent
  alias would reinterpret an old `-s 15000` as a stratification file path. The docs must show the new form
  and say plainly that the old one now fails.
- **New:** `-md`/`--max-dist`; `-st`/`--stratification <manifest.tsv>` (D5 §4.2);
  `--stringency {lm,am,gm,pm}` and `--require-phase`/`--no-require-phase`
  ([`D4`](./D4_ga4gh-compatibility.md) §6.6).
- **Changed:** `-b`/`--bed` accepts gzipped and bgzipped BED files (D5 §5.1).
- **Still hidden but parsed:** `-x`, `-o`, `-e`, `-mr`/`--max-retries` are accepted by `parse_args` but
  commented out of `print_usage`. **Decision still needed** for these four: document them as advanced flags
  or remove the parsing — undocumented-but-live is the one state not to ship. `-mr`'s fate is coupled to
  PR #92 (D8 §2.2); if single-pass merges, the flag may go the way of `-i`.
- **Every default value** on this page is the one [`D8`](./D8_vcfdist-v3-release.md) §4.3 settles.

### 2.4 Output-file delta to document

| Change | Detail |
| :-- | :-- |
| Renamed | `parameters.txt` → `parameters.tsv` |
| Added | `runtime.tsv`, `genotype-errors.tsv` |
| Removed | `distance-summary.tsv`, `distance.tsv`, `edits.tsv`, `superclusters.tsv`, `orig-query.vcf`, `orig-truth.vcf` |
| Removed (decided) | `query.vcf` / `truth.vcf` — [#96](https://github.com/TimD1/vcfdist/issues/96), closed before rc1 per D8 G7, so `summary.vcf` is the only VCF vcfdist emits |
| Split | `phasing-summary.tsv` → `phasing-variants-summary.tsv` (per-variant, stratified) + `phasing-blocks-summary.tsv` (genome-level contiguity), D5 §7.1 |
| New columns | `STRATUM` leads every aggregate TSV, `*` first (D5 §7.1); `STRATA` set column on `query.tsv`/`truth.tsv` (D5 §7.2); `FP_GT`, `FP_AL`, `QUERY_UNK`, `FRAC_NA`, TiTv and het/hom ratios in the P/R summaries (D4 §7.1); `ASSESSED_PAIRS`, `PHASED_HET_VARIANTS`, `TOTAL_HET_VARIANTS`, `PHASED_HET_FRACTION`, `PHASED_REGION_FRACTION` in the phasing files (D6 §5.3) |
| Changed semantics | **Counting is per site, not per haplotype copy** (D4 §4) — the one change most likely to be read as a bug by a v2 user comparing numbers, so it needs its own callout on `06-Precision-and-Recall`, not a table row |
| Changed format | Switch/flip rates are bare-float **fractions** with a terminating newline, replacing the `%.6f%%` percent-string of `phasing-summary.tsv` (D5 §7.4). Same column names, values 100× apart — call this out explicitly, since a reader comparing a v2 run will otherwise read it as a hundredfold improvement |
| Changed semantics | Switch/flip rates use `ASSESSED_PAIRS`, not every merged query variant (D6 §5.1) — document both denominators explicitly, including that NG50 normalizes by contig length while `PHASED_REGION_FRACTION` uses evaluated BED bases (D6 §5.2) |
| Changed format | `summary.vcf`: one record per variant (het-alt still split), per-variant ploidy, no `ploidy=` on `##contig`, `Number=.` per-haplotype fields, `FORMAT/SC` → `INFO/BS`, `FORMAT/BS` → `FORMAT/PBS`, `BD=N` retention with `VCFDIST_*` FILTER tags, `/` separators for unphased genotypes |
| Fixed | `query.tsv`/`truth.tsv` `LOCATION` no longer emits a trailing space, [#95](https://github.com/TimD1/vcfdist/issues/95) |

### 2.5 README changes

Refresh the version string (currently `v2.5.0`; source is `3.0.0-b0`), the demo command and its expected
`PRECISION-RECALL SUMMARY` block (regenerate — the numbers will have moved), the feature bullets
(stratification, unphased query input), and add a link to the sandbox.bio tutorial once published.

### 2.6 GA4GH export caveats to document

The new `GA4GH-Compatibility` page has to state three things a user cannot infer from the output, all
established in [`D4`](./D4_ga4gh-compatibility.md):

1. **The native TSVs are vcfdist's numbers** (D4 §1.2). `summary.vcf` is conformant so that an existing
   GA4GH pipeline *can* consume it, but counts derived that way are a compatibility export rather than the
   vcfdist result, because `BD` is a hard label and fractional `BC` credit cannot survive it.
2. **`BK=pm` must be read as a match** by any consumer applying the published Method #3 `BK`→`BD` table
   (D4 §6.4.1). `pm` is a fourth-method extension and appears on the majority of TP records at the default
   stringency.
3. **A preserved non-`PASS` input `FILTER` is demoted by a quantifier** (D3 §4.6, D4 §3.3): filtered TPs
   become FNs and filtered FPs become Ns. vcfdist preserves the caller's `FILTER` deliberately, so a user
   who accepts a non-`PASS` filter via `--filter` and then quantifies will see those calls demoted. Also
   worth showing here: the `bcftools view -i 'FORMAT/BD=="FP"'` one-liners that replace vcfeval's separated
   `tp`/`fp`/`fn` VCFs (D4 §7.2).

### 2.7 `demo/pr_plot.py`

The fix is small — filter `precision-recall.tsv` rows to `STRATUM == "*"`, treating an absent column as
`*`. Without it the script silently aggregates all 182 strata into one curve, the same failure mode
[`D5`](./D5_multi-bed-stratification.md) §10 documents for the MultiQC module. It belongs here because the
demo output it plots is regenerated in this deliverable anyway (§2.5), as is `demo/results/`.

---

## 3. Workstream B — sandbox.bio tutorial

### 3.1 Division of labor

Per the author-provides-outline model the SOW describes (collaboration with sandbox.bio creator Robert
Aboukhalil, who proposed the tutorial): **we** finalize the outline and produce the demo dataset;
**upstream** writes the step Markdown; **we** review drafts. Tutorials live at
`src/content/<tutorial>/` with `steps/*.md`, `data/`, and a `config.js`, and commands are surfaced via
`<Execute command={...} />` components.

### 3.2 WebAssembly build — **out of scope**

A WASM build **is** required: sandbox.bio executes tools client-side via
[biowasm](https://github.com/biowasm/biowasm), with no server-side terminal option, and vcfdist is not
currently a biowasm tool. **Robert owns this**, along with the rest of the sandbox.bio side — biowasm is his
repository — so it is **not vcfdist-side scope** and is not costed here.

Recorded only so the hand-off is informed — the friction points a porter will hit, and the one that
could land on our side:

- **`std::thread` is constructed unconditionally**, even at `-t 1`: `main.cpp:70,87` and
  `dist.cpp:777,803`, with `main.cpp:69` looping over `HAPS × contigs` regardless of `max_threads`.
  Emscripten needs `-pthread` for that (SharedArrayBuffer + COOP/COEP headers, and biowasm's htslib is
  built `USE_PTHREADS=0`, so it would need a pthread-enabled rebuild). The cheap alternative is an inline
  path when `max_threads == 1`. Biowasm carries per-tool `patches/<tag>.patch` files for exactly this, so
  it can stay entirely upstream — but **if Robert would rather not carry the patch, a `-t 1` serial path
  upstreamed into vcfdist is the one plausible ask on us**, and it is worth having in vcfdist regardless
  for deterministic debugging.
- `Makefile` hardcodes `-g -pg -O1`; gprof is unsupported under Emscripten, so `CXXFLAGS` must be
  overridden wholesale (as `tools/ViralConsensus/compile.sh` does).
- Defaults `max_threads = 64` and `max_ram = 64` GB exceed wasm32's 4 GB address space — flags only
  (`-t 1 -r 1`), no code change.

The precedent is good: **ViralConsensus** (C++ + htslib) is already a biowasm tool with a ~10-line
`compile.sh` calling `emmake make CXX=em++` against biowasm's prebuilt htslib, and `htslib`/`bcftools`/
`samtools` are all existing modules. vcfdist's htslib surface is narrow (`vcf.h`, `kseq.h`) across ~7.9k
lines and 8 translation units.

### 3.3 Demo dataset

The existing `demo/` inputs are far too large for the browser (6.8 MB FASTA, 530 KB query VCF, against a
<500 KB-per-file guideline — ideally <100 KB). Build a new self-contained set, checked in under
`demo/tutorial/` and mirrored into the tutorial's `data/`:

- **`ref.fa`** — a single short synthetic contig (target ~20–50 kb), not a GRCh38 slice, so it stays small
  and license-free.
- **`truth.vcf`, `query.vcf`** — phased, single-sample, a few dozen variants, hand-built so each result
  bucket is reachable in a printed table: exact TPs, a genotype error, a partial-credit complex variant
  where representation differs but sequence matches (the case `hap.py`/`vcfeval` miscount — this is the
  point of the tutorial), an FP, an FN, and a phase flip plus a switch.
- **`regions.bed`** — one or two intervals, including a variant deliberately left outside to demonstrate
  `-b`.

Reproducibility: generate via a checked-in script with a fixed seed rather than by hand-editing, so the
dataset can be regenerated when output formats shift.

**This dataset has three consumers, not one.** That raises the bar on the generator script slightly and is
worth knowing before it is written:

1. **The sandbox.bio tutorial** — the driver above, and the source of the size guideline.
2. **The `snakemake-wrappers` test fixtures** ([`D11`](./D11_pipeline-ecosystem-integration.md) §7.3), which
   need exactly this: a tiny, phased, license-free set where every result bucket is non-empty, checked in
   next to the wrapper.
3. **`nf-core/test-datasets`** (D11 §5.6), for the nf-core module's `nf-test` cases. No phased truth/query
   VCF pair exists in that repository today, and an unphased pair would give a green test asserting
   near-empty results — worse than no test.

Two implications for the generator: the truth VCF must be *phased* (vcfdist still requires phased truth,
[`D6`](./D6_vcfdist-v3-unphased-eval.md) §9.1), and its output must be stable enough to be
snapshot-compared by `nf-test` and `compare_results_with_expected`, which makes a fixed seed a requirement
rather than good practice.

### 3.4 Outline shape

Install/`--help` → the input triple (query, truth, reference) and why phasing is required → run on the
demo data, read `precision-recall-summary.tsv` → the complex-variant step (show the same variant written
two ways; contrast with exact-match tools) → partial credit and `-ct` → phasing errors from
`switchflips.tsv` → restrict with `-b` → where to go next (wiki, MultiQC). Quiz/exercise points: predict
the TP/FP/FN counts before running; change `-ct` and explain the shift.

---

## 4. Verification

- Every command in the README and wiki is executed against the v3 binary, and every quoted output block
  is regenerated from a real run rather than edited by hand.
- The v3 flag list in `02-Parameters-and-Usage` is diffed against `vcfdist --help` output.
- The v3 output-file list in `09-Outputs` is diffed against an `ls` of a real run's `-p` directory.
- The demo dataset is run end-to-end and each intended result bucket is confirmed non-empty in
  `query.tsv`/`truth.tsv`/`switchflips.tsv`.
- Total dataset size is confirmed under the sandbox.bio per-file guideline.

---

## 5. Notes

- **Documenting follows the feature work.** `06-Precision-and-Recall` and the new `GA4GH-Compatibility`
  page describe D3/D4's output; `10-Variant-Stratification` describes D5's; `07-Phasing-Analysis` and the
  new `Unphased-Evaluation` page describe D6's. The deletions and the settled deltas (§2.3, §2.4) can be
  written at any point.
- **`02-Parameters-and-Usage` is regenerated after** [`D8`](./D8_vcfdist-v3-release.md) §4.3 settles the
  defaults, since `print_usage` prints every one of them on that page.
- **The WASM port is Robert's** (§3.2). The outline and dataset proceed in parallel; publication waits on
  the port. The only piece that could come back to us is a `-t 1` serial path, if he would rather not carry
  a biowasm patch.
- **[`D8`](./D8_vcfdist-v3-release.md) snapshots the finished wiki** into `docs/v3.0.0/` at release.

---

## 6. References

- SOW: `D0_vcfdist-v3-SOW.md` (#9); design doc: [`D1`](./D1_vcfdist-v3-design-doc.md) §5.9
- Tracking issue [#52](https://github.com/TimD1/vcfdist/issues/52)
- [sandbox.bio](https://sandbox.bio) — contributing guide, `sandbox-bio/sandbox.bio` README
- [biowasm](https://github.com/biowasm/biowasm) — WebAssembly modules for genomics
