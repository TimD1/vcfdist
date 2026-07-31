# D11: Pipeline Ecosystem Integration — nf-core & snakemake-wrappers

- **Status:** Design — complete, pending review. Decisions recorded in §11.
- **Issue:** [#54](https://github.com/TimD1/vcfdist/issues/54) (D11: Pipeline Ecosystem Integration)
- **Branch:** none yet
- **Area:** three upstream repositories — `nf-core/modules`, `nf-core/variantbenchmarking`,
  `snakemake/snakemake-wrappers`. No vcfdist source changes are in scope.
- **Builds on:** [`D8`](./D8_vcfdist-v3-release.md) (stable tag, bioconda, Docker),
  [`D5`](./D5_multi-bed-stratification.md) (`-sc`/`-st`, widened TSVs),
  [`D3`](./D3_retain-info-format-fields.md) / [`D4`](./D4_ga4gh-compatibility.md) (`summary.vcf` schema),
  [`D9`](./D9_docs-and-tutorial.md) §3.3 (demo dataset, reused as upstream test data)

---

## 1. Background

The SOW's thesis (#11) is ecosystem reach: `hap.py` and `vcfeval` are embedded in the pipelines
people actually run, and vcfdist is not. The current upstream state, verified rather than assumed:

| Fact | Evidence |
| ---- | -------- |
| No vcfdist module exists in `nf-core/modules` | `gh api "search/code?q=vcfdist+repo:nf-core/modules"` → `total_count: 0` |
| No vcfdist wrapper exists in `snakemake-wrappers` | `bio/` contains `hap.py`, `whatshap`, `vcftools`, `vcf2maf`; no `vcfdist` |
| `nf-core/variantbenchmarking` wraps the incumbents | modules `happy`, `rtgtools`, `truvari`, `svanalyzer`, `wittyer`; `params.method = 'truvari,svanalyzer,happy,rtgtools,wittyer,sompy,bndeval'` |
| bioconda ships vcfdist, one major behind | `recipes/vcfdist/meta.yaml` → version `2.6.4`, `run_exports: pin_subpackage(max_pin="x")`, recipe test is `vcfdist --version` only |

The closest precedents to follow are `modules/nf-core/happy/happy` (nf-core) and `bio/hap.py/hap.py`
(snakemake-wrappers). Both were read in full; §5 and §7 follow their current shape, not the shape
described in older third-party tutorials.

Two structural facts constrain everything below:

1. **Both ecosystems consume a *released, packaged* artifact.** nf-core modules pin a
   biocontainers image plus a bioconda `environment.yml`; snakemake wrappers pin a conda
   environment and a `snakedeploy`-generated lockfile. Neither can reference a locally-built binary.
2. **vcfdist's CLI and output filenames are the wrapper interface.** Every flag rename and every
   renamed output file is a breaking change to an artifact under someone else's review control, which
   is what D8 §3's output-format freeze protects.

---

## 2. Goals / Non-Goals

**Goals**

- An `nf-core/modules` module (`VCFDIST`) that runs vcfdist v3 with named outputs for every file it
  writes, passing `nf-core modules lint` and `nf-test` on all three container profiles.
- vcfdist offered as a first-class `method` in `nf-core/variantbenchmarking`, alongside
  `happy`/`rtgtools`, with its **native** metrics reaching the pipeline's report layer.
- A `snakemake-wrappers` entry (`bio/vcfdist`) meeting the repository's current contribution
  requirements, including a lockfile and example rules that double as rendered documentation.
- A recorded contract (§4) stating exactly which parts of the v3 CLI and output surface the two
  wrappers depend on, so a later change to any of them is visibly a breaking change.

**Non-Goals**

- **Any vcfdist source change.** §8 records observations that surfaced while writing this design;
  they are recorded here only, not filed, and none gates this deliverable.
- **Any GA4GH-quantifier route as the reporting path.** [`D4`](./D4_ga4gh-compatibility.md) makes
  `summary.vcf` conformant so an existing GA4GH pipeline *can* consume it, but T2 wires the native
  TSVs: a GA4GH export re-denominates the counts, drops fractional `BC` credit (`BD` is a hard label),
  and loses phasing errors entirely, since the intermediate VCF has no representation for them at all
  — reporting vcfdist as a slightly-worse `hap.py`. The export stays available to users; it is not a
  path this deliverable builds on.
- **MultiQC.** Owned by [`D10`](./D10_vcfdist-v3-multiqc-module.md).
- **The structural-variant benchmarking path** in variantbenchmarking (§6.5).
- **Somatic** (`sompy`-equivalent) integration — vcfdist is germline (D0 #14).

---

## 3. Upstream preconditions

Authoring starts once v3.0.0 is tagged and bioconda/biocontainers ship it. Submitting against the 2.6.4
package bioconda carries today would document a CLI that changes in three ways — `-s` renamed to `-sc`
with no alias retained (D5 §4.1), `-st` added (D5 §4.2), and output files renamed, split, and added
(D3/D4, D5 §7.1) — so the first review round would ship a `meta.yml` documenting flags that error on the
next release.

One consequence worth stating: D0's upstream-review latency — "weeks, controlled by those communities,
not by Fulcrum" — falls *after* the launch announcement, so the blog post (D8 §8) claims bioconda and
Docker availability, not nf-core.

Each precondition is observable:

| # | Precondition | Check |
| - | ------------ | ----- |
| P1 | v3.0.0 stable tag published | D8 §2.1 checklist complete |
| P2 | CLI and output surface frozen | D8 §3 freeze in effect since rc1 |
| P3 | bioconda recipe updated to 3.0.0 and merged | BiocondaBot PR verified per D8 §7 — the recipe's only test is `vcfdist --version`, so it passes even if flags changed; the flag review is manual |
| P4 | biocontainers image + exact build string available | `quay.io/biocontainers/vcfdist:3.0.0--<build>` resolvable; needed verbatim in `main.nf` |
| P5 | Phased truth/query test data available upstream | §5.6 — the only precondition that is itself an upstream PR |
| P6 | `bio.tools` entry for vcfdist, or an empty `identifier` accepted by lint | §5.2 |

---

## 4. The vcfdist contract both wrappers depend on

Written against the frozen v3 surface (current code plus the D3–D6 changes their design docs
commit to). Anything in this section changing after publication is a breaking change to two
upstream repositories.

### 4.1 Invocation shape

```
vcfdist <query.vcf> <truth.vcf> <ref.fasta> [options]
```

The three positional arguments must come **first**: `parse_args` inspects `argv[1..3]` and, if any
begins with `-`, warns `"Optional arguments should be provided AFTER mandatory arguments"` and
prints usage instead of running (`globals.cpp`). Both wrappers therefore build the command
positionals-first, and neither may interpolate a user-supplied `extra`/`args` string ahead of them.

### 4.2 Input formats

| Input | Accepted | Notes |
| ----- | -------- | ----- |
| query / truth VCF | plain, gzip, bgzip | opened with `bcf_open`; read sequentially |
| VCF index (`.tbi`/`.csi`) | **not used** | no random access; must not be a module input |
| reference FASTA | **uncompressed only** | `fasta.h` is `KSEQ_INIT(int, read)` over a raw `fopen`/`fileno`, so a `.fa.gz` is parsed as binary rather than rejected (§8.1) |
| reference index (`.fai`) | **not used** | vcfdist reads and uppercases the whole FASTA itself |
| `-b` BED | plain; `.bed.gz` after D5 §5.1 | |
| `-st` manifest | TSV; stratum BEDs may be `.bed.gz` | relative paths resolve against **the manifest's own directory** (D5 §4.2) — see §5.4 |

Truth input must be **locally phased**. D6 relaxes the phasing requirement for the query callset only;
unphased truth heterozygotes are still discarded (D6 §3, §9.1). This matters most of anything in this
section, because the failure mode is a successful run with a collapsed truth denominator, not an error
(§6.4, §8.3).

### 4.3 Output files at the freeze

`-p/--prefix` is a **path prefix**, not a directory: files are written as `<prefix><name>`. With
`-p out/run.` the run writes `out/run.precision-recall-summary.tsv`.

| File | Post-D5/D3/D4 status |
| ---- | -------------------- |
| `precision-recall-summary.tsv` | gains leading `STRATUM` column (D5 §7.1) |
| `precision-recall.tsv` | gains leading `STRATUM` column |
| `genotype-errors.tsv` | gains leading `STRATUM` column |
| `phasing-variants-summary.tsv` | **new** — per-variant half of `phasing-summary.tsv`, with `STRATUM` |
| `phasing-blocks-summary.tsv` | **new** — genome-level contiguity half, no `STRATUM` |
| `switchflips.tsv`, `phase-blocks.tsv`, `query.tsv`, `truth.tsv` | retained; `query.tsv`/`truth.tsv` gain a `STRATA` column |
| `summary.vcf` | rewritten by D3/D4: one record per variant (het-alt still split), per-variant ploidy, `Number=.` per-haplotype fields, INFO/FORMAT pass-through, `BD=N` retention, `INFO/BS`, `FORMAT/BS`→`PBS` |
| `parameters.tsv`, `runtime.tsv` | retained |
| phasing rate columns | bare-float **fractions** with a terminating newline, replacing `%.6f%%` (D5 §7.4) — see §8.2 |
| `query.vcf`, `truth.vcf` | **removed** by [#96](https://github.com/TimD1/vcfdist/issues/96) |
| `distance*.tsv`, `edits.tsv`, `orig-*.vcf` | already removed in v3 |

Exact headers, needed for §6.3's mapping:

```
precision-recall-summary.tsv:  STRATUM VAR_TYPE THRESHOLD MIN_QUAL TRUTH_TP QUERY_TP TRUTH_FN
                               QUERY_FP PREC RECALL F1_SCORE F1_QSCORE
                               FP_GT FP_AL QUERY_UNK FRAC_NA                        (D4 §7.1)
                               TRUTH_TITV QUERY_TITV TRUTH_HET_HOM QUERY_HET_HOM    (D4 §7.1)
genotype-errors.tsv:           STRATUM VAR_TYPE ALLELE_COUNT_0_TO_1 … ALLELE_COUNT_2_TO_2
```

`THRESHOLD` takes `NONE` (all calls) or `BEST` (the max-F1 quality threshold, per stratum after D5).

**`-n/--no-output-files` is unsupported in both wrappers.** It clears the `g.write` flag the file
writers check (e.g. `main.cpp`, `print.cpp`), so passing it through `ext.args`/`extra` makes
every declared output missing. Documented as unsupported rather than defended against in code.

### 4.4 Streams, exit status, version

- `INFO`, `WARN`, and `ERROR` all write to **stderr** (`defs.h`); `ERROR` then calls `std::exit(1)`,
  so a failed run is a nonzero exit and both ecosystems' default failure handling is correct.
- **stdout** carries the human-readable summary table and nothing a wrapper should parse.
- `vcfdist --version` prints `vcfdist v<VERSION>` to stdout. Use the **long form**: `-v` is
  `--verbosity` in the normal argument path and only reaches `--version` in the `argc < 4` branch
  (`globals.cpp`). This overloading is accepted as-is; the wrappers simply use `--version`.
- Output is deterministic across `-t` values, so byte-comparison of the metric TSVs is a valid test
  assertion (§9).

---

## 5. T1 — `nf-core/modules` module

### 5.1 Layout

Generated with `nf-core modules create vcfdist`, then edited:

```
modules/nf-core/vcfdist/
├── environment.yml
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    └── main.nf.test.snap
```

Process name `VCFDIST`, no subcommand subdirectory — vcfdist is a single command, unlike
`happy/{happy,prepy,sompy,ftxpy}`.

`environment.yml`:

```yaml
channels: [conda-forge, bioconda]
dependencies:
  - bioconda::vcfdist=3.0.0
```

Container pinned to the biocontainers build resolved at P4, with the `depot.galaxyproject.org`
singularity URL, exactly as `happy/happy/main.nf` does.

### 5.2 Interface

```groovy
input:
tuple val(meta),  path(query_vcf), path(truth_vcf), path(regions_bed)
tuple val(meta2), path(fasta)                  // uncompressed; §4.2
tuple val(meta3), path(stratification_tsv)
tuple val(meta4), path(stratification_beds)
```

Four channels, deliberately mirroring `HAPPY_HAPPY`'s shape so the two are interchangeable inside a
benchmarking subworkflow. No `fasta_fai` and no VCF index channel — §4.2 shows neither is read, and
declaring unused inputs invites callers to stage files for nothing.

Outputs: one named `emit` per file in §4.3, each `tuple val(meta), path('*<suffix>')`, plus the
versions topic:

```groovy
tuple val("${task.process}"), val('vcfdist'),
      eval("vcfdist --version | sed 's/vcfdist v//'"), topic: versions, emit: versions_vcfdist
```

`eval()` is used rather than the hardcoded `val('0.3.15')` string `happy/happy` carries, because
vcfdist does report its own version (§4.4) — the reason for hap.py's hardcoding does not apply.

Script body:

```groovy
def args   = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
def regions = regions_bed        ? "--bed ${regions_bed}"            : ''
def strat   = stratification_tsv ? "--stratification ${stratification_tsv}" : ''
"""
vcfdist \\
    ${query_vcf} ${truth_vcf} ${fasta} \\
    ${regions} ${strat} ${args} \\
    --max-threads ${task.cpus} \\
    --max-ram ${task.memory.toGiga()} \\
    --prefix ${prefix}.
"""
```

Three constraints this encodes, each documented in `meta.yml`:

- **`task.ext.prefix` must not contain `/`.** `-p` is a path prefix (§4.3) and a separator would
  write outputs outside the task directory, where the output globs do not reach.
- **`--max-ram` is derived from `task.memory`, not left to the default 64 GB**, which would let
  vcfdist plan around memory the executor has not granted.
- `stratification_beds` is staged but never named on the command line; it exists so the manifest's
  referenced files are present in the work directory (§5.4), which is also how `HAPPY_HAPPY`
  handles `--stratification`.

`label 'process_high'`: vcfdist defaults to `-t 64 -r 64`, and
[`D7`](./D7_vcfdist-v3-benchmarking-design.md) §7.1 benchmarks it on a 64-vCPU / 256 GiB AWS
`m6i.16xlarge`. Callers can lower
via `withName` config.

A `stub` block touches every declared output.

`meta.yml` needs `identifier:` — either `biotools:vcfdist` once registered, or empty (P6). Worth
registering regardless; it is a free discoverability win.

### 5.3 Snapshot strategy

Two outputs are not byte-stable and must not be content-snapshotted:

- `parameters.tsv` records the full command line, which contains the absolute Nextflow work
  directory (a per-run hash).
- `runtime.tsv` records wall-clock timings.

Snapshot their **names** only (`file(process.out.parameters[0][1]).name`), the pattern
`happy/happy`'s test already uses for `metrics_json` and `runinfo`. Everything else is
content-snapshotted; §4.4 makes that safe.

### 5.4 The stratification-manifest staging problem

D5 §4.2 resolves relative manifest paths against **the manifest's own directory**, so that GIAB's
shipped `v3.1-GRCh38-all-stratifications.tsv` works unmodified. GIAB's entries carry
subdirectories:

```
refseq_cds	FunctionalRegions/GRCh38_refseq_cds.bed.gz
```

Nextflow stages inputs **flat** into the work directory. The manifest lands beside the BEDs, so
manifest-directory resolution and working-directory resolution coincide — but only for a **flat**
manifest whose second column is a bare filename. A GIAB manifest staged this way resolves to
`FunctionalRegions/…`, which does not exist, and D5 §4.2 makes that an `ERROR` naming both the
manifest-relative and resolved absolute path — a loud failure, which is the right outcome.

**Decision: require a flat manifest in the module, and document it.** The alternatives — rewriting
the manifest inside the module, or staging with `stageAs` to recreate the directory tree — either
add a hidden transformation to a module whose job is to run one binary, or push a fragile staging
expression onto every caller. Flattening a GIAB manifest is a one-line `awk` in the calling
pipeline, and §7.4 notes snakemake needs no flattening at all.

### 5.5 Tests

`tests/main.nf.test`, tagged `modules`, `modules_nfcore`, `vcfdist`:

1. **`homo_sapiens phased vcf bed`** — query, truth, `-b` BED, reference; empty stratification
   channels. Asserts `process.success`, snapshots the metric TSVs and `summary.vcf`, names only for
   `parameters.tsv`/`runtime.tsv`.
2. **`homo_sapiens phased vcf stratified`** — adds a two-stratum flat manifest, asserting the
   per-stratum rows appear and that the `*` row is unchanged from test 1. This is D5 §9.2's central
   invariant, re-asserted at the module boundary because it is exactly what a pipeline user would
   otherwise silently get wrong.
3. **`homo_sapiens phased vcf bed - stub`** — required by nf-core.

### 5.6 Test data — the one real prerequisite (P5)

vcfdist requires phased truth (§4.2). The `nf-core/test-datasets` `modules` branch has no phased
truth/query VCF pair: the GATK HaplotypeCaller call sets that `HAPPY_HAPPY`'s test uses are
unphased, and a search of the branch tree for phasing-related paths returns phased **BAMs**
(`nanopore/bam/test.sorted.phased.bam`, `pacbio/bam/…phased.bam`) and 1000G resource VCFs, not a
truth/query pair. Running vcfdist against unphased truth would produce a green test asserting
near-empty results — the worst possible test.

**Reuse [`D9`](./D9_docs-and-tutorial.md) §3.3's demo dataset.** It is already specified as a small
(<100 KB/file), self-contained, *phased*, synthetic set built by a checked-in fixed-seed script,
deliberately exercising exact TPs, a genotype error, a partial-credit complex variant, an FP, an FN,
and a flip plus a switch. That is a better module fixture than any real-data slice: tiny, license-
free, regenerable, and every result bucket is non-empty. It requires its own PR to
`nf-core/test-datasets` (a synthetic contig is unusual there and may draw discussion), which is why
it starts first (§10).

The same dataset serves T3 directly from the vcfdist repo, with no upstream dependency.

---

## 6. T2 — `nf-core/variantbenchmarking`

### 6.1 Where vcfdist plugs in

`subworkflows/local/small_germline_benchmark/main.nf` dispatches on `params.method`, one `if` block
per tool, mixing two channels: `summary_reports` and `tagged_variants`. vcfdist adds a third block
next to `rtgtools` and `happy`:

```groovy
if (params.method.contains('vcfdist')) {
    VCFDIST(input_ch.map { meta, tvcf, _tidx, truth, _bidx, regions, _targets ->
                              [ meta, tvcf, truth, regions ] },
            fasta, stratification_tsv, stratification_bed)

    summary_reports = summary_reports.mix(
        VCFDIST.out.pr_summary
               .map { _meta, file -> tuple([vartype: params.variant_type,
                                            benchmark_tool: "vcfdist"], file) }
               .groupTuple())
}
```

Also required: `vcfdist` added to `params.method` in `nextflow.config` and to the corresponding
`nextflow_schema.json` enum, and a `conf/tests/test_vcfdist.config` profile mirroring
`test_happy.config`.

### 6.2 What is *not* wired in v1

- **`tagged_variants`.** The pipeline expects per-decision VCFs (`FN`, `FP`, `TP_base`, `TP_comp`),
  which `rtgtools` emits natively as four files. vcfdist emits one two-sample `TRUTH`/`QUERY`
  `summary.vcf`; splitting it into four would mean a `bcftools view -i 'FORMAT/BD="FP"'` fan-out
  whose semantics differ per sample column and, after D4 §4's per-site counting, per site rather than
  per haplotype.
  The channel is a `mix`, so contributing nothing is legal. Deferred, documented.
- **Phasing metrics.** `phasing-variants-summary.tsv` / `phasing-blocks-summary.tsv` are vcfdist's
  differentiator and have no column in the pipeline's canonical summary schema (§6.3). They are
  published as files; surfacing them belongs to [`D10`](./D10_vcfdist-v3-multiqc-module.md), whose
  report has a section for exactly this. Proposing a schema extension in the same PR that
  introduces the tool would enlarge the review surface for the least-reviewable part.
- **Stratification.** The pipeline passes GIAB manifests to `HAPPY_HAPPY`; wiring the same channels
  through `VCFDIST` needs the §5.4 flat-manifest step. Deferred to a follow-up so the first PR is
  one tool, one path.

### 6.3 `bin/merge_reports.py` — a `vcfdist` reader

`MERGE_REPORTS` is tool-aware (`merge_reports.py $inputs -b $meta.benchmark_tool …`), and each
reader normalizes onto one canonical frame. From `get_happy_results`, the canonical columns are:

```
Tool File Caller Type Filter TP_base TP_comp FN TP_Total FP UNK FP_gt FP_al
Recall Precision Frac_NA F1 [TiTv / het_hom ratios]
```

`get_vcfdist_results()` reads `*precision-recall-summary.tsv` with `csv.DictReader`-equivalent
`pandas.read_csv(sep='\t')` and maps:

| Canonical | vcfdist source |
| --------- | -------------- |
| `Type` | `VAR_TYPE` (`SNP`, `INDEL`, `SV`, `ALL`) |
| `Filter` | `THRESHOLD` — `NONE` → `ALL`, `BEST` → `BEST` (a vcfdist-specific value, documented) |
| `TP_base` / `TP_comp` | `TRUTH_TP` / `QUERY_TP` |
| `FN` / `FP` | `TRUTH_FN` / `QUERY_FP` |
| `TP_Total` | `QUERY_TP + QUERY_FP` |
| `Recall` / `Precision` / `F1` | `RECALL` / `PREC` / `F1_SCORE` |
| `UNK` / `Frac_NA` | `QUERY_UNK` / `FRAC_NA`, emitted natively by D4 §7.1 — no derivation needed |
| `FP_gt` / `FP_al` | `FP_GT` / `FP_AL`, also native per D4 §7.1. Note these do **not** sum to `QUERY_FP`: an FP that matched nothing carries `BK=.` and belongs to neither sub-column, exactly as in hap.py's own output, so a reader must not treat the shortfall as a bug |

Three details that will otherwise silently produce wrong rows:

- **Filter to `STRATUM == "*"`.** After D5 every aggregate TSV carries a leading `STRATUM` column;
  a reader that ignores it aggregates 182 strata into one dataset and reports inflated counts with
  no error. D5 §10 flags the same trap for the D10 MultiQC module — it is the standard way to get this
  wrong.
- **Counting is per site.** D4 §4 makes a homozygous variant contribute 1 to TP/FP/FN rather than 2,
  matching what `happy` and `rtgtools` report, so vcfdist's rows are directly comparable to theirs in
  the merged table. This is what makes a shared summary schema meaningful at all; before D4 the same
  table would have put differently-denominated numbers side by side. Worth stating in the PR
  description, since a reviewer comparing vcfdist against `happy` on the same callset is the first
  person who would notice if it were not true.
- **`Caller` is positional.** `merge_reports.py` derives it as `basename.split('.')[2]`, so
  `ext.prefix` must be `vcfdist.<id>.<caller>`; vcfdist appends `precision-recall-summary.tsv`
  after the trailing dot, giving `vcfdist.test.gatk.precision-recall-summary.tsv` → token 2 =
  `gatk`. Getting this wrong mislabels every row rather than failing.

Plus `assets/datavzrd/vcfdist.datavzrd.template.yaml`, which the pipeline requires per tool
(`REPORT_BENCHMARK_STATISTICS` resolves `assets/datavzrd/${meta.id}.datavzrd.template.yaml` with
`checkIfExists: true`, so its absence is a hard failure).

### 6.4 Truth-phasing validation — the item most likely to bite

The pipeline preprocesses the truth VCF (`prepare_vcfs_truth`: normalization, deduplication,
multiallelic splitting) before benchmarking. If any of those steps drops `|` separators, vcfdist
discards the truth heterozygotes (§4.2) and reports a collapsed truth denominator from a run that
exits 0.

**Validation step before the PR:** run the pipeline's own truth-preparation path on a phased truth
VCF and confirm the output retains phased genotypes; if it does not, the vcfdist block must either
document a phased-truth requirement on `params` or skip preparation for its input. §8.3 records the
vcfdist-side guard that would make this self-diagnosing.

### 6.5 SVs out of scope

vcfdist evaluates SVs (`-sv`, `-l`) and `sv_benchmark` is where Truvari/SVanalyzer live, so a second
integration point exists. It is deferred: the practical SV ceiling is 1–10 kb (D0 #12/#13, D1 §4),
below what SV benchmarking users expect, and offering vcfdist there before #12 lands would invite a
poor first impression on the axis where it is currently weakest.

---

## 7. T3 — `snakemake-wrappers`

### 7.1 Layout

Per the repository's current `docs/contributing.rst` (not the older `bio/hap.py` layout):

```
bio/vcfdist/
├── environment.yaml
├── environment.linux-64.pin.txt      # snakedeploy pin-conda-envs environment.yaml
├── meta.yaml
├── wrapper.py
└── test/
    ├── Snakefile
    └── <fixtures: ref.fa, truth.vcf, query.vcf, regions.bed, strat.tsv, strata/>
```

Plus a `test_wrappers.py::test_vcfdist` entry **per example rule** in `test/Snakefile` — the
repository requires one-to-one coverage.

`environment.yaml` pins **major.minor only**, per the guide's semantic-versioning rule, with
`nodefaults` last:

```yaml
channels: [conda-forge, bioconda, nodefaults]
dependencies:
  - vcfdist =3.0
```

`meta.yaml` requires `name`, `description`, `url`, `authors`, `input`, `output`, and optionally
`params` and `notes`. `notes` is where the §4.2 constraints belong: uncompressed reference,
phased truth, `-n` unsupported.

### 7.2 `wrapper.py`

```python
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

bed = snakemake.input.get("bed", "")
if bed:
    bed = f"--bed {bed}"

strat = snakemake.input.get("strat", "")
if strat:
    strat = f"--stratification {strat}"

mem_gb = snakemake.resources.get("mem_mb", 0) / 1024
max_ram = f"--max-ram {mem_gb:.2f}" if mem_gb else ""

prefix = snakemake.params.get("prefix") or common_prefix(snakemake.output)

shell(
    "vcfdist"
    " {snakemake.input.query} {snakemake.input.truth} {snakemake.input.ref}"
    " {bed} {strat}"
    " --max-threads {snakemake.threads} {max_ram}"
    " --prefix {prefix}"
    " {extra}"
    " {log}"
)
```

(`shell()` formats against the calling scope, so every interpolated value is bound to a local
first — the `hap.py` wrapper's idiom — rather than computed inside the format string.)

Four decisions:

- **`stdout=True, stderr=True`.** vcfdist's diagnostics are on stderr and its summary table on
  stdout (§4.4); the guide requires both go to `log:`.
- **Prefix derived from the outputs, `params.prefix` as override.** The wrapper must "deal with
  arbitrary `input:`/`output:` paths" and "automatically infer command line arguments wherever
  possible". `common_prefix()` takes the declared outputs, strips the known vcfdist suffixes from
  §4.3, and raises if they do not share one prefix — which is the only real failure mode and a clear
  one. It stays local to the wrapper; if a second vcfdist wrapper is ever added it moves to
  `snakemake-wrapper-utils`, per the guide's rule on repeated functionality.
  `hap.py`'s wrapper requires `params.prefix` unconditionally; deriving it is friendlier
  and keeps the outputs the single source of truth. (This is the same suffix-stripping rule D10 §2
  uses to collapse one run's files onto one sample.)
- **`--max-ram` from `resources.mem_mb`** when set, else the vcfdist default — same reasoning as
  §5.2.
- **`extra` after the derived flags**, so a user can override anything except the positionals
  (§4.1).

Formatted with `black`.

### 7.3 `test/Snakefile`

Two rules, `snakefmt`-formatted, snake_case, doubling as the rendered documentation:

- `benchmark_variant_calls` — the minimal case: query, truth, reference, `-b` BED, the metric TSVs
  as outputs, `threads: 4`.
- `benchmark_variant_calls_stratified` — adds `strat` (manifest + stratum BEDs) and shows
  per-stratum output.

Fixtures are the D9 §3.3 demo dataset (§5.6), copied in. `compare_results_with_expected` pins
`precision-recall-summary.tsv` against a checked-in expected file — valid because §4.4 establishes
determinism — and deliberately excludes `parameters.tsv` (absolute paths) and `runtime.tsv`
(timings).

### 7.4 One asymmetry worth documenting

Snakemake passes **real filesystem paths**, so a GIAB `-st` manifest with subdirectory entries
resolves correctly with no preprocessing — the opposite of the flat-manifest requirement Nextflow's
staging imposes (§5.4). The `notes` field says so, since a user moving between the two ecosystems
will otherwise hit the difference blind.

---

## 8. Recorded vcfdist-side observations

Surfaced while writing this design. **Recorded here only — not filed, and none gates D11.** Each is
a place where a wrapper must work around behavior rather than rely on it.

### 8.1 Gzipped reference FASTA is mis-parsed, not rejected

`fasta.h`'s `KSEQ_INIT(int, read)` reads through raw POSIX `read` on a `fopen`ed handle, so a
bgzipped or gzipped reference yields no usable contigs rather than an error, and the run fails later
with a confusing contig-lookup message. D5 §5.1 sets the precedent for the fix — it replaces
`bed.cpp`'s `ifstream` with htslib's `hts_open`/`hts_getline` for exactly this reason — and the
cheap alternative is a magic-byte check with a clear `ERROR`. Both wrappers document uncompressed
FASTA instead (§4.2).

### 8.2 `phasing-summary.tsv` rate format

`SWITCH_ERROR_RATE` and `FLIP_ERROR_RATE` are written as `%.6f%%` (`phase.cpp:810`) — **strings with a
trailing percent sign** rather than floats — and the single data row carries no trailing newline.
[`D5`](./D5_multi-bed-stratification.md) §7.4's split writers emit **bare-float fractions** with a
terminating newline instead. The column names do not change, so the same column reads `0.020930%` before
and `0.00020930` after — a 100× difference that consumers discriminate on the legacy `%` sigil.
[`D10`](./D10_vcfdist-v3-multiqc-module.md) carries that logic; D11 is unaffected in v1, since §6.2 defers
phasing metrics.

### 8.3 Unphased truth degrades silently

Discarded truth heterozygotes are reported only as an aggregate `WARN` among nine other drop
counters. In an automated pipeline nobody reads that line, and the visible result is a successful
run with a small truth denominator. A guard that reports the *fraction* of truth heterozygotes
dropped for missing phase — `WARN` above some level, or `ERROR` at the extreme — would turn §6.4's
manual validation step into a self-diagnosing one. Analogous in spirit to D5 §6.4's decision to
promote the VCF sort order from an unchecked assumption to an enforced precondition.

---

## 9. Verification

| Target | Checks |
| ------ | ------ |
| T1 | `nf-core modules lint vcfdist`; `nf-test test modules/nf-core/vcfdist/… --profile docker`, `conda`, `singularity`; snapshots reviewed for the §5.3 exclusions; `prettier`/editorconfig clean |
| T2 | `nf-core pipelines lint`; `nf-test` on `small_germline_benchmark`; a full `-profile test_vcfdist,docker` run producing a datavzrd report with populated vcfdist rows; the §6.4 truth-phasing check performed and its result recorded |
| T3 | `pytest test_wrappers.py -k vcfdist`; `black wrapper.py`; `snakefmt test/Snakefile`; snakemake lint; docs render inspected (the `test/Snakefile` *is* the published example) |
| Cross | The same demo dataset produces matching `precision-recall-summary.tsv` values under T1, T3, and a bare CLI run — three invocation paths, one number |

---

## 10. Order of work

1. **Demo dataset** — finalize the D9 §3.3 set, the fixture for both T1 and T3, confirming it stays under
   the size guidance while keeping every result bucket non-empty.
2. **`nf-core/test-datasets` PR** (P5) — the longest-running item, since a synthetic contig is unusual
   there and may draw discussion.
3. **T3, `snakemake-wrappers`** — carries its own fixtures, so it is the cheapest place to discover CLI
   ergonomics problems, and likely the first to merge.
4. **T1, `nf-core/modules`** — uses the test data from step 2 and the biocontainers build string from P4.
5. **T2, `nf-core/variantbenchmarking`** — wraps the T1 module.
6. **Follow-ups, each its own PR:** stratification wiring in T2, `tagged_variants`, and the SV path once
   #12 lands.

Per D0, the upstream PRs are "submitted and shepherded," not merged on a date.

---

## 11. Decision record

| # | Decision | Rationale |
| - | -------- | --------- |
| 1 | Author against the tagged v3.0.0 release, not bioconda's 2.6.4 | 2.6.4 would document `-s`, omit `-st`, and name output files D3/D4/D5 rename — a guaranteed second review round on someone else's queue (§3) |
| 2 | Native TSVs as T2's reporting path; the GA4GH export documented as interop only | A GA4GH consumer counts hard integer `BD`, losing fractional `BC` credit, and VCF-I cannot carry phasing errors; routing through one would report vcfdist as a slightly-worse `hap.py` (§2) |
| 3 | Module interface mirrors `HAPPY_HAPPY`'s channel shape | Makes the two drop-in interchangeable inside a benchmarking subworkflow, which is the actual use case (§5.2) |
| 4 | No `fasta_fai` or VCF-index inputs | Neither is read (§4.2); declaring them would have callers stage files for nothing |
| 5 | Uncompressed reference documented rather than worked around | A wrapper cannot fix a mis-parse; §8.1 records the source-side option |
| 6 | Flat `-st` manifest required in T1; no in-module rewriting | Keeps the module a thin binary invocation; D5 §4.2 already makes the mismatch a loud ERROR, and flattening is a one-line `awk` upstream (§5.4) |
| 7 | `parameters.tsv`/`runtime.tsv` snapshotted by name only | Absolute work-dir paths and timings; same pattern `happy/happy` uses for `runinfo` (§5.3) |
| 8 | Reuse the D9 demo dataset as upstream test data | Tiny, phased, license-free, regenerable, and every result bucket non-empty — better than any real-data slice, and no phased pair exists upstream (§5.6) |
| 9 | Wrapper derives its prefix from declared outputs | The guide requires inferring arguments where possible; keeps `output:` the single source of truth, unlike `hap.py`'s mandatory `params.prefix` (§7.2) |
| 10 | `-n/--no-output-files` unsupported in both wrappers | It clears the flag the file writers check, so all declared outputs vanish (§4.3) |
| 11 | `tagged_variants`, phasing metrics, stratification, SVs deferred in T2 | One tool, one path in the first PR; phasing has no canonical column and belongs to D10 (§6.2) |
| 12 | T3 before T1 | Self-contained fixtures, and it surfaces CLI ergonomics problems before they are baked into two nf-core repos (§10) |
| 13 | `-v` overloading accepted; wrappers use `--version` | Long form is unambiguous in every argument path; no source change warranted (§4.4) |
| 14 | Observations recorded, not filed | The §8 items are wrapper-visible but none is D11 work; recording keeps them available to D8/D9/D10 |

---

## 12. References

- SOW: `D0_vcfdist-v3-SOW.md` #11; design doc:
  [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md) §5.11
- Sibling designs: [`D5`](./D5_multi-bed-stratification.md) (`-sc`/`-st`, `STRATUM` columns),
  [`D3`](./D3_retain-info-format-fields.md), [`D4`](./D4_ga4gh-compatibility.md) (GA4GH export),
  [`D8`](./D8_vcfdist-v3-release.md) (tag, freeze, bioconda),
  [`D9`](./D9_docs-and-tutorial.md) §3.3 (demo dataset),
  [`D10`](./D10_vcfdist-v3-multiqc-module.md) (MultiQC)
- [nf-core modules documentation](https://nf-co.re/docs/contributing/modules) and
  [`modules/nf-core/happy/happy`](https://github.com/nf-core/modules/tree/master/modules/nf-core/happy/happy)
- [`nf-core/variantbenchmarking`](https://github.com/nf-core/variantbenchmarking) —
  `subworkflows/local/small_germline_benchmark`, `bin/merge_reports.py`, `assets/datavzrd/`
- [snakemake-wrappers contributing guide](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html)
  and [`bio/hap.py/hap.py`](https://github.com/snakemake/snakemake-wrappers/tree/master/bio/hap.py/hap.py)
- [bioconda `recipes/vcfdist`](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/vcfdist)
- Krusche et al. 2019, *Best practices for benchmarking germline small-variant calls in human
  genomes*, Nat. Biotechnol. — GA4GH benchmarking standard
