# vcfdist v3.0.0 — Benchmarking vs. Prior Work (Design / Plan)

|  |  |
| :-- | :-- |
| **Version:** | 1 (Draft) |
| **Authored By:** | Tim Dunn |
| **Status:** | Draft — plan of record for SOW Deliverable #7 |
| **Companion docs:** | [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md) §5.7, `D0_vcfdist-v3-SOW.md` #7 |

---

## 1. Purpose & Scope

This document is the design for **Deliverable #7** of the v3.0.0 SOW:
a reproducible harness that benchmarks the current **vcfdist v3** (`dev` branch) against the most recent
**vcfdist v2** release and the two established germline benchmarking tools, **vcfeval** and **hap.py
(`xcmp` engine)**, regenerating the accuracy and runtime/RAM figures from the prior publications.

The harness produces four things:

1. **FNR / FDR bar charts** across **vcfdist v3 (dev) / vcfdist v2.6.4 / vcfeval / hap.py (xcmp)**.
   Truvari is not included.
2. **Precision-recall curves with all tools on the same plot**, plus **runtime and peak-RAM** figures.
3. **Per-variant concordance** (§6.5) — the SOW's "safe to displace" evidence: cross-tabulate
   vcfdist v3's TP/FP/FN decision against each other tool's decision per variant, so residual
   disagreements can be counted and characterized rather than inferred from aggregate rates.
4. **Phasing correctness** (§6.6) — switch/flip metrics validated against `whatshap compare` and
   checked for v2→v3 regression. It covers the half of vcfdist's headline claim that the P/R figures do
   not; rationale in §11.

The truth set, query call-sets, region filtering, and variant size classes match the prior publications,
for continuity with the published figures.

### Out of scope

- Truvari (removed from the comparison).
- New truth sets or query call-sets beyond those already used in v2.
- The MultiQC module / downstream toolkit (**SOW #10**) — this harness emits plain figures + TSVs only.
- **Per-stratum (multi-BED) reporting** — the harness runs a single evaluation BED; stratified reporting
  arrives with #5 and is a follow-on comparison, not part of this deliverable (§10).
- **Unphased-input evaluation** — a second run with query phase stripped becomes possible only once #6
  lands; scoped as the companion run described in §10, not as part of the current figures.

---

## 2. Tools Compared

| Tool | Version | Source | Notes |
| :-- | :-- | :-- | :-- |
| **vcfdist v3** | `dev` branch (this repo) | Built from local `src/` | The tool under test; not yet on bioconda. |
| **vcfdist v2** | **v2.6.4** (latest release) | bioconda (`vcfdist=2.6.4`) | Regression baseline: v3 must not regress vs. v2. |
| **vcfeval** | RTG Tools (pin in `config.yml`) | bioconda (`rtg-tools`) | `--all-records --ref-overlap --vcf-score-field=QUAL`. |
| **hap.py** | `pkrusche/hap.py` image tag pinned in `config.yml` | **Docker** (not bioconda) | `--engine=xcmp`. Currently disabled in the harness; runs natively on the `x86-64` benchmarking host (§7.1). |

**Run host: an on-demand AWS `x86-64` instance.** All published figures are generated on a
**purpose-provisioned on-demand AWS EC2 `x86-64` (linux-64) instance**, not on the Fulcrum M3 Max named in
the SOW's Materials section. The M3 Max remains the **development host** — unit tests, DAG dry-runs, and the
chr20 `smoke` run, all of which are minutes-scale and platform-agnostic. Rationale, in the order that
decided it:

1. **hap.py runtime becomes comparable.** hap.py is a legacy **Python-2**, `linux-64`-only tool. On Apple
   silicon it can only run as an emulated `linux/amd64` container, so its wall-clock would be an artifact of
   emulation rather than a measurement — making the runtime figure misleading exactly where the SOW wants a
   head-to-head. On native `x86-64` all four tools run on the same ISA, uncontended, on one machine.
2. **Reproducibility.** A pinned instance type + AMI + region is a specification a third party can re-create;
   "Tim's laptop" is not. It also removes the last platform-specific shim in the harness (see `time` below).
3. **Comparability with the published v2 figures.** The prior harness ran on a 64-thread Linux server with
   `-t 64 -r 64`; a 64-vCPU instance restores that operating point (§7.1).

**hap.py execution — Docker, now unemulated.** hap.py still runs from the pinned `pkrusche/hap.py` **Docker
image** rather than the bioconda package, so the exact hap.py build is captured by an image digest:

- `rule eval_happy` shells out to `docker run --rm` with the data directory and reference bind-mounted
  read-only and the results directory writable, passing container-side paths. The image tag/digest is pinned
  in `config.yml` alongside the other tool versions (§8).
- On the `x86-64` host the image runs **natively — no emulation penalty** — so hap.py's wall-clock is a real
  measurement and belongs on the runtime figure without a caveat. The containerization itself is worth a
  footnote, nothing more.
- **Peak-RAM still needs a container-aware mechanism.** `/usr/bin/time -v docker run ...` measures the Docker
  *client*, not the container, so its RSS is meaningless. Capture container RSS from cgroup v2 accounting
  inside the container (`/sys/fs/cgroup/memory.peak`, read by the container command after hap.py exits) or
  from `docker stats`; if neither lands cleanly, omit hap.py from the RAM figure and say so on it.
- *Alternative, if container RAM accounting proves fiddly:* the already-scaffolded `happy` bioconda pixi env
  now installs fine on the `linux-64` host, and `/usr/bin/time -v` would then measure hap.py directly. That
  trades an image digest for a conda solve as the version record — acceptable, but Docker is preferred for
  reproducibility.

Until that rule is written, the harness ships hap.py **wired but disabled**: `rule eval_happy`
(`evaluate.smk`) and `rule parse_happy` (`parse.smk`) are commented out, and `happy` is removed from
`ALL_TOOLS` (`plot.smk`) and the plotters' `TOOLS` lists — so current runs compare **three** tools
(vcfdist v3, vcfdist v2.6.4, vcfeval). Re-enabling it is tracked as §12 item 1.

**Credit-threshold pin (apples-to-apples vcfdist comparison).** vcfdist's `--credit-threshold` (`-ct`) —
the minimum partial credit at which a variant counts as TP — has a **different default in v2 and v3**
(v3 `src/globals.h:50` sets `0.98`; the v2.6.4 default is `0.70`). The harness sets it **explicitly and
identically for both** (`params.credit_threshold: 0.98`), so v2-vs-v3 differences reflect engine changes
rather than a default change. Any change to the v3 default must be mirrored here or the regression
comparison is invalid.

**Other sourcing notes:**

- **vcfdist v3** is built from the local checkout (`src/vcfdist`, path in `params.vcfdist_v3_bin`), so the
  harness depends on a prior build step, not a package.
- **vcfeval** needs a reference **SDF** (`rtg format`), built once in the prep stage (now pre-staged, §5).
- **GNU `time -v`** is required for the resource metrics; macOS `/usr/bin/time` is BSD and lacks `-v`, so
  the harness invokes the conda `time` package via `command time -v`.

---

## 3. Datasets

**Exact versions pinned** in `config.yml` for reproducibility.

| Role | Name | Version | Sample |
| :-- | :-- | :-- | :-- |
| **Truth** | T2T-Q100 | v0.9 | HG002 |
| Query | HPRC (hifiasm-dipcall) | v1 | HG002 |
| Query | PAV (Q100-PAV) | v4 | HG002 |
| Query | GIAB-TR (hifiasm-GIAB-TR) | v4.20 | HG002 |
| Reference | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` | — | GRCh38 |

The truth set's benchmark BED (`bench.bed`) is used as both the region filter and the evaluation region for
every tool. All evaluation is over HG002 (`--bSample`/`--cSample HG002` where applicable). A chr20-only
subset of that BED (`bench-chr20.bed`) drives the fast end-to-end `smoke` run (§9).

---

## 4. Directory Layout

New sibling to the existing analysis directory, structured on the Fulcrum
[tertiary-analysis-template](https://github.com/fulcrumgenomics/tertiary-analysis-template) (pixi + Snakemake):

```
analysis-v3/
├── data/                        # pre-prepared inputs, shared (gitignored): refs/ + <id>-<ver>/split/
└── vs_prior_work/
    ├── pixi.toml                # envs: default (analysis), bench (analysis+tools), happy (hap.py, py2)
    ├── pixi.lock                # pinned environment
    ├── README.md                # setup + run instructions, incl. how to re-enable hap.py
    ├── workflows/
    │   ├── config/
    │   │   ├── config.yml       # datasets, versions, URLs, tool params, threads, data/cache paths
    │   │   └── config_schema.yml
    │   ├── download.smk         # fetch pinned inputs → local cache (retained; not run, §5)
    │   ├── benchmark.smk        # top-level: evaluate → parse → plot
    │   └── rules/
    │       ├── prep.smk         # normalize/split VCFs, build SDF (retained; not run, §5)
    │       ├── evaluate.smk     # one rule per tool × dataset
    │       ├── parse.smk        # tool outputs → unified TSV shards
    │       └── plot.smk         # concatenate shards → figures (+ confusion matrices)
    ├── scripts/
    │   ├── parse_common.py      # size classes, GT/allele counting, shared QUAL sweep
    │   ├── parse_vcfdist.py     # v2 + v3 (same TSV format)
    │   ├── parse_vcfeval.py
    │   ├── parse_happy.py
    │   ├── aggregate.py         # parser dispatch + `time -v` log → TSV shards
    │   ├── confusion_matrix.py  # NEW: per-variant concordance vs. vcfdist v3 (§6.5)
    │   ├── parse_phasing.py     # PLANNED: vcfdist + whatshap phasing metrics (§6.6)
    │   ├── plot_fnr_fdr.py      # FNR/FDR bar charts
    │   ├── plot_pr_curves.py    # NEW: all tools on one plot
    │   ├── plot_runtime_ram.py  # NEW: runtime + peak RAM
    │   └── plot_phasing.py      # PLANNED: switch/flip counts + whatshap confusion (§6.6)
    ├── tests/                   # pytest unit tests for every parser/plotter, tiny text fixtures
    ├── results/                 # genome-wide run (gitignored)
    │   ├── {tool}/{dataset}/... # raw per-tool eval outputs + `time -v` logs
    │   ├── parsed/*.tsv         # unified counts, quality sweeps, runtime, confusion
    │   └── img/*.pdf            # figures
    └── results-chr20/           # same tree, chr20-only `smoke` run
```

`pixi run -e bench download`, `... benchmark`, and `... smoke` task aliases wrap the Snakemake entry points,
matching the template's `pixi run` convention; `pixi run test` runs the Python unit tests in the
cross-platform `default` env.

---

## 5. Data Acquisition & Caching

**As built, the benchmark consumes pre-prepared inputs** and skips download/prep:

- `config.yml` `data_dir` (default `../data`) points at an already-downloaded,
  already-filtered data tree: `<data_dir>/<id>-<version>/split/<id>.most.vcf.gz` for truth and each query
  call-set, `<truth-id>-<version>/split/bench.bed` for the evaluation regions, and
  `<data_dir>/refs/<reference.name>` (plus `.fai` and `.sdf`). `benchmark.smk` reads these directly.
- `download.smk` and `rules/prep.smk` are **retained for from-scratch reproducibility** but not part of the
  current run; `benchmark.smk`'s `include: "rules/prep.smk"` is commented out. Re-running them requires the
  `TODO:` download URLs in `config.yml` to be filled in first (§12).
- Downloads, when run, are idempotent: a rule whose cached output already exists is skipped. The existing
  reference is **moved** into the cache rather than re-downloaded.

*Reproducibility caveat to close before publication:* because the current figures come from a hand-staged
data tree rather than a URL-pinned download, the provenance of the inputs is recorded only by version label.
Filling in the URLs (§12 item 2) is what makes the harness independently re-runnable.

**The AWS run host forces this closed.** A freshly-provisioned instance has no local `data/` tree, so
the hand-staged tree cannot be the input for the published run — the inputs must either be regenerated from
pinned URLs or pulled from S3 on instance start. See §7.3, which is now the operative section for input
staging.

---

## 6. Workflow Stages

### 6.1 prep (`rules/prep.smk`)

Input filtering:

- Filter each query and truth VCF to the evaluated set: keep **SNPs** and indels/SVs with
  **`|ILEN| < 1000`**, and **drop `INV`** records.
- bgzip + tabix each filtered VCF.
- Build the RTG **SDF** from the reference (`rtg format`) once, shared by all vcfeval runs.

These filters are **already applied** in the staged `*.most.vcf.gz` inputs and the staged `refs/*.sdf`, so
the rules do not run in the current pipeline (§5).

### 6.2 evaluate (`rules/evaluate.smk`)

One rule per **tool × dataset** (3 enabled tools × 3 datasets = 9 eval jobs; 12 with hap.py re-enabled, plus
3 `whatshap compare` jobs once §6.6 lands). Every invocation is wrapped in `command time -v` and its stderr
captured to a `*.time.log` for the runtime/RAM figures. The pipeline runs serially (`-j1`) so runtime and RAM
are measured without contention.

- **vcfdist v3** — local `src/vcfdist`, `--bed bench.bed -t {threads} -r {max_ram_gb} -l 1000 -ct 0.98`
  plus the **v3-only** flags `-md {max_dist} -mr {max_retries}`, `-p {outdir}/vcfdist-v3/{ds}.`.
- **vcfdist v2.6.4** — bioconda binary, same arguments **minus `-md`/`-mr`** (v2.6.4 does not accept them
  and errors out).
- **vcfeval** — `rtg vcfeval -b truth -c query -t ref.sdf --bed-regions bench.bed
  --evaluation-regions bench.bed --threads {threads} --ref-overlap --all-records --vcf-score-field=QUAL`.
- **hap.py** (currently disabled, §2) — the same
  `hap.py truth query -r ref.fa -f bench.bed --engine=xcmp -o /out/{ds}` command line, executed inside
  `docker run --rm pkrusche/hap.py:<pinned tag>` with the data/reference paths bind-mounted read-only and the
  results directory writable. On the `x86-64` host its wall-clock is a real, comparable measurement; only its
  **RSS** needs container-aware capture, since `time -v` around `docker run` sees the client (§2).
- **whatshap** (planned, §6.6) — `whatshap compare --names {ds},{truth} --tsv-pairwise ... --switch-error-bed ...`
  per dataset; a phasing oracle, not an evaluation tool, so it is absent from the P/R and resource figures.

Threads, RAM cap, and the credit threshold come from `config.yml`. The committed values (`threads: 10`,
`max_ram_gb: 32`) are sized for the M3 Max development host and **must be re-tuned to the AWS instance before
the published run** (§7.2); `credit_threshold: 0.98` is host-independent and stays pinned for both vcfdist
versions (§2).

### 6.3 parse (`rules/parse.smk` + `scripts/parse_*.py` + `scripts/aggregate.py`)

Each tool's output is parsed into **unified tidy TSVs** — one of single-operating-point counts (for the
bar charts), one of quality-swept precision/recall points (for the curves), and one of runtime/RAM. Because
the tools do not all emit an SV-specific ROC, two methods are used, and the difference is documented on the
figures:

- **vcfdist (v2 & v3)** — **native curve.** `precision-recall-summary.tsv` gives the operating-point counts;
  `precision-recall.tsv` gives the full per-`MIN_QUAL` sweep already split into `SNP/INDEL/SV/ALL` with
  `TRUTH_TP/TRUTH_FN/QUERY_TP/QUERY_FP`. This preserves vcfdist's fractional/partial-credit scoring. v2 and
  v3 share the format, so one parser serves both.
- **vcfeval** — **computed from classified VCFs.** vcfeval emits no SV-specific ROC, so the parser reads the
  per-variant `tp.vcf` / `fp.vcf` (query side, carrying QUAL) and `tp-baseline.vcf` / `fn.vcf` (truth side),
  classes each record with the shared 50 bp rule, and sweeps QUAL to produce a full curve per size class.
- **hap.py** — **computed from the annotated VCF**, same method: query-side `BD=TP`/`BD=FP` and truth-side
  `BD=TP`/`BD=FN` records, classed and swept identically. *Quality field:* hap.py carries its ROC quality in
  the per-sample **`QQ`** FORMAT field and frequently leaves `QUAL` as `.`, so the parser prefers `QQ` and
  falls back to `QUAL` — to be confirmed against real `--engine=xcmp` output, since a wrong choice here
  silently flattens hap.py's curve to a single point.

The QUAL sweep is a single shared function so vcfeval and hap.py are computed identically; per threshold `t`,
`precision = TPq(>=t) / (TPq(>=t) + FP(>=t))` and `recall = TPq(>=t) / (all truth in class)`. The
**methodological difference** — vcfdist's native partial-credit curve vs. the hard-label QUAL sweep for
vcfeval/hap.py — is noted on the PR-curve figures and in the README.

`aggregate.py` dispatches to the right parser per `--tool`, stamps `tool`/`dataset` onto every row, and
writes per-(tool, dataset) shards that `plot.smk` concatenates into `results/parsed/`:

| File | Columns |
| :-- | :-- |
| `counts.tsv` | `tool, dataset, size_class, tp_query, tp_truth, fp, fn` (all-records operating point) |
| `pr_curve.tsv` | `tool, dataset, size_class, min_qual, precision, recall` |
| `runtime.tsv` | `tool, dataset, wall_seconds, max_rss_kb` |
| `confusion_{other}.tsv` | `dataset, other_tool, size_class, vcfdist-v3_label, {other}_label, count` (§6.5) |
| `phasing.tsv` | `tool, dataset, switch_errors, flip_errors, switch_rate, flip_rate, phased_frac` (§6.6) |
| `phasing_confusion.tsv` | `dataset, category, vcfdist_label, whatshap_label, count` (§6.6) |

The **size class** boundary is:

| Class | Rule |
| :-- | :-- |
| SNP | `len(ref) == 1 and len(alt) == 1` |
| INDEL | size difference `> 0` and `< 50 bp` |
| **SV** | size difference **`>= 50 bp`** (up to the 1000 bp prep cap) |

(The old code already used a 50 bp indel/SV split; it is stated explicitly here and applied uniformly across
all tools' parsers.)

### 6.4 plot (`rules/plot.smk` + `scripts/plot_*.py`)

1. **FNR / FDR bar charts** — one panel per size class (SNP / INDEL / SV), datasets on the x-axis, Q-score axis as before.
   *Terminology:* the requested **FPR** is reported as **FDR** = `FP / (FP + TP)` = `1 − precision`, the
   conventional germline-benchmarking false-positive axis; **FNR** = `FN / (FN + TP)` = `1 − recall`. This
   mapping is stated on the figures and in the README.
2. **Precision-recall curves** — **all tools overlaid on one plot**, one plot per
   dataset × size class, drawn from each tool's quality sweep. Consistent per-tool colors across all figures.
3. **Runtime & peak RAM** — parsed from the `time -v` logs (wall-clock + maximum resident set
   size), one grouped bar chart per metric across tools × datasets. The instance type is stated on the figure
   (§7.1); hap.py is footnoted as containerized, and its RAM bar comes from container accounting or is marked
   not captured (§2).
4. **Concordance matrices** — one heat-map per other tool (§6.5).
5. **Phasing figures** (planned) — switch/flip counts and rates per tool × dataset, plus the
   vcfdist-vs-whatshap flip/switch confusion matrices (§6.6).

### 6.5 concordance / confusion matrices (`scripts/confusion_matrix.py`)

This is the SOW's **"safe to displace" evidence**: SOW #7 requires validating that vcfdist's TP/FP/FN
classifications agree with the incumbent on simple variants, and characterizing the residual disagreements
(expected around complex variants, where vcfdist is *designed* to differ).

- **Reference tool:** vcfdist v3. One matrix per *other* tool (`OTHER_TOOLS = ALL_TOOLS − vcfdist-v3`), so
  the same machinery yields both the **incumbent** comparison (vcfeval today, hap.py once enabled) and the
  **v2→v3 regression** comparison.
- **Join key:** standard-VCF `(CHROM, POS, REF, ALT)`, read from each tool's own classified output —
  vcfdist's `summary.vcf` (`BD` per `TRUTH`/`QUERY` sample), vcfeval's `tp`/`fp`/`fn` VCFs, hap.py's
  annotated two-sample VCF. Multi-allelic records are split per ALT; per-haplotype duplicate rows are
  collapsed with TP winning.
- **Labels:** `{TP, FP, FN, N}`, where **`N` = the variant is absent from that tool's evaluated set** — the
  representation-difference bucket (e.g. a repeat indel the two tools left-align differently). The `N` row
  and column are what make representation disagreement countable rather than invisible.
- **Output:** `results/parsed/confusion_{other}.tsv` (long form, per dataset × size class) and
  `results/img/confusion_{other}.pdf`. The diagonal is agreement; off-diagonal cells are the disagreements
  to bucket and explain in the report.
- **Interpretation in the report:** SNP/INDEL off-diagonal mass is the trust-building number (should be
  small); SV and `N` mass is where vcfdist's credit-based decisions legitimately diverge, and each bucket
  gets a worked example. The chr20 smoke run already shows the expected shape — e.g. HPRC SNPs at
  78,429 TP/TP against 20 `TP/FP` and 91 `TP/N` cells.

*Caveat:* the join is exact-coordinate, so a variant the two tools represent at different positions lands in
the `N` bucket rather than being recognized as the same call. That is deliberate (it is the quantity being
measured), but it means `N` mass is an **upper bound** on representation disagreement, not a count of
genuinely missing calls. #3's `BD=N` retention makes vcfdist's own not-assessed set explicit, which
tightens this.

### 6.6 Phasing correctness (`scripts/parse_phasing.py`, `scripts/plot_phasing.py`)

Rationale in §11. vcfeval and hap.py report no phasing metrics, so this stage has two comparison axes
instead of four tools:

- **v2 vs. v3 — the regression signal.** v3 moved phase decisions from per-supercluster to **per-variant**
  (`D1` §3.2), so switch/flip counts are not guaranteed comparable across versions. Both versions already run
  on every dataset, so their `phasing-summary.tsv` outputs are free to collect.
- **vcfdist vs. `whatshap compare` — the independent oracle.** `whatshap compare --names {ds},{truth}
  --tsv-pairwise {ds}.whatshap.tsv --switch-error-bed {ds}.switches.bed` scores the query's phasing against
  the T2T-Q100 truth phasing by a completely different method (block-wise switch/flip counting over shared
  heterozygous sites, no alignment or credit model), which is what makes agreement meaningful.

**Stages.** Add `whatshap` to the `bench` pixi feature and a `rule eval_whatshap` per dataset; then:

- `parse_phasing.py` → `results/parsed/phasing.tsv`
  (`tool, dataset, switch_errors, flip_errors, switch_rate, flip_rate, phased_frac`), reading vcfdist v2's and
  v3's `phasing-summary.tsv` and whatshap's pairwise TSV (switch/flip arrive as a single `switches/flips`
  field summed across blocks; parse it **by header name** against the pinned whatshap release rather than by
  position, since the column order has changed across versions). `phased_frac` is vcfdist-only and blank for whatshap.
- `plot_phasing.py` → `results/img/phasing.pdf` (grouped switch/flip bars per dataset, v2 vs. v3 vs. whatshap)
  and `results/img/phasing_confusion_{ds}.pdf` — flip and switch confusion matrices joining vcfdist's
  `switchflips.tsv` events against
  whatshap's `switches.bed` by overlap and cross-tabulating category labels
  (`NONE / SWITCH / FLIP`, plus vcfdist's `UNKNOWN` states) into `phasing_confusion.tsv`.
- Unit tests over tiny inline fixtures for both, matching the existing parser/plotter test pattern.

**What this is expected to show, and the gates it feeds:** switch/flip counts within a few percent of v2 (a
large divergence is a per-variant-phasing regression to investigate, not a new result), and vcfdist-vs-whatshap
agreement comparable to the published v2 matrices. It also establishes the **pre-#6 phasing baseline** the
baseline §10 describes, since `D6` redefines these metrics over the phaseable subset and adds
`phased_frac`.

**Known limits, stated on the figures:** whatshap scores only sites both call-sets genotype heterozygous, so
its denominator differs from vcfdist's; and both v2's and v3's phase-group construction requires adjacent
variants (`D1` §4), so neither tool's switch/flip counts are absolute truth — the comparison is agreement
between two independent methods, not accuracy against a gold standard.

---

## 7. Execution Environment — On-Demand AWS `x86-64`

All published figures come from a single run on a purpose-provisioned **on-demand AWS EC2 `x86-64`
instance** (§2). The instance is created for the run and terminated after the results are pulled off it;
nothing about the harness assumes a long-lived machine.

### 7.1 Instance specification

| Property | Value | Why |
| :-- | :-- | :-- |
| Platform | **`linux-64` (`x86-64`)** | Native hap.py — no emulation, so its runtime is a measurement (§2); also removes the macOS-only shims. |
| Instance type (primary) | **`m6i.16xlarge`** — 64 vCPU / 256 GiB | Restores the published v2 harness's 64-thread operating point, so runtimes are comparable to the prior figures. |
| Instance type (economy) | `m6i.8xlarge` — 32 vCPU / 128 GiB | Sufficient for correctness and cross-tool ordering; runtimes then are *not* comparable to the published 64-thread numbers. |
| Purchase model | **On-demand — not spot** | A spot interruption mid-run invalidates every timing measurement in the run. |
| Storage | ~200 GB `gp3` EBS | Reference (~2.9 GB) + RTG SDF + filtered VCFs + four tools' outputs, with headroom. |
| Region | Pinned in `config.yml` | Part of the reproducibility record; also where the S3 staging bucket must live (§7.3). |

**Recorded into `results/` for every run** (alongside the tool versions and `dev` commit hash, §8):
instance type, AMI id, region, kernel version, EBS volume type/size, and the hap.py Docker image **digest**
(not just its tag). Without the instance type, the runtime figure is uninterpretable.

### 7.2 Configuration changes this implies

- **`params.threads: 10` → 64** (or 32 on the economy instance) and **`params.max_ram_gb: 32` → ~200**
  (~100 on the economy instance). The current values were sized for the 48 GB M3 Max.
- **`max_ram_gb` is a parameter of the measurement, not just a guardrail** — vcfdist's `-r`/`--max-ram`
  feeds its supercluster RAM scheduling, so it changes *what the tool does*. It must be identical for v2 and
  v3 (as `credit_threshold` already is, §2) and reported with the figures. Re-tuning it also means the
  chr20 `smoke` numbers are not comparable to the genome-wide run — fine, since smoke exists to check the
  DAG, not to produce figures.
- **`-j1` (serial) stays**, so each tool's runtime and RAM are uncontended even though the instance is large.
- **`command time -v` (conda `time`) can stay as-is.** It is now redundant — the instance has GNU
  `/usr/bin/time` — but keeping it means the same rules still run on the macOS development host for smoke
  runs. One shim, two platforms.
- **`pixi`**: `bench` (`analysis` + `tools`) is the run environment on `linux-64`; `default` stays
  cross-platform for tests and dry-runs. The `happy` feature becomes optional (Docker path, §2) but installs
  natively here if the bioconda fallback is taken.
- **Docker** must be installed on the instance (or an AMI with it preinstalled selected) for the hap.py rule.

### 7.3 Data staging

An ephemeral instance has **no local `data/` tree**, so the hand-staged inputs the current run reads
(§5) stops being a valid input. Two paths, and the deliverable needs both:

- **(a) Regenerate on the instance** — fill in the `TODO:` download URLs so `download.smk` + `prep.smk`
  rebuild the filtered VCFs and the SDF from public sources. This is the *documented* reproducibility story
  and the reason the retained-but-unused stages exist. It moves URL pinning from a nice-to-have to a
  a prerequisite for a from-scratch run (§12 item 2).
- **(b) Mirror the prepared tree to S3** — `aws s3 sync` the existing `data/` tree (with the S3 prefix
  recorded in `config.yml`) and pull it on instance start. Practically faster and avoids re-running prep and
  `rtg format` on every instance, at the cost of trusting a private artifact.

**Recommendation: do both** — (a) as the provenance record a third party can follow, (b) as the mechanism the
actual runs use, with the S3 objects' checksums recorded so the two are demonstrably the same inputs.

**Before terminating:** sync `results/` (parsed TSVs, figures, `*.time.log`, tool-version and instance
records) off the instance. Raw per-tool eval output is large and regenerable; the parsed tables and logs are
not, and they are what the report is written from.

### 7.4 Cost

`m6i.16xlarge` on-demand runs on the order of **$3/hour** in `us-east-1` (verify current pricing before
provisioning — this figure is not pinned). A full serial run of four tools × three datasets plus staging is
expected in the **several-hours** range, so **roughly $25–75 per full run**, and few full runs are needed.
That sits comfortably inside the SOW's **"Materials — Compute / AWS — < $500"** line. Terminate the instance
when done: EBS volumes bill while they exist, even for a stopped instance.

> **SOW amendment needed.** `D0`'s Materials *prose* states there are no
> anticipated material costs and that "benchmarking can be run on existing hardware (Fulcrum Apple M3 Max
> MacBook Pro with 48GB RAM)". That is no longer the plan. Its cost table already anticipates the spend
> (the `Materials — Compute / AWS — < $500` row), so the fix is a one-line prose amendment naming the
> on-demand `x86-64` instance and why. Not edited here.

---

## 8. Reproducibility

- All dataset versions, URLs, tool versions, and parameters live in `workflows/config/config.yml`, validated
  against `config_schema.yml`.
- `pixi.lock` pins the exact environment; `pixi run` builds it on first use.
- Tool versions are recorded into `results/` (each tool's `--version`) alongside the figures.
- The git commit hash is stamped into the run outputs (template convention) so every figure is traceable to
  a code state. Because the tool under test is a **branch build**, the `dev` commit hash is the only record
  of what "vcfdist v3" meant for a given figure — it must appear in the report.
- Every parser and plotter has pytest unit tests over tiny inline text fixtures (`pixi run test`), so the
  Python core is verifiable without the multi-hour pipeline.
- **The execution environment is part of the record, not an implementation detail** — instance type, AMI,
  region, kernel, EBS configuration, and the hap.py image digest are captured per run (§7.1). Runtime and RAM
  figures are meaningless without them, and cross-run comparison is only valid within one instance type.
- **So is the operating point.** Every vcfdist parameter that differs from the shipped default is recorded
  with the figure, and the published figures' operating point must match the defaults the release ships
  (§12 item 11, [`D8`](./D8_vcfdist-v3-release.md) §4.3). `credit_threshold` already sets this precedent by
  being pinned and stated; `-md`, `-mr`, and `-sc` are held to the same standard.
- **The counting convention is part of the record too.** After #4 a homozygous variant counts once, not
  twice (§10). Any figure comparing v3 to v2 or to a published number states which convention each side
  used.

---

## 9. Current State

The harness is **built and unit-tested**, and has completed an end-to-end **chr20 smoke run**. Its
committed configuration:

| Setting | Value |
| :-- | :-- |
| **hap.py** | Wired and tested, run from the pinned Docker image on the `x86-64` host (§2, §7) |
| **Inputs** | Pre-staged `data_dir` tree; `download.smk` and `rules/prep.smk` retained for from-scratch reproducibility (§5) |
| **`credit_threshold`** | Pinned to `0.98` for both v2 and v3, since their defaults differ (§2) |
| **`-md`/`-mr`** | Passed to v3 only; v2.6.4 rejects them |
| **`threads` / `max_ram_gb`** | Re-tuned to the AWS instance for the published run — 64 / ~200 on `m6i.16xlarge` (§7.2) |
| **`-j1` (serial)** | Retained, so each tool's runtime and RAM are uncontended |
| **`command time -v`** | The conda `time` package, so the same rules run on macOS for smoke checks (§7.2) |
| **`smoke` task** | chr20 only, into a separate `results-chr20/` tree — a minutes-scale end-to-end check |
| **`bench` pixi env** | Installs on `osx-arm64` and `linux-64`. **All published figures come from the `linux-64` AWS host** (§7); the extra platform is a development convenience, and runtime/RAM numbers from a macOS smoke run are never reported |
| **`confusion_matrix.py`** | Implements the concordance requirement (§6.5), sharded per (tool, dataset) by `aggregate.py` |

**To finish the deliverable:** (a) the hap.py Docker eval rule, including container-side RAM capture (§2);
(b) the phasing stage (§6.6) — whatshap, two scripts, tests; (c) input staging for an ephemeral host — URL
pinning plus the S3 mirror (§7.3); (d) instance provisioning and the config re-tune (§7.1–7.2); (e) the
genome-wide run on that instance; (f) the written report with regenerated figures and the
disagreement-bucket narrative.

---

## 10. What This Harness Measures for Other Deliverables

Three changes move the numbers this harness reports, and each is measured and attributed separately rather
than as one aggregate delta:

| Change | Effect on the numbers | What the harness does |
| :-- | :-- | :-- |
| **#4 per-site counting** ([`D4`](./D4_ga4gh-compatibility.md) §4) | A homozygous variant contributes 1 unit instead of 2, so **every count and every derived rate changes**. Largest of the three, and it changes what a v2-vs-v3 comparison means, since v2.6.4 counts per haplotype. | Capture the baseline, re-run after, report the delta. For the v2 comparison, state explicitly that hom variants are weighted differently in the two versions — convention-mismatched for those sites, and not presentable as like-for-like. |
| **#6 clustering restructure** ([`D6`](./D6_vcfdist-v3-unphased-eval.md) §4.2) | Moves the haplotype merge ahead of clustering, shifting supercluster boundaries and therefore P/R **even for phased input**. | Capture a baseline, re-run after, and quantify the shift. |
| **PR #92 single-pass** (D8 §2.2) | Replaces the retry-based FN-dropping loop. | Register `vcfdist-v3-singlepass` as its own tool id, so one genome-wide run yields baseline and comparison together. |

Per-deliverable couplings:

- **#3 record shape + `BD=N`** ([`D3`](./D3_retain-info-format-fields.md)) — rewrites `summary.vcf`, which
  `confusion_matrix.py` parses: a homozygous variant is **one** record rather than two, so the parser's
  per-haplotype duplicate-row collapsing (§6.5) must not assume two rows. `BD=N` retention makes vcfdist's
  not-assessed set explicit and tightens the `N` bucket (§6.5 caveat).
- **#4 GA4GH compatibility** ([`D4`](./D4_ga4gh-compatibility.md)) — beyond the counting change above,
  renames `FORMAT/SC` → `INFO/BS` and `FORMAT/BS` → `FORMAT/PBS`, and makes `BD`/`BK` single-valued per
  site. `confusion_matrix.py` reads `BD` per sample, so it is re-validated against the new shape.
- **#5 Multi-BED stratification** ([`D5`](./D5_multi-bed-stratification.md)) — adds a leading `STRATUM`
  column to every aggregate TSV and splits `phasing-summary.tsv` into `phasing-variants-summary.tsv` and
  `phasing-blocks-summary.tsv`. The parsers filter to `STRATUM == "*"` and read the new filenames. A
  stratified cross-tool comparison is a follow-on, not part of these figures.
- **#6 Unphased evaluation** — redefines the switch/flip denominator as `ASSESSED_PAIRS` and adds
  phased-fraction columns (D6 §5.3), so the phasing baseline (§6.6) is captured on the old denominator and
  the change noted rather than diffed. A second run with query phase stripped is the companion comparison —
  the only way to put vcfdist head-to-head with vcfeval and hap.py on the input *they* accept.
- **#8 Release / announcement** — the accuracy and concordance figures feed the blog post, and the parameter
  sweep feeds the default-value decision in [`D8`](./D8_vcfdist-v3-release.md) §4.3 (G12). See §12 item 8.
- **#12 Faster alignment** — the runtime/RAM figures are the before-baseline for any alignment speedup. If
  #12 proceeds, that comparison becomes a configured tool id rather than a side artifact.

---

## 11. Phasing-Correctness Benchmarking — Rationale

Phasing validation is in scope for #7; the stage design is §6.6. Three reasons it belongs here:

1. **Phasing is the stated differentiator.** The SOW's executive summary sells vcfdist on uniquely reporting
   flip and switch errors, and the P/R figures validate only the other half of that claim.
2. **v3 changed how phasing works.** Per-supercluster phase decisions became **per-variant** (D1 §3.2), so
   v3's switch/flip counts are not guaranteed comparable to v2's — exactly the kind of regression #7 exists
   to catch.
3. **#6 perturbs it again.** [`D6`](./D6_vcfdist-v3-unphased-eval.md) §5 redefines the phasing metrics over
   the *phaseable subset* and adds a reported phased fraction, so a baseline is what makes "no phasing
   regression" assertable at release.

**Scope** — one increment; the stage design is §6.6. **Cost estimate: 2–4 hours.**

**Deliberately excluded:**

- **A third-party phasing benchmark beyond whatshap** (e.g. `vcfeval --squash-ploidy` tricks or
  `whatshap compare` against a second oracle) — one independent method is enough to catch a regression, and
  no incumbent benchmarking tool reports switch/flip at all.
- **Worked flip examples** — useful for the paper narrative, not for a regression check; revisit if the
  report needs illustrations.
- **Unphased-input phasing behavior** — belongs to #6's own validation, over the metrics `D6` defines.

---

## 12. Open Items

1. **hap.py execution.** Write the `docker run` eval rule with bind mounts and container-side paths, pin
   the image **digest**, and solve container RAM capture (`time -v` around `docker run` measures the
   client, not the container) — or annotate hap.py's RAM bar as not captured (§2).
2. **Input staging for an ephemeral host.** Truth, query, and reference URLs are still `TODO:` in
   `config.yml`, and a fresh instance has no local input tree. Pin the URLs so the inputs regenerate on
   the instance, and mirror the prepared tree to S3 (§7.3).
3. **Genome-wide run.** Only the chr20 `smoke` run has completed end-to-end. All published figures come
   from a single run on one documented instance.
4. **vcfeval version pin.** RTG Tools version to pin in `config.yml`.
5. **Instance confirmation.** Confirm the type (64 vCPU `m6i.16xlarge` preferred, matching the published
   v2 harness's `-t 64`), and verify current on-demand pricing against the SOW's `< $500` compute line
   (§7.1, §7.4). Absolute runtimes are comparable to the published v2 figures only insofar as the vCPU
   count and thread settings match — state the instance type on the runtime figure either way.
6. **`credit_threshold` drift.** If the v3 default moves off 0.98, `config.yml` is updated to match, or
   the v2-vs-v3 comparison silently changes meaning (§2).
7. **whatshap TSV column names.** Check the pairwise-TSV header against the pinned whatshap release; the
   column order has changed across versions (§6.6).
8. **The harness runs vcfdist v3 at non-default `-md 1000 -mr 10`.** `workflows/config/config.yml` sets
   `max_dist: 1000` and `max_retries: 10`, while v3's defaults are `max_dist = 100` and
   `max_retries = 0` (`src/globals.h`). Left as-is, the figures would describe an operating point a user
   does not get by default while the release claims accuracy on them, and PR #92's decision rule — the PR
   *replaces* the retry loop — would be evaluated against a retry-enabled baseline the default never
   exercises. Resolution: run at the defaults, or **sweep** `-ct`, `-md`, `-mr`, and `-sc` and let the
   sweep be the evidence for [`D8`](./D8_vcfdist-v3-release.md) §4.3's default-value decision (G12). The
   sweep is the better use of one AWS run, since it produces the release evidence *and* the figures.
   Either way, the published operating point and the shipped defaults agree, and the operating point is
   stated on every figure.
9. **AWS operational risk.** The run is a single long serial job on a machine that bills by the hour: an
   interrupted run wastes the spend, and results left on a terminated instance are gone. Mitigated by the
   S3 input mirror (§7.3), syncing `results/` before termination, and running the chr20 `smoke` check
   **on the instance** before launching the genome-wide run.
10. **SOW Materials prose.** `D0`'s Materials section names the Fulcrum M3 Max
    as sufficient; the published figures come from the AWS instance instead. Its cost table already
    carries the `Compute / AWS — < $500` row, so this is a one-line prose amendment (§7.4).

---

## 13. References

- Companion design doc: [`D1_vcfdist-v3-design-doc.md`](./D1_vcfdist-v3-design-doc.md) §5.7
- Companion SOW: `D0_vcfdist-v3-SOW.md` Deliverable #7
- Unphased-evaluation design (release-gate coupling, §10): [`D6_vcfdist-v3-unphased-eval.md`](./D6_vcfdist-v3-unphased-eval.md)
- Martin et al. 2016, *WhatsHap: fast and accurate read-based phasing* (bioRxiv) — `whatshap compare`, the
  phasing oracle in §6.6
- hap.py Docker image: https://hub.docker.com/r/pkrusche/hap.py
- AWS EC2 instance types (`m6i` family) — https://aws.amazon.com/ec2/instance-types/m6i/
- AWS EC2 on-demand pricing (verify before provisioning, §7.4) — https://aws.amazon.com/ec2/pricing/on-demand/
- Fulcrum tertiary-analysis-template — https://github.com/fulcrumgenomics/tertiary-analysis-template
- Krusche et al. 2019, *Best practices for benchmarking germline small-variant calls*, Nat. Biotechnol.
- Dunn & Narayanasamy 2023, *vcfdist*, Nat. Commun.
