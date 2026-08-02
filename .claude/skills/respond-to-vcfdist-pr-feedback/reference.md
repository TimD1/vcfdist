# Reference — A/B run and comment template

Mechanics for `respond-to-vcfdist-pr-feedback`. The decisions — whether to run an A/B, which
tier, and the gates around pushing — are in `SKILL.md`.

## The A/B run

Both sides use the committed chr20 fixtures in `tests/integration/data/`: a chr20-only
reference, truth and query VCFs, and a BED, all four tracked in git. No external data, works in
a fresh clone, and the same inputs `tests/integration/test-integration.yml` drives.

**PR side**, from the branch worktree after step 3 has built `src/vcfdist`:

```bash
D=tests/integration/data
./src/vcfdist $D/query_chr20.vcf.gz $D/truth_chr20.vcf.gz \
  $D/GCA_000001405.15_GRCh38_no_alt_analysis_set_chr20.fasta \
  -b $D/chr20.bed -p out/pr.
```

**Target side**, in a separate worktree so both binaries exist at once:

```bash
git worktree add /tmp/vcfdist-base origin/<base>
cd /tmp/vcfdist-base/src && make && cd ..
D=tests/integration/data                      # re-set: a new shell has no $D
./src/vcfdist $D/query_chr20.vcf.gz $D/truth_chr20.vcf.gz \
  $D/GCA_000001405.15_GRCh38_no_alt_analysis_set_chr20.fasta \
  -b $D/chr20.bed -p out/base.
```

Target output changes only when the target moves — cache it keyed on
`git rev-parse origin/<base>` rather than rebuilding each iteration. Then:

```bash
diff out/base.precision-recall-summary.tsv out/pr.precision-recall-summary.tsv
```

**Escalation tier**, only for accuracy-critical changes (alignment, clustering, credit
assignment, phasing) where chr20 fixtures may not reach the affected path:

```bash
cd analysis-v3/vs_prior_work && pixi run -e bench smoke   # writes results-chr20/
```

Slow, serial, and needs the 2.9 GB genome-wide reference in `analysis-v3/data/`. Name the tier
you ran either way.

## Comment template

One visible summary line, then every section in a collapsed `<details>` block so the PR
conversation stays scannable. **A blank line after `</summary>` is required** or GitHub renders
the markdown inside as literal text.

````markdown
**<model> 🤖:** <one line: what was addressed, and the headline count impact>

<details>
<summary><b>Responses to review comments</b> (5 addressed, 1 deferred)</summary>

- `src/dist.cpp:190` — renamed `relax()` to `<new>`. Resolved.
- `tests/…:15` — moved the helper as suggested. Resolved.
- `src/dist.cpp:165` — **not done**: <one line why>. Left unresolved, replied in thread.

</details>

<details>
<summary><b>Count impact</b> — no change on the chr20 fixtures</summary>

Tier: chr20 fixtures (`tests/integration/data`), `chr20.bed`.
`precision-recall-summary.tsv` is byte-identical between base and PR.

</details>

<details>
<summary><b>Proposed root causes</b> (3 variants reclassified)</summary>

<!-- only when counts moved; omit this block entirely otherwise -->

</details>
````

Section contents are unchanged from the rules below; only the wrapping is new.

**1. Summary line.** Begins `**<model> 🤖:**`, required by the repo's GitHub-write hooks.
**REQUIRED SUB-SKILL:** `github-ai-authorship` owns the format. Editing the PR *body* also needs
the verbatim `> [!NOTE]` authorship block, which that skill supplies.

**2. Responses to review comments.** One line per restated item: what changed, and for anything
not done, why not. Say which threads you resolved — the summary is where a reader reconciles the
list against the resolved threads.

**3. Count impact.** Always present. Name the tier you ran and what it showed; when output was
identical that sentence is the whole section, and the `<summary>` line should say so.

When counts moved, give the table — one row per variant type at both the `NONE` and `BEST`
thresholds, read from `precision-recall-summary.tsv`. Include raw counts: a delta without
`TRUTH_TP` / `TRUTH_FN` / `QUERY_FP` cannot be audited.

| VAR_TYPE | THRESHOLD | PREC (base → PR) | RECALL (base → PR) | F1 (base → PR) | ΔF1 | TRUTH_TP | TRUTH_FN | QUERY_FP |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |

State which fixtures and BED produced it, and whether each move was expected.

**4. Proposed root causes** — only when counts moved. Diff `query.tsv` / `truth.tsv` between the
runs, bucket the variants that changed classification, one bucket per paragraph with a count and
**at least one `CONTIG:POS REF>ALT` example**. Write each as a hypothesis. If zero variants
changed classification, say so in the count-impact block — an absent statement reads as an
omission.


## Troubleshooting

**A test errors on a missing reference.** `tests/README.md` says several tests need a
genome-wide FASTA in `data/refs/`, which ships **empty**. The committed integration test
actually uses the chr20 reference in `tests/integration/data/`, so that instruction is the
likely culprit rather than a real missing dependency.
