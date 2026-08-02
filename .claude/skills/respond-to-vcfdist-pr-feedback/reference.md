# Reference — A/B run and comment template

Mechanics for `respond-to-vcfdist-pr-feedback`. Open this when you have reached the A/B run or
are about to write the comment. The decisions — whether to run an A/B at all, which tier, and
the gates around pushing — live in `SKILL.md` and are not repeated here.

## The A/B run

Both sides use the committed chr20 fixtures in `tests/integration/data/`: a chr20-only
reference, truth and query VCFs, and a BED, all four tracked in git. The comparison therefore
needs no external data and works in a fresh clone — these are the same inputs
`tests/integration/test-integration.yml` already drives.

### PR side

From the PR branch's worktree, after step 3 of the sequence has built `src/vcfdist`:

```bash
D=tests/integration/data
./src/vcfdist $D/query_chr20.vcf.gz $D/truth_chr20.vcf.gz \
  $D/GCA_000001405.15_GRCh38_no_alt_analysis_set_chr20.fasta \
  -b $D/chr20.bed -p out/pr.
```

### Target-branch side

Build the target branch in a **separate worktree** so both binaries exist at once, then run the
identical command with `-p out/base.`:

```bash
git worktree add /tmp/vcfdist-base origin/<base>
cd /tmp/vcfdist-base/src && make
cd /tmp/vcfdist-base
D=tests/integration/data                      # re-set: a new shell has no $D
./src/vcfdist $D/query_chr20.vcf.gz $D/truth_chr20.vcf.gz \
  $D/GCA_000001405.15_GRCh38_no_alt_analysis_set_chr20.fasta \
  -b $D/chr20.bed -p out/base.
```

The worktree carries its own copy of the fixtures, so `$D` resolves there. Collect both result
sets in one place before diffing.

Target-branch output changes only when the target moves, so cache it keyed on
`git rev-parse origin/<base>` rather than rebuilding every iteration.

### Compare

```bash
diff out/base.precision-recall-summary.tsv out/pr.precision-recall-summary.tsv
```

- **Identical** → the count-impact section is one sentence. A table of zero-deltas is noise.
- **Different** → build the table and the root-cause buckets below.

### Escalation tier

Only for accuracy-critical changes — alignment, clustering, credit assignment, phasing — where
the chr20 fixtures may not exercise the affected path:

```bash
cd analysis-v3/vs_prior_work && pixi run -e bench smoke   # writes results-chr20/
```

Slow, serial (`-j1`), and needs the 2.9 GB genome-wide reference in `analysis-v3/data/`. Name
the tier you ran in the comment either way.

### Never quote runtime or RAM

`src/Makefile` defaults to `CXXFLAGS = -g -pg -O1`, with `-O3` commented out. Accuracy is
unaffected; timings from this build are meaningless.

## Comment template

Four parts, in this order. Print the full body in your reply before posting it.

### 1. Prefix

First line begins `**<model> 🤖:**`, required by the repo's GitHub-write hooks. **REQUIRED
SUB-SKILL:** `github-ai-authorship` owns the current format. If you also edit the PR *body*, it
needs the verbatim `> [!NOTE]` authorship block, which that skill supplies.

### 2. Per-item responses

One line per item enumerated in step 1 of the sequence: what changed, the commit, and for
anything not done, why not.

### 3. Count impact

Always present, even when it is one sentence. Name the tier you ran and what it showed.

When counts moved, give the table — PR vs target branch, one row per variant type, at both the
`NONE` and `BEST` thresholds, read from `precision-recall-summary.tsv`. Include the raw counts:
a delta without `TRUTH_TP` / `TRUTH_FN` / `QUERY_FP` cannot be audited.

| VAR_TYPE | THRESHOLD | PREC (base → PR) | RECALL (base → PR) | F1 (base → PR) | ΔF1 | TRUTH_TP | TRUTH_FN | QUERY_FP |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |

State which fixtures and BED produced it, and whether each move was expected from the change.

### 4. Root causes — only when counts moved

Diff `query.tsv` / `truth.tsv` between the two runs, bucket the variants that changed
classification, and give one bucket per paragraph with a count and **at least one concrete
`CONTIG:POS REF>ALT` example per bucket**.

Label the section **Proposed root causes** and write each as a hypothesis. You are inferring
mechanism from output tables, that inference can be wrong, and the comment is public. A cited
coordinate a reviewer can check beats a confident explanation they cannot.

If zero variants changed classification, say so explicitly — it is a meaningful result, and an
absent section reads as an omission.
