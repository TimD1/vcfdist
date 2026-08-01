---
name: respond-to-vcfdist-pr-feedback
description: Use when a vcfdist pull request authored by TimD1 or TimD1-bot has received blocking review feedback (a review in CHANGES_REQUESTED state, or a maintainer asking for changes) and that feedback needs to be addressed, verified, and pushed. Does not apply to pull requests from any other author, including outside contributors. Triggers include "address the PR feedback", "the reviewer requested changes", "respond to that review".
---

# Respond to blocking PR feedback

Address a `CHANGES_REQUESTED` review on a vcfdist PR: fix, verify locally, re-review, measure
whether the change moved any counts, push, and comment.

## Scope gate — run this first

One query answers all three gates:

```bash
gh pr view <N> --json author,isCrossRepository,headRefName,baseRefName,reviews,url \
  --jq '{author: .author.login, fork: .isCrossRepository, base: .baseRefName,
         head: .headRefName, reviews: [.reviews[] | {author: .author.login, state: .state}]}'
```

**Gate 1 — ownership.** `.author.login` must be exactly `TimD1` or `TimD1-bot`. Any other
author: **stop immediately and report that the PR is out of scope.** Do not read the feedback,
do not check out the branch, do not comment.

vcfdist is a public repository that takes outside contributions. Pushing commits to a
contributor's PR — even correct ones, even when they asked for help — silently rewrites someone
else's work under their name. If a contributor's PR genuinely needs changes, the maintainer
reviews it and the contributor pushes; that is not this skill's job. "Allow edits by
maintainers" makes this *possible* on fork PRs, which is exactly why the gate is explicit
rather than relying on permissions to stop you.

**Gate 2 — blocking feedback.** Some review must be in `CHANGES_REQUESTED` state, or a
maintainer must have explicitly asked for changes. If not, **stop and report that.** Do not act
on `COMMENTED` reviews, nitpicks, or your own reading of what a comment implies. A single
unrequested force-push onto a PR costs more trust than a day of waiting.

**Gate 3 — not already handled.** Look for an existing comment from you that post-dates the
review. If one exists, stop; re-running duplicates commits and comments.

## Hard rules

- **Never write to a PR authored by anyone but `TimD1` or `TimD1-bot`** — no commits, no
  pushes, no comments. Being *able* to push is not authorization. There is no exception for
  "the fix is trivial", "the contributor asked", or "the PR is stalled".
- **Never merge.** Not `gh pr merge`, not a local merge into `dev` or `master`. Tim
  squash-merges every branch himself. (`gh pr merge` and `gh pr review --approve` are
  hook-denied anyway — treat the denial as correct, not as an obstacle to route around.)
- **Never force-push** unless the PR branch contains only your own commits. Check first:
  `git log --format='%an' origin/<base>..HEAD | sort -u`. If any name but yours appears, append
  commits instead.
- **Work in the branch's existing worktree** under `.claude/worktrees/`, per the repo's
  CLAUDE.md. Do not create a second worktree for a branch that already has one.
- **Do not push before `pytest -vv` passes.**

## Sequence

**1. Read the feedback and restate it.** Enumerate each requested change as a discrete item
before editing anything. Items you decide not to act on must be listed with a reason — silent
omission is what makes a reviewer re-review from scratch.

**2. Fix.** Follow the repo's documentation conventions from CLAUDE.md (`.h` files get
`@brief` only; `.cpp` files get full parameter/return/throws blocks; 100-character separators).

**3. Test.**

```bash
cd src && make            # produces src/vcfdist
cd ../tests && pytest -vv # unit (GoogleTest) + integration (pytest-workflow)
```

Tests must pass before proceeding. If one fails for a reason unrelated to your change, say so
explicitly with the output rather than proceeding quietly. (`tests/README.md` claims several
tests need a genome-wide FASTA in `data/refs/`, which ships **empty** — but the committed
integration test uses the chr20 reference in `tests/integration/data/`. If a test errors for a
missing reference, that stale instruction is the likely cause.)

## Deciding whether to re-evaluate counts

Most review feedback cannot change a single count. Re-evaluating anyway wastes time and pads
the comment with a table of zeroes. Pick the tier from what the diff actually touches:

```bash
git diff --stat origin/<base>...HEAD
git diff -w origin/<base>...HEAD -- src/   # -w: ignore whitespace-only churn
```

| What the diff touches | Do this |
| :-- | :-- |
| Nothing under `src/` (docs, tests, CI, comments in other trees) | **No evaluation.** State in the comment that no evaluated code changed |
| `src/` but only comments, docstrings, whitespace, or identifier renames | **Confirm, don't assume** — run the A/B below and expect identical output. A "pure rename" that moves a number is exactly the bug worth catching |
| Any change to executable logic in `src/` | **Full evaluation** — run the A/B and report the table |

Judging "this is just a nitpick" by eye is not the gate; the A/B is. It is cheap enough that
the only case for skipping it is that no `src/` file changed at all.

### The A/B run

Both sides use the committed chr20 fixtures in `tests/integration/data/` — chr20-only
reference, truth and query VCFs, and BED, all tracked in git, so this works in a fresh clone
with no external data:

```bash
D=tests/integration/data
./src/vcfdist $D/query_chr20.vcf.gz $D/truth_chr20.vcf.gz \
  $D/GCA_000001405.15_GRCh38_no_alt_analysis_set_chr20.fasta \
  -b $D/chr20.bed -p out/pr.
```

Build the target branch's binary in a separate worktree and run it the same way with
`-p out/base.`, so both result sets exist at once. Then:

```bash
diff out/base.precision-recall-summary.tsv out/pr.precision-recall-summary.tsv
```

**Identical → say so and stop there.** "No change to any count on the chr20 fixtures" is the
result; a table of zero-deltas is noise. **Different → build the table and the root causes.**

Target-branch output only changes when the target moves, so cache it keyed on
`git rev-parse origin/<base>` instead of rebuilding each iteration.

**Escalate to the full harness only when the change is accuracy-critical** — alignment,
clustering, credit assignment, or phasing logic — where chr20 fixtures may not exercise the
affected path. That is `analysis-v3/vs_prior_work` (`pixi run -e bench smoke`), which is slow,
serial, and needs the 2.9 GB genome-wide reference in `analysis-v3/data/`. Say in the comment
which tier you ran.

**Never quote runtime or RAM** from these runs. `src/Makefile` defaults to
`CXXFLAGS = -g -pg -O1` with `-O3` commented out — accuracy is unaffected, timings are
meaningless.

## Independent review

Before pushing, get fresh eyes on the diff — your own second read is not review.
**REQUIRED SUB-SKILL:** use whichever reviewer-persona review skill is available in the
environment, matched to whoever reviews this repo; for changes to alignment, clustering, or
variant logic also use `reviewing-bioinformatics-code`. Apply findings you agree with; record
why for those you reject. **REQUIRED SUB-SKILL:** `superpowers:receiving-code-review` governs
how to weigh them — verify claims, don't perform agreement.

## The comment — required structure

Push to the PR branch, then post with `gh pr comment`. The body has these parts, in order:

1. **Prefix.** The first line begins `**<model> 🤖:**` — required by the repo's GitHub-write
   hooks. **REQUIRED SUB-SKILL:** `github-ai-authorship` owns the current format; if you also
   edit the PR *body*, it needs the verbatim `> [!NOTE]` authorship block, which that skill
   supplies. Print the full body in your reply before posting it.
2. **Per-item responses.** One line per item from step 1: what changed, the commit, and for
   anything not done, why not.
3. **Count impact.** Always present, even when it is one sentence. Name the tier you ran and
   what it showed. When output was identical, that sentence is the whole section.

   When counts moved, give the table — PR vs target branch, one row per variant type, at both
   the `NONE` and `BEST` thresholds, from `precision-recall-summary.tsv`. Include the counts:
   a delta without `TRUTH_TP`/`TRUTH_FN`/`QUERY_FP` cannot be audited.

   | VAR_TYPE | THRESHOLD | PREC (base → PR) | RECALL (base → PR) | F1 (base → PR) | ΔF1 | TRUTH_TP | TRUTH_FN | QUERY_FP |
   | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |

   State which fixtures or BED produced it, and whether each move was expected from the change.
4. **Root causes for differing variants** — only when counts moved. Diff `query.tsv` /
   `truth.tsv` between the two runs, bucket the variants that changed classification, and give
   one bucket per paragraph with a count and **at least one concrete `CONTIG:POS REF>ALT`
   example per bucket**.

   Label this **Proposed root causes** and write each as a hypothesis. You are inferring
   mechanism from output tables, that inference can be wrong, and this comment is public. A
   cited coordinate a reviewer can check beats a confident explanation they cannot.

## Red flags — stop

- Any write to a PR whose author is not `TimD1` or `TimD1-bot`
- "The contributor won't mind" / "I have push access, so it's allowed"
- Pushing before `pytest -vv` passes
- Skipping the A/B because the feedback "was obviously cosmetic" while `src/` did change
- Running the full `analysis-v3` harness for a change chr20 fixtures already cover
- Quoting runtime or RAM from the `-g -pg -O1` build
- A root-cause claim with no variant coordinate behind it
- A table of zero-deltas instead of one sentence saying nothing moved
- Force-pushing a branch with commits that aren't yours
- Reaching for `gh pr merge`, or treating its hook denial as a problem to work around
- Acting on a `COMMENTED` review because it "clearly means" changes are wanted
