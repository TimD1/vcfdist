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

**Gate 2 — an actual request to act.** Exactly one of these must hold:

1. **Rejected** — a review in `CHANGES_REQUESTED` state authored by `TimD1`.
2. **Handed off** — the PR is assigned to `TimD1-bot`, *and* `TimD1` is the most recent voice
   on it. "Most recent voice" spans issue comments **and** reviews: a `COMMENTED` review
   carries real feedback but never appears in `.comments`, so checking comments alone reads
   the bot as the last voice and misses the handoff.

If neither holds, **stop and report that.** In particular, a `COMMENTED` review from `TimD1`
with no assignment is *not* a trigger — reading intent out of prose is what produces
force-pushes over nitpicks. Assignment is the signal; the comment is the content.

A rejection authored by `TimD1-bot` never triggers anything. Otherwise reviewing a PR would
dispatch a run to fix it, and the two would ping-pong.

**Gate 3 — not already handled.** Look for an existing comment from you that post-dates the
triggering review or comment. If one exists, stop; re-running duplicates commits and comments.

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
- **Push only the exact commit you tested.** `pytest -vv` must pass, and `HEAD` at push time
  must be the SHA that passed, with a clean worktree. Any edit after the test run — however
  small, whatever its source — sends you back to step 2 to fix and re-test.

## Sequence

**1. Read the feedback and restate it.** Enumerate each requested change as a discrete item
before editing anything. Items you decide not to act on must be listed with a reason — silent
omission is what makes a reviewer re-review from scratch.

**2. Fix.** Follow the repo's documentation conventions from CLAUDE.md (`.h` files get
`@brief` only; `.cpp` files get full parameter/return/throws blocks; 100-character separators).

**3. Commit, build, test.** Commit first, so the tested state has an identity:

```bash
git add <explicit paths> && git commit      # never `git add -A` or `-u`
cd src && make                              # produces src/vcfdist
cd ../tests && pytest -vv                   # unit (GoogleTest) + integration (pytest-workflow)
git rev-parse HEAD                          # ← the TESTED SHA. Write it down.
```

**Every test must pass.** If any fails, stop: report the failing test with its output, leave
the PR assigned to `TimD1-bot`, and push nothing. This holds even when the failure looks
unrelated to your change or reproduces on the target branch — triaging pre-existing breakage is
a decision for the maintainer, not a reason to continue.

(`tests/README.md` claims several tests need a genome-wide FASTA in `data/refs/`, which ships
**empty** — but the committed integration test uses the chr20 reference in
`tests/integration/data/`. If a test errors for a missing reference, that stale instruction is
the likely cause.)

### The tested tree is the pushed tree

**Code never changes after the tests pass. If it does, you are back at step 2.**

Every later stage — independent review, the A/B count comparison, anything you notice while
writing the comment — can produce an edit. When one does, that is not a patch to apply on the
way out: return to step 2, fix, and run step 3 again from the top. The tested SHA is replaced,
not amended.

This is mechanically checkable, so check it rather than trusting your memory of what you
touched. Immediately before pushing:

```bash
git rev-parse HEAD        # must equal the TESTED SHA from step 3
git status --porcelain    # must be empty — no unstaged or untracked changes
```

If either check fails, **do not push**. Go to step 2 and start the loop again.

A stale binary is the quiet version of this failure: editing `src/` and running the A/B without
`make` compares the *old* build, so the numbers in your comment describe code that is not on
the branch. Re-running step 3 rebuilds, which is why the fix is to loop rather than to patch.

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

The A/B builds the PR and target-branch binaries and runs both over the committed chr20
fixtures in `tests/integration/data/`. Commands, baseline caching, and the escalation tier for
accuracy-critical changes: see `reference.md` in this directory.

Two results, two obligations. **Identical** → the count-impact section is one sentence; a table
of zero-deltas is noise. **Different** → report the table and the root-cause buckets.

**Never quote runtime or RAM** from these runs. `src/Makefile` defaults to
`CXXFLAGS = -g -pg -O1` with `-O3` commented out — accuracy is unaffected, timings are
meaningless.

## Independent review

Before pushing, get fresh eyes on the diff — your own second read is not review.
**REQUIRED SUB-SKILL:** use whichever reviewer-persona review skill is available in the
environment, matched to whoever reviews this repo; for changes to alignment, clustering, or
variant logic also use `reviewing-bioinformatics-code`. Record why for findings you reject.
**REQUIRED SUB-SKILL:** `superpowers:receiving-code-review` governs how to weigh them — verify
claims, don't perform agreement.

**Findings you accept send you back to step 2.** This stage exists to produce edits, so it is
the most likely origin of an untested tree. Applying a review finding here and pushing without
re-running step 3 is the exact failure the tested-tree rule prevents.

## Push and comment

Run the two tested-tree checks above. Only then push to the PR branch, and post with
`gh pr comment`.

The body has four required parts — prefix, per-item responses, count impact, root causes — with
the full template, table columns, and authorship requirements in `reference.md`. Print the body
in your reply before posting it.

Two properties of that template are decisions rather than formatting, so they are stated here
too: the count-impact section is **always present**, collapsing to one sentence when nothing
moved; and root causes are labelled **proposed**, each carrying a concrete
`CONTIG:POS REF>ALT` a reviewer can check. You are inferring mechanism from output tables, that
inference can be wrong, and the comment is public.

## Hand the PR back

**Last step, and not optional.** Once the push has landed and the comment is posted:

```bash
gh pr edit <N> --add-assignee TimD1 --remove-assignee TimD1-bot
```

This is the visible signal that the ball is back in the reviewer's court — an assignment they
can see in their PR list, rather than a comment they have to notice. It also clears the handoff
trigger: while the PR stays assigned to `TimD1-bot`, it remains a candidate for another run.

Assignment carries no body text, so it needs no authorship note — that requirement applies to
`--body` edits. Do this even when the run ended without changes (identical A/B output, or
nothing to fix): reassign and say so in the comment. A PR left assigned to the bot reads as
still in progress.

**If the run ends early or fails**, leave the assignment on `TimD1-bot` and say what stopped
you. That is what makes a stalled run visible instead of silently abandoned.

## Red flags — stop

- Any write to a PR whose author is not `TimD1` or `TimD1-bot`
- "The contributor won't mind" / "I have push access, so it's allowed"
- Acting on a `COMMENTED` review with no assignment to `TimD1-bot`
- Finishing the work but leaving the PR assigned to `TimD1-bot`
- Pushing when `HEAD` is not the tested SHA, or `git status --porcelain` is non-empty
- "This edit is too small to re-test" / "I'll just fix it on the way out"
- Applying a review finding and pushing without returning to step 2
- Running the A/B after editing `src/` without rebuilding — that compares the old binary
- Continuing past a test failure because it looks unrelated or pre-existing
- Skipping the A/B because the feedback "was obviously cosmetic" while `src/` did change
- Running the full `analysis-v3` harness for a change chr20 fixtures already cover
- Quoting runtime or RAM from the `-g -pg -O1` build
- A root-cause claim with no variant coordinate behind it
- A table of zero-deltas instead of one sentence saying nothing moved
- Force-pushing a branch with commits that aren't yours
- Reaching for `gh pr merge`, or treating its hook denial as a problem to work around
- Acting on a `COMMENTED` review because it "clearly means" changes are wanted
