---
name: respond-to-vcfdist-pr-feedback
description: Use when a vcfdist pull request authored by TimD1 or TimD1-bot has received blocking review feedback (a review in CHANGES_REQUESTED state, or a maintainer asking for changes) and that feedback needs to be addressed, verified, and pushed. Does not apply to pull requests from any other author, including outside contributors. Triggers include "address the PR feedback", "the reviewer requested changes", "respond to that review".
---

# Respond to blocking PR feedback

Fix, verify locally, re-review, measure whether counts moved, push, comment, hand back.

## Scope gate — run this first

One query answers all three gates:

```bash
gh pr view <N> --json author,assignees,reviews,comments,baseRefName,headRefName,url
```

**Gate 1 — ownership.** `.author.login` must be exactly `TimD1` or `TimD1-bot`. Anything else:
stop, report the PR as out of scope, and touch nothing — no checkout, no comment. vcfdist takes
outside contributions, and "allow edits by maintainers" means you *can* push to a contributor's
fork branch. That is precisely why this is an explicit gate and not left to permissions.

**Gate 2 — a request to act.** One of:

1. **Rejected** — a `CHANGES_REQUESTED` review authored by `TimD1`.
2. **Handed off** — assigned to `TimD1-bot`, and `TimD1` is the most recent voice, counting
   issue comments **and** reviews. A `COMMENTED` review never appears in `.comments`, so
   checking comments alone reads the bot as the last voice and misses the handoff.

A `COMMENTED` review with no assignment is **not** a trigger: assignment is the signal, the
comment is the content. A rejection authored by `TimD1-bot` triggers nothing — otherwise
reviewing a PR would dispatch a run to fix it and the two would ping-pong.

**Gate 3 — not already handled.** An existing comment from you post-dating the trigger means
stop; re-running duplicates commits and comments.

## Hard rules

- **Never write to a PR authored by anyone but `TimD1` or `TimD1-bot`** — no commits, no
  pushes, no comments. Being *able* to push is not authorization.
- **Never merge.** Tim squash-merges every branch himself. `gh pr merge` and
  `gh pr review --approve` are hook-denied; the denial is correct, not an obstacle.
- **Never force-push** unless every commit is yours:
  `git log --format='%an' origin/<base>..HEAD | sort -u`.
- **Never work in the main checkout.** Resolve the PR's worktree first — see step 0.
- **Push only the exact commit you tested** — see below.

## Sequence

**0. Get into the PR's worktree.** You are almost certainly starting in the main checkout,
which sits on whatever branch was last used and often holds unrelated uncommitted work. Resolve
the right tree before touching anything:

```bash
BR=$(gh pr view <N> --repo TimD1/vcfdist --json headRefName --jq .headRefName)
WT=$(git worktree list --porcelain \
     | awk -v b="refs/heads/$BR" '/^worktree /{p=$2} $0=="branch "b{print p}')
if [ -z "$WT" ]; then WT=".claude/worktrees/$BR"; git worktree add "$WT" "$BR"; fi
cd "$WT"
git rev-parse --abbrev-ref HEAD   # must equal $BR
git status --porcelain            # must be empty
```

**A dirty worktree at this point is a stop, not a cleanup.** Those changes are someone else's
work in progress — do not commit, stash, or discard them. Report what you found, leave the PR
assigned to `TimD1-bot`, and end the run.

**1. Restate the feedback** as discrete items before editing anything. Items you decide not to
act on are listed with a reason; silent omission forces a re-review from scratch.

**2. Fix.** Follow CLAUDE.md's documentation conventions: `.h` gets `@brief` only, `.cpp` gets
full parameter/return/throws blocks, separators are 100 characters.

**3. Commit, build, test.** Commit first, so the tested state has an identity:

```bash
git add <explicit paths> && git commit   # never `git add -A` or `-u`
cd src && make
cd ../tests && pytest -vv
git rev-parse HEAD                       # ← the TESTED SHA
```

**Every test must pass.** On any failure: report it with its output, leave the PR assigned to
`TimD1-bot`, and push nothing. This holds when the failure looks unrelated to your change or
reproduces on the target branch — triaging pre-existing breakage is the maintainer's call.

### The tested tree is the pushed tree

**Code never changes after the tests pass. If it does, you are back at step 2** — not patching
on the way out. Review findings, A/B fallout, anything noticed while writing the comment: all
loop back, and the tested SHA is replaced rather than amended.

Check it rather than trusting memory. Immediately before pushing:

```bash
git rev-parse HEAD        # must equal the TESTED SHA
git status --porcelain    # must be empty
```

Either check failing means do not push. The quiet version of this failure is a stale binary:
editing `src/` and running the A/B without `make` compares the old build, so the comment
reports numbers for code that is not on the branch.

## Whether to re-evaluate counts

```bash
git diff -w origin/<base>...HEAD -- src/   # -w: ignore whitespace-only churn
```

| Diff touches | Action |
| :-- | :-- |
| Nothing under `src/` | No evaluation; say so in the comment |
| `src/`, but only comments, whitespace, or renames | Run the A/B anyway, expecting identical output — a rename that moves a number is the bug worth catching |
| Executable logic in `src/` | Run the A/B and report the table |

Commands, baseline caching, and the escalation tier for accuracy-critical changes are in
`reference.md`. **Identical output** → count impact is one sentence; a zero-delta table is
noise. **Different** → table plus root causes.

**Never quote runtime or RAM.** `src/Makefile` defaults to `CXXFLAGS = -g -pg -O1` with `-O3`
commented out; accuracy is unaffected, timings are meaningless.

## Independent review

Fresh eyes on the diff before pushing — your own second read is not review. **REQUIRED
SUB-SKILL:** whichever reviewer-persona review skill the environment provides, plus
`reviewing-bioinformatics-code` for changes to alignment, clustering, or variant logic.
**REQUIRED SUB-SKILL:** `superpowers:receiving-code-review` — verify claims, don't perform
agreement. Record why for findings you reject; **findings you accept send you back to step 2.**

## Push, comment, hand back

Run both tested-tree checks, push, then `gh pr comment`. The template is in `reference.md`;
print the body in your reply before posting.

Two parts of it are decisions rather than formatting: count impact is **always present**,
collapsing to one sentence when nothing moved; and root causes are **proposed** hypotheses,
each carrying a `CONTIG:POS REF>ALT` a reviewer can check. You are inferring mechanism from
output tables, that inference can be wrong, and the comment is public.

Then hand it back — required, and it clears the handoff trigger:

```bash
gh pr edit <N> --add-assignee TimD1 --remove-assignee TimD1-bot
```

Assignment carries no body text, so it needs no authorship note. Do this even when nothing
changed; a PR left on the bot reads as still in progress. **If the run fails or stops early**,
leave it assigned to `TimD1-bot` and say what stopped you — that is what makes a stall visible
instead of silently abandoned.

## Red flags — stop

- "The contributor won't mind" / "I have push access, so it's allowed"
- Running `make`, `pytest`, or `git commit` before step 0 has confirmed the branch
- Committing, stashing, or discarding uncommitted work you did not create
- "This edit is too small to re-test" / "I'll just fix it on the way out"
- "The feedback was obviously cosmetic" — while `src/` did change
- "This `COMMENTED` review clearly means changes are wanted" — with no assignment
- Continuing past a test failure because it looks unrelated or pre-existing
- Running the A/B after editing `src/` without rebuilding
- Running the `analysis-v3` harness for a change the chr20 fixtures already cover
- Quoting runtime or RAM from the `-g -pg -O1` build
- A root-cause claim with no variant coordinate behind it
- A zero-delta table instead of one sentence saying nothing moved
- Finishing the work but leaving the PR assigned to `TimD1-bot`
- Treating a `gh pr merge` denial as an obstacle to route around
