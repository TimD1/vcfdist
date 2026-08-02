---
name: respond-to-vcfdist-pr-feedback
description: Use when a vcfdist pull request authored by TimD1 or TimD1-bot has been assigned to TimD1-bot, which is the hand-off signal that its review feedback should be addressed, verified, and pushed. Does not apply to pull requests from any other author, nor to PRs merely carrying review comments or a CHANGES_REQUESTED review without that assignment. Triggers include "address the PR feedback", "respond to that review".
---

# Respond to blocking PR feedback

Fix, verify locally, re-review, measure whether counts moved, push, comment, hand back.

## Scope gate — run this first

One query answers both gates:

```bash
gh pr view <N> --json author,assignees,reviews,comments,baseRefName,headRefName,url
```

**Gate 1 — ownership.** `.author.login` must be exactly `TimD1` or `TimD1-bot`. Anything else:
stop, report the PR as out of scope, and touch nothing — no checkout, no comment. vcfdist takes
outside contributions, and "allow edits by maintainers" means you *can* push to a contributor's
fork branch. That is precisely why this is an explicit gate and not left to permissions.

**Gate 2 — assigned to `TimD1-bot`.** That assignment is the *only* trigger. A
`CHANGES_REQUESTED` review is not one. Nor is a comment, however blocking its content.

**Expect the assignee list to be empty.** When the watcher dispatched you it already un-assigned
the bot — that is how a failed run is stopped from looping. An empty list is therefore the
normal case, not a failed gate. Confirm the assignment happened and was recent:

```bash
gh api repos/TimD1/vcfdist/issues/<N>/timeline --paginate \
  --jq '.[] | select(.event=="assigned" and .assignee.login=="TimD1-bot")
            | "\(.created_at) by \(.actor.login)"' | tail -1
```

Gate 2 passes if `TimD1` assigned the bot within roughly the last hour. If the bot is *still*
an assignee, you were invoked by hand — un-assign it at step 0a as usual. If no such event
exists at all, stop: nothing triggered this.

## Hard rules

- **Never write to a PR authored by anyone but `TimD1` or `TimD1-bot`** — no commits, no
  pushes, no comments. Being *able* to push is not authorization.
- **Never merge.** Tim squash-merges every branch himself. `gh pr merge` and
  `gh pr review --approve` are hook-denied; the denial is correct, not an obstacle.
- **Force-push only after the rebase, only with `--force-with-lease`, and only when every
  commit is yours**: `git log --format='%an' origin/<base>..HEAD | sort -u`. If any other name
  appears, stop — do not rebase, do not push. Plain `--force` is never correct here; the lease
  is what catches the remote having moved since you fetched.
- **Never work in the main checkout.** Resolve the PR's worktree first — see step 0.
- **Push only the exact commit you tested** — see below.

## Sequence

**0a. Claim the trigger.** Before anything else, consume it:

```bash
gh pr edit <N> --repo TimD1/vcfdist --remove-assignee TimD1-bot
```

The watcher normally does this before dispatching; do it yourself when invoked by hand. Claiming
first is what stops a failed run from being re-dispatched forever — the trigger is spent whether
or not the run succeeds, and **re-assigning is how a human retries**. Never re-add the bot as an
assignee to "keep your place".

**0b. Get into the PR's worktree.** You are almost certainly starting in the main checkout,
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

**Name a worktree you create after the branch, exactly** — `.claude/worktrees/$BR`, never a
shortened alias. The directory name is the only thing that tells a `git worktree list` reader
which tree belongs to which PR. Older trees may carry abbreviated names (`single-pass` for
`91_D12_td_single-pass-prec-recall`); use one where it already exists rather than moving a
directory someone may have a shell open in, and note the mismatch in the comment.

**Only the PR's worktree has to be clean.** The main checkout you started in is very often
dirty with unrelated work; that is normal and blocks nothing, because you never operate there.
Run the status check *after* `cd "$WT"`, never before.

**A dirty PR worktree is a stop, not a cleanup.** It means someone is mid-edit on this very
branch. Do not commit, stash, or discard their changes, and do not sidestep with a second
worktree: git refuses the same branch twice, and the `--force` and detached-`origin/<branch>`
routes both push over work you cannot see. Report what you found, assign `TimD1` back, and end
the run.

**0c. Rebase onto `dev`.** Always, before any new work:

```bash
git fetch origin dev
git log --oneline origin/<base>..HEAD    # anything here is unpushed pre-existing work
git rebase origin/dev
```

Two checks before you start editing. If the rebase reports **conflicts**, `git rebase --abort`,
stop, and report them — resolving someone else's conflicts unattended is beyond this run's
remit. And if that `git log` shows commits **not yet on the remote**, your eventual push will
carry them along: name them in the comment rather than shipping them silently.

Rebasing rewrites the branch, so the final push becomes `--force-with-lease`. That is expected
and does not contradict the never-force-push rule, which is about not overwriting *other
people's* commits — confirm that first, per the hard rules.

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

**Every test must pass.** On any failure: push nothing, report it with its output, and hand the
PR back to `TimD1` with that explanation. This holds when the failure looks unrelated to your change or
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

## Close out the review threads

After pushing, give every inline thread a disposition. List them:

```bash
gh api graphql -f query='
{ repository(owner:"TimD1",name:"vcfdist"){ pullRequest(number:<N>){
    reviewThreads(first:100){ nodes{ id isResolved path line
      comments(first:10){ nodes{ author{login} body } } } } } } }'
```

**Fully addressed** → resolve it; no reply needed, the diff is the answer:

```bash
gh api graphql -f query='mutation($t:ID!){ resolveReviewThread(input:{threadId:$t}){ thread{ isResolved } } }' -f t=<threadId>
```

**Not fully addressed** — partially done, deferred, or declined → reply, and leave it
**unresolved**:

```bash
gh api graphql -f query='mutation($t:ID!,$b:String!){ addPullRequestReviewThreadReply(input:{pullRequestReviewThreadId:$t,body:$b}){ comment{ url } } }' -f t=<threadId> -f b='**<model> 🤖:** ...'
```

One or two sentences: what you did, or why you did not. The summary comment carries detail; the
reply exists so someone scanning the diff sees the disposition in place. Replies take the
`**<model> 🤖:**` prefix like any other comment.

**Never resolve a thread you did not address.** Unresolved is the reviewer's queue — resolving
to tidy the PR erases the only record that something was left undone.

## Push, comment, hand back

Run both tested-tree checks, then push with `--force-with-lease` (step 0c rebased the branch),
then `gh pr comment`. The template is in `reference.md`;
print the body in your reply before posting.

Two parts of it are decisions rather than formatting: count impact is **always present**,
collapsing to one sentence when nothing moved; and root causes are **proposed** hypotheses,
each carrying a `CONTIG:POS REF>ALT` a reviewer can check. You are inferring mechanism from
output tables, that inference can be wrong, and the comment is public.

Then hand it back — required, and the last action of every run:

```bash
gh pr edit <N> --add-assignee TimD1
```

Assignment carries no body text, so it needs no authorship note. Do this even when nothing
changed, and **also when the run fails or stops early** — in that case say plainly what stopped
you. The bot un-assigned itself at step 0a, so a run that dies without reaching here leaves the
PR with no assignee at all: that silence is the failure signal, and it cannot re-trigger.

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
- Ending a run without assigning `TimD1` back — including when it failed
- Re-adding `TimD1-bot` as an assignee for any reason
- Resolving a review thread you did not fully address
- Starting work without rebasing onto `dev`, or resolving rebase conflicts unattended
- Reaching for plain `git push --force` instead of `--force-with-lease`
- Creating a worktree named anything other than the branch
- Treating a `gh pr merge` denial as an obstacle to route around
