# CircleCI to GitHub Actions migration guide

Notes from migrating this repo's CI off CircleCI onto GitHub Actions. Split into two parts on purpose:

- **Part 1** is a standalone guide to setting up GitHub Actions CI — useful even in a brand-new repo with no CircleCI history at all.
- **Part 2** is specific to retiring an existing CircleCI integration once Part 1's workflow is proven. Skip it entirely if there's no CircleCI to remove.

If you're migrating an existing repo, do Part 1 first, get it green, *then* do Part 2 — don't remove the old CI before the new one is proven and required checks are updated.

## Part 1: Setting up GitHub Actions CI

### 1. Gather requirements
- If migrating from CircleCI: read `.circleci/config.yml` for the base image/executor, install steps, lint/test commands, resource class (memory/CPU), and any orbs (e.g. the Slack orb).
- Also check CircleCI's project dashboard (web UI) for things not in the YAML: **scheduled pipelines** (nightly/cron triggers configured outside `config.yml`), contexts, and org-level secrets — these won't show up anywhere in git history.
- If this is a new repo: just decide what should run (lint, tests, build) and on what triggers (push, PR, nightly).

### 2. Write the workflow file
- Add `.github/workflows/ci.yml` (or similar) on a feature branch.
- Prefer official/maintained actions over reimplementing shell steps (`actions/checkout`, `conda-incubator/setup-miniconda`, etc.) — but pin versions and expect to need to bump them later (see "Lessons learned" below).
- Set triggers: `on: push` / `on: pull_request`, plus `schedule:` if you want nightly/cron builds (see step 5).
- If adding a Slack (or similar) notification step, look up the current action's actual input names rather than assuming — see "Lessons learned."

### 3. Prove it on a branch first
- Push the workflow file on a feature branch and open a PR. GitHub Actions runs workflow files from PR branches automatically — no need to merge first.
- Confirm the run is actually green before doing anything else (updating required checks, removing old CI, etc.).

### 4. Set required status checks / branch protection
- In `Settings → Branches` (or `Settings → Rules → Rulesets`), edit each protected branch's rule and add the new job under "Require status checks to pass before merging."
- The job only becomes selectable **after it has run at least once** against that branch (via a PR targeting it) — so this step has to come after step 3, not before.

### 5. Nightly/scheduled builds, if needed
- GitHub Actions only evaluates `schedule:` triggers from the workflow file **on the repository's default branch** — a `schedule:` block committed only to a non-default branch will never fire on its own.
- To get a nightly run against a non-default branch (e.g. `dev` when `main` is default): add a **separate workflow file that lives on the default branch**, triggered by `schedule:` (+ `workflow_dispatch:` for manual reruns), whose job explicitly does `actions/checkout@vN` with `ref: <other-branch>`.
- This duplicates steps between the two workflow files — acceptable for a single job; consider a reusable workflow (`workflow_call`) if you have several.

### 6. Add the status badge
- README badge: `https://github.com/<org>/<repo>/actions/workflows/<file>.yml/badge.svg?branch=<branch>`, linking to `https://github.com/<org>/<repo>/actions/workflows/<file>.yml`.

### Web UI actions checklist (Part 1)
- [ ] Create a Slack **incoming webhook** app (api.slack.com/apps → Incoming Webhooks) if you want Slack notifications — this is a distinct mechanism from Slack API tokens/orbs used elsewhere.
- [ ] Add the webhook URL (or any other secret) as a GitHub repo secret (`Settings → Secrets and variables → Actions`).
- [ ] Update branch protection / ruleset required status checks for every protected branch, after the new check has run at least once against it.
- [ ] Open/merge PRs from the web UI if you don't have an authenticated `gh` CLI session available.

## Part 2: Removing a CircleCI integration

Only do this once Part 1's workflow is proven green and required checks have been updated to point at it (Part 1, steps 3-4).

### 1. Update required checks before touching CircleCI
- Confirm (again) that the new GitHub Actions check is required and the old CircleCI check is not, on every branch CircleCI used to protect.
- **Order matters**: if you disable/remove CircleCI before this, the old required check will show "Expected — waiting for status to be reported" forever and permanently block merges.

### 2. Decommission CircleCI itself
- Stop/unfollow the project in CircleCI's dashboard (or delete it) to kill its webhook — deleting just the config file still leaves it trying (and failing) to build, posting error statuses on every push.
- Then remove `.circleci/config.yml` from the repo.

### 3. Update docs and badges
- Swap any CircleCI status badge in the README for the GitHub Actions one (Part 1, step 6).
- Grep the repo for the old CI's name (e.g. `circleci`) — badges are easy to remember, prose references in contributing/testing docs are not.

### 4. Handle divergent branches deliberately
- If `main` is behind your active development branch, do **not** merge a dev-based feature branch directly into `main` — it'll drag in every unrelated commit dev has that main doesn't.
- Instead: merge into the active branch normally, then `git cherry-pick` the specific CI-related commit(s) onto a fresh branch based on `main`, and open a separate PR.
- Check `git merge-base HEAD origin/main` vs `origin/dev` (or whatever branches are relevant) up front — don't assume the branch you're on is based on the branch you think it is.

### Web UI actions checklist (Part 2)
- [ ] Stop/unfollow the project in CircleCI's dashboard (or delete it) to remove its webhook trigger.
- [ ] Double check a PR/branch is *actually* merged (check the target branch's latest commit) before deleting it — a branch can look "handled" in the UI while the PR is still open.

## Lessons learned in this migration

- **Third-party actions drift out of Node-runtime compatibility.** GitHub deprecated Node 20 runners in favor of Node 24 mid-2026; actions pinned to older major versions (`actions/checkout@v4`, `conda-incubator/setup-miniconda@v3`, `slackapi/slack-github-action@v1.x`) broke or warned until bumped to the versions that ship Node 24 support (`@v5`, `@v4`, `@v3.x` respectively). Expect to need to re-check/bump action versions periodically, not just at initial setup.
- **Don't trust a web-fetched summary of an action's inputs — check the actual error or `action.yml`.** `slackapi/slack-github-action` v3's real input names are `webhook` and `webhook-type: incoming-webhook`; an initial guess of `webhook-url` (from a web search summary) failed at runtime with a clear list of valid inputs in the error log. When an action rejects an input, the error message itself is the most reliable source.
- **`schedule:` only runs from the default branch's copy of the workflow file** (Part 1, step 5). Easy to add a schedule block everywhere for consistency and assume it "just works" on every branch; it doesn't.
- **Required status checks must be updated in the right order** relative to decommissioning the old CI (Part 2, steps 1-2) — otherwise merges get permanently blocked on a check that will never report again.
- **Verify merge state directly (`git log <branch>`, compare commit hashes) rather than relying on "I merged it" / "I deleted it" from memory** — during this migration a branch was assumed merged-and-safe-to-delete when the PR was actually still open, and separately a "deleted" branch turned out to still exist on the remote. A quick `git fetch` + `git log --oneline origin/<branch>` (using an explicit refspec if a plain `fetch` seems to return stale results) settles it in seconds.
- **When cherry-picking a fix across diverged branches, re-verify the target file's exact wording** before reusing a diff — a commit that removes redundant text on one branch may not apply cleanly (or correctly) to another branch's slightly different copy of the same file.

## Reusable checklists

**Setting up GitHub Actions (new repo, or alongside existing CI):**
- [ ] Gather requirements (existing CI config and its web dashboard settings, or a fresh decision on what should run)
- [ ] Write workflow file on a feature branch
- [ ] Get it green on a PR
- [ ] Set required status checks on protected branches
- [ ] Add nightly/scheduled builds if needed, remembering the default-branch-only rule
- [ ] Add README badge

**Removing CircleCI (only after the above is done and required):**
- [ ] Re-confirm required checks point at GitHub Actions, not CircleCI
- [ ] Stop/unfollow the CircleCI project
- [ ] Delete `.circleci/config.yml`
- [ ] Update README badge(s) and any prose docs mentioning CircleCI
- [ ] Handle divergent branches with targeted cherry-picks, not blanket merges
- [ ] Re-verify actual merge/branch state via `git log`/`git fetch` before deleting anything
