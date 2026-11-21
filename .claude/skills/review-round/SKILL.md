---
name: review-round
description: Review an implementer's branch of qespace_cef against PLAN.md, CLAUDE.md and the accepted CONVs, independently of the implementer's own report, and report to the maintainer with numbers. Use when an implementer started by implementation-round reports done, or when the maintainer asks to review a branch.
---

# Review round

You are the advisor (CLAUDE.md section 10). You review; you do not fix. Every
correction goes back to the implementer that wrote the branch (SendMessage to
its session), or to the maintainer if the plan itself was wrong. Your
checkout stays on `main` and unchanged.

## 1. Read what was done

In the package root, read-only:

- `git log --format='--- %h %s%n%b' main..<branch>`: one logical change per
  commit, Conventional Commits (CLAUDE.md 8), `!` and `BREAKING CHANGE:`
  where the API or a computed result changes (8.2), the body saying why,
  and each commit green on its own (8.3).
- `git diff --stat main...<branch>`: only the files PLAN.md, *Now*, lists.
  Any other file is a finding.
- `git diff main...<branch>`, read in full. Check against:
  - PLAN.md *Now*: each step is done, and nothing beyond it;
  - the accepted CONVs the step touches: the code follows them; every newly
    accepted CONV has tests named `test_conv_xxxx_*` (CLAUDE.md 6.6);
  - CLAUDE.md 5 (code) and 6 (tests): tests assert results with expected
    values that do not come from the code under test; no test is loosened
    or deleted to pass (6.4); a fix comes with a test that fails without it
    (6.3);
  - CLAUDE.md 7: ARCHITECTURE.md, README.md, CHANGELOG.md and AUDIT.md
    follow the change.

## 2. Run the checks yourself

Wait until the implementer is idle, so no two uv installs run at once
(CLAUDE.md 4.4). Export the branch to the scratchpad, which leaves the
repository untouched:

    git archive <branch> | tar -x -C <scratchpad>/rev

Then, in that folder, one after another:

1. `uv run --no-project --isolated --with-editable . --with pytest pytest -q -p no:cacheprovider`
2. the doctests: the same with `--doctest-modules qespace_cef --ignore=qespace_cef/symbolic`
3. `uv run --no-project --isolated --with ruff ruff check . --exclude examples`
   and `ruff format --check . --exclude examples`
4. the step-specific checks of PLAN.md against your own data, e.g. the
   snapshot implementation-round saved from `main` (bit-identity of a
   refactor), not the implementer's.

Compare the numbers with `main` (passed, xfailed): every change must be
explained by the plan, e.g. new tests, an xfail that turns green.

## 3. Report

To the maintainer, in this order:

- one sentence: clean, or what is wrong;
- a table: commits, scope, the step-specific checks, pytest, doctests, ruff,
  each with its numbers;
- what the branch changes, in a few bullets;
- findings, each with the rule it breaks; observations that break no rule
  labelled as such;
- what was not checked;
- whether `main` can fast-forward (`git merge-base --is-ancestor main
  <branch>`) and the command to merge:

      git merge --ff-only <branch>

The maintainer merges, or tells you to. After the merge, replace PLAN.md's
*Now* with the next step (or say there is none) and remove what the branch
fixed from AUDIT.md, if the implementer has not.
