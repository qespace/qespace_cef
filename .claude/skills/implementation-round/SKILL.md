---
name: implementation-round
description: Start an implementation round of qespace_cef from PLAN.md. Checks that the plan is ready, starts a fresh implementer session in its own worktree with the standard brief, opens it in split view, and hands over to review-round when it reports. Use only when the maintainer says to start the round ("start", "commit and start", "start a worker"), never on your own initiative.
---

# Implementation round

You are the advisor (CLAUDE.md section 10). This skill starts one implementer
for the step in PLAN.md, section *Now*, and nothing else. You do not edit
package code yourself.

## 1. Preflight

Run these read-only checks in the package root and stop at the first failure,
telling the maintainer what is missing:

1. `git status --short --branch`: on `main`, tree clean. PLAN.md, and every
   decision it relies on, is committed.
2. PLAN.md, *Now*, has the four parts the brief relies on: a branch name, a
   **Files** list, numbered **Steps** (commits) and **Tests**. A step that
   needs a decision still marked `OPEN:` or a CONV still in draft is not
   ready, nor is a plan that contradicts an accepted ADR or CONV (check the
   ones it cites and the ADRs on the API it touches): report it instead of
   starting.
3. No branch for this step exists yet (`git branch --list`), and no
   implementer you started is still running (the linked sessions).
4. If the plan promises an unchanged result (a refactor, "no computed result
   changes"), save a snapshot of those results now, from `main`, to the
   scratchpad, so review-round can compare against your own copy.

## 2. Start the implementer

By default, send the brief with SendMessage to the existing worker "CEF
implementer" (find it with `list_sessions` linked; it must be idle, its
worktree clean), telling it to cut the plan's branch from `main` in its
worktree; skip step 3 if its split is already open. Start a new one only if
it is gone or the maintainer asks: before the call, write a bold title line
and one sentence on what it will do. Then call `start_session` with:

- `initiation`: `"user_asked"`; `context`: `"fresh"` (a `"fork"` only when the
  maintainer asks for it, or the step depends on discussion not written into
  PLAN.md or `decisions/`: say so);
- `use_worktree`: `true`; `branch_name`: the slug of the branch in PLAN.md;
- `name`: a short handle; `title`: "Implement <branch> (<findings, CONVs>)";
- `wake_when_done`: `true`;
- `prompt`: the brief below, with the placeholders filled in;
- `background`: one paragraph on the package, the advisor/implementer split,
  and why this step exists (the finding and the CONV it serves).

The brief:

> Implement the section "Now: `<branch>`" of `PLAN.md`, following `CLAUDE.md`
> (read it first, section 10 in particular, then the decisions in
> `decisions/` that PLAN.md cites, then `ARCHITECTURE.md` and `DEV.md`).
>
> You are the implementer, in your own git worktree of `qespace_cef`, on your
> own branch (its name, given by the app, stands in for `<branch>`; don't
> rename it). The maintainer has approved this plan: commit on your branch
> without asking (CLAUDE.md 1.3). Never commit on `main`, never merge, never
> push, never rewrite history.
>
> Scope: touch only the files PLAN.md lists: `<files>`. If anything else seems
> needed, stop and report instead. Decide no convention and no `OPEN:` item.
> `<step-specific checks, e.g. the bit-identity snapshot>`
>
> Environment (CLAUDE.md 4): `uv` only. Don't create or change a shared
> `.venv`: use `uv run --no-project --isolated --with-editable . --with
> pytest pytest -q` and the same with `ruff`. Run uv installs one after
> another. Give `uv run python` a
> script file or `-c`, never an empty `python -` on stdin.
>
> When done, stop and report (CLAUDE.md 1.3.5): `git log --oneline main..HEAD`,
> the branch name, pytest / ruff numbers, the results of the
> step-specific checks, what was not checked, and the command for the diff.
> Then end your turn; the advisor reviews the branch.

## 3. Show it

Right after the session starts, call `open_session_in` with its session id and
`target: "split"`. If this session is not on screen and the split does not
happen, tell the maintainer where to find the worker (sidebar, right-click,
Open in → Split view).

## 4. While it runs

End your turn with one line saying you will review when it reports. Do not
poll. Do not run uv installs of your own while the implementer runs
(CLAUDE.md 4.4).

## 5. When it reports

Say in one line that it reported, then go back to the discussion. Follow the
`review-round` skill only when the maintainer says to review.
