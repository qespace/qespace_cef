# qespace_cef: instructions for Claude Code

`qespace_cef` is scientific code for the problem of a magnetic ion in a crystal
electric field (CEF): building the CEF Hamiltonian, diagonalising it, and
computing the observables that follow from it. Results go into the maintainer's
research, so **correctness and traceability come before speed and features.**

Read `ARCHITECTURE.md` (module map, conventions, test tree) and `DEV.md`
(commands, dependencies, releases) before planning anything. When they and this
file disagree, this file wins; flag the disagreement to the maintainer.

Lines marked `OPEN:` are not decided yet. Do not decide them yourself; ask.

## 1. Nothing reaches `main` unreviewed

No change reaches `main` before the maintainer has read its diff. The same
holds for subagents and worktrees.

### 1.1 Before you start

1. The working tree must be clean, so the diff shows only your change. If it
   isn't, stop and ask; don't commit other work to make room.
2. **Estimate the work package** and say it in one line: how many logical
   changes (commits, section 8), which modules, what you will test. One change
   → 1.2. Two or more changes that belong together (e.g. a core change and its
   test, then the API using it, then docs) → **suggest a branch** (1.3), with
   its name, and wait for the maintainer's OK. The maintainer may also ask for
   a branch, or for single changes, at any time.

### 1.2 A single change

1. Write the change to the files. Do **not** stage or commit it.
2. Summarise it in chat: what changed, why, how it was checked (numbers), what
   was not checked. Print the proposed commit title and body (section 8).
3. The maintainer reviews the diff, may edit it, and approves by committing or
   by telling you to commit, with the chosen message. Approval covers that
   diff only; anything you change afterwards needs a new review.
4. One unreviewed change at a time: don't start the next one while the last is
   uncommitted.

### 1.3 A branch

1. Only after the maintainer's OK to the suggested branch, or on their request.
2. One branch per work package, named `<type>/<topic>` (types as in section 8,
   e.g. `feat/wybourne-conversion`), cut from an up-to-date `main`.
3. Commit on the branch without asking, one logical change per commit,
   following section 8, each green on its own (`uv run pytest`).
4. Stop and report after the work package is done or after about 5 commits,
   whichever comes first. If the work turns out larger than estimated, stop
   and say so before going on.
5. The report: the commits (`git log --oneline main..<branch>`), test numbers,
   what was not checked, and the command for the diff
   (`git diff main...<branch>`).
6. Corrections after review are new commits on the branch. The maintainer
   merges, or tells you to merge a branch they have reviewed. Never merge into
   `main` or commit on `main` yourself otherwise.

### 1.4 Always

1. Stage explicit paths, never `git add -A` or `git add .`.
2. Never push unless asked. Never tag, amend, rebase, reset, force or otherwise
   rewrite history. Delete a branch only after it is merged and the maintainer
   agrees.

## 2. Working with the maintainer

1. The maintainer is a physicist, strong in Python. Explain physics choices in
   physics terms, code choices briefly.
2. Read-only requests ("review", "advise", "plan", "read only") change no files.
3. Plan first for anything touching more than one module, the public API, or the
   physics: show the plan (goal, steps, files, risks, test plan) and wait for OK.
4. Ask when a choice changes the outcome: conventions, units, scope, API. Give
   every question a short ID (`Q1`, `Q2`, ...) so the maintainer can answer by
   ID.
5. Report outcomes with numbers: tests passed / failed / xfailed, deviations
   from references and the tolerance used, versions. Say what was not verified.
6. Where you may look and write:
   1. Scratch work (temporary tests, plans, renders) goes to a place of your
      choice, never into tracked files.
   2. Other projects in the qespace folder only on the maintainer's written
      request naming the project and what to look for (e.g. "see
      `qespace_crystallography` for how the lattice is implemented"). Read them;
      don't edit them.
   3. Never read, list or search anything outside the qespace folder. If
      something isn't there, ask where it is.

## 3. The physics is the maintainer's decision

1. Physics conventions (signs, units, normalisation, operator basis, frames)
   are the maintainer's to set. They are recorded in the table below and in
   `ARCHITECTURE.md`. Don't pick one. When code, a reference code or the
   literature disagree, record the case as a strict `xfail` with the reason
   and ask.
2. **Never invent a value silently.** No hard-coded physical constants, CEF
   parameters or literature values without a source:
   1. physical constants live in one module, with their source (e.g. CODATA
      year);
   2. literature parameters carry units and a reference (DOI) next to the
      value;
   3. a value that is missing is reported as missing, not guessed.
3. Every quantity in the API and in docstrings states its unit.

| Topic | Convention |
|---|---|
| Energy unit | `OPEN:` (meV?) |
| CEF parameters | `OPEN:` Stevens `B_k^q` / Wybourne `B^k_q`, and the conversion between them |
| Stevens operator normalisation and tables | `OPEN:` |
| Quantisation axis and local frame | `OPEN:` and how it relates to the crystal frame |
| Basis | `OPEN:` ground `J` multiplet only, or intermediate coupling |
| Zeeman term | `OPEN:` sign, `g_J` vs full `L + 2S` |
| Cross sections | `OPEN:` normalisation, form factor, units |
| Degenerate levels | Eigenvectors within a degenerate subspace are basis dependent: compare projectors or sums over the subspace, never single eigenvectors |

`qespace_cef` is a package of the qespace suite. Suite rules reach this
repository as entries in `CONVENTIONS.md` (accepted CONVs, with the test that
enforces them) or as verbatim copies of suite ADRs, which keep their ID and say
where they came from. Those take precedence over this table; a copy is never
edited here.

## 4. Environment

1. Windows. Managed with `uv` only: no `pip`, `conda`, `venv`, `poetry`.
   `python` is not on the PATH; use `uv run`.
2. **Ask before adding any dependency**, runtime or dev, with one sentence on
   why numpy, the standard library or existing code can't do it.
3. Don't change the maintainer's `.venv` for experiments: use
   `UV_PROJECT_ENVIRONMENT=<short temp path>` or
   `uv run --no-project --isolated ...`. Keep temp paths short; deep paths
   break numpy DLL loading.
4. Run `uv` installs one after another, never in parallel (they race in the
   uv cache on Windows).
5. `uv.lock` is tracked. After editing dependencies run `uv lock`; change the
   version only with `uv version`.

## 5. Code

1. The core is computational and needs only `numpy` (`OPEN:` `scipy`?).
   Plotting and other optional packages are extras, imported lazily inside
   functions; `import qespace_cef` never imports them.
2. Prefer the smallest change that works: plain functions, explicit code, no
   layers or abstractions before a feature needs them. When one is added,
   update `ARCHITECTURE.md` in the same change.
3. Reformatting the maintainer's hand-written code is welcome. Do it as its own
   change (`style:`), not mixed into a change of behaviour.
4. Every public module, class and function has a NumPy-style docstring: what it
   does or returns, units, conventions, and array shapes.
5. Validate array arguments (shape, hermiticity) at the public boundary, not
   deep inside computations.
6. Numerical checks are relative to the scale of the checked object and use
   tolerances from one config place, read at call time. No ad hoc absolute
   thresholds.
7. numpy 1 and 2 must both work:
   1. arrays hashed via `tobytes()` get an explicit dtype (`np.int64`);
   2. convert numpy scalars with `float()` before `Fraction`;
   3. keep numpy scalar reprs (`np.float64(1.0)`) out of strings and labels;
   4. no APIs removed in numpy 2 (`np.float_`, `np.product`, `np.in1d`, ...).
8. Pure computations don't print or log. Warnings and logging happen at the
   edges (file I/O, fits, calls into outside libraries).

## 6. Tests

1. `tests/` mirrors the package: `qespace_cef/<pkg>/<module>.py` is tested in
   `tests/<pkg>/test_<module>.py`. Move and rename tests with their modules.
2. Test behaviour through the public entry point and assert a result, not just
   that the code ran. Prefer `parametrize` or tables over copied tests.
3. A bug fix comes with a test that fails without the fix.
4. **Never loosen or delete a test to make it pass.** Fix the code, or record
   the bug as a strict `xfail` with the reason.
5. Check physics against independent results: analytic limits (free ion,
   cubic CEF, Kramers degeneracy, sum rules, traces) and reference codes
   (`OPEN:` which ones, e.g. McPhase, PyCrystalField, Mantid CrystalField).
6. Regression reference data are regenerated only after a deliberate change of
   results, and the summary and commit body say so, with the size of the change.
7. Before finishing: `uv run pytest` (`OPEN:` plus `ruff`, and a numpy 1 / 2
   matrix with `nox`).

## 7. Documentation

Keep in sync with the code, in the same change: the module map and test tree in
`ARCHITECTURE.md`, the commands and rules in `DEV.md`, the usage in `README.md`
(copy-pasteable commands). Record significant decisions in a short ADR
(context, decision, consequences; under 20 lines).

## 8. Commit messages (section 1.2 and 1.3)

1. [Conventional Commits 1.0](https://www.conventionalcommits.org):
   `type(scope): description`, lower-case imperative, no full stop, at most 72
   characters. A body after a blank line says why. Footers last.
2. Types: `feat`, `fix`, `docs`, `test`, `refactor`, `perf`, `style`, `build`,
   `ci`, `chore`. A change that breaks the API or changes computed results gets
   `!` and a `BREAKING CHANGE:` footer saying what changed.
3. One logical change per commit, each commit green on its own (e.g. core
   change and its test, then the API using it, then docs).
4. Staging: see 1.4.

## 9. Definition of done

1. Tests added or updated, `uv run pytest` green (numbers reported).
2. No new dependency without approval; no `OPEN:` decided without the
   maintainer.
3. Docs updated (section 7).
4. A summary for the maintainer: what changed, how it was checked, what was
   not, the proposed commit message, or the branch report (1.3.5) when a
   branch was used. Nothing on `main` until the maintainer has reviewed it
   (section 1).
