# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `fix/cefion`, P5, P6, P8 and P11

**Goal.** Four bugs that need no convention, each fixed with its strict
`xfail` turning green:

- `CEFion.diagonalize` can be called twice;
- `sortWithE=False` works;
- `str(Ion)` shows a half-integer J;
- `magnetization_exchange` stops when the mean field doesn't converge.

No computed result of a correct call changes.

**Steps.** Branch `fix/cefion`, 4 commits. Each removes its `# AUDIT`
marker, Warnings entry, `xfail` marker and AUDIT row, and adds a CHANGELOG
line.

1. `fix(cefion): rotate J from the free-ion operators (P5)`.
   - `diagonalize` rotates `J_x(J)`, `J_y(J)`, `J_z(J)`, built anew, instead
     of the current `Jx`, `Jy`, `Jz`, which are already in the eigenbasis
     after a first call. The operators were rotated twice; now a repeated
     call gives the same state.
   - `test_second_diagonalize_changes_nothing` turns green.
   - The first call is unchanged: in the scratch run the rotated operators
     match today's to 0.0.
2. `fix(cefion): keep the order of eigh with sortWithE=False (P6)`.
   - The indices are `np.arange(2J + 1)` (an int) instead of
     `np.arange(J)` (a float, of the wrong length).
   - `test_diagonalize_without_sorting` turns green. Its marker has
     `raises=IndexError`, so it can't turn green any other way.
   - `eigh` already returns ascending energies, so the order is the same as
     with sorting; the docstring says so.
3. `fix(ion): show a half-integer J as a fraction (P8)`.
   - `str(Ion)` prints J with `Fraction(float(J))` (CLAUDE.md 5.7.2):
     `'Ce: J = 5/2, gJ = 0.86'`, and `'Ho: J = 8, ...'` stays as it is.
   - `test_str_shows_half_integer_J` turns green (Ce, Dy).
4. `fix(observables): stop the mean-field iteration (P11)`.
   - At most `max_iterations` steps (a new keyword argument, default
     1000), then `RuntimeError`.
   - The convergence test becomes `|M_k+1 − M_k| ≤ precision · max(|M_k+1|,
     |M_k|)`, so M = 0 converges at once, without the 0/0.
   - The placeholder `residual = 42` goes.
   - Both P11 tests turn green: the stand-in that never converges raises
     `RuntimeError` before its 10000-call guard, and the zero field gives
     exactly (0, 0, 0) without a warning.
   - The same commit updates the test tree in `ARCHITECTURE.md` and AUDIT
     section 7.

**Files.** `qespace_cef/cefion.py`, `qespace_cef/ion.py`,
`qespace_cef/observables.py`, `tests/test_cefion.py`, `tests/test_ion.py`,
`tests/test_observables.py`, `AUDIT.md`, `CHANGELOG.md`,
`ARCHITECTURE.md`.

**Checked in advance** (scratch script on `efb4a0b`, not a test):
- P5: rotating the free-ion operators gives the first call's operators
  exactly.
- P6: the indices are `[0 … 5]` for Ce.
- P8: Ce 5/2, Dy 15/2, Ho 8, Pr 4.
- P11, with the new loop: zero field gives (0, 0, 0) with warnings as
  errors; Ce³⁺ (tetragonal test values) at 5 K, 1 T along z, λ = 1 T/μB
  converges to M_z = 0.15001 μB.

**Risks.**
- **The cap.** `max_iterations = 1000` is my choice of default, not a
  physics value; CLAUDE.md 5.6 asks for tolerances from one config place
  (A11). `precision` stays a keyword argument as today. Both move to the
  config with A11.
- **The error type.** `RuntimeError` on no convergence is a new error the
  caller can see. A10 (`ValueError` at the boundary) is about invalid
  arguments, not this.
- **Not in this step:**
  - `CEFion.J` is still a list before `diagonalize` and an array after
    (A2, v2);
  - P7 (degeneracy grouping) waits for A11 and CONV-0008;
  - the mean field still uses `magnetization`, one `CEFion` per step (A3).

**Tests.** `uv run pytest`: the five `xfail`s of P5, P6, P8 and P11 (two
of them) turn green; 36 → 31 xfailed. Then `ruff check`,
`ruff format --check` and `uv run nox`.

## In parallel: decisions only you can take

- The `OPEN:` parts of the drafts CONV-0003 to CONV-0008 (AUDIT section 6.2,
  Q11), and a draft for the g-tensor (C7). CONV-0003 (units) should agree
  with CONV-0009: χ in μB/T per ion.
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## Next

Most remaining fixes wait for a convention (AUDIT section 6.2, Q11):
- P3 for CONV-0003;
- P4 and P17 for C7 and CONV-0008;
- P9 and P13 for CONV-0007;
- P7 for A11 and CONV-0008.

The others need no CONV:
- P14 and P15 need sources: the ion table and one CODATA set;
- P16 needs your answer to Q12.

`refactor/api-v2` waits for Q29 and Q31.
