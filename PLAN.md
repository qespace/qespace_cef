# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `refactor/j-matrices`, one construction of J (A5, CONV-0005)

**Goal.** J_x, J_y, J_z, J₊ and J₋ are built in one place, the module
functions at the top of `cefmatrices.py`, and every other matrix comes from
them (CONV-0005). Today `StevensBase.__init__` repeats their formulas line
for line (AUDIT A5). The two agree bit for bit, so no result is wrong, but a
change to one would not reach the other. No computed result changes.

**Steps.** Branch `refactor/j-matrices`, 3 commits.

1. `refactor(cefmatrices): build StevensBase from the module's J matrices
   (A5)`.
   - `StevensBase.__init__` takes `Jz`, `Jp`, `Jm`, `Jx`, `Jy` from `J_z`,
     `J_plus`, `J_minus`, `J_x`, `J_y` instead of its own formulas. The
     attribute names stay.
   - `J_plus` is written with the m of the lower state,
     `√(J(J+1) − m(m+1))` for m = J−1, ..., −J, as CONV-0005 states it. It
     currently writes the same values through the symmetry m → −m−1.
   - Check: all 784 operator matrices (49 operators, J = 1/2 … 8) are
     bit-identical to a snapshot taken from `ea26fff`, and the existing
     tests pass, including the 432 CONV-0001 matrices.
2. `test(cefmatrices): check the J matrices against CONV-0005` (accepted;
   CLAUDE.md 6.6).
   - `test_conv_0005_*`, from its Consequences:
     - `J_z = diag(J, …, −J)`;
     - the Condon–Shortley elements of `J₊`, alone on its first
       superdiagonal;
     - `[J_x, J_y] = i J_z` and `J² = J(J+1)`;
     - `StevensBase(J)` holds these matrices;
   - for J = 1/2 … 8.
3. `docs: make the J matrices the root of the basis convention`.
   - The `cefmatrices` module docstring cites CONV-0005 and CONV-0001,
     instead of "need to be cross-checked against some reliable book" (AUDIT
     6.1).
   - `ARCHITECTURE.md`: the module map and the test tree.
   - `AUDIT.md`: A5 narrows to the sympy copy in `symbolic/`, and the 6.1
     row goes.
   - `CHANGELOG.md`.

**Files.** `qespace_cef/cefmatrices.py`, `tests/test_cefmatrices.py`,
`ARCHITECTURE.md`, `AUDIT.md`, `CHANGELOG.md`.

**Risks.**
- **None in the results:** the matrices are bit-identical by the check in
  step 1. If the rewritten `J_plus` differs in the last bit for some J, that
  is reported, not hidden: a difference of one rounding step would change
  no test, but it would change the snapshot comparison.
- **Not in this step:**
  - `symbolic/cefmatrices.py`, a third, sympy copy, stays; it is
    experimental, not imported and not tested (A5, A9);
  - generating the Stevens operators from one formula (A6) is later;
  - the McPhase test reference (`tests/reference/mcphase.py`) keeps its own
    J matrices on purpose, as an independent check.

**Tests.** `uv run pytest`: the existing tests pass unchanged, and after
step 2 the `test_conv_0005_*` tests pass too. Then `ruff check`,
`ruff format --check` and `uv run nox`.

## In parallel: decisions only you can take

- The `OPEN:` parts of the drafts CONV-0004 to CONV-0008 and CONV-0011
  (AUDIT section 6.2, Q11), and a draft for the g-tensor (C7).
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## Next

Most remaining fixes wait for a convention (AUDIT section 6.2, Q11):
- P3 for CONV-0003 (its zero of U and F is CONV-0010, accepted);
- P4 and P17 for C7 and CONV-0008;
- P9 and P13 for CONV-0007;
- P7 for A11 and CONV-0008.

The others need no CONV:
- P14 and P15 need sources: the ion table and one CODATA set;
- P16 needs your answer to Q12.

`refactor/api-v2` waits for Q29 and Q31.
