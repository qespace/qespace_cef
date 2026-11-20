# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `test/conv-0001`, the test of the accepted convention

**Goal.** CONV-0001 is enforced by tests (CLAUDE.md 6.6):

- every even-rank Stevens operator equals the McPhase definition copied into
  CONV-0001, signs included;
- a parameter `Blm` multiplies `O_lm`, and `B4m3` multiplies the sine-type
  `O_4,−3`;
- comparisons with PyCrystalField record its sign differences instead of
  accepting ± (AUDIT P18).

Tests only: no code under `qespace_cef/` changes, and no result changes.

**Steps.** Branch `test/conv-0001`, 3 commits:

1. `test(cefmatrices): check the Stevens operators against CONV-0001`.
   - A new reference `tests/reference/mcphase.py` builds `O_lm` from `J±`,
     `J_z` and the `P_lm` table of CONV-0001. The polynomials are typed from
     the CONV, not from the package.
   - `test_conv_0001_operators` compares it with `StevensBase` for l = 2, 4, 6,
     every m and J = 1/2 … 8. That is 27 operators × 16 J = 432 matrices.
     The tolerance is 1e-10 of the largest element, as in the module's
     `_close`.
2. `test(cefpars): check the parameter labels against CONV-0001`.
   - `test_conv_0001_labels`: for each even-rank label (`B20` … `B6m6`),
     `CEFion(ion, CEFpars({label: x})).hamiltonian` equals `x · O_lm` from the
     reference.
   - The matrices themselves are compared, not the energies. The spectrum
     cannot see the sign of a single `B_lm` with m ≠ 0: a rotation about z by
     π/m maps `O_lm` to `−O_lm` and leaves the spectrum unchanged.
3. `test(cefmatrices): compare with PyCrystalField sign by sign`.
   - `test_matrices_match_pycrystalfield` drops the ±.
   - The 16 operators that differ in sign become a strict `xfail` naming
     CONV-0001: O10, O21, O2m1, O30, O32, O3m2, O41, O4m1, O43, O4m3, O61, O6m1,
     O63, O6m3, O65, O6m5. All of them have l + m odd (run on `eae002d`).
   - The same commit removes P18 from AUDIT, updates the test tree in
     `ARCHITECTURE.md`, and drops "not a test yet" from `ROADMAP.md`.

**Files.** `tests/reference/mcphase.py` (new), `tests/test_cefmatrices.py`,
`tests/test_cefpars.py`, `ARCHITECTURE.md`, `AUDIT.md`, `ROADMAP.md`.

**Risks.**
- The hand check of 2026-11-20 found no difference in 601 matrices, so steps
  1 and 2 should pass. If one fails, the case becomes a strict `xfail` and I
  ask you; the code is not changed on this branch.
- Odd l are not in CONV-0001's table, so they stay checked only against the
  old implementation and PyCrystalField (Q12).
- The label `B43S` for McPhase files (Q33) is not part of this step.
- CONV-0001 has no `enforced_by` line, and an accepted CONV is never edited
  (CLAUDE.md 7). The test names the CONV instead (Q34).

**Tests.** `uv run pytest`: 432 + 27 new checks pass. The PyCrystalField
test passes for 21 operators and is a strict `xfail` for 16. Also
`uv run ruff check .`, `uv run ruff format --check .` and `uv run nox`.

## In parallel: decisions only you can take

- The `OPEN:` parts of the drafts CONV-0003 to CONV-0008 (AUDIT section 6.2,
  Q11), and a draft for the g-tensor (C7).
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## Next

`fix/susceptibility` (P1, P2, P10, A12): it needs no convention, and each fix
already has its strict `xfail`. `refactor/api-v2` can be planned once Q29 and
Q31 are answered, and built once the drafts it touches (CONV-0003, CONV-0006,
CONV-0008) are accepted.
