# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0004). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `fix/pre-v2` (the independent review, CONV-0003, CONV-0013, CONV-0015)

**Goal.** v1 gives correct results before v2 is built against it (ADR-0005:
v2's tests compare with v1). The findings of the independent review
`AUDIT_GPT.md` that change results, reproduced on `9b49c3d`, are fixed in
one run. Physical results change only where they were wrong.

**Steps.** Branch `fix/pre-v2`, one commit each:

1. `fix(observables)!: susceptibility at zero field (H1, CONV-0013)`.
   Both methods use the zero-field eigensystem of the ion, whatever its
   field: χ is not a function of B. Test: free Ce³⁺ solved at 10 T gives
   the Curie value with both methods.
2. `fix(observables)!: thermodynamics in meV and k_B per ion (P3, H2, H3,
   M8)`. Boltzmann factors from `E_n − E_min`; U and F in meV per ion from
   the zero of CONV-0014, S and the heat capacity in k_B per ion
   (CONV-0003); the heat-capacity variance about the thermal mean. The P3
   oracles get k_B (`S = −(dF/dT)/k_B`, `C = (dU/dT)/k_B`), stated in the
   commit body as a correction of the oracle, and their xfails turn green.
   Tests: S and C unchanged under `B00`; finite results at 0.1 K, S → ln 2
   for a ground doublet.
3. `fix(observables): the g-tensor needs exactly one ground doublet (H4)`.
   A `ValueError` instead of the `assert` when the ground level is not
   twofold (A10). P4 and P17 stay for C7.
4. `fix(cefion): pair Kramers partners only for time-even Hamiltonians
   (H5, L2)`. `summary` pairs by time reversal only in zero field and
   without odd-rank parameters; every printed state lies inside its
   level; `precision=0` keeps coefficients at the cutoff. `CEFpars` warns
   (`UserWarning`) when odd-rank parameters are given (Q12: allowed, with
   a warning).
5. `fix(cefion)!: degenerate levels by the spread of the energies (M5,
   CONV-0015)`. The tolerance `rtol · max(max E − min E, E_floor)`; the
   `test_conv_0008_*` tests become `test_conv_0015_*`, plus `B00`
   invariance of the grouping and of χ.
6. `fix: validate public inputs (M3, M6, L1)`. A field of three
   components; parameter labels `B…` only; finite energies in the level
   grouping; J a non-negative multiple of ½; finite crystal Q;
   temperatures > 0 for populations; a scalar temperature in
   `scattering_function`; `max_iterations ≥ 1`. Each a `ValueError` with a
   test.
7. `docs: describe the fixes`: `ARCHITECTURE.md`, `README.md` where usage
   changes, `CHANGELOG.md`, `AUDIT.md` (P3, P16 and the review's findings
   fixed here leave section 10).

**Files.** `qespace_cef/observables.py`, `qespace_cef/cefion.py`,
`qespace_cef/cefpars.py`, `qespace_cef/cefmatrices.py`,
`tests/test_observables.py`, `tests/test_cefion.py`,
`tests/test_cefpars.py`, `tests/test_cefmatrices.py`, `ARCHITECTURE.md`,
`README.md`, `CHANGELOG.md`, `AUDIT.md`.

**Not in this step.** P21 (the Zeeman signs) and M4 (mutable inputs) go
with v2 (ADR-0005); H6 (the point-group table) waits for the
crystallography symmetry API; P4, P17 for C7; M7 (the fitting sketch).

**Tests.** pytest (536 passed, 28 xfailed on `9b49c3d`; the P3 xfails turn
green), doctests, ruff. No nox.

## Then: decisions before v2 (advisor and maintainer)

1. **The CONVs that name code** (ADR-0004): CONV-0007 and CONV-0011, their
   values unchanged. CONV-0003 is done, by editorial edits.
2. **Plan `refactor/api-v2`:** the branches, each green, batched results
   checked against v1 (ADR-0005, ADR-0007), with P21, M4 and the other
   AUDIT items v2 rewrites.

## After that

1. **`refactor/api-v2`**, planned in *Then*, step 2.
2. **The g-tensor (C7):** a draft CONV on the CONV-0015 grouping and
   CONV-0011's basis, returning G or the principal values and axes. Then P4
   and P17.

Other fixes wait for a source: P14 (sources of the ion table).
