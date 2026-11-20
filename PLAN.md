# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `fix/susceptibility`, P1, P2, P10, A12 and CONV-0009

**Goal.** `observables.susceptibility` follows CONV-0009 (accepted):

- it returns the tensor `χ_αβ` in μB/T per ion, of shape `(nT, 3, 3)`, with
  α the moment component and β the field component;
- the formula is Van Vleck's:
  `gJ² μB Σ_nm W_nm Re(<n|J_α|m><m|J_β|n>)`, where
  `W_nm = p_n/k_B T` for `|E_m − E_n| < constants.numerical_zero` and
  `(p_n − p_m)/(E_m − E_n)` for every other pair, in both orders;
- both methods return it: `'perturbation'` from the formula,
  `'magnetization'` as `M_α(ε e_β)/ε`;
- it prints nothing.

The argument `Hfield_direction` goes. `χ_n̂ = n̂·χ·n̂` and `Tr χ / 3` are
the caller's.

**Steps.** Branch `fix/susceptibility`, 5 commits. Each commit removes the
`# AUDIT` markers, Warnings entries, `xfail` markers and AUDIT rows of what
it fixes, and adds a CHANGELOG line.

1. `fix(observables)!: count both orderings in the van Vleck term (P1)`.
   - The sum over `E_m > E_n` becomes a sum over every non-degenerate pair.
     Each ordering gives the same amount, so the code's sum was half.
   - Still on today's API, so the physics change can be reviewed on its own.
   - `test_susceptibility_perturbation_along_z` turns green.
   - Results change: along z, χ × μB goes from 0.12704 / 0.02680 / 0.00847
     to 0.13130 / 0.02981 / 0.00999 μB/T at 5 / 20 / 100 K, which is +3 %,
     +11 % and +18 %.
2. `fix(observables): print nothing in susceptibility (A12)`.
   - The `print(jumps)` goes, with the Notes line about it.
   - A `capsys` test checks that `'perturbation'` writes nothing to stdout.
     It fails before the fix.
3. `feat(observables)!: return the susceptibility tensor in μB/T (CONV-0009)`.
   - Both methods return `(nT, 3, 3)` in μB/T, and `Hfield_direction` is
     removed.
   - `'perturbation'` builds the tensor from `Re(<n|J_α|m><m|J_β|n>)`. With
     no direction to project, P2 is gone by construction.
   - `'magnetization'` takes three finite differences, `M(ε e_β)/ε`, one per
     field axis, without the norm, so P10 is gone too.
   - The existing tests move to the tensor with the same assertions:
     - along z, `χ_zz` against the exact formula;
     - for P2, `χ_xx = χ_yy`, `χ_xy = 0`, and n̂·χ·n̂ for (1, 1, 0) equal to
       that for (1, 0, 0);
     - for P10, n̂·χ·n̂ of `'magnetization'` for (1, 0, 1) against the exact
       formula.
   - Their `xfail` markers go.
   - Results change twice. `'perturbation'` was in μB²/meV and is now ×μB
     (0.05788 meV/T). Off the axes, the values change as P2 and P10 describe:
     (1, 1, 0) is −29 %, and (1, 0, 1) at 20 K is −8.5 %.
4. `test(observables): check the susceptibility against CONV-0009`.
   - The tests of its Consequences, named `test_conv_0009_*` as its
     `enforced_by` says:
     - the free-ion Curie law;
     - the high-T limit of `Tr χ / 3`;
     - `'perturbation'` against `'magnetization'`;
     - χ symmetric and positive semi-definite;
     - invariance under a unitary rotation inside a degenerate level;
     - tetragonal `χ_xx = χ_yy`, `χ_xy = 0`.
   - The expected values come from formulas in the tests.
5. `docs: describe the susceptibility tensor`.
   - `ARCHITECTURE.md` (the module map and test tree) and the README table.
   - AUDIT: section 7 and the CONV-0009 row in 6.2.

**Files.** `qespace_cef/observables.py`, `tests/test_observables.py`,
`ARCHITECTURE.md`, `README.md`, `AUDIT.md`, `CHANGELOG.md`.

**Checked in advance** (scratch script on `423827b`, not a test). The fixed
formulas, with the package's constants and degeneracy rule, agree with the
tests' exact formula:

- to at most 1.8·10⁻⁷ relative along z, x, (1, 1, 0) and (1, 0, 1);
- to 3.4·10⁻⁷ for `'magnetization'`.

The tests allow 1e-6 and 1e-5. The difference comes from the constants: the
tests use CODATA 2018, the package `meV2K = 11.60452` (P15). The CONV-0009
formula matches −∂²F/∂B² to ≤ 2.4·10⁻⁶, which is finite-difference noise,
for Ce³⁺ of low symmetry and cubic Pr³⁺. It is positive definite in all four
cases, and Tr χ/3 → Curie (0.99852 at 100 K, 1.00000 at 10⁴ K).

**Risks.**
- **The tolerance.** CONV-0009 fixes the degeneracy rule to the absolute
  `constants.numerical_zero` (10⁻¹² meV). CLAUDE.md 5.6 asks for relative
  tolerances from one config place; the CONV takes precedence. A11 and
  CONV-0008 (draft) stay open for the general case.
- **The constants.** μB and k_B come from `constants` (P15). The unit change
  in step 3 uses `constants.uB`, so it carries P15's mixed CODATA vintages.
- **P5.** `cefion.J` is used as it is. After a second `diagonalize()` it is
  rotated twice (P5), and χ with it. That belongs to `fix/cefion`.
- **Speed.** `'magnetization'` builds three `CEFion` per temperature (A3), so
  it gets 3× slower. The perturbation method still loops over T.
  Vectorising waits for v2 (ADR-0003).
- **Size.** 5 commits is the limit of CLAUDE.md 1.3.4. If step 3 grows, I
  stop after it and report.

**Tests.** `uv run pytest`:
- the three susceptibility `xfail`s are green (P1 in step 1, P2 and P10 in
  step 3);
- the A12 test and the `test_conv_0009_*` tests pass, and each new test of a
  fix fails without that fix;
- then `ruff check`, `ruff format --check` and `uv run nox`.

## In parallel: decisions only you can take

- The `OPEN:` parts of the drafts CONV-0003 to CONV-0008 (AUDIT section 6.2,
  Q11), and a draft for the g-tensor (C7). CONV-0003 (units) should agree
  with CONV-0009: χ in μB/T per ion.
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## Next

`fix/cefion` (P5, P6, P8, P11), which needs no convention. P7 waits for A11
and CONV-0008.
