# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `fix/van-vleck`, the van Vleck weight without cancellation (P20)

**Goal.** `susceptibility(method='perturbation')` evaluates the van Vleck
weight of CONV-0009 without cancellation, as CONV-0008 (accepted) asks of
every formula that is continuous through a degeneracy. Today it computes
`(p_n − p_m)/(E_m − E_n)` as a plain difference, which loses digits for a
pair just above the degeneracy tolerance (AUDIT P20). No result changes
beyond rounding, except for such pairs, which become accurate.

**Steps.** Branch `fix/van-vleck`, 1 commit:

1. `fix(observables): evaluate the van Vleck weight without cancellation
   (P20)`.
   - In the perturbation branch of `susceptibility`, the weight of a pair
     with `E_m > E_n` (not degenerate) becomes
     `p_n · (−expm1(−(E_m − E_n)/k_B T)) / (E_m − E_n)`. `expm1` is
     evaluated only on those entries, where its argument is never positive,
     so nothing overflows. The pairs `E_m < E_n` take the same value
     mirrored, `W_mn = W_nm`, since W is symmetric. Degenerate pairs keep
     the Curie weight `p_n / k_B T`.
   - The degeneracy rule stays `constants.numerical_zero`, as CONV-0009
     fixes it; moving to the CONV-0008 tolerance is a later step.
   - New tests in `tests/test_observables.py`:
     - `test_susceptibility_continuous_through_degeneracy`: Ce³⁺
       (`B20=0.5, B40=0.01, B44=0.1` meV) at 5 K, `chi_zz` with a field of
       10⁻⁹, 10⁻⁸ and 10⁻⁷ T along x, which splits the ground doublet just
       above the tolerance, against `chi_zz` in zero field, to 10⁻¹²
       relative. It fails without the fix: 6.9·10⁻⁷ at 10⁻⁹ T, measured on
       `2446607`.
     - `test_susceptibility_low_temperature_finite`: the same ion at
       0.05 K with numpy warnings turned into errors gives a finite χ.
   - `AUDIT.md` (P20 goes) and `CHANGELOG.md` in the same commit.

**Files.** `qespace_cef/observables.py`, `tests/test_observables.py`,
`AUDIT.md`, `CHANGELOG.md`.

**Checked in advance** (scratch run on `2446607`, not a test): the mirrored
`expm1` form gives `chi_zz` rel. changes of 1.8·10⁻¹⁵, 4.4·10⁻¹⁶ and
1.9·10⁻¹⁵ at 10⁻⁹, 10⁻⁸ and 10⁻⁷ T, against 6.9·10⁻⁷, 1.3·10⁻⁸ and
1.4·10⁻⁹ for the plain difference, and stays finite at 0.05 K.

**Risks.**
- The existing susceptibility and CONV-0009 tests must pass unchanged: the
  new form equals the old one up to rounding for every pair that is not
  nearly degenerate.
- Evaluating `expm1` on the whole matrix and masking afterwards would
  overflow for the reversed pairs and warn; the second new test catches it.

**Tests.** `uv run pytest`: two new tests pass, all others unchanged
(447 passed, 30 xfailed on `d1ee227` before the change). Then
`ruff check`, `ruff format --check` and `uv run nox`.

## In parallel: decisions only you can take

- The `OPEN:` parts of the drafts CONV-0004, CONV-0006 and CONV-0007 (AUDIT
  section 6.2, Q11), and a draft for the g-tensor (C7).
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## Next

1. **A CONV that supersedes CONV-0009's degeneracy rule** with the tolerance
   of CONV-0008 (accepted). A short decision for you to accept, then:
2. **`feat/degenerate-levels`**, built on CONV-0008 and CONV-0011:
   - one config place for `rtol` and `E_floor` (part of A11);
   - one function that groups the energies into levels;
   - `susceptibility`, `get_degeneracies` (P7) and `summary` use it;
   - `summary` prints the canonical basis of CONV-0011;
   - the `test_conv_0008_*` and `test_conv_0011_*` tests.
3. **The g-tensor (C7):** a draft CONV, once 2 exists. It identifies the
   ground doublet by CONV-0008 and returns G, or the principal values and
   axes. Then P4 and P17.

Other fixes wait for a convention or a source: P3 (CONV-0003 units of the
thermodynamics, Q44), P9 and P13 (CONV-0007), P14 and P15 (sources), P16
(Q12). `refactor/api-v2` waits for Q29 and Q31.
