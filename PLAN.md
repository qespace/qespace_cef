# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `feat/degenerate-levels` (CONV-0008, CONV-0009, CONV-0011, P7)

**Goal.** One rule decides which energies form a level, CONV-0008's
`|E_m − E_n| ≤ rtol · max(max|E|, E_floor)`, from one config place, and
every consumer uses it. Today `susceptibility` and `gtensor` use the absolute
`constants.numerical_zero`, against CONV-0009 (which defers to CONV-0008);
`get_degeneracies` truncates decimals (P7); `summary` has its own threshold
and prints `eigh`'s basis instead of CONV-0011's.

**Steps.** Branch `feat/degenerate-levels`, 5 commits:

1. `feat: one config place and the level grouping of CONV-0008`.
   `qespace_cef/config.py` with `rtol = 1e-10` and `E_floor = 1.0` (meV),
   read at call time; `cefion.degenerate_levels(energies)`, the levels as
   index slices of the sorted energies. Tests `test_conv_0008_*`: every exact
   Kramers doublet of 200 low-symmetry Hamiltonians is one level; the free
   ion is one level; grouping unchanged under a random Hermitian perturbation
   of size ε `max|E|`.
2. `fix(observables)!: degenerate pairs by CONV-0008 in susceptibility and
   gtensor`. Both use the grouping; `numerical_zero` leaves `constants.py`.
   The continuity test is re-chosen so the split falls below and above the
   new tolerance, χ continuous in both. Docstrings and comments cite
   CONV-0008.
3. `fix(cefion)!: group get_degeneracies by CONV-0008 (P7)`. The
   `precision` argument goes; each level is reported with its mean energy.
   The P7 xfail turns green.
4. `feat(cefion): print the canonical basis in summary (CONV-0011)`.
   `summary` groups by step 1 (`_PRINT_DEGENERACY_ATOL` goes) and prints the
   basis of CONV-0011: J_z inside each level, then J_x, J_y where J_z is
   degenerate there (by `rtol · J`); Kramers partners by `Θ` for half-integer
   J in zero field; per-state phase otherwise. `eigenvectors` is unchanged.
   Tests `test_conv_0011_*`: the printed basis is the same after a random
   unitary rotation inside each level; `Θ|−> = |+>` for Kramers pairs;
   `<J_z>` ascending in each level.
5. `docs: describe the degenerate levels`: `ARCHITECTURE.md` (module map,
   tolerance row, test tree), `AUDIT.md` (P7, the CONV-0008 part of A11).
   CHANGELOG lines go with each commit.

**Files.** `qespace_cef/config.py` (new), `qespace_cef/__init__.py`,
`qespace_cef/cefion.py`, `qespace_cef/observables.py`,
`qespace_cef/constants.py`, `tests/test_cefion.py`,
`tests/test_observables.py`, `ARCHITECTURE.md`, `AUDIT.md`, `CHANGELOG.md`.

**Choices for your review.** `get_degeneracies` loses `precision` and
reports mean energies; the J_z degeneracy inside a level uses `rtol · J`
(CONV-0011 names no tolerance).

**Risks.** χ changes only for pairs split between 10⁻¹² meV and the new
tolerance, where the Curie and van Vleck weights agree to O(ΔE/k_BT); the
review measures it. `gtensor` keeps its other problems (P4, P17, A10).

**Tests.** pytest (451 + new, 29 xfailed after P7), doctests, ruff. No nox.

## In parallel: decisions only you can take

- The `OPEN:` parts of the draft CONV-0006 (AUDIT
  section 6.2, Q11), and a draft for the g-tensor (C7).
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## After that

1. **The g-tensor (C7):** a draft CONV on the CONV-0008 grouping and
   CONV-0011's basis, returning G or the principal values and axes. Then P4
   and P17.

Other fixes wait for a convention or a source: P3 (CONV-0003 units of the
thermodynamics, Q44), P14 (sources of the ion
table), P16 (Q12). `refactor/api-v2` waits for Q29 and Q31.
