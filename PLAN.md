# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `fix/constants`, CODATA 2022 (P15, Q43)

**Goal.** `constants.py` takes its values from one CODATA set, 2022 (your
answer to Q43), and cites it, as CONV-0003 asks ("constants from one module
that cites its CODATA set"). The unit conversions are derived from those
values instead of being typed as separate numbers (AUDIT P15). Computed
results move by a few 10⁻⁷ relative, the difference between today's mixed
1986/2010 values and CODATA 2022.

**Values** (read from NIST on 2026-11-21; each page says "Source: 2022
CODATA"):

| Quantity | Value | Kind |
|---|---|---|
| elementary charge e | 1.602 176 634 × 10⁻¹⁹ C | exact (SI 2019) |
| Planck constant h | 6.626 070 15 × 10⁻³⁴ J s | exact (SI 2019) |
| speed of light c | 299 792 458 m s⁻¹ | exact (SI 2019) |
| Boltzmann constant k | 1.380 649 × 10⁻²³ J K⁻¹ | exact (SI 2019) |
| Bohr magneton μB | 9.274 010 0657(29) × 10⁻²⁴ J T⁻¹ | measured |

Checks against NIST's own derived values: μB = 5.788 381 7982(18) × 10⁻⁵
eV T⁻¹, and 1 m⁻¹ ≙ 1.239 841 984… × 10⁻⁶ eV (exact).

Reference, in the module docstring with the link:
P. J. Mohr, D. B. Newell, B. N. Taylor, and E. Tiesinga, "CODATA
recommended values of the fundamental physical constants: 2022", Rev. Mod.
Phys. 97, 025002 (2025), doi:10.1103/RevModPhys.97.025002; companion article
J. Phys. Chem. Ref. Data 54, 033105 (2025); values at
https://physics.nist.gov/cuu/Constants/index.html.

**Steps.** Branch `fix/constants`, 2 commits.

1. `fix(constants)!: use CODATA 2022 and derive the conversions (P15)`.
   - `constants.py`:
     - a module docstring with the reference and the link above;
     - the five values of the table, each with its unit and "exact" or its
       uncertainty;
     - the conversions the package uses derived from them:
       - `K2meV = 1e3 · k / e`, and `meV2K = 1 / K2meV`;
       - `uB = 1e3 · μB / e` (meV/T);
       - `invcm2meV = 1e3 · 100 · h c / e`;
     - `numerical_zero` unchanged, since CONV-0009 names it.
   - The names nothing uses go: `pi` (use `np.pi` where needed), the neutron
     and electron-radius values (`_gn`, `_hq`, `_gamn`, `_r0`, `R0`),
     `kB_SI`, `NA`, `C1`–`C3`.
   - The tests keep their own reference values, written from the table
     above, not imported from the module (CLAUDE.md 3.2.1, 6.2):
     - `KB` and `MUB` in `tests/test_observables.py`;
     - `K_IN_MEV` and `INVCM_IN_MEV` in `tests/test_cefpars.py`;
     - `MUB_CODATA` in `tests/test_cefion.py`;
     - each with the CODATA 2022 source.
   - New `tests/test_constants.py`:
     - `K2meV · meV2K == 1` to rounding;
     - `K2meV`, `uB` and `invcm2meV` equal the values from the table,
       computed in the test, to 10⁻¹⁵ relative;
     - `uB` agrees with NIST's 5.788 381 7982 × 10⁻⁵ eV/T within its
       uncertainty.

     These tests fail today: the conversions aren't reciprocal, and `uB` is
     4.8·10⁻⁷ off.
   - `BREAKING CHANGE:` results move by a few 10⁻⁷ relative. The unused
     names are removed from `qespace_cef.constants`.
2. `docs: describe the constants`.
   - `ARCHITECTURE.md`: the `constants.py` row and the test tree.
   - `AUDIT.md`: P15 and Q43 go.
   - `CHANGELOG.md`.

**Files.** `qespace_cef/constants.py`, `tests/test_constants.py` (new),
`tests/test_observables.py`, `tests/test_cefpars.py`, `tests/test_cefion.py`,
`ARCHITECTURE.md`, `AUDIT.md`, `CHANGELOG.md`.

**Risks.**
- **The size of the change.** Every result that uses `uB`, `meV2K` or
  `invcm2meV` moves. The review measures it against a snapshot from `main`
  (energies in a field, χ, M) and states it in the report and the commit
  body (CLAUDE.md 6.7). Expected: about 5·10⁻⁷ for μB, about 2·10⁻⁷ for
  `meV2K`.
- **Tests with tight tolerances**, at 10⁻⁶ against CODATA 2018 values,
  become tighter-fitting once module and tests share CODATA 2022: 2018 and
  2022 differ only in μB, by about 10⁻⁹. No tolerance is loosened.
- **Removing names** is an API change in a public module. Nothing in the
  package, the tests or the docs uses them (grepped on `1bfac76`).

**Tests.** `uv run pytest`: `test_constants.py` passes, everything else
unchanged in count (451 passed, 30 xfailed on `053650e`). Then the
doctests, `ruff check` and `ruff format --check`. No nox: the step changes
no dependency (CLAUDE.md 6.8).

## In parallel: decisions only you can take

- The `OPEN:` parts of the drafts CONV-0004, CONV-0006 and CONV-0007 (AUDIT
  section 6.2, Q11), and a draft for the g-tensor (C7).
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## After that

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
thermodynamics, Q44), P9 and P13 (CONV-0007), P14 (sources of the ion
table), P16 (Q12). `refactor/api-v2` waits for Q29 and Q31.
