# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `feat/neutron-spectra` (CONV-0007, CONV-0012, P9)

**Goal.** `scattering_function`, `cross_section_crystal` and `cross_section_powder` replace
`spectral_weight`, returning excitations as CONV-0012 defines them: the
scattering function without neutron factors, and the cross section per ion
in absolute units with the full projector.

**Steps.** Branch `feat/neutron-spectra`, 4 commits:

1. `feat(constants): neutron and electron constants of CODATA 2022`.
   `gamma_n = 1.91304276` (|μ_n/μ_N|, CODATA 2022 −1.913 042 76(45)) and
   `r0 = 2.8179403205e-15` m (2.817 940 3205(13)), with source, and
   `(γr₀)²` in barn (≈ 0.2906).
2. `feat(observables): scattering_function (CONV-0012)`.
   `scattering_function(cefion, temperature)` returns energies (n²,) in meV
   and scattering weights (n², 3, 3), complex,
   `W_if^{αβ} = p_i <i|J_α|f><f|J_β|i>`. Tests `test_conv_0012_*`: sum rule
   `Σ_k Tr W_k = J(J+1)`; detailed balance; summed weights per level pair
   unchanged under a unitary rotation inside each level; each `W_k`
   Hermitian positive semidefinite.
3. `feat(observables)!: cross sections replace spectral_weight (CONV-0007)`.
   One kernel `_cross_section_unprojected(cefion, temperature, Qnorm, P, Ei,
   u2)` applies `(γr₀)²`, `(g_J/2)² F(|Q|)²`, `e^{−|Q|²u2/3}`, `k_f/k_i`
   and `Σ P_αβ W^{αβ}`; NaN where `ħω > E_i`. Two wrappers only build P
   and check Q's shape: `cross_section_crystal(cefion, temperature, Q, Ei,
   u2=0)` with Q `(..., 3)` in Å⁻¹ (local frame, CONV-0004) and
   `P = δ − Q̂Q̂`; `cross_section_powder(...)` with `|Q|` `(...,)` and
   `P = (2/3) δ`. Both return energies and cross-section weights (n²,) in
   barn sr⁻¹ per ion; `Ei` (meV) is required. `spectral_weight` goes. Tests
   `test_conv_0007_*`: absolute sum rule
   `Σ_k w_k / (k_f/k_i)_k = (γr₀)² (g_J/2)² (2/3) J(J+1)` for Q → 0,
   `u2 = 0`; NaN exactly where `ħω > E_i`; powder equals the crystal
   averaged over Q directions; full projector against the diagonal for Q
   along an axis. The P9 xfail turns green; `test_performance.py` calls
   `cross_section_powder`.
4. `docs: describe the neutron spectra`: `README.md`, `ARCHITECTURE.md`,
   `AUDIT.md` (P9, the CONV-0007 and CONV-0012 rows). CHANGELOG lines
   go with each commit.

**Files.** `qespace_cef/constants.py`, `qespace_cef/observables.py`,
`qespace_cef/__init__.py`, `tests/test_constants.py`,
`tests/test_observables.py`, `tests/test_performance.py`, `README.md`,
`ARCHITECTURE.md`, `AUDIT.md`, `CHANGELOG.md`.

**Decided in discussion.** One kernel, two wrappers (Q2); `Ei` required,
fixed-E_f kinematics later (Q3a); `k_f/k_i` taken at ε_k, exact for δ-peaks
and slightly off across a broadened line, moved to the grid later (Q4a).

**Risks.** The cross section changes for every caller: absolute units and
`(g_J/2)²` scale it by `(γr₀)² (g_J/2)²`, and P9 changes single-crystal
spectra off the axes. The review measures both against `spectral_weight` on
`main`.

**Tests.** pytest (496 + new, 29 xfailed after P9), doctests, ruff. No nox.

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
