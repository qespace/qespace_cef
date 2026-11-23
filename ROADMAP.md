# Roadmap

Ideas for features, written for people. Nothing here is decided or planned:
an idea becomes work when it moves into [`PLAN.md`](PLAN.md), and a rule only
when it becomes an ADR or a CONV in [`decisions/`](decisions/). Known
problems with the current code are in [`AUDIT.md`](AUDIT.md).

## Version 2: the broadcasting API

Decided in ADR-0005 (the objects), ADR-0006 (the names) and ADR-0007
(broadcasting and array layout); still to build.

- `import qespace_cef as cef`; `cef.Ion`, `cef.Blm`, `cef.Hamiltonian`.
- `Blm` holds the parameters by name, with values of any loop shape
  `S_p + (K,)`; a grid of parameter sets is one object.
- `Hamiltonian(ion, blm, field)` with the field `(..., 3)`: one batched
  diagonalisation for the whole loop shape. Once solved it holds the energies,
  the eigenvectors and J in the eigenbasis, so the observables don't rotate
  operators themselves.
- Observables take T and Q that broadcast against the loop shape: a grid with
  `None` axes, or a measurement log paired point by point.
- Open: when a `Hamiltonian` is solved and whether it is immutable;
  observables as functions or methods (AUDIT Q29, Q31).
- Memory: J in the eigenbasis is `S + (3, n, n)`, 139 MB per 10⁴ points at
  J = 8; compute only what an observable needs, or work in chunks.
- The print of the level scheme (`summary`, in v1 today): its grouping of
  levels follows CONV-0008 once accepted, and the basis it prints inside a
  degenerate level follows CONV-0011: J_z diagonalised inside the level,
  Kramers partners by time reversal.

## Fitting

Sketched in `examples/fit_cef_mock.py` (a sketch, written on the
rejected bare-array objects).

- One `model(blm)` for every data set (peak energies and intensities, χ(T),
  M(H)), each with its own trailing data axis, so the same function serves a
  single parameter set and a whole grid.
- A vectorised scan over a parameter grid for start values (9261 sets in 6 s
  for Ce³⁺ with three parameters), then a local fit with lmfit from the best
  point. lmfit would be an optional extra, not a dependency of the core.
- Constraints between parameters as lmfit expressions, e.g. cubic
  `B44 = 5 B40`.
- An analytic Jacobian: `∂E_n/∂B_lm = <n|O_lm|n>` (Hellmann–Feynman) from the
  same eigensystem; the derivatives of χ and M need more work.
- Whole spectra with a resolution function instead of peak lists.
- Several sites or ions in one fit.

## Neutron spectra

- `convolve_spectrum(hamiltonian, E, Q, resolution)`, the user-facing
  spectrum on an energy grid E (meV), built from the cross-section weights.
- `k_f/k_i` at the grid energy instead of at each ε_k: exact for a broadened
  line, where the factor varies across it. Moves `Ei` to
  `convolve_spectrum`, with NaN on grid points above E_i; needs a draft CONV
  superseding CONV-0007 and CONV-0014 (the weights without `k_f/k_i`).
- Fixed-E_f kinematics (triple-axis, indirect geometry),
  `k_f/k_i = √(E_f/(E_f + ħω))`: exactly one of `Ei=` or `Ef=`; needs a CONV
  superseding CONV-0007.
- A side-by-side table of the CONV-0007 formula against the functions,
  showing which factor each layer applies itself and which it takes from
  below: `convolve_spectrum` → `cross_section_powder` / `_crystal` →
  `_cross_section_unprojected` → `scattering_function`.

## Physics

- Intermediate coupling: the full fⁿ basis, J-mixing (Sm³⁺, Eu³⁺). Today the
  basis is the ground J multiplet (CONV-0005).
- Other parameter types (`Alm`, `Llm`, Wybourne `Dlm`, ...) converted with
  McPhase's Tables 5 and 6 (copied in CONV-0001); today `Blm` only.
- Rotation of the parameters between the crystal frame and the local frame
  (CONV-0004).
- Stevens operators generated from one closed formula, with the hand-written
  table kept as a test fixture (AUDIT A6).
- The g-tensor of a non-Kramers ground doublet.

## Symmetry

- Point-group support: the allowed parameters of a site and the relations
  between them (cubic `B44 = 5 B40`, ...), in the setting of CONV-0004
  (AUDIT A7).
- Exact degeneracies by construction. A q-fold axis along z couples only m
  and m ± q, so H splits into blocks that are diagonalised separately. For a
  Kramers ion, the partner of |a⟩ is Θ|a⟩, built by time reversal, so the
  doublets are exactly paired, not paired up to rounding (CONV-0008,
  CONV-0011).
- Levels labelled by the irreducible representations of the site group
  (Γ6, Γ7, Γ8, ...), computed with the rotations `exp(−iθ n̂·J)` of the
  group's operations on the multiplet; half-integer J needs the double
  group, and time reversal gives the Kramers pairing.
  - A level is degenerate because its irrep is, not because two energies
    agree within a tolerance; only accidental degeneracies still need
    CONV-0008.
  - The partner functions of each irrep, with a standard phase
    convention, give the basis inside a level for Kramers and non-Kramers
    ions alike. This replaces CONV-0011, which leaves the degenerate
    levels of integer-J ions unpaired until then (AUDIT Q46).
  - Selection rules: transitions that vanish by symmetry come out exactly
    zero, not numerically small.

## Crystal structure and point charges

- Take a `Crystal` from `qespace_crystallography`: the lattice, the site of
  the magnetic ion and its neighbours. It gives the site symmetry (above),
  the local frame and its rotation from the crystal frame (methods in the
  CONV-0004 draft at `d19e90e`), and
  fields and Q given in crystal coordinates. Whether `qespace_cef` may depend
  on `qespace_crystallography` is a question for the suite's layer rules.
- Point-charge model: the parameters from the charges of the neighbours,
  `A_lm ∝ Σ_j q_j Z_lm(r̂_j) / r_j^(l+1)`, then `B_lm = A_lm <r^l> θ_l` and
  the other types of CONV-0001 (Table 5). It needs the radial integrals
  `<r^l>` and, if used, shielding factors, each with its source, and it gives
  starting values for a fit rather than final parameters.

## Validation against reference codes

Cross-checks with McPhase, PyCrystalField and Mantid CrystalField: the same
ion and parameters, compared level by level and for the observables. Each code
has its own conventions and bugs, so each disagreement is recorded as a strict
`xfail` with the reason until the maintainer decides which side is right. The
codes would live in a test dependency group, not in the core.

On 2026-11-20 the Stevens operators of the package were compared by hand with
McPhase's table (manual node132): 601 matrices, l = 1–6 (except l = 5, m = 0),
J = 1/2–8, no difference. The even ranks are now tests (`test_conv_0001_*`);
against PyCrystalField the operators with l + m odd have the opposite sign,
recorded as a strict `xfail`.
