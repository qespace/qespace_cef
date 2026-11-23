---
id: CONV-0014
title: Energies, excitation energies and spectral weights
status: accepted
owner: qespace_cef
supersedes: CONV-0010, CONV-0012
superseded_by:
---

## Statement

**Energies.** The energies `E_n` are the eigenvalues of the Hamiltonian as
written, `ℋ = Σ_lm B_lm O_lm + ℋ_Z` (CONV-0006), in meV, not shifted. Every
`O_lm` with l ≥ 1 and the Zeeman term are traceless, so with `B00 = 0` the
energies average to zero: the zero of energy is the centre of gravity of
the ground J multiplet, the free-ion energy, in every field. `B00` (`O00`
is the identity) shifts all energies by `B00`. The internal energy U and
the free energy F are measured from the same zero. Populations, entropy,
heat capacity, magnetisation, susceptibility and spectral weights do not
depend on it.

**Excitation energies.** For every ordered pair of eigenstates (i, f),
initial i and final f, the excitation energy is `ε_if = E_f − E_i`, in meV;
`ε_if > 0` is neutron energy loss (CONV-0007). They form a matrix with
element (i, f) (CONV-0002): `ε_ii = 0` and `ε_fi = −ε_if`. The excitations
from the ground state, `Δ_n = E_n − E_0`, are its row of the lowest state.

**Spectra as excitations.** A spectrum computed here is a sum of δ-peaks
and is given in that discrete form, never as a continuous function of
energy:

    S^{αβ}(ħω) = Σ_{i,f} W_if^{αβ} δ(ħω − ε_if)

- the **scattering weights** `W_if^{αβ} = p_i <i|J_α|f><f|J_β|i>`, with
  `p_i` the populations: dimensionless, a 3×3 matrix per pair;
- the **cross-section weights** `w_if`, the coefficients of the δ-peaks in
  the cross section of CONV-0007: barn sr⁻¹ per ion. The differential
  cross section dσ/dΩ of an energy window is the sum of the cross-section
  weights of the pairs whose ε_if lie inside it.

There is one entry per ordered pair, elastic ones (i = f, and pairs inside
a level) included. A pair names two eigenstates; which are degenerate
follows CONV-0008. Inside a degenerate level single entries depend on the
basis, so only sums over the entries at one energy are meaningful.

A continuous spectrum needs a line shape L (resolution, lifetime),
`S(ħω) = Σ_if W_if L(ħω − ε_if)`. L belongs to the instrument and the
sample and is not part of this convention.

## Rationale

This supersedes CONV-0010 and CONV-0012. CONV-0010 gave the name
`excitation_energies` to `E_n − E_0`, while the excitation energies of the
spectra are all the pairs `E_f − E_i`; one quantity per name (ADR-0006)
needs one definition. CONV-0012 gave the excitations as an unordered set
without level information, a misunderstanding of the maintainer's intent
(AUDIT A16): the pairs keep their states, so a weight traces back to its
levels. Both also named code or array shapes, which are ADR-0006's and
ADR-0007's (ADR-0004). The zero of energy and the weights are unchanged.

The zero above is fixed by the Hamiltonian itself, the same for every ion,
parameter set and field, and needs no choice. Shifting the ground state to
0 is a different zero at every field and every point of a parameter grid:

- In a field, `E_0(B) ≡ 0` hides the Zeeman energy of the ground state, so
  an `E(B)` diagram is wrong, and F computed from shifted energies breaks
  `M = −∂F/∂B`.
- On a grid of fields or parameters, energies shifted point by point cannot
  be compared or differentiated across the grid.
- Hellmann–Feynman, `∂E_n/∂B_lm = <n|O_lm|n>`, holds for the eigenvalues.
  For shifted energies it gains `−<0|O_lm|0>`. An analytic Jacobian for
  fits rests on it.

Level schemes and neutron spectra are read from the ground state; that is
the row `Δ_n` of the excitation energies, not a shift of the energies.

For discrete excitations the pairs (ε_if, W_if) are the scattering function
exactly; a continuous S(ħω) would need a line shape, which belongs to the
instrument and the sample's lifetimes, not to the crystal field. One entry
per pair of states gives every Hamiltonian the same number of entries, so
spectra broadcast over parameter and Q–E grids without grouping levels.

Checked on `c94cb46` (CONV-0010): `|Tr O_lm| = 0` for every implemented
l ≥ 1 operator at J = 5/2 and 8, and `O00` is the identity. Ce³⁺ with
`B20 = 0.5, B40 = 0.01, B44 = 0.1` meV has eigenvalues −3.584, −2.800 and
+6.384 meV (each twice), with mean 0 to 10⁻¹². Free Ho³⁺ at 2 T along z:
`E_0 = −gJ μB J B` = −1.15768 meV exactly, mean 0. Low-symmetry Ce³⁺ in a
field: `<n|O_40|n>` equals the finite difference `∂E_n/∂B40` to
1.4·10⁻⁷ meV/meV (largest value 60.5); from shifted energies the difference
is up to 23.5.

## Consequences

- Tests, each naming CONV-0014, comparing sums over equal energies or
  spectra after a line shape, never single entries inside a level:
  - the energies are the eigenvalues of the Hamiltonian, and with
    `B00 = 0` they average to zero, with and without a field;
  - for the free ion in a field B along z, `E_0 = −gJ μB J |B|`;
  - `ε_if = E_f − E_i`, `ε_ii = 0`, `ε_fi = −ε_if`, and the row of the
    ground state starts at 0;
  - `∂E_n/∂B_lm = <n|O_lm|n>` for a non-degenerate level, against a
    finite difference;
  - the sum rule `Σ_{i,f} Σ_α W_if^{αα} = J(J+1)`;
  - detailed balance `W_fi^{αβ} = e^{−ε_if/k_BT} W_if^{βα}`;
  - the summed weights at each energy are unchanged under a unitary
    rotation inside each level.
- The tests of CONV-0010 and CONV-0012 carry over under this ID.
