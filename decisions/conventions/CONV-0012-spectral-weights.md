---
id: CONV-0012
title: Spectral weights: spectra as excitations, not functions
status: superseded
owner: qespace_cef
supersedes:
superseded_by: CONV-0014
---

> Superseded by CONV-0014 on 2026-11-23: the weights are unchanged; each
> pair (i, f) keeps its states, and the shapes moved to ADR-0007.

## Statement

A spectrum computed by the package is a sum of δ-peaks and is returned in
that discrete form, never as a continuous function of energy:

    S^{αβ}(ħω) = Σ_k W_k^{αβ} δ(ħω − ε_k)

- the **excitation energies** `ε_k`, in meV;
- the **scattering weights** `W_k`, the coefficients of the δ-peaks in the
  scattering function: dimensionless, a 3×3 matrix per excitation;
- the **cross-section weights** `w_k`, the coefficients of the δ-peaks in
  the cross section (CONV-0007): barn sr⁻¹ per ion. The differential
  cross section dσ/dΩ of an energy window is the sum of the
  cross-section weights of the excitations inside it.

The excitations are an unordered set of n² entries, n = 2J + 1, one per pair
of states (i, f), elastic ones included, with `ε_if = E_f − E_i` and
`W_if^{αβ} = p_i <i|J_α|f><f|J_β|i>`. They carry no level information;
entries may share an energy or have zero weight. Only sums over the entries
at one energy are meaningful (CONV-0008). Weights between specific levels
are a separate calculation.

A continuous spectrum is the caller's: a line shape L (resolution, lifetime)
replaces each δ, `S(ħω) = Σ_k W_k L(ħω − ε_k)`. The package does not choose
L.

## Rationale

For discrete excitations the pair (ε_k, W_k) is the scattering function,
exactly; a continuous S(ħω) would need a line shape, which belongs to the
instrument and the sample's lifetimes, not to the CEF. One entry per pair of
states gives every Hamiltonian the same n², so spectra broadcast over
parameter and Q–E grids (ADR-0003) without grouping levels. The spin-wave
code returns the same pair, so the convention may become a suite one.

## Consequences

- Functions return `(energies, weights)`: energies of shape (n²,), scattering
  weights (n², 3, 3), cross-section weights (n²,), after the loop shape of ADR-0003.
- Tests compare sums over equal energies, or spectra after a line shape,
  never single entries:
  - the sum rule `Σ_k Σ_α W_k^{αα} = J(J+1)`;
  - detailed balance `W_fi^{αβ} = e^{−ε_if/k_BT} W_if^{βα}`;
  - the summed weights are unchanged under a unitary rotation inside each
    level.
