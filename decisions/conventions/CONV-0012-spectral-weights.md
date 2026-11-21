---
id: CONV-0012
title: Spectral weights: spectra as excitations, not functions
status: draft
owner: qespace_cef
supersedes:
---

## Statement

A spectrum computed by the package is a sum of δ-peaks and is returned in
that discrete form, never as a continuous function of energy:

    S^{αβ}(ħω) = Σ_k W_k^{αβ} δ(ħω − ε_k)

- the **excitation energies** `ε_k`, in meV;
- the **spectral weights** `W_k`, the coefficients of the δ-peaks in the
  scattering function: dimensionless, a 3×3 matrix per excitation;
- the **intensities** `I_k`, the coefficients of the δ-peaks in the cross
  section (CONV-0007): barn sr⁻¹ per ion.

An excitation is a pair of levels (a, b) of CONV-0008, elastic pairs a = b
included, with `ε_ab = E_b − E_a` and the weight summed over the states of
both levels:

    W_ab^{αβ} = Σ_{i∈a, f∈b} p_i <i|J_α|f><f|J_β|i>

A continuous spectrum is the caller's: a line shape L (resolution, lifetime)
replaces each δ, `S(ħω) = Σ_k W_k L(ħω − ε_k)`. The package does not choose
L.

## Rationale

For discrete excitations the pair (ε_k, W_k) is the scattering function,
exactly; a continuous S(ħω) would need a line shape, which belongs to the
instrument and the sample's lifetimes, not to the CEF. The sum rule, detailed
balance and basis invariance are statements about the weights. Summing over
the states of each level makes the weights independent of the basis inside
degenerate levels (CONV-0008), which single-state transitions are not. The
spin-wave code returns the same pair, so the convention may become a suite
one.

## Consequences

- Functions return `(energies, weights)`: `energies` of shape (K,), spectral
  weights (K, 3, 3), intensities (K,).
- Tests:
  - the sum rule `Σ_k Σ_α W_k^{αα} = J(J+1)`;
  - detailed balance `W_ba^{αβ} = e^{−ε_ab/k_BT} W_ab^{βα}`;
  - the weights are unchanged under a unitary rotation inside each level.
