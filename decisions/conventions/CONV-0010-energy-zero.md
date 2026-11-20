---
id: CONV-0010
title: Zero of energy and excitation energies
status: accepted
owner: qespace_cef
---

## Statement

The energies `E_n` are the eigenvalues of the Hamiltonian as written,
`H = Σ_lm B_lm O_lm + H_Z`, in meV, not shifted. Every `O_lm` with l ≥ 1 and
the Zeeman term are traceless, so with `B00 = 0` the energies average to zero:
the zero of energy is the centre of gravity of the ground J multiplet, the
free-ion energy, in every field. `B00` (`O00` is the identity) shifts all
energies by `B00`.

The excitation energies `Δ_n = E_n − E_0`, with `E_0` the lowest energy, are a
separate quantity, `excitation_energies` in the API, next to `energies`. They
are what level schemes and neutron spectra show.

The internal energy U and the free energy F are measured from the same zero
as `E_n`. Populations, entropy, heat capacity, magnetisation, susceptibility
and spectral weights do not depend on the zero.

The Hamiltonian (v1 `CEFion`) has no `str()` of its own. A print of the
level scheme is a separate feature, and what it shows is decided with it.

## Rationale

The zero above is fixed by the Hamiltonian itself, the same for every ion,
parameter set and field, and needs no choice. Shifting the ground state to 0
is a different zero at every field and every point of a parameter grid:

- In a field, `E_0(B) ≡ 0` hides the Zeeman energy of the ground state, so an
  `E(B)` diagram is wrong, and F computed from shifted energies breaks
  `M = −∂F/∂B`.
- On a grid of fields or parameters (broadcasting), energies shifted point by
  point cannot be compared or differentiated across the grid.
- Hellmann–Feynman, `∂E_n/∂B_lm = <n|O_lm|n>`, holds for the eigenvalues. For
  shifted energies it gains `−<0|O_lm|0>`. An analytic Jacobian for fits
  rests on it.
- One name with one meaning: an option that shifts `energies` makes every
  consumer (observables, thermodynamics, fits) depend on how it was set.

Level schemes and neutron spectra are read from the ground state. That is
served by the excitation energies, a separate, named quantity, rather than by
changing what `energies` means.

Checked on `c94cb46`: `|Tr O_lm| = 0` for every implemented l ≥ 1 operator at
J = 5/2 and 8, and `O00` is the identity. Ce³⁺ with `B20 = 0.5, B40 = 0.01,
B44 = 0.1` meV has eigenvalues −3.584, −2.800 and +6.384 meV (each twice),
with mean 0 to 10⁻¹². Free Ho³⁺ at 2 T along z: `E_0 = −gJ μB J B` =
−1.15768 meV exactly, mean 0. Low-symmetry Ce³⁺ in a field: `<n|O_40|n>`
equals the finite difference `∂E_n/∂B40` to 1.4·10⁻⁷ meV/meV (largest value
60.5). From the shifted energies, the difference is up to 23.5.

## Consequences

- `energies` equals the eigenvalues of `hamiltonian`. No option shifts them
  (today: `diagonalize(shiftToZero=True)` shifts them by default).
- `excitation_energies` is provided next to `energies` and starts at 0.
- `CEFion.__str__` goes; `get_degeneracies` groups `energies`.
- `thermodynamics` gives U and F from this zero (with the fix of its units).
- Tests:
  - `energies` are the eigenvalues of `hamiltonian`, and with `B00 = 0` they
    average to zero, with and without a field;
  - for the free ion in a field B along z, `E_0 = −gJ μB J |B|`;
  - `excitation_energies` is `energies − energies[0]` and starts at 0;
  - `∂E_n/∂B_lm = <n|O_lm|n>` for a non-degenerate level, against a finite
    difference.
