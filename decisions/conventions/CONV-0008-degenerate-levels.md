---
id: CONV-0008
title: Degenerate levels
status: draft
owner: qespace_cef
---

## Statement

Two energies belong to the same level when

    |E_m − E_n| ≤ rtol · max(max_k |E_k|, E_floor)

with `E_k` the energies of CONV-0010, and `rtol` and `E_floor` read at call
time from the package's one config place. Levels are formed from the sorted
energies, split wherever two consecutive energies differ by more than that.

Results do not depend on the basis inside a level: they use projectors onto
the level, traces and sums over it, never single eigenvectors. A result with
a limit as two levels merge comes from a formula that is continuous through
the degeneracy and evaluated without cancellation; the tolerance then only
picks between equal branches of it. Grouping into levels is for counting and
reporting: degeneracies, level lists, per-level intensities, the ground
doublet of the g-tensor.

How a basis inside a level is shown is CONV-0011.

`OPEN:` the defaults, proposed `rtol = 1e-10` and `E_floor = 1` meV (AUDIT
Q47).

## Rationale

The solver returns one orthonormal basis of a degenerate level out of
infinitely many, chosen by rounding, so it can change between numpy
versions or platforms. Anything computed from one vector, e.g. `<a|J_z|a>` in
a Kramers doublet, can take any value between the eigenvalues of J_z in the
doublet; the projector and sums over the level are fixed. The g-tensor is an
example: g depends on the basis, `G = g gᵀ` does not.

The tolerance follows the accuracy of the symmetric eigensolver,
`|λ̂ − λ| ≲ p(n) ε ‖H‖` (LAPACK); with the unshifted energies of CONV-0010,
`max_k |E_k|` is that size in every field, and `E_floor` sets it when H = 0.
In 1200 low-symmetry Kramers Hamiltonians `eigh` split the exact doublets by
at most 13 ε `max|E|` (measured on `c94cb46`), so `rtol = 1e-10` is 10⁵ above
rounding and far below a physical splitting.

Continuous formulas (divided differences, as in linear response) avoid
depending on the threshold, but as a plain difference they cancel: the van
Vleck weight `(p_n − p_m)/(E_m − E_n)` loses 5·10⁻⁵ relative at a splitting
of 10⁻¹² meV, while the form with `expm1` stays exact.

## Consequences

- The package gets one config place for `rtol` and `E_floor` (AUDIT A11).
- Level grouping uses the tolerance, never rounding to decimals.
- CONV-0009 fixes the degeneracy of the susceptibility to
  `constants.numerical_zero`; once this CONV is accepted, a new CONV
  supersedes CONV-0009 to use this rule.
- Tests:
  - every exact Kramers doublet of low-symmetry Hamiltonians is one level;
  - grouping is unchanged when H gets a random Hermitian perturbation of
    size ε `max|E|`;
  - results are unchanged under a unitary rotation inside each level;
  - a result with a limit at the degeneracy (χ) is continuous as a splitting
    shrinks through the tolerance.
