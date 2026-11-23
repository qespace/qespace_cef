---
id: CONV-0015
title: Degenerate levels
status: draft
owner: qespace_cef
supersedes: CONV-0008
superseded_by:
---

## Statement

Two energies belong to the same level when

    |E_m − E_n| ≤ rtol · max(max_k E_k − min_k E_k, E_floor)

with `E_k` the energies of CONV-0014, and the defaults `rtol = 1e-10` and
`E_floor = 1` meV. Levels are formed from the sorted energies, split
wherever two consecutive energies differ by more than that. The scale is
the spread of the energies, so the grouping doesn't change when all
energies shift together (B00, CONV-0014).

Results do not depend on the basis inside a level: they use projectors onto
the level, traces and sums over it, never single eigenvectors. A result with
a limit as two levels merge comes from a formula that is continuous through
the degeneracy and evaluated without cancellation; the tolerance then only
picks between equal branches of it. Grouping into levels is for counting and
reporting: degeneracies, level lists, per-level intensities, the ground
doublet of the g-tensor.

How a basis inside a level is shown is CONV-0011.

## Rationale

This supersedes CONV-0008, whose scale `max_k |E_k|` includes the zero of
energy. With `B20 = 5·10⁻¹¹` meV, Ce³⁺ has three doublets; adding
`B00 = 10` meV merged them into one sixfold level and changed χ, although
every gap was unchanged (independent review `AUDIT_GPT.md`, M5, reproduced
on `9b49c3d`). CONV-0014 makes every observable independent of the zero,
so the grouping must be too. With `B00 = 0` the energies average to zero
and the spread lies between `max|E|` and `2 max|E|`, so the margin below is
unchanged.

The solver returns one orthonormal basis of a degenerate level out of
infinitely many, chosen by rounding, so it can change between numpy
versions or platforms. Anything computed from one vector, e.g. `<a|J_z|a>` in
a Kramers doublet, can take any value between the eigenvalues of J_z in the
doublet; the projector and sums over the level are fixed. The g-tensor is an
example: g depends on the basis, `G = g gᵀ` does not.

The tolerance follows the accuracy of the symmetric eigensolver,
`|λ̂ − λ| ≲ p(n) ε ‖H‖` (LAPACK). Removing the mean of the energies (the
trace of H) before solving keeps that error independent of the zero as
well. In 1200 low-symmetry Kramers Hamiltonians `eigh` split the exact
doublets by at most 13 ε `max|E|` (measured on `c94cb46`), so
`rtol = 1e-10` is 10⁵ above rounding and far below a physical splitting.

Continuous formulas (divided differences, as in linear response) avoid
depending on the threshold, but as a plain difference they cancel: the van
Vleck weight `(p_n − p_m)/(E_m − E_n)` loses 5·10⁻⁵ relative at a splitting
of 10⁻¹² meV, while the form with `expm1` stays exact.

## Consequences

- Tests, each naming CONV-0015:
  - every exact Kramers doublet of low-symmetry Hamiltonians is one level;
  - grouping is unchanged when H gets a random Hermitian perturbation of
    size ε times the spread;
  - grouping and χ are unchanged when `B00` shifts all energies;
  - results are unchanged under a unitary rotation inside each level;
  - a result with a limit at the degeneracy (χ) is continuous as a splitting
    shrinks through the tolerance.
- The tests of CONV-0008 carry over under this ID.
