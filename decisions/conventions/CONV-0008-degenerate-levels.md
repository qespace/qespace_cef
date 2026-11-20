---
id: CONV-0008
title: Degenerate levels
status: draft
owner: qespace_cef
---

## Statement

Eigenvectors within a degenerate level are basis dependent: results and tests
use quantities that don't depend on that basis (projectors onto the level,
traces, sums over the level), never single eigenvectors.

Two energies belong to the same level when

    |E_m − E_n| ≤ rtol · max(max_k |E_k|, E_floor)

with `E_k` the energies of CONV-0010 and `rtol`, `E_floor` read at call time
from the package's one config place. The proposed defaults are `rtol =
1e-10` and `E_floor = 1` meV; `E_floor` sets the scale when H = 0 (the free
ion). Levels are formed from the sorted energies, splitting wherever two
consecutive energies differ by more than that tolerance.

A result that has a limit as two levels merge is computed from a formula
that is continuous through the degeneracy and evaluated without
cancellation (e.g. the van Vleck weight of CONV-0009). The tolerance then
only chooses between equal branches of the formula. Grouping into levels is
for counting and reporting: degeneracies, level lists, per-level
intensities, the ground doublet of the g-tensor.

Where a basis inside a level is shown (the eigenvectors of `summary`), it is
the canonical basis below. Numbers the package returns stay basis invariant.

1. Diagonalise J_z (in the frame of the J operators, CONV-0004) inside the
   level, and order its eigenstates by `<J_z>`, ascending.
2. Kramers ion (half-integer J) in zero field: every level is made of
   Kramers pairs. For each eigenstate |μ> with `<J_z> < 0`, fix its phase so
   that its largest component is real and positive; its partner is `Θ|μ>`,
   with time reversal `Θ|m> = (−1)^(J−m) |−m>` followed by complex
   conjugation, and has `<J_z> > 0`.
3. Every other case (integer J, or a field): each eigenstate gets its own
   phase, largest component real and positive.

`OPEN:` the operator when J_z is itself degenerate inside a level, and the
phase of a state with `<J_z> = 0`, which `Θ` maps to itself (integer J).
`OPEN:` the defaults `rtol = 1e-10`, `E_floor = 1` meV.

## Rationale

If |a> and |b> have the same energy, so does every combination
`α|a> + β|b>`. The solver returns one orthonormal pair out of infinitely many,
chosen by rounding, so it can change between numpy versions or platforms.
Anything computed from one vector of the pair, e.g. `<a|J_z|a>` in a Kramers
doublet, can then take any value between −j and +j, where ±j are the
eigenvalues of J_z within the doublet. What is fixed is the level as a
whole: the projector `P = |a><a| + |b><b|`, and sums
such as `Σ_{i,j ∈ level} |<i|J_α|j>|²`. The g-tensor is an example: g depends
on the basis inside the doublet, `G = g gᵀ` does not. In a small field the
degeneracy lifts and the eigenvectors become well defined, but they jump as
the field direction changes.

The tolerance follows the accuracy of the symmetric eigensolver,
`|λ̂ − λ| ≲ p(n) ε ‖H‖` (LAPACK), so it is relative to the size of H, and
with the unshifted energies of CONV-0010, `max_k |E_k|` is that size, the
same in every field. Measured on `c94cb46`: in 1200 low-symmetry Kramers
Hamiltonians (Ce, Nd, Sm, Dy, Er, Yb), `eigh` split the exact doublets by at
most 13 ε `max|E|` ≈ 3·10⁻¹⁵ `max|E|`. `rtol = 1e-10` is 10⁵ above that and
far below any physical splitting (10⁻⁹ meV on a 10 meV scale).

Formulas continuous through the degeneracy are the standard way of linear
response and matrix-function derivatives (divided differences, the
Daleckii–Krein formula) to avoid depending on a threshold. Evaluated as a
plain difference they cancel: the van Vleck weight `(p_n − p_m)/(E_m − E_n)`
loses 5·10⁻⁵ relative at `E_m − E_n = 10⁻¹²` meV and 10⁻³ at 10⁻¹⁴ meV
(k_B T = 0.4 meV), while `p_n (1 − e^{−(E_m − E_n)/k_B T})/(E_m − E_n)`,
with `expm1`, stays exact.

Grouping consecutive sorted energies can chain a slow drift into one level;
with a tolerance this close to rounding it does not happen in practice.

The canonical basis is degenerate perturbation theory for an infinitesimal
field along z: diagonalising J_z inside the level gives the B → 0 limit of
the states in a field along z, without a fake field in H. Diagonalising
leaves one phase per state free, not one common phase: a common phase
cancels in every `<i|O|j>`, while the relative phase of two states rotates
the off-diagonal elements. For Ce³⁺ (`B20 = 0.5, B40 = 0.01, B44 = 0.1` meV)
the ground doublet has `<−|J_x|+> = −0.602` and `<−|J_y|+> = +0.602 i` with
`|+> = Θ|−>`; an extra relative phase `e^{iπ/3}` turns them into
`−0.301 − 0.521 i` and `−0.521 + 0.301 i`, with the same modulus.
Taking the partner from time reversal fixes the relative phase physically:
the pair then transforms like a spin ½. Two different, equally valid outputs
of `eigh` for that doublet gave the identical canonical states (scratch run,
2026-11-21).

## Consequences

- Observables sum over each degenerate level; tests compare projectors or sums,
  or rotate the basis inside the level and check that the result is unchanged.
- Level grouping uses the tolerance, never rounding to a number of decimals.
- The package gets one config place for `rtol` and `E_floor`, read at call
  time (AUDIT A11).
- CONV-0009 fixes the degeneracy of the susceptibility to the absolute
  `constants.numerical_zero`. Once this CONV is accepted, a new CONV
  supersedes CONV-0009 so that it uses this rule.
- `summary` groups levels with this tolerance and prints the canonical
  basis.
- Tests:
  - the canonical basis is the same after a random unitary rotation inside
    each level, and `Θ|−> = |+>` for every Kramers pair;
  - grouping is unchanged when H gets a random Hermitian perturbation of
    size ε `max|E|`;
  - every exact Kramers doublet of low-symmetry Hamiltonians is one level;
  - results are unchanged under a unitary rotation inside each level;
  - a result with a limit at the degeneracy (χ) is continuous as a
    splitting shrinks through the tolerance.
