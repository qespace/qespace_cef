---
id: CONV-0011
title: Canonical basis inside a degenerate level
status: accepted
owner: qespace_cef
---

## Statement

Where a basis inside a level (CONV-0008) is shown, e.g. the eigenvectors
printed by `summary`, it is this one. Numbers the package returns stay basis
invariant.

1. Diagonalise J_z, in the frame of the J operators (CONV-0004), inside the
   level, and order its eigenstates by `<J_z>`, ascending. Where J_z is
   itself degenerate inside the level (e.g. a doublet with `<J_z> = 0` for
   both states), diagonalise J_x inside that subspace, then J_y.
2. Kramers ion (half-integer J) in zero field: every level is made of
   Kramers pairs. For each eigenstate |μ> with `<J_z> < 0`, fix its phase so
   that its largest component is real and positive; its partner is `Θ|μ>`,
   with time reversal `Θ|m> = (−1)^(J−m) |−m>` followed by complex
   conjugation, and has `<J_z> > 0`.
3. Every other case (integer J, or a field): each eigenstate gets its own
   phase, largest component real and positive. This includes the degenerate
   levels of integer-J ions, which the site symmetry makes degenerate
   (non-Kramers doublets, triplets in a cubic field): their states are not
   paired by time reversal.

Once the site symmetry is implemented, a basis from the irreducible
representations of the site group replaces rules 1–3, by a new CONV that
supersedes this one.

## Rationale

Diagonalising J_z inside a level is degenerate perturbation theory for an
infinitesimal field along z: it gives the B → 0 limit of the states in that
field, without a fake field in H.

It leaves one phase per state free, not one common phase. A common phase
cancels in every `<i|O|j>`; the relative phase of two states rotates the
off-diagonal elements. For Ce³⁺ (`B20 = 0.5, B40 = 0.01, B44 = 0.1` meV) the
ground doublet has `<−|J_x|+> = −0.602` and `<−|J_y|+> = +0.602 i` with
`|+> = Θ|−>`. An extra relative phase `e^{iπ/3}` turns them into
`−0.301 − 0.521 i` and `−0.521 + 0.301 i`: same modulus, different phases.
Taking the partner from time reversal fixes the relative phase physically,
and the pair then transforms like a spin ½.

Two different, equally valid outputs of `eigh` for that doublet gave
identical canonical states (scratch run, 2026-11-21).

Integer-J ions have degenerate levels only where the site symmetry makes
them: Pr³⁺ in a tetragonal field (`B20 = 0.3, B40 = 0.01, B44 = 0.05` meV)
has two doublets, which `B22 = 0.1` splits; Ho³⁺ in a cubic field has
doublets and triplets. A Kramers doublet survives any CEF: Ce³⁺ with the same
`B22` stays three doublets. Inside both Pr doublets, `Θ` maps the state with
`<J_z> < 0` onto the one with `<J_z> > 0` (overlap 1.0), and J_x and J_y have
no matrix elements: they are Ising doublets (scratch run, 2026-11-21).
Pairing them is possible, but whether a level is a symmetry-enforced
degeneracy rests on the tolerance of CONV-0008. The irreducible
representations settle it exactly: the dimension of the irrep, not two
energies within a tolerance, makes a level degenerate, and the partner
functions of the irrep, with a standard phase convention, give the basis.
Pairing integer-J levels is therefore left to that basis.

## Consequences

- `summary` prints the canonical basis, with levels grouped by CONV-0008.
- The basis depends on the frame of the J operators (CONV-0004). The basis
  from the irreducible representations of the site group (`ROADMAP.md`,
  Symmetry) supersedes this CONV when it comes.
- Tests:
  - the canonical basis is the same after a random unitary rotation inside
    each level;
  - `Θ|−> = |+>` for every Kramers pair, and `<J_z>` is ascending inside each
    level.
