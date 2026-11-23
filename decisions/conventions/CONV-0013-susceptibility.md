---
id: CONV-0013
title: Single-ion magnetic susceptibility
status: accepted
owner: qespace_cef
supersedes: CONV-0009
superseded_by:
---

## Statement

The single-ion susceptibility is the linear response of the moment to the
field at zero field, `χ_αβ = ∂<m_α>/∂B_β` at B = 0, with the moment
operator `m = −gJ μB J` and the Zeeman term `−m·B` of CONV-0006.
Equivalently, `χ_αβ = −∂²F/∂B_α∂B_β` at B = 0, with `F = −k_B T ln Z`.
For the eigenstates |n> of `H_CEF`, with energies `E_n` and populations
`p_n = exp(−E_n/k_B T)/Z`:

    χ_αβ = gJ² μB² Σ_n Σ_m  W_nm  Re(<n|J_α|m><m|J_β|n>)

    W_nm = p_n / (k_B T)               if E_m = E_n   (Curie)
    W_nm = (p_n − p_m) / (E_m − E_n)   if E_m ≠ E_n   (van Vleck)

Both sums run over all 2J + 1 states, so each pair (n, m) appears in both
orders. The susceptibility along a unit vector n̂ is
`χ_n̂ = n̂·χ·n̂ = gJ² μB² Σ_nm W_nm |<n|n̂·J|m>|²`, and the powder average is
`Tr χ / 3`.

Two levels count as degenerate according to the rules of CONV-0008.

χ is in μB/T per ion (CONV-0003): α is the component of the moment, β that
of the field. With energies in meV and μB in meV/T it is
`gJ² μB Σ_nm W_nm Re(<n|J_α|m><m|J_β|n>)`. `χ_αβ` is real and symmetric.

## Rationale

This supersedes CONV-0009 for two reasons. Its moment operator,
`m = +gJ μB J`, is the moment-parallel convention, while CONV-0006 takes
the electron's sign, `m = −gJ μB J`. χ is quadratic in that sign, so every
value of CONV-0009 stays. And CONV-0009 said how the code returns χ (array
shapes, what the caller computes), which is an ADR's matter (ADR-0004).

This is Van Vleck's formula. In second-order perturbation theory in the
Zeeman term, the first-order shifts within a degenerate level give the Curie
term, and the second-order shifts,
`gJ² μB² Σ_{E_m ≠ E_n} |<n|B·J|m>|² / (E_n − E_m)`, give the van Vleck
term. Weighting them with `p_n` and making the sum symmetric in (n, m) gives
`W_nm`. The Curie term has no `−<m>²` part, because `<J> = 0` at
B = 0: `H_CEF` is invariant under time reversal.

Jensen and Mackintosh give the same formula for the susceptibility of ⟨J⟩
to an effective field in energy units, `χ°_αβ = ∂<J_α>/∂h_β` with the
Zeeman term `−J·h` (Eq. 2.1.18, below). With CONV-0006, `h = −gJ μB B` and
`m = −gJ μB J`, so `χ_αβ = gJ² μB² χ°_αβ`. Their states p, p′ are our n, m.
Their van Vleck weight `(n_p − n_p′)/(E_p′ − E_p)` and Curie term `β n_p`
are our `W_nm`. They add `−β <J_α><J_β>` for a state with `<J> ≠ 0`, as in
a field or a mean field; at B = 0 without exchange it is zero, so it is not
in the Statement.

The literature often writes the van Vleck term over the pairs `E_m > E_n`
only, with a factor 2. That is the same sum; leaving out the factor halves
the van Vleck term.

Both terms are basis independent inside a degenerate level: the Curie term is
`Tr(P J_n̂ P J_n̂)` over that level, with P its projector. `W_nm ≥ 0` for every
pair, so χ is positive semi-definite. The formula needs only the zero-field
eigensystem. A finite difference `n̂·<m>(ε n̂)/ε` gives the same value as
ε → 0 and checks it.

Scope: one ion, the ground J multiplet, no exchange, no diamagnetic core term,
no demagnetisation. The van Vleck term of the excited J multiplets is
missing; it matters for Sm³⁺ and Eu³⁺.

References:
- J. H. Van Vleck, "On Dielectric Constants and Magnetic Susceptibilities in
  the New Quantum Mechanics Part I. A General Proof of the Langevin-Debye
  Formula", Phys. Rev. 29, 727 (1927), doi:10.1103/PhysRev.29.727.
- J. H. Van Vleck, "On Dielectric Constants and Magnetic Susceptibilities in
  the New Quantum Mechanics Part III—Application to Dia- and Paramagnetism",
  Phys. Rev. 31, 587 (1928), doi:10.1103/PhysRev.31.587.
- J. H. Van Vleck, *The Theory of Electric and Magnetic Susceptibilities*
  (Clarendon Press, Oxford, 1932).
- J. Jensen and A. R. Mackintosh, *Rare Earth Magnetism: Structures and
  Excitations* (Clarendon Press, Oxford, 1991), Section 2.1, Eq. (2.1.18),
  p. 76.

The maintainer checked the formula on 2026-11-20 against Van Vleck's 1932
book and against Eq. (2.1.18) of Jensen and Mackintosh. The two papers
exist (checked on the APS pages the same day), but their equations were not
compared.

## Consequences

- Tests in μB/T per ion, as for CONV-0009, each naming CONV-0013:
  - the free ion (`H_CEF = 0`) gives the Curie law,
    `χ_αβ = δ_αβ gJ² μB² J(J+1) / (3 k_B T)`;
  - at high T, `Tr χ / 3` approaches that Curie law;
  - the formula and the finite difference with the moment of CONV-0006
    agree, which checks the signs of the Zeeman term and the moment
    together;
  - χ is symmetric and positive semi-definite;
  - χ doesn't change under a unitary rotation inside a degenerate level;
  - in a tetragonal field, `χ_xx = χ_yy` and `χ_xy = 0`.
- Every pair (n, m) with `E_m ≠ E_n` contributes, in both orders.
- A susceptibility with exchange (mean field) or from more than one multiplet
  is a different quantity and needs its own CONV.
- How the code returns χ (shape, loop shape) belongs in an ADR.
