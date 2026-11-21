---
id: CONV-0001
title: CEF parameters and Stevens operators as in McPhase
status: accepted
owner: qespace_cef
supersedes:
---

## Statement

We follow the representation of crystal field parameters from McPhase:
https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node130.html
The default representation is Stevens, with parameter labelling `Blm`:
`H_CEF = Σ_lm B_lm O_lm(J)`, with `O_lm` the Stevens operators of the ground
J multiplet as tabulated by McPhase (below) and `B_lm` in the package's energy
unit. Parameters with m < 0 (the sine-type operators) are labelled `B4m3`.

As of today the package works with `Blm` and `O_lm` only. Any future
implementation of other parameter types (the other rows of Table 5) will be a
matter of further review.

## Definitions (from the McPhase manual)

Copied from the McPhase manual by M. Rotter (pages dated 2013-09-19):
node130 (conventions, Tables 5 and 6), node131 (tesseral harmonics), node132
(Stevens operators). Read on 2026-11-20.

**Stevens operators** (node132), with `X = J(J+1)` and `J± = J_x ± i J_y`.
For m > 0, with a polynomial `P_lm(J_z, X)`:

    O_l0      = P_l0
    O_lm      =  1/4 [ (J+^m + J-^m) P_lm + P_lm (J+^m + J-^m) ]
    O_l,-m    = -i/4 [ (J+^m - J-^m) P_lm + P_lm (J+^m - J-^m) ]

McPhase writes `O_l,-m` as `O^s_lm`. For m = l, `P_ll = 1`, so
`O_ll = (J+^l + J-^l)/2`. The polynomials for the even ranks:

| l | m | `P_lm` |
|---|---|---|
| 2 | 0 | `3 J_z² − X` |
| 2 | 1 | `J_z` |
| 2 | 2 | `1` |
| 4 | 0 | `35 J_z⁴ − (30X − 25) J_z² + 3X² − 6X` |
| 4 | 1 | `7 J_z³ − (3X + 1) J_z` |
| 4 | 2 | `7 J_z² − X − 5` |
| 4 | 3 | `J_z` |
| 4 | 4 | `1` |
| 6 | 0 | `231 J_z⁶ − (315X − 735) J_z⁴ + (105X² − 525X + 294) J_z² − 5X³ + 40X² − 60X` |
| 6 | 1 | `33 J_z⁵ − (30X − 15) J_z³ + (5X² − 10X + 12) J_z` |
| 6 | 2 | `33 J_z⁴ − (18X + 123) J_z² + X² + 10X + 102` |
| 6 | 3 | `11 J_z³ − (3X + 59) J_z` |
| 6 | 4 | `11 J_z² − X − 38` |
| 6 | 5 | `J_z` |
| 6 | 6 | `1` |

**Hamiltonians** (node130). Stevens:
`H = Σ_lm A_lm <r^l> <J||θ_l||J> O_lm(J) = Σ_lm B_l^m O_lm(J)`, with the
Stevens factors `θ_l = α_J, β_J, γ_J` for l = 2, 4, 6. Wybourne:
`H = Σ_lm D_l^m C_lm`, with `C_lm` the tensor operators of
`C_lm(θ, φ) = √(4π/(2l+1)) Y_lm(θ, φ)`. McPhase uses their Hermitian
combinations `T_l0 = C_l0`, `T_l,±|m| = √±1 [C_l,−|m| ± (−1)^|m| C_l,|m|]`,
`H = Σ_lm L_l^m T_lm`: real parameters with the Wybourne normalisation.

**Table 5: parameter notations** (node130). `half(m) = 1/2` for m ≠ 0 and
`1` for m = 0. `lm` is replaced by numbers, e.g. `A20`.

| Literature | Name in McPhase | Relation to `L_l^m` | Normalisation |
|---|---|---|---|
| `A_l^m` in Hutchings 1964, Kassman 1970, Abragam 1970 | `Alm` | `λ_lm L_l^m / <r^l>` | Stevens |
| `B_l^m` in Hutchings 1964, Jensen 1991; `B^l_m` in Mulak 2000, FOCUS, AMOS | `Blm` | `λ_lm L_l^m <J‖θ_l‖J>` (rare earths), `λ_lm L_l^m <L‖θ_l‖L>` (transition metals) | Stevens |
| `V_l^m` in Elliott 1953 | | `λ_lm L_l^m` | Stevens |
| | `Wlm` | `half(m) λ_lm L_l^m / <r^l>` | Stevens |
| | `Vlm` | `half(m) λ_lm L_l^m <J‖θ_l‖J>` | Stevens |
| | `Llm` | `L_l^m` | Wybourne |
| `B^l_0`, `B^l_m(c)`, `B^l_m(s)` in SPECTRE | | `L_l^0`, `L_l^m`, `−L_l^−m` (m > 0) | Wybourne |
| `B^l_m` in Wybourne 1965, Kassman 1970; `B_lm` in Mulak 2000; `Ak()`, `CAk()` in XTLS | `Dlm` | `L_l^0` (m = 0); `(−1)^m (L_l^m − i L_l^−m)` (m > 0); `L_l^−m + i L_l^m` (m < 0) | Wybourne |
| `B^l_m` in Newman 2000 | | `(−1)^m L_l^m` | Wybourne |

**Table 6: `λ_lm`**, the ratios of Stevens to real Wybourne parameters,
`λ_lm ≡ B_l^m / (θ_l L_l^m) = A_lm <r^l> / L_l^m`, the same for ±m:

| l \ \|m\| | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| 0 | 1 | | | | | | |
| 2 | 1/2 | √6 | √6/2 | | | | |
| 4 | 1/8 | √5/2 | √10/4 | √35/2 | √70/8 | | |
| 6 | 1/16 | √42/8 | √105/16 | √105/8 | 3√14/16 | 3√77/8 | √231/16 |

## Rationale

To pay respects to McPhase and base somewhere historically. McPhase is a
long-standing reference implementation, and its manual documents both the
Stevens operators and the conversions to the other parameter types, so a
published parameter set can be checked against it. Stevens operator
equivalents are the natural form inside one J multiplet and are what the
neutron-scattering codes (McPhase, Mantid, PyCrystalField) and most of the
rare-earth magnetism literature use. The Wybourne form on spherical tensors is
independent of the multiplet and comparable between ions; it is needed for
intermediate coupling and common in optical spectroscopy. The normalisation of
`O_lm` fixes the meaning of every `B_lm`: another table gives other numbers, and
a sign difference in one operator flips the sign of its parameter.

References: K. W. H. Stevens, Proc. Phys. Soc. A 65, 209 (1952),
doi:10.1088/0370-1298/65/3/308; M. T. Hutchings, Solid State Phys. 16, 227
(1964), doi:10.1016/S0081-1947(08)60517-2.

## Consequences

- `B_lm` contains the Stevens factor θ_l of the ion, so a parameter set belongs
  to one ion and one multiplet and is not transferable between ions.
- A test compares every operator with the table above, signs included. A
  reference code that differs by a sign has the opposite sign of that
  parameter; a comparison records the case instead of accepting ±.
- Another parameter type, if accepted, is converted once at the boundary with
  Tables 5 and 6, e.g. `B_lm = λ_lm θ_l L_l^m`.
- McPhase input files write the sine-type parameters with an `S` suffix
  (`B22S` for `B_2^−2`, manual node123), so reading them needs a name mapping.
