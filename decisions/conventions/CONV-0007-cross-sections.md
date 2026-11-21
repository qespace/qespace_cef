---
id: CONV-0007
title: Neutron scattering function and cross section
status: draft
owner: qespace_cef
supersedes:
---

## Statement

The **scattering function** of one ion, with eigenstates |i>, |f>, energies
`E_i`, `E_f` (CONV-0010) and populations `p_i`:

    S^{αβ}(ħω) = Σ_{i,f} p_i <i|J_α|f><f|J_β|i> δ(ħω − (E_f − E_i))

returned as transitions: the energies `E_f − E_i` in meV, positive for
neutron energy loss, and a 3×3 weight per transition. It holds no neutron
factor.

The **cross section** per ion is built from it:

    d²σ/dΩdE = (γr₀)² (k_f/k_i) [½ g_J F(Q)]² e^{−Q²⟨u²⟩/3}
               Σ_{αβ} (δ_αβ − Q̂_α Q̂_β) S^{αβ}(ħω)

with the dipole form factor `F(Q)`, `k_f/k_i = √(1 − ħω/E_i)` and the
mean-square displacement `⟨u²⟩` in Å² (`u2`). A powder averages the
projector to `(2/3) δ_αβ`. Q is in the local frame of CONV-0004.

`OPEN:` intensities in arbitrary units (`(γr₀)² = 1`) or absolute
(barn sr⁻¹ meV⁻¹ per ion). `OPEN:` whether the transitions with
`E_f = E_i` (elastic, `ħω = 0`) are returned.

## Rationale

This is Jensen and Mackintosh, sec. 4.2 and (4.2.2a), and Raymond (1.4), for
a single ion (`REFERENCES.md`). Splitting the scattering function from the
cross section keeps the CEF physics apart from the neutron and instrument
factors. Against the code before this CONV:

| Factor | References | Code | Comment |
|---|---|---|---|
| prefactor | `(γr₀)²` | 1 | arbitrary units |
| ion factor | `[½ g F(Q)]²`; Raymond `f²(Q)` (spin only, g = 2) | `F(Q)²` | `(g_J/2)²` was missing |
| geometry | `Σ (δ_αβ − Q̂_αQ̂_β) S^{αβ}` | diagonal only (single crystal) | cross terms missing (AUDIT P9); powder `2/3 Σ_α` right |
| `k_f/k_i` | `k′/k` | `√(1 − ΔE/E_i)` | same |
| energy sign | `δ(ħω + E_i − E_f)`: `ħω = E_f − E_i` > 0 is neutron energy loss | `E_f − E_i` | same |
| Debye–Waller | `e^{−2W}`, `W ≈ κ²⟨u²⟩/6` | `exp(−B_iso Q²/3)` | same if `B_iso` is `⟨u²⟩`; renamed `u2` |

The transition sum gives the total scattering function, static and
dynamic parts of Jensen and Mackintosh (4.2.3a–c) together. Weights of
single transitions inside a degenerate level depend on the basis; sums over
a level do not (CONV-0008).

## Consequences

- `scattering_function` and `cross_section` replace `spectral_weight`.
- Tests: detailed balance, `S^{αβ}(−ħω) = e^{−ħω/k_BT} S^{βα}(ħω)`; the sum
  rule `Σ_{i,f} p_i Σ_α |<i|J_α|f>|² = J(J+1)`; the powder average equals
  the projector averaged over Q directions; the full projector against the
  diagonal for Q along an axis.
