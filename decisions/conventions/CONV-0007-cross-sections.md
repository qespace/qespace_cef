---
id: CONV-0007
title: Neutron scattering function and cross section
status: accepted
owner: qespace_cef
supersedes:
---

## Statement

The **scattering function** of one ion, with eigenstates |i>, |f>, energies
`E_i`, `E_f` (CONV-0010) and populations `p_i`:

    S^{αβ}(ħω) = Σ_{i,f} p_i <i|J_α|f><f|J_β|i> δ(ħω − (E_f − E_i))

It holds no neutron factor and is returned as excitations with scattering
weights (CONV-0012), elastic ones included; `E_f − E_i` > 0 is neutron
energy loss.

The **cross section** per ion is built from it:

    d²σ/dΩdE = (γr₀)² (k_f/k_i) [½ g_J F(Q)]² e^{−Q²⟨u²⟩/3}
               Σ_{αβ} (δ_αβ − Q̂_α Q̂_β) S^{αβ}(ħω)

with the dipole form factor `F(Q)`, `k_f/k_i = √(1 − ħω/E_i)` and the
mean-square displacement `⟨u²⟩` in Å² (`u2`). A powder averages the
projector to `(2/3) δ_αβ`. Q is in the local frame of CONV-0004. It is
returned as excitations with cross-section weights (CONV-0012), in barn sr⁻¹ per ion
(the cross section itself is in barn sr⁻¹ meV⁻¹ per ion). γ and r₀ come from CODATA 2022
in `constants.py`.

## Rationale

This is Jensen and Mackintosh, sec. 4.2 and (4.2.2a), and Raymond (1.4), for
a single ion (`REFERENCES.md`). Splitting the scattering function from the
cross section keeps the CEF physics apart from the neutron and instrument
factors. Against the code before this CONV:

| Factor | References | Code | Comment |
|---|---|---|---|
| prefactor | `(γr₀)²` | 1 | arbitrary units; now absolute |
| ion factor | `[½ g F(Q)]²`; Raymond `f²(Q)` (spin only, g = 2) | `F(Q)²` | `(g_J/2)²` was missing |
| geometry | `Σ (δ_αβ − Q̂_αQ̂_β) S^{αβ}` | diagonal only (single crystal) | cross terms missing (AUDIT P9); powder `2/3 Σ_α` right |
| `k_f/k_i` | `k′/k` | `√(1 − ΔE/E_i)` | same |
| energy sign | `δ(ħω + E_i − E_f)`: `ħω = E_f − E_i` > 0 is neutron energy loss | `E_f − E_i` | same |
| Debye–Waller | `e^{−2W}`, `W ≈ κ²⟨u²⟩/6` | `exp(−B_iso Q²/3)` | same if `B_iso` is `⟨u²⟩`; renamed `u2` |

The transition sum gives the total scattering function, static and
dynamic parts of Jensen and Mackintosh (4.2.3a–c) together.

## Consequences

- `scattering_function` and `cross_section` replace `spectral_weight`.
- Tests: detailed balance, `S^{αβ}(−ħω) = e^{−ħω/k_BT} S^{βα}(ħω)`; the sum
  rule `Σ_{i,f} p_i Σ_α |<i|J_α|f>|² = J(J+1)`; the powder average equals
  the projector averaged over Q directions; the full projector against the
  diagonal for Q along an axis.
