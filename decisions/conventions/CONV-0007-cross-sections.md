---
id: CONV-0007
title: Neutron cross section: normalisation, form factor and units
status: draft
owner: qespace_cef
enforced_by:
since:
supersedes:
---

## Statement

`OPEN:` the normalisation, form factor and units of the neutron cross section.

As implemented:
- intensities in arbitrary units, `N (γ r0)² = 1`;
- the dipole form factor `F(Q) = <j0> + (2 − gJ)/gJ <j2>` with the ILL
  coefficients (for Ce those of Ce²⁺);
- the Debye–Waller factor `exp(−B_iso Q²/3)`, so `B_iso` is <u²> in Å², not
  the crystallographic B;
- the factor `k_f/k_i = √(1 − ΔE/Ei)`, with ΔE > 0 the neutron energy loss;
- the powder average `⅔ Σ_α |<i|J_α|f>|²`; for a single crystal the factor
  `Σ_α (1 − Q̂_α²) |<i|J_α|f>|²`;
- transition matrices `[i, f]`, rows the initial states, `E_f − E_i`.

`OPEN:` arbitrary or absolute units (barn sr⁻¹ meV⁻¹ per ion), the
single-crystal factor (the full projector `δ_αβ − Q̂_α Q̂_β` has cross terms
that this one drops), and the sign of the energy transfer.

## Rationale

Intensities are compared with measured spectra, so every factor that changes
with Q, E or the transition has to be stated. The energy-transfer sign and the
units are shared with other qespace packages and may come from the suite.

## Consequences

- Each factor is a separate, tested function, so it can be switched off or
  replaced.
- Sum rules (the total weight over all transitions) test the normalisation.
