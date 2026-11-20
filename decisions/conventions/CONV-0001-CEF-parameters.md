---
id: CONV-0001
title: CEF parameters in the Stevens representation of McPhase
status: accepted
owner: qespace_cef
enforced_by:
since:
supersedes:
---

## Statement

We follow the representation of crystal field parameters from McPhase:
https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node130.html
The default representation is Stevens, with parameter labelling `Blm`:
`H_CEF = Σ_lm B_lm O_lm(J)`, with `O_lm` the Stevens operators of the ground
J multiplet and `B_lm` in the package's energy unit.
Parameters with m < 0 (the sine-type operators) are labelled `B4m3`.

Any future implementation of other base will be a matter of further review.

## Rationale

To pay respects to McPhase and base somewhere historically. McPhase is a
long-standing reference implementation, and its manual documents both the
Stevens parameters and their conversion to Wybourne, so a published parameter
set can be checked against it. Stevens operator equivalents are the natural
form inside one J multiplet and are what the neutron-scattering codes
(McPhase, Mantid, PyCrystalField) and most of the rare-earth magnetism
literature use. The alternative, Wybourne `B^k_q` on spherical tensors, is
independent of the multiplet and comparable between ions; it is needed for
intermediate coupling and common in optical spectroscopy.

References: K. W. H. Stevens, Proc. Phys. Soc. A 65, 209 (1952),
doi:10.1088/0370-1298/65/3/308; M. T. Hutchings, Solid State Phys. 16, 227
(1964), doi:10.1016/S0081-1947(08)60517-2.

## Consequences

- `B_lm` contains the Stevens factor θ_l of the ion, so a parameter set belongs
  to one ion and one multiplet and is not transferable between ions.
- A Wybourne input, if accepted, is converted once at the boundary with
  `B_lm = λ_lm θ_l B^l_m`; the λ table needs a source and a test.
