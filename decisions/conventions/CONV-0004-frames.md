---
id: CONV-0004
title: Quantisation axis, local frame and crystal frame
status: draft
owner: qespace_cef
supersedes:
---

## Statement

`OPEN:` the frame of the J operators and how it relates to the crystal frame.

As implemented: J_x, J_y, J_z and the CEF parameters are in one Cartesian
frame, with z the quantisation axis; the field and Q are 3-vectors in that
frame. Nothing connects it to the crystal axes.

Options:
1. A local frame with z along the main symmetry axis of the site (the
   setting of the point-group tables); the user gives field and Q in it.
2. As 1, plus a rotation from the crystal frame to the local frame; field and
   Q are given in the crystal frame and rotated at the boundary.
3. McPhase's CEF axes (manual node130): x ∥ a, y ∥ b, z ∥ c for orthogonal
   lattices; otherwise y ∥ b, z ∥ a × b, x perpendicular to y and z.

`OPEN:` the setting of x for each point group, and whether parameters are
rotated between frames.

## Rationale

The allowed parameters of a point group, and the cubic relations such as
B44 = 5 B40, hold only in the setting that defines them. A field "along c" is
a statement in the crystal frame, so a fit to magnetisation or to a single
crystal needs the relation between the two frames. The canonical basis of
CONV-0011 diagonalises J_z, so it depends on this frame too.

Ways to find a local frame, for the discussion (not decided):

1. From the site symmetry, given the point-group operations of the site as
   rotation matrices in the crystal frame (e.g. from a `Crystal` of
   `qespace_crystallography`): z is the axis of the highest-order proper
   rotation, the eigenvector of R with eigenvalue 1; x is a perpendicular
   twofold axis or mirror normal, as the setting of the group says.
2. From the parameters, when there is no symmetry to use (e.g. a point-charge
   sum in the crystal frame): the rank-2 part of H is a symmetric traceless
   tensor, `Q_αβ = Tr(H {J_α, J_β}/2)`, whose principal axes give z (as the
   electric-field-gradient axes in NMR and Mössbauer, ordered
   `|Q_zz| ≥ |Q_yy| ≥ |Q_xx|`). About z, `B44` and `B4−4` turn like
   `cos 4φ` and `sin 4φ`, so x is the angle with `B4−4 = 0` and `B44 > 0`.
   A cubic site has no rank-2 part and needs rank 4 or method 1.
3. Rotating parameters needs no tables of Wigner matrices:
   `H′ = U H U†` with `U = exp(−iθ n̂·J)`, and the parameters are read back
   by `B_k = Tr(O_k H′) / Tr(O_k O_k)`, since the Stevens operators are
   orthogonal under the trace.

Checked in a scratch run (2026-11-21): the Stevens operators of J = 5/2 are
orthogonal to 1.4·10⁻¹⁴. Ce³⁺ with `B20 = 0.5, B40 = 0.01, B44 = 0.1` meV,
tilted by an arbitrary rotation, has 14 nonzero parameters; the rank-2
principal axes (eigenvalues −28, −28, 56) and the rank-4 angle bring back
`B20 = 0.5, B40 = 0.01, B44 = 0.1` exactly, with the same energies.

## Consequences

- Each point-group table states its setting.
- If parameters can be rotated, the rotation acts within each rank l (Wigner
  D-matrices, or numerically as in method 3), and a test checks that the
  levels don't change.
