---
id: CONV-0004
title: Quantisation axis, local frame and crystal frame
status: draft
owner: qespace_cef
enforced_by:
since:
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

`OPEN:` the setting of x for each point group, and whether parameters are
rotated between frames.

## Rationale

The allowed parameters of a point group, and the cubic relations such as
B44 = 5 B40, hold only in the setting that defines them. A field "along c" is
a statement in the crystal frame, so a fit to magnetisation or to a single
crystal needs the relation between the two frames.

## Consequences

- Each point-group table states its setting.
- If parameters can be rotated, the rotation acts within each rank l (Wigner
  D-matrices), and a test checks that the levels don't change.
