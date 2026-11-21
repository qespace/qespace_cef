---
id: CONV-0004
title: The local frame
status: accepted
owner: qespace_cef
supersedes:
---

## Statement

The J matrices (CONV-0005), the parameters `B_lm`, the magnetic field and Q
are all in one local Cartesian frame, with z the quantisation axis. The
caller gives the field and Q in that frame. The package does not relate it to
the crystal frame.

## Rationale

The parameters fix the frame: the same ion in a rotated frame has other
`B_lm` and the same energies, so a field or Q is meaningful only in the frame
of the parameters. Relating the local frame to the crystal frame needs the
site symmetry or the crystal structure, which the package does not have yet
(`ROADMAP.md`); a later CONV adds that rotation. The canonical basis of
CONV-0011 depends on this frame.

## Consequences

- Docstrings state that the field and Q are in the frame of the parameters.
- Once point groups are supported, each table states its setting.
