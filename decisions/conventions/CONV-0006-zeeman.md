---
id: CONV-0006
title: Zeeman term and magnetic moment
status: draft
owner: qespace_cef
enforced_by:
since:
supersedes:
---

## Statement

`OPEN:` the sign of the Zeeman term and of the moment, and `gJ J` or `L + 2S`.

As implemented: `H_Z = −μB gJ J·B`, B in T, and the moment `m = +gJ μB <J>`,
so the moment is parallel to J.

Options:
1. As implemented: the moment parallel to J, the usual choice in magnetism
   codes.
2. The electron's sign: `μ = −gJ μB J`, so `H_Z = +gJ μB J·B`.

## Rationale

Both give the same energies, the same M(B) and the same intensities; they
differ in the sign of <J> relative to the moment, which shows up wherever <J>
itself is reported or combined with other operators. Inside the ground J
multiplet `L + 2S = gJ J` (projection theorem), so `gJ J` is exact there;
between multiplets `L + 2S` is needed.

## Consequences

- One sign everywhere: Zeeman term, moment, magnetisation, mean field.
- A test checks that the moment points along a field much larger than the
  CEF and approaches `gJ J μB` there.
