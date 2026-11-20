---
id: CONV-0005
title: Basis: the ground J multiplet and the order of states
status: draft
owner: qespace_cef
---

## Statement

The basis is the ground J multiplet of the free ion (LS coupling, Hund's
rules), 2J + 1 states, with the Landé `gJ`; there is no intermediate coupling.

As implemented: the states are ordered `|J>, |J−1>, ..., |−J>`: the first
vector component is m = +J.

`OPEN:` whether this order is the convention or only an implementation detail.

## Rationale

For most rare-earth ions the spin–orbit splitting is much larger than the CEF,
so J is a good quantum number and the ground multiplet is enough. It fails
when an excited multiplet is close, as for Sm³⁺ and Eu³⁺, where the CEF mixes
multiplets (J-mixing). The order of states matters for anyone who reads an
eigenvector component by index.

## Consequences

- Intermediate coupling would need the full fⁿ basis and Wybourne parameters
  on spherical tensors. It may come later; a new CONV then supersedes this one.
- Eigenvectors in the API are in this order; a test checks `<J|J_z|J> = J` for
  the first state.
