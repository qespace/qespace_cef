---
id: CONV-0002
title: Stevens operators, their normalisation and tables
status: draft
owner: qespace_cef
enforced_by:
since:
supersedes:
---

## Statement

`OPEN:` which table defines the Stevens operators `O_lm(J)`, l = 0 to 6, and
with it their normalisation.

As implemented: the tables of McPhase and EasySpin; operators with m < 0 are
the sine-type (imaginary-part) combinations, written `O4m3`; odd l are
implemented and accepted as parameters.

Options:
1. The McPhase tables, the same source as the CEF parameters.
2. PyCrystalField's tables (agree with today's up to a sign for some operators).
3. Generated from spherical tensor operators with a stated normalisation, the
   hand-written tables kept as a test fixture.

`OPEN:` whether odd l are accepted as CEF parameters.

## Rationale

The normalisation of `O_lm` fixes the meaning of every `B_lm`: the same
physics with another table gives other numbers, and a sign difference in one
operator flips the sign of its parameter. Inside one J multiplet the odd-l
terms vanish by parity, so they are building blocks (e.g. for operators
between multiplets), not CEF parameters.

## Consequences

- A test compares every operator with the chosen table, including the sign.
- Where another code differs by a sign, the parameter has the opposite sign
  there; a comparison test records the case instead of accepting ±.
