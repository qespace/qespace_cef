---
id: ADR-0004
title: Rules, descriptions and tracking: how the maintainer works with agents
status: accepted
date: 2026-11-22
supersedes: ADR-0001
---

## Context

ADR-0001 ranked the rules `CLAUDE.md` < ADR < CONV and froze accepted
decisions. CONVs then took on implementation (function names, return
values, array shapes), and through the precedence a CONV overrode an ADR on
a question of design: CONV-0012's flat `(n²,)` against ADR-0003's
`S + (n, n)`. Freezing everything also made small corrections as costly as
changes of meaning.

## Options

1. Keep ADR-0001.
2. A chain with CONV on top, but narrower CONVs.
3. ADRs and CONVs side by side above `CLAUDE.md`, each with its own domain.

## Decision

Option 3. The project's documents are of three kinds.

1. **Rules** bind every change.
   - `CLAUDE.md`: how to work (review, git, environment, code, tests). No
     design and no physics decision. It tells the agent to read `decisions/`
     before implementing.
   - CONVs in `decisions/conventions/`: what a quantity or rule of the
     physics is: definition, units, signs, normalisation, frames, bases.
     Not how the code exposes it: no function, argument, return type or
     array shape.
   - ADRs in `decisions/adr/`: design choices, among them how the code
     exposes the quantities of the CONVs (names, arguments, shapes). An ADR
     cites the CONVs it builds on and defines no convention.

   Precedence: `CLAUDE.md` < ADR and `CLAUDE.md` < CONV. ADRs and CONVs
   don't overlap, so neither ranks above the other; an overlap or a
   contradiction between them is a mistake in one of them, raised with the
   maintainer. One file per decision, from the `_template.md` in its folder.
   Only the maintainer accepts one. An accepted one is edited only
   editorially (typos, links, formatting, wording that keeps the meaning);
   any other change is a new one that supersedes it. Decisions taken from
   the suite or elsewhere are dropped in with their ID and source.
2. **Descriptions** say what is: `ARCHITECTURE.md`, `DEV.md`, `README.md`,
   the docstrings and `docs/` (the Sphinx sources, which keep the name tools
   expect). When one disagrees with a rule, the description is fixed.
3. **Tracking** says what to do next and changes all the time: `AUDIT.md`
   (findings and ideas from internal or external reviews), `PLAN.md` (the
   next step) and `ROADMAP.md` (feature ideas, written for people). Rules
   don't depend on them; what a rule needs from them becomes an ADR or a
   CONV.

## Consequences

An agent reads `CLAUDE.md`, then the ADRs and CONVs that touch its task,
then the descriptions, and takes its task from `PLAN.md`. Every change to a
rule is checked against the domains above. Accepted CONVs that name
functions, returns or shapes (CONV-0003, 0007, 0009, 0010, 0011, 0012) are
superseded by a CONV without them and an ADR with them. Tracking can be rewritten or discarded
without touching a rule. Decisions dropped in from elsewhere keep their
IDs, which can clash with ours. `OPEN:` how to number them.
