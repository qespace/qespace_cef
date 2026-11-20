---
id: ADR-0001
title: Rules, descriptions and tracking: how the maintainer works with agents
status: accepted
date: 2026-11-20
supersedes:
---

## Context

The maintainer works with agents through files. Rules, descriptions of the code
and lists of work were mixed: decisions sat in review notes and chat, physics in
the agent instructions. An agent needs to know which files bind it and which
win when they disagree.

## Decision

The project's documents are of three kinds.

1. **Rules** bind every change. Each kind is self-contained and does not cite
   the others:
   - `CLAUDE.md`: how to work (review, git, environment, code, tests). No
     physics. It tells the agent to read `decisions/` before implementing.
   - ADRs in `decisions/adr/`: design choices.
   - CONVs in `decisions/conventions/`: physics conventions (units, signs,
     normalisation, frames, bases).

   Precedence: `CLAUDE.md` < ADR < CONV. One file per decision, from the
   `_template.md` in its folder. Only the maintainer accepts one; an accepted
   one is never edited, a new one supersedes it. Decisions taken from the suite
   or elsewhere are dropped in with their ID and source.
2. **Descriptions** say what is: `ARCHITECTURE.md`, `DEV.md`, `README.md`,
   the docstrings and `docs/` (the Sphinx sources, which keep the name tools
   expect). When one disagrees with a rule, the description is fixed.
3. **Tracking** says what to do next and changes all the time: `AUDIT.md`
   (findings and ideas from internal or external reviews) and `PLAN.md` (the
   next implementation step). Rules never cite them; what a rule needs from
   them becomes an ADR or a CONV.

## Consequences

An agent reads `CLAUDE.md`, then the decisions that touch its task, then the
descriptions, and takes its task from `PLAN.md`. Tracking can be rewritten or
discarded without touching a rule. Decisions dropped in from elsewhere keep
their IDs, which can clash with ours; how to number them is not settled.
