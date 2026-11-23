# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: decisions before v2 (advisor and maintainer, no worker)

**Goal.** The rules `refactor/api-v2` builds on are settled and keep to
their domains (ADR-0004), so v2 is planned once and not reworked. The
naming audit is done (ADR-0005, ADR-0006, `NAMING.md`). No package code
changes.

**Steps.**

1. **CONV-0002, vectors in matrices:** columns or rows for every matrix
   of vectors (Q77), then accepted. If arrays need their own rule, an ADR
   next to ADR-0003.
2. **A16:** a CONV superseding CONV-0012 and CONV-0010 for the physics of
   energies, excitation energies `E_f − E_i` and the weights, as matrices
   over (i, f) by CONV-0002; an ADR for the shapes the code returns.
3. **The CONVs that name code** (ADR-0004): superseding CONV-0003 (with
   Q42, Q44), CONV-0007 and CONV-0011, their values unchanged.
4. **Plan `refactor/api-v2`:** the branches, each green, batched results
   checked against v1 (ADR-0003, ADR-0005), with P21 and the other AUDIT
   items v2 rewrites. That plan becomes the next *Now*, for the worker.

**Questions.** Q77 (CONV-0002); Q42 and Q44 (CONV-0003).

**Files.** `decisions/`, `PLAN.md`, `AUDIT.md`.

## After that

1. **`refactor/api-v2`**, planned in step 4.
2. **The g-tensor (C7):** a draft CONV on the CONV-0008 grouping and
   CONV-0011's basis, returning G or the principal values and axes. Then P4
   and P17.

Other fixes wait for a source: P14 (sources of the ion table), P16 (Q12).
