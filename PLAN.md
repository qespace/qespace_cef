# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001).

## Now: `docs/conventions`, deciding the CONV drafts

**Goal.** The physics conventions are decided as CONVs, so that the fixes and
the v2 API can be built on them.

**Steps.**
1. The maintainer reviews the drafts CONV-0001 to CONV-0008 in
   `decisions/conventions/` and decides their `OPEN:` parts.
2. A draft for the g-tensor (AUDIT C7), which the old table had no row for.
3. Each accepted CONV gets the test named in its `enforced_by`.

**Files.** `decisions/conventions/`, the tests for the accepted CONVs.

**Risks.** A convention that differs from the code turns its test into a
strict `xfail` until the fix lands; the fixes are separate steps.

**Tests.** `uv run pytest`: the new convention tests pass or are strict
`xfail` with the reason.

## Next

Not chosen yet. The candidates in AUDIT section 7: `fix/susceptibility`
(needs no convention) and `refactor/api-v2` (needs ADR-0002's open parts and
the accepted CONVs).
