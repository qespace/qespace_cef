# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `docs/conventions`, deciding the conventions

**Goal.** The physics conventions are decided as CONVs and each accepted one
has a test, so that the fixes and the v2 API can be built on them.

**Steps.**
1. A test for CONV-0001 (accepted): every Stevens operator against the
   McPhase table copied into it, signs included. Checked by hand on
   2026-11-20 (601 matrices, no difference), so the test should pass.
2. The maintainer decides the `OPEN:` parts of the drafts CONV-0003 to
   CONV-0008, and the labelling question AUDIT Q33.
3. A draft for the g-tensor (AUDIT C7), which has no CONV yet.
4. Each newly accepted CONV gets a test that names its ID.

**Files.** `decisions/conventions/`, `tests/test_cefmatrices.py` or a new
test module for the conventions.

**Risks.** A convention that differs from the code turns its test into a
strict `xfail` until the fix lands; the fixes are separate steps.

**Tests.** `uv run pytest`: the new convention tests pass or are strict
`xfail` with the reason.

## Next

Not chosen yet. The candidates in AUDIT section 7: `fix/susceptibility`
(needs no convention) and `refactor/api-v2` (needs ADR-0002's open parts and
the accepted CONVs).
