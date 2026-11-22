# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: the naming audit (advisor and maintainer, no worker)

**Goal.** Every public name of v1, and every word the decisions use for a
quantity, with its meaning, unit and shape, and the v2 name it should get;
then the decisions `refactor/api-v2` needs, so it can be planned. No
package code changes.

**Steps.**

1. **Inventory** (advisor, read-only): classes, functions, methods,
   attributes, arguments, returns, modules and constants of v1, and the
   terms of the CONVs, ADRs and docstrings (energies, excitation energies,
   levels, transitions, scattering and cross-section weights, ...). Each
   with where it is, what it means, its unit, its shape in v1 and under
   ADR-0003, and the proposed v2 name with a reason.
2. **Collisions and synonyms**: one name for two quantities (`excitation
   energies` is `E_n − E_0` in CONV-0010 and `ε_k` in CONV-0012), and one
   quantity under two names (AUDIT Q2).
3. **The naming rule** (Q26), drafted from the table.
4. **Your decisions**, drafted by the advisor: Q26, Q29, Q31, Q50, Q51 and
   A16 (one CONV superseding CONV-0010 and CONV-0012: energies `S + (n,)`,
   excitation energies `S + (n, n)`, weights `S + (n, n, 3, 3)` and
   `S + (n, n)`; the shift to the ground state goes to the print).
5. **Plan `refactor/api-v2`**: the branches, each green, batched results
   checked against a loop over v1 (ADR-0003).

**Questions.**

- **Q26 (Q2).** The naming rule for everything but the classes
  (ADR-0002): snake_case with physics symbols as in the formulas (which
  ones?); for this package or for the suite?
- **Q29 (ADR-0002).** Is a `Hamiltonian` immutable and solved on first
  access (a, recommended), or solved in place by `H.solve()` (b)? A16
  relies on (a).
- **Q31 (ADR-0002).** Observables as functions, `cef.magnetization(H, T)`
  (recommended), or as methods, `H.magnetization(T)`?
- **Q50.** Where the audit lives: (a) the rule as a draft ADR-0005 and the
  old → new table as a tracking file next to `AUDIT.md` (recommended); (b)
  everything in AUDIT section 4, as Q2.
- **Q51.** v2 next to v1, in new modules, v1 removed when v2 is complete
  (recommended: `main` stays usable, batched results compare against v1);
  or v1 rewritten in place?

**Files.** `PLAN.md`, `AUDIT.md`, and by Q50 a draft ADR-0005 and the
table file.

## In parallel: decisions only you can take

- The `OPEN:` parts of the draft CONV-0006 (AUDIT
  section 6.2, Q11), and a draft for the g-tensor (C7).

## After that

1. **`refactor/api-v2`**, planned in step 5.
2. **The g-tensor (C7):** a draft CONV on the CONV-0008 grouping and
   CONV-0011's basis, returning G or the principal values and axes. Then P4
   and P17.

Other fixes wait for a convention or a source: P3 (CONV-0003 units of the
thermodynamics, Q44), P14 (sources of the ion
table), P16 (Q12).
