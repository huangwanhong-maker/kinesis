# Answers

Status: **draft** (2026-11-13). Module `kinesis::kernel`.

## 1. Concepts

Every tool of kinesis answers in one vocabulary. A system is known at a level of access: by a presentation, a law one can read; by an oracle, something one can drive and read; or by a record, finite data seen through an observation. The levels are ordered: a presentation gives an oracle, and an oracle gives records. A tool states the access it needs, and one that is asked for more than a system gives refuses with the reason; it never answers with a made-up value.

An answer drawn from a record is not one value but the set of values of a feature Φ over every history Γ the record D admits under stated assumptions A, $I_\Phi(D; A) = \{\Phi(\Gamma) : \Gamma \in R(D; A)\}$. It is certified when the set is a single value, ambiguous when it has several, and inconsistent when no history satisfies the assumptions. A certificate keeps the assumptions it was issued under and how it stands to the one before it: retained by refinement, transported across a change of representation, withdrawn, or issued again. Theory: [docs/theory.md](../docs/theory.md).

## 2. Content

- `answer.hpp`:
  - `Access{record, oracle, presentation}`, ordered, with `name(Access)`; `require_access(have, need, tool)`, which refuses when `have` is below `need`, naming what the tool needs and what it was given.
  - `Refusal`: a request that cannot be answered, with the reason; thrown by tools and reported by their callers.
  - `Identified<T>{values, complete, assumptions, reason}`: the values of a feature over every admissible history, kept distinct and sorted; `add(v)`, `certified()`, `ambiguous()`, `inconsistent()`, `value()` (refused with the reason unless certified), `refines(earlier)` (every value of this set is a value of the earlier one).
  - `Disposition{issued, retained, transported, withdrawn, reissued, ambiguous, inconsistent}` with `name(Disposition)`; `disposition(before, now)`, how a new set stands to the one issued before.

*Remark.* `complete` is false when the list is known to be a part of the set (an enumeration that was cut short). A set that is not complete certifies nothing, and when it is empty it is not inconsistent: nothing was found, which is not that nothing is there.

*Caveat.* `disposition` compares two sets under the same feature and representation, so it never answers `transported`: that is a caller's to say, after a change of representation.

## 3. Checks

17: the order of access (a system known by more than a tool needs is given to it; a tool asked for more than a system gives refuses, naming what it needs and what it was given); identified values kept distinct and sorted; one value certified, several ambiguous, none inconsistent; a set known to be a part certifying nothing, and not inconsistent when empty; the value of a set that is not certified refused with the reason; a consistent refinement only contracting the set (Huang 2026, Proposition 7.1); and the dispositions of a certificate against the one before it, in the cases of the draft's Table 2: issued, retained by refinement, withdrawn, issued again, ambiguous, inconsistent.

## 4. Dependencies

None among the modules of kinesis or the sibling libraries. The checks use calx's `verification`.

## 5. References

Huang (2026), "Rethinking trajectory-based credibility under symbol inflation" (a draft): Eqs. 19–20 for the identified values, §7 for how a certificate stands to the one before it, Proposition 7.1 for refinement, Table 2 for the cases checked.
