# The Language of Dynamical Systems

Status: **draft** (2026-11-13). Module `kinesis::language`.

## 1. Concepts

A script is a list of systems. A system is a declaration in parts, or a derivation from another system. The parts say what a dynamical system is made of: `state`, `param`, `noise`, `unknown` and `input` give the carrier, its parameters, its sources of noise, what is not yet known and the signals it takes from outside; `change` gives the generators and how each acts (equations of a flow, a map or a stochastic differential, an action, a Hamiltonian with its jumps, rules, events with a guard and resets); `structure` says how changes compose into histories and what a history weighs; `constraint` says what every change respects; `difference` says what is told apart; `part` and `connect` make a system of others.

The grammar is written in BNF over bytes and read by morphe's reader; every script is parsed by Earley's parser. A text that is not a script is refused with the place as far as which it is right and what may follow there, from the parser's chart. A text that is a script becomes declarations that hold what was written, in the order and under the names of the grammar, with every expression as a term in normal form. Writing is the inverse on what it covers, so that a system a tool induces is a word of the language. Theory: [docs/theory.md](../docs/theory.md).

## 2. Content

- `grammar.hpp`: `grammar_text`, the grammar: `script`, `system` (a `body` or a `derivation` by `expand`, `reverse`, `lift`, `vary`, `legendre` or `unravel`), the items `state`, `param`, `noise`, `unknown`, `input`, `part`, `connect`, `change`, `structure`, `constraint`, `difference`, and the expressions `expr`.
- `declaration.hpp`: `Place` (line and column, from 1); `ScriptError`, a `Refusal` with its place, of syntax or of meaning; `Domain`, `Declarator`, `Space`, `Ranging`, `States`, `Binding`, `Noise`, `Unknown`, `Input`, `Equation`, `Named`, `Rule`, `Assignment`, `Event`, `Change`, `Structure`, `Constraint`, `Observation`, `Difference`, `Part`, `Wire`; `System` with `derived()`; `Script` with `find(name)`.
- `reader.hpp`: `read(text)`, the `Script` a text is, or a `ScriptError`; `read_term(text)`, one expression as a term; `write(system)` and `write(script)`, the text the language writes.

*Remark.* Blocks are in braces. A product is written with `*`, since juxtaposition would be ambiguous with application. A line may continue before a binary operator, and no item begins with one. A power binds tighter than a negation, a negation tighter than a product, a product tighter than a sum; a power associates to the right. The draft is ASCII only.

*Remark.* A change to the grammar changes the reader, the declarations and the checks together.

*Caveat.* The declarations decide nothing: which kind of system a declaration is, and whether its names resolve, is decided by `kinesis::systems`. A derivation is read as its operator and its source, and is not carried out here.

*Caveat.* The language reads more than the library runs. Of the derivations (`vary`, `legendre`, `expand`, `reverse`, `lift`, `unravel`) only `vary` is carried out; a field weighted by its action (`histories weight exp(-S)`), with or without a background, is read and not run; an unknown function is read and not filled. Each is answered where it is asked with a note that says so (the theory note, §4 and §11).

## 3. Checks

40: the grammar compiles in morphe's reader, and the three tests the language was set (a generative grammar, a group field theory with no background, a diffusion model with noise of a stated law) are read, each part landing where it was written; precedence and association, on expressions whose value decides them (−2² is −4, powers associating to the right, a − b − c, a / b / c, a / b * c); terms built at random, written and read back to the same term; scripts written from their declarations reading back to declarations that write the same; events, inputs, parts and connections, named observations, unknowns and derivations read; and refusals with their place: a product by juxtaposition, an unclosed block, a system declared twice. A part's state in a term may carry its primes (`a.q'`, the grammar's `qvar`): it is read and written back.

## 4. Dependencies

`kinesis::terms` (and through it `kinesis::kernel`); morphe's `formal` (`bnf.hpp`: grammars written in BNF, Earley's parser and its diagnostics). The checks use calx's `verification`.

## 5. References

The idea of 2026-11-12, §9.1, for the choices of the grammar. The headers cite no other source.
