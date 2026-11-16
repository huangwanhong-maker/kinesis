# CLAUDE.md

This file guides Claude Code in the kinesis repository.

kinesis is one of the foundation libraries of [Noetarium, the Science and Humanity Infrastructure](https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure), which includes it as the git submodule `foundations/kinesis`. It is a project in its own right and must make sense without the infrastructure.

## 1. Purpose and Scope

kinesis holds dynamical systems as objects of study and of inference:

- a language in which a system is declared by its change, structure, constraint and difference, whatever its kind: a flow, a map, a system with events, a stochastic differential, a system given by its action, an open quantum system, a rewriting system, a system made of systems;
- what a declared law does over time: runs, exponents, equilibria, diagrams, basins, landscapes, signs of criticality, the lift to operators;
- what a record identifies: braids of trajectories, read and certified under stated assumptions;
- laws recovered from records: templates filled, sparse discovery, candidates ranked, a catalogue of models; balance laws of a density on a graph of cells, run and learned from a few snapshots (decision 0005);
- later, statistical fields that run, the reversal of stochastic differentials, and oracles.

It builds on four sibling libraries:

- **calx**, for linear algebra, random streams and the report of its checks;
- **morphe**, for grammars written in BNF with Earley's parser, and for braid groups;
- **physis**, for the Runge–Kutta step, spectra of exponents and the Lindblad generator;
- **tyche**, for descriptive statistics and red-noise surrogates.

Out of scope:

| Topic | Where it belongs |
|---|---|
| The laws of physics as generators: mechanics, fields, quantum dynamics, kinetics | physis |
| Integrators, linear algebra, random streams, signal processing | calx |
| Groups, grammars and parsers, braid groups and their algorithms | morphe |
| Distributions, estimation, time series, samplers of laws known up to normalisation | tyche |
| Search over spaces too large to enumerate | heuresis |
| What a system is *of*: a neuron, an epidemic, an ocean, a market, with its measured values | the domains of the infrastructure |
| Drawing systems | the infrastructure's `visualization/` and its works |

The boundary is [decision 0002](development/decisions/0002-scope-boundaries-and-dependencies.md): physis writes generators; kinesis studies what a generator of any kind does over time, and recovers generators from records.

## 2. Invariants

1. **Dependencies.** kinesis depends only on calx, morphe, physis and tyche among project libraries, and never on peira, nomos, eikon, heuresis, syntaxis, kategoria or the infrastructure's domains, edges or applications.
2. **No particular objects.** Nothing in kinesis names an object in nature. A model of the catalogue is a form of law with the values in common use; measured values belong to the domain that owns the object.
3. **Pure.** No I/O and no global state. Randomness comes only from streams keyed by what they are for, so that a path depends on its key alone.
4. **Access is stated.** A tool states the level at which it needs a system (presentation, oracle, record) and refuses with a reason when it is not given. No tool answers with a made-up value.
5. **A record identifies a set.** An answer drawn from a record is the set of values over the histories it admits, with the assumptions it was drawn under. Nothing is reported as certified that is one of several admissible values.
6. **Limits are said.** An estimate says it is one; an indicator says what misleads it; a search says what it did not reach; what the language reads and the library does not run is answered with a note that says so.
7. **Exact where it can be.** Derivatives, variations and conversions between calculi are of terms. Braid words, their reductions and their coordinates are integers.
8. **Truth known.** Checks are held to closed forms, published values, exact solutions or laws the check itself declares. Reference values from other software come from running it, never from memory.
9. **Reuse.** What a sibling holds is used through its public headers, never written again here.
10. **Cited sources.** Every implemented method cites its source; section and equation numbers are given only where certain.

## 3. Design Philosophy

kinesis shares its governing principle with its sibling repositories: **modular, extensible, revisable**.

- **Status vocabulary:**
  - `planned`: named, no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, with a pointer to its replacement.
- **One declaration, many readings.** A system is declared once, and each tool reads the declaration: it is run, differentiated, varied, written out, drawn in LaTeX and fitted from the same terms.
- **A result is a system again.** What a derivation, a composition or a fit returns is a declaration of the language, so that every tool applies to it.
- **A new kind of system is a new form.** It adds a clause to the grammar, a case to `form_of` and a header in `systems`; the tools that can serve it state so.

## 4. Repository Layout

```
kinesis/
├── CMakeLists.txt   the project: one header-only target kinesis::<module> per module (kinesis_module); the siblings added once
├── kernel/ · terms/ · language/ · systems/ · analysis/ · braids/ · induction/ · coupling/ · catalogue/ · balance/
├── docs/            the theory note
└── development/     decisions · ideas · release notes
```

Each module directory holds:
- `README.md`: concepts first, then the content, then the checks;
- `include/kinesis/<module>/`;
- `tests/check_<module>.cpp`, a check program on calx's `verification`, registered with CTest by `kinesis_module` in the top-level `CMakeLists.txt`.

The modules and their order are in [decision 0003](development/decisions/0003-the-first-modules.md).

## 5. Working Procedures

### 5.1 Decisions and Ideas

- Decisions are numbered records in `development/decisions/`. Records are never rewritten; a later record supersedes an earlier one.
- Ideas are files named `development/ideas/YYYY-MM-DD-slug.md`, with status `seed`, `promoted` or `abandoned`.
- **Release notes.** At each commit, add to `development/release-notes/YYYY-MM-DD.md` (one file per day, one short section per commit) a concise summary in plain prose of what is new, what one can now do and what works better. It is not a list of changes (the git log keeps those), and it is not written between commits.

### 5.2 Change Discipline

- Before restructuring or changing an invariant (§2), propose the change and wait for the user's agreement.
- Keep each change to one concept.
- A change to the language changes its grammar (`language/grammar.hpp`), the declaration it is read into, the reader and the writer, `form_of`, and their checks together. What the grammar reads and nothing runs is said in the theory note's §4 and §11.
- A model is added to the catalogue with its source and with a check of what is known of it.

### 5.3 Shared Conventions

calx, morphe, physis, tyche, peira, nomos, eikon, heuresis, syntaxis, kategoria, kinesis and the infrastructure keep these conventions identical:

- writing style;
- status vocabulary;
- decision, idea and release-note conventions.

Change them in every repository together.

### 5.4 Working Inside the Infrastructure

When you edit kinesis inside the infrastructure's `foundations/kinesis`:

1. commit in kinesis first;
2. then commit the updated pointer in the infrastructure.

When a change needs a sibling to change as well, commit them in dependency order: calx and morphe first, then physis, then tyche, then kinesis.

## 6. Conventions

- **The prime.** `x'` is a derivative when time is real and the next value when time is in the naturals. An equation of order k makes the components q, q', …, q⁽ᵏ⁻¹⁾, named so.
- **Products are written.** `sigma * (y - x)`; a name followed by a bracket is an application.
- **Parts.** In a composite a part's state is written `a.q` and its velocity `a.q'`; in the system the parts merge into they are `a_q` and `a_q'`. A coupling is wires (`connect`), or, between parts given by actions, an action of the composite's own, never both.
- **Stochastic differentials** are read in Itô's calculus unless the declaration states Stratonovich's.
- **Events.** An event of a flow is a crossing into its relation; an event of a map is a level.
- **Braids.** σ_i is the counterclockwise half-exchange of the strands at positions i and i + 1, the left one passing below; a word is read from the left, as time runs. Strands and positions are counted from 1 in words and from 0 in code.
- **Quantum operators.** For a spin, level 0 is up and `sm` lowers up to down. A product of operators on one factor is written `ordered(…)`, since terms commute.
- **Time and units.** kinesis carries no units: a declaration's numbers are in whatever units its author states.
- **Writing style.**
  - Headings are declarative, numbered and academic, never questions.
  - Concepts come before structure.
  - Label asides with *Rationale.*, *Remark.* or *Caveat.*
  - British spelling; English throughout.

## 7. Technology

- **Language: C++.** This follows ADR 0003 of the infrastructure.
- **Provisional:** C++20 without extensions, header-only targets `kinesis::<module>`; CMake 3.20 or later. Checks are programs on calx's `verification`, registered with CTest. On the development machine the compiler is GCC in WSL.
- **Dependency resolution.** Inside the infrastructure, kinesis uses the targets of `foundations/calx`, `foundations/morphe`, `foundations/physis` and `foundations/tyche`, adding each that no larger project has added. Never nest submodules.
- **Learning.** Methods that need no unknown function are written here in C++. Where an unknown function is to be learnt, a work may use torch on the CPU beside the library; kinesis itself takes no such dependency.
