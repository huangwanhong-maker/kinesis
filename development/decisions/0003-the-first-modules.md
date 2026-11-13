# 0003. The First Modules, and the Language

- **Status:** Accepted, provisionally
- **Date:** 2026-11-13
- **Decided by:** Claude, on the user's instruction to implement what was discussed, for the infrastructure's ADR 0103; the language's parts and its first tests are the user's

## 1. Context

The user set the first test of the language: that it can declare a generative grammar, a background-independent dynamics (a group field theory), and a random model (a diffusion model). They turned an early sketch away as "mainly for PDE" and asked that a system may be given by its action; that their four primitives (Change, Structure, Constraint, Difference) be considered in the design; and that the work cover simulation, topological features, criticality, landscapes, heterogeneous systems connected and run together, and models trained and matched from data.

## 2. Decision

Nine modules, each a header-only target `kinesis::<module>`, in this order:

| Module | Depends on | Holds |
|---|---|---|
| `kernel` | — | `Access`, `Refusal`, `Identified`, `Disposition` |
| `terms` | kernel | terms in normal form, derivatives, text and LaTeX; compiled programs |
| `language` | terms, morphe's `formal` | the grammar; `System` and `Script`; `read` and `write` |
| `systems` | language, physis's `dynamics` and `quantum`, calx | `form_of`; `Dynamics`, `simulate`; `Stochastic`; `vary`; `Rewriting`; `Quantum` |
| `analysis` | systems, tyche | exponents, equilibria, sections, diagrams, basins; criticality; landscapes; the transfer operator |
| `braids` | kernel, morphe's `algebra` | `read`, `certify`, windings, entropy |
| `induction` | systems | `fit`, sparse discovery, stochastic templates, `match`, operators, a forecaster |
| `coupling` | systems | `merged`, `Cosimulation` |
| `catalogue` | language | models as declarations; candidates for a record |

Decisions fixed with them:

- **The four primitives are the parts of a declaration:** `change`, `structure`, `constraint`, `difference`, after the carrier (`state`, `param`, `noise`, `unknown`, `input`). The background geometry of a field is a constraint.
- **The grammar is a text** in morphe's BNF over bytes, read by Earley's parser: blocks in braces, products written with `*`, ASCII only. An error is answered with its line, its column and what may follow.
- **A law is a term.** Derivatives are symbolic and exact. A system given by an action is varied into a flow as terms.
- **The prime** is a derivative in real time and the next value in the naturals.
- **A result is a declaration again:** what `vary`, `merged` and `fit` return is written out and read back.
- **A template is a declaration with unknowns** (`?a`), and matching, discovery and fitting are one mechanism.
- **Noise is a declared source** with a keyed stream; the calculus is stated, and Itô's is the default.
- **Events** of a flow are crossings found by bisection; of a map, levels.
- **A composite** of flows is merged exactly; parts of different forms are co-simulated by the explicit scheme with inputs held.
- **Braids:** σ_i is the counterclockwise half-exchange, the left strand below; a reading is on a line, with ties at a sample broken one way everywhere.

## 3. Alternatives Considered

- **A language of equations only.** The first sketch. Turned away by the user: a field theory is given by its action, a grammar by its rules, an open quantum system by a Hamiltonian and jumps.
- **A library of classes with no language.** A system built in code cannot be kept in a project, shown to a person, or returned by a fit as something to read. The declaration is the interface.
- **An existing language** (Modelica, SBML, the notations of the scientific machine-learning libraries). The survey found none that declares a grammar, an action and a stochastic differential in one notation, and none whose answers from data carry their assumptions.
- **Automatic differentiation in place of terms.** calx's jets differentiate a program; they do not give an equation to show, to write out or to vary. Terms do, and the systems declared are small.
- **A parser written by hand.** morphe's reader exists and diagnoses exactly; genomics' notations and kategoria's scripts set the precedent.
- **Indentation in place of braces.** A context-free grammar over bytes cannot read it.

## 4. Consequences

- 397 checks at founding: kernel 17, terms 71, language 38, systems 64, analysis 58, braids 44, induction 37, coupling 21, catalogue 47.
- **The three tests.** A generative grammar is declared and run: its words, the histories of each, its samples. A random model is declared and run, with its noise as a source and its calculus stated; its reversal, which is the diffusion model's generation, is read and not carried out. A group field theory is declared with no background, and read; it is not run, and its expansion into rules is not carried out. The first is passed; the second and third are passed as declarations and not yet as computations, and the library says so where it is asked.
- **Read and not run** (theory note, §11): the derivations other than `vary`, statistical fields, unknown functions. Each is answered with a note.
- Systems are small: terms are evaluated by a stack program, a state of an open quantum system has dimension at most 64, and Euler–Lagrange equations are solved as terms for at most three coordinates.
- Found while building, and now part of the design: a number that multiplies a sum is distributed only inside a sum, or a term does not read back as it was written; letters that appear only as context in the rules are not rewritten letters; a crossing that falls exactly on a sample needs one rule at both ends of every interval; Milstein's condition is on every driven component, whichever source is declared first; a forecaster with delays has modes of its own, which too small a ridge lets grow.
