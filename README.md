# kinesis

*κίνησις, change: Aristotle's word for motion in the widest sense, of place, of quality, of quantity.* kinesis holds dynamical systems as objects of study and of inference: a language in which a system is declared, whatever its kind; what a declared law does over time; what a record of a system identifies; and laws recovered from records. It is one of the foundation libraries of [Noetarium, the Science and Humanity Infrastructure](https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure), but it stands on its own (with calx, morphe, physis and tyche beside it).

## 1. Concepts

The theory, with what each method rests on, is in [docs/theory.md](docs/theory.md).

- **Four primitives.** A dynamical system is given by its Change (the generators and how each acts), its Structure (how changes compose into histories, and what a history weighs), its Constraint (what every change respects) and its Difference (what is told apart). These are the user's primitives for generative relations, and they are the four parts of a declaration (theory §1, §4).
- **One language for every kind.** A flow, a map, a system with events, a stochastic differential, a system given by its action, an open quantum system, a rewriting system and a system made of systems are declared in one notation, read by a grammar that says where a script goes wrong. A law need not be an equation: an action is varied, exactly, into the equations of its stationary paths.
- **Terms.** Every expression is a term with exact symbolic derivatives. The Jacobian of a flow, the Euler–Lagrange equations of an action and the correction from Stratonovich's calculus to Itô's are derived from the terms, not approximated.
- **Levels of access.** A system is known by a presentation, by an oracle or by a record. A tool states what it needs, and refuses with a reason when a system does not give it.
- **What a record identifies.** An answer drawn from a record is the set of values over every history the record admits under stated assumptions. It is certified when the set is one value, and is otherwise reported as it is. Braids of trajectories are the first case: a braid read from samples is certified only under a bound on the speed of the points, and the certificate says how far it is from failing (theory §3, §7).
- **A model is a template filled.** A declaration with unknowns is a template, and a law is matched, discovered or fitted by filling them. A fitted model is a system of the language again, and every tool applies to it.
- **The lift to operators.** A change lifted to the algebra of observables is a linear operator, whatever the system. Koopman's operator, a transfer operator on a partition, the generator of an open quantum system and the Burau matrix of a braid are that lift (theory §2).

A method belongs here when it studies what a law does over time, or recovers a law from a record, and names no object in nature. physis writes the laws of physics as generators; kinesis studies generators of any kind ([decision 0002](development/decisions/0002-scope-boundaries-and-dependencies.md)).

## 2. Modules

| Module | Holds | Status |
|---|---|---|
| [kernel](kernel/README.md) | What every tool answers with: the levels of access, a refusal with its reason, the values a record identifies, and how a certificate stands to the one before it | draft |
| [terms](terms/README.md) | Terms in a normal form, their exact derivatives, their text and LaTeX; programs compiled from them | draft |
| [language](language/README.md) | The grammar of declarations, in BNF read by morphe; a script as it was written; the reader and the writer | draft |
| [systems](systems/README.md) | The form of a declaration; flows and maps with events; stochastic differentials (Itô and Stratonovich; Euler–Maruyama and Milstein); actions varied into flows; rewriting systems and their histories; open quantum systems | draft |
| [analysis](analysis/README.md) | Lyapunov's exponents with the error of a finite time, equilibria, sections, bifurcation diagrams, basins and basin stability; the averages of what a system observes, with their errors, and the power law by which they scale; signs of criticality (slowing down, fluctuation analysis, avalanches, the branching ratio); landscapes (potential, stationary landscape, Kramers's time); the transfer operator on a partition | draft |
| [braids](braids/README.md) | The braid of a record of planar trajectories; what the record certifies under a bound on speed and an error of observation; the identified windings of a pair; topological entropy between Burau's bound and an estimate | draft |
| [induction](induction/README.md) | Templates filled from records: least squares and Levenberg–Marquardt on any declaration, unknowns inside actions among them; sparse discovery and its weak form; stochastic templates; candidates ranked on a part kept aside; linear operators; a forecaster with no law | draft |
| [coupling](coupling/README.md) | Systems made of systems: flows merged into one declaration, exactly; parts given by actions joined by an interaction action, with what the whole keeps; parts of different forms co-simulated | draft |
| [catalogue](catalogue/README.md) | Thirty models from the literature as declarations, with their sources, a shell model of turbulence among them; plain models as templates for a record, in every order of its components | draft |
| [balance](balance/README.md) | Balance laws of a density on a graph of cells (diffusion, crowding, drift; growth), run, and learned from a few snapshots: the integral form, shooting, a negative diffusion set aside; forecasts judged as a whole and as a pattern; the modes of a field from a few snapshots | draft |
| [prisms](prisms/README.md) | What a record of positions identifies of a path under a bound on speed (Hägerstrand's space-time prism): the deepest point of discs; consistency with a witness history or the pair that cannot be joined; the slice at an instant; when a region can be reached; presence in a region over every admissible history, an alibi certified ([decision 0006](development/decisions/0006-prisms.md)) | draft |
| [scales](scales/README.md) | Systems across scales: quantities on a hierarchy of units (cells within districts within a city), carried between levels by named operators (sums and reductions up, a parent's value down, shares, logits and z-scores among siblings, the graph Laplacian, the production-constrained inflow) and advanced by laws each at its own pace (a flow, a map, a fast equilibrium); every input with how it is known, what is missing refused by name; ensembles over what is declared and its screening; unknowns fitted to records at any level, with what the record does not identify said; perfect aggregation of linear laws ([decision 0007](development/decisions/0007-systems-across-scales.md)) | draft |
| fields | Statistical fields that run: histories weighted by an action on a declared background; the expansion of a field theory into rules | planned |
| reversal | The time reversal of a stochastic differential; scores, and learning them | planned |
| oracles | Systems known only by driving and reading: reservoirs, simulators and devices as declarations | planned |

## 3. Checks

`node development/tools/run-cpp.mjs foundations/kinesis foundations/kinesis` from the infrastructure's root builds and runs the 562 checks. Each holds against an answer known beforehand: closed forms (sums of exponents equal to a constant divergence; a map conjugate to a shift; the period of a firing neuron; the stationary variance of a process; Kramers's time; counted words of a grammar; the two invariants of a cascade and the balance of its fluxes); the windows of a space-time prism and a disc's deepest point; a hierarchy's sums, allocations and ratios of sums, the orders of its multirate and fast–slow schemes, Harris and Wilson's equilibria; published values (Lorenz's largest exponent; the conditional exponents of Pecora and Carroll; entropies of braids as logarithms of algebraic numbers); exact solutions on the same Brownian path; derivatives against central differences; and laws that the checks declare and then recover from their own records. The numbers of the draft on which the braids rest are reproduced.

## 4. Dependencies

calx (linear algebra, integrators, random streams, `verification`), morphe (the BNF reader and Earley's parser, grammars, braid groups, Laurent polynomials), physis (the Runge–Kutta step and spectra of exponents, the Lindblad generator), tyche (descriptive statistics, red-noise surrogates, Latin hypercubes and Morris's screening, the normal quantile). Header-only C++20.

## 5. Licence

MIT; see [LICENSE](LICENSE) and [CONTRIBUTORS.txt](CONTRIBUTORS.txt).
