# 0004. Averages, Coupling by an Action, and a Cascade of Turbulence

- **Status:** Accepted, provisionally
- **Date:** 2026-11-13
- **Decided by:** the user ("try a fluid dynamics and investigate turbalence in the lab"; "several integrated heteogenous dynamics examples … They co-evoluated"; "I think such utilize a structured field theory based method and utilize action S_{int} will be easier? I not sure") and Claude; the infrastructure's ADR 0103

## 1. Context

On the day the library was founded the user asked for two case studies: turbulence in a fluid, and systems of different kinds (a fluid, its chemistry, a species; the behaviours of two animals) joined and evolving together. They asked whether a coupling is easier to give as an interaction term of an action.

Three things were missing. A system that neither settles nor repeats is described by its statistics, and the library had runs and exponents but no averages. Parts given by actions could be wired at inputs they do not have, and could not be joined as actions are joined, by adding a term. And the catalogue held no model of turbulence.

## 2. Decision

1. **Averages** (`analysis/statistics.hpp`). The means of a system's declared observations along one long path, each with the standard error of its block means and with its second moment over its mean squared; and the power law through means set against a scale, by least squares on logarithms, with the error of the slope. A system that declares nothing to observe is refused.
2. **Coupling by an action** (`coupling/composite.hpp`). A composite whose parts are each given by an action may give its coupling as an action of its own, with `histories stationary`: the whole is the system whose Lagrangian is the sum of the parts' and the interaction's, and is varied as one. In a term, a part's state is `a.q` and its velocity `a.q'` (the grammar's `qvar`).
3. **The answer to the question asked.** An interaction action is the form of a coupling that is reciprocal: each part acts on the other through the one term, and what the whole action is invariant under is kept (Noether). A wire is the form of a coupling that need not be: `b.u = k * (a.q - b.q)` moves b and leaves a alone. Both stay in the language, and each says which it is. Couplings between living things are mostly wires.
4. **Recorded, not built:** the response-field action of Martin, Siggia and Rose, in which every wired composite of flows, reciprocal or not and with noise or without, is the stationarity of one action in doubled variables. It is the field-theoretic form the user's question points to for dissipative and one-way couplings.
5. **A shell model of turbulence** (`catalogue/turbulence.hpp`). The Sabra model on fourteen shells, written out as a declaration of 28 real components with what is observed of it: each shell's energy, the whole energy, the second quadratic quantity, the power put in and dissipated, and the flux of energy through two shells. It is a model of the cascade and is said not to be a discretisation of a fluid.

## 3. Alternatives Considered

- **A discretised fluid** (the Navier–Stokes equations on a grid or in Fourier modes). That is the turbulence of a fluid itself. It needs fields that run and thousands of components; the library reads fields and does not run them, and its terms are for small systems. A shell model is what a library of declared small systems can hold honestly, and it keeps the cascade, the exact balance of fluxes and intermittency.
- **Indexed families in the language** (`state u[1..N]`), so that fourteen shells are one line. Wanted, and a change to the grammar, the scope and every tool. The declaration is written out by a function meanwhile, which a person can read line by line.
- **Every coupling as an action.** It would give one form. It cannot state a one-way coupling without doubling the variables, and the doubled form is not what a person writes for "fish graze plankton".
- **Wires only.** They state every coupling, and say nothing of what a coupling keeps: two pendulums and a spring would be two wires that happen to be equal and opposite, with the energy of the whole nowhere declared.
- **Averages over an ensemble of paths.** Right for systems that are not ergodic. One long path with block errors is the usual estimate for a chaotic attractor, and its caveat is stated.

## 4. Consequences

- 433 checks: analysis 68 (the block mean and its error on a series with known block means, a power law recovered exactly and with scatter, cos² over whole turns, the logistic map's mean and second moment), coupling 32 (two oscillators joined by an action against their normal modes, the energy of the whole kept, two free masses keeping their total momentum where a wire one way does not, four refusals), catalogue 60 (the shell model's two invariants as identities of its field, the balance of power, the flux through two shells, the exponent of its spectrum between 2/3 and 0.9, the growth of flatness, a positive exponent), language 40 (a part's velocity read and written back).
- **Found while building.** The shell model's nonlinear term keeps its two quadratic quantities to rounding as an identity of the field, which holds every sign and coefficient to the model's; along a path the same is true only to the integrator's error, so the identity is checked on the field. A path started from zero stays on the forced shell, an invariant subspace: a start must move more than one shell. With a larger viscosity the last shell is stiff and the step must shorten, or the run leaves the finite numbers, which the answer says.
- Parts joined by an action are solved as terms for at most three coordinates in all, as any action is.
- The catalogue holds thirty models.
- Not built: the response-field action; indexed families; fields that run, on which a fluid, a concentration it carries and animals moving through both would be coupled as fields.
