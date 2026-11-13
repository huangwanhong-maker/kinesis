# Systems of Systems

Status: **draft** (2026-11-13). Module `kinesis::coupling`.

## 1. Concepts

A composite declares parts, each a system of the script, and wires: `connect { a.I = 0.5 * (b.v - a.v) }` gives the input I of the part a as a term of the parts' states (written part.state), the composite's parameters and the time. A part's input that no wire gives keeps its own value. There are two ways of running a composite, and each says what it is.

Merged: when every part is a flow, or an action through its stationary paths, the composite is itself a flow, and its declaration is built, with each part's states renamed part_state, each input replaced by the term wired to it, and the parts' parameters inlined as numbers. The result is a system of the language; no error is added by the composition, and every tool applies to the whole. Co-simulated: parts of different forms have no common presentation to merge into, and are advanced side by side. At each communication step H the wires are evaluated on the states of all parts, each part's inputs are held at those values, and each part is advanced over H by its own method. Theory: [docs/theory.md](../docs/theory.md).

## 2. Content

- `composite.hpp`:
  - `merged(script, composite, as)`: the flow a composite of flows is, as a declaration (a `System`), with the parts' equations, events and observations in the merged names.
  - `Cosimulation`: its `parts` (each a `Part{name, form, system, offset, dimension}`, the system a `Dynamics` or a `Stochastic`), the `names` of every component (part.state), `dimension()`, `index(qualified)`; `Options{t0, t1, communication, step, every, seed, event_tolerance}`; `run(x, o)`, a `systems::Trajectory`. A flow is advanced by Runge–Kutta steps with its events, a map by one iteration, a stochastic differential by Euler–Maruyama steps from its own streams.
  - `cosimulation(script, composite)`: the co-simulation of a composite of flows, maps and stochastic systems.

*Remark.* The order of the parts does not matter: the wires are all read before any part moves. A map iterates once in each step H, so H is the time between its iterations.

*Caveat.* Only flows are merged into one declaration; a part of another form is refused there, with the advice to co-simulate. A part is a declared system: one that is itself derived or composed is refused.

*Caveat.* Co-simulation is the explicit (Jacobi) scheme with inputs held constant. Its error is of the first order in H where a wired signal varies within a step; it is exact where the signals are constant between steps, as when a map drives a flow. An input of a stochastic part is not wired yet, and parts of forms other than flows, maps and stochastic systems are refused.
- Joined by an action, in the same header: `joined_by_action(composite)`; `joined_action(script, composite, as)`, the declaration by an action whose Lagrangian is the sum of the parts' and of the composite's interaction, its coordinates named part_coordinate; `merged` varies it when a composite gives an action. In the interaction a part's state is `a.q` and its velocity `a.q'`.

*Remark.* A coupling given by an action is reciprocal and the whole keeps what its Lagrangian is invariant under; a wire need not. A wired composite is the stationarity of the response-field action of Martin, Siggia and Rose (1973), which is recorded in the header and not built.

## 3. Checks

32: merged, a response driven by Lorenz's x (Pecora and Carroll 1990) as one flow of five components, the response falling onto the drive and the composite's exponents holding the drive's three and the response's two conditional ones, −1.81 and −1.86 (their Table I, the x drive); the merged declaration, written out, reading back with the same field; parameters set on a part being those of the part, and an input that is not wired keeping its own value; events of the parts as events of the whole, two integrators that fire, each with its period. Co-simulated, two flows side by side against their merged flow, the error falling by 2 when the communication step is halved, and the order of the parts not mattering; a map that drives a flow exactly, a counter n and x' = n giving x(KH) = H K (K − 1)/2; an Ornstein–Uhlenbeck process through a filter of rate a with the variance σ²/(2θ) · a/(a + θ). Refusals with their reasons: a part the script does not declare, a wire to what is not an input, a map merged.

## 4. Dependencies

`kinesis::systems` (`dynamics.hpp`, `lagrangian.hpp`, `simulate.hpp`, `stochastic.hpp`), and through it `kinesis::language`, `kinesis::terms` and `kinesis::kernel`. The checks use `kinesis::analysis` (Lyapunov's exponents of the merged flow) and calx's `verification`.

## 5. References

Kübler and Schiehlen (2000), Math. Comput. Model. Dyn. Syst. 6, 93–113; Gomes et al. (2018), ACM Comput. Surv. 51, 49. In the checks: Pecora and Carroll (1990), Phys. Rev. Lett. 64, 821–824.
