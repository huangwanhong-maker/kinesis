# Systems

Status: **draft** (2026-11-13). Module `kinesis::systems`.

## 1. Concepts

A declaration means something before it runs: its form, and the names its terms may use. The form is read from what its `change` holds and how its `structure` reads it: rules make a rewriting system, whose histories are derivations; a Hamiltonian or jumps an open quantum system; an action with stationary histories a Lagrangian system; an action with weighted histories a statistical field; equations with a source of noise a stochastic system; equations in discrete time a map; equations otherwise a flow; parts a composite; a derivation a system derived from another. A scope resolves the names of a declaration once, and a name that a term uses and the scope does not hold is an error of the script, reported with its place.

Each form is then made into a system that runs. A flow or a map has its field and its exact Jacobian from the terms; a stochastic differential is written as it is on paper and read as a deterministic map of state and noise, each source with a stream of its own; an action is varied, exactly, into the equations of its stationary paths, which are a declaration of a flow again; rules give the graph of their histories; a Hamiltonian and jumps give physis's Lindblad generator. Theory: [docs/theory.md](../docs/theory.md).

## 2. Content

- `model.hpp`: `State`; `Form{empty, flow, map, stochastic, rules, lagrangian, field, quantum, composite, derived}` with `name(Form)`, `discrete_time(s)`, `form_of(s)`; `Source`, `InitialLaw`; `Scope` with `is_state`, `param_index`, `declares`, `resolve` (functions inlined, `pi` and `e` given their values unless declared) and `constant`; `scope(s)`; `component(sc, t, i, prefix)`, the i-th component of a term over vectors.
- `dynamics.hpp`: `Dynamics`, a flow x' = f(t, x; p, u) or a map, with `field`, `jacobian`, `jacobian_vector`, `jacobian_terms` (the entries that are not zero), `observe`, `conserved`, `guard`, `reset`, `set_param`, `param`, `compile`, `value`, and its `Event`s; `dynamics_of(s, sc, equations)`, `dynamics(s)`. An equation of order k makes the components q, q', …, q⁽ᵏ⁻¹⁾; a state in real^n is n components.
- `simulate.hpp`: `Occurrence`, `Trajectory` (with `column(i)`), `Options`; `advance` (one step of a flow, through its events), `iterate` (one iteration of a map), `simulate(d, x, o)`, `initial(d, values)`. A flow is stepped by physis's classical Runge–Kutta method at a fixed step.
- `stochastic.hpp`: `Stochastic` with `drift`, `diffusion`, `step` and `diagonal`; `stochastic(s)`, which reads the coefficients off the right side of `dX = a*dt + b*dW` and converts a Stratonovich equation once to Itô's; `draw(source, stream)` for the laws normal, uniform, exponential, poisson, cauchy and levy; `StochasticOptions`, `StochasticTrajectory`, `simulate(d, x, o, increments)` by Euler–Maruyama or Milstein.
- `lagrangian.hpp`: `Lagrangian{name, time, coordinates, lagrangian, momenta, mass, force, acceleration, energy, autonomous}`; `vary(s)`, the Euler–Lagrange equations M q'' = F as terms; `varied(s, as)`, the flow of the stationary paths as a declaration, with the energy declared conserved when the Lagrangian does not hold the time.
- `rewriting.hpp`: `Word`; `Rewriting` with `steps(w, leftmost)`, `accepts(w)`, `terminal`, its `grammar` as morphe's, `context_free` and `monotone`; `rewriting(s)`; `Histories`, `Bounds`, `explore(r, b)`; `Sentence`, `language(r, h)`; `sample(r, stream, max_steps)`.
- `quantum.hpp`: `Quantum` with `factors`, `generator`, `basis(levels)`, `expectations(rho)`; `quantum(s)`, from `state rho : density(2, 12)` and terms of the operators `sx`, `sy`, `sz`, `sm`, `sp`, `a`, `adag`, `num`, `id`; `QuantumTrajectory`, `evolve(q, rho0, t1, step)` by the exact propagator of one step.

*Caveat.* An event of a flow is a crossing, not a level: a guard that is zero or beyond before a step does not fire. Its time is found by bisection, and is known to `event_tolerance` and to the integrator's error. Events of a map are levels.

*Caveat.* Milstein's scheme is given for diagonal noise and refused otherwise, since the general scheme needs iterated integrals that are not sampled here. `vary` solves for the accelerations as terms for at most three coordinates; a Lagrangian whose M is singular as a term is refused.

*Caveat.* Histories are explored within bounds, and nothing is reported as the whole language unless the exploration was complete; membership is decided for context-free rules only. A product of two operators on one factor is refused, since terms are commutative: its order is written `ordered(adag, adag, a, a)`. The state of a quantum system has dimension at most 64.

*Caveat.* `form_of` names two forms that nothing here runs: `field`, a system whose histories are weighted by its action, and `derived`, for the derivations other than `vary`. A tool given one refuses with the reason.

## 3. Checks

64: Lorenz's field and its Jacobian as terms, 8 entries that are not zero, equal to central differences; an equation of the second order as two components, the harmonic oscillator's error falling by 16 when the step is halved and its energy kept to that order; parameters with arguments inlined; the logistic map at r = 4 following x_n = sin²(2ⁿθ) (Ulam and von Neumann's conjugacy); a ball dropped from h₀ striking at √(2h₀/g) and then at intervals 2eᵏv/g, and a leaky integrator firing with the period τ ln(I/(I − 1)); the coefficients of a stochastic differential read off its right side, one that is not linear refused; dX = σX∘dW becoming Itô's with the drift σ²X/2; on geometric Brownian motion against the exact solution on the same Brownian path, the strong error of Euler–Maruyama falling as h^½ and Milstein's as h; the Ornstein–Uhlenbeck variance σ²/(2θ), a path depending on (seed, path) alone; an autoregression's variance 1/(1 − a²); the pendulum's action giving q'' = −(g/l) sin q as a term, the derived system reading back and keeping its energy, the double pendulum keeping its energy, a central force its angular momentum, a degenerate action refused; a spin and an oscillator at resonance with populations cos²(gt) and sin²(gt), a driven, decaying atom reaching Ω²/(γ² + 2Ω²), `ordered(adag, adag, a, a)` equal to num² − num; aⁿbⁿ within a bound with one history each and P(aⁿbⁿ) = 0.7ⁿ⁻¹·0.3, an ambiguous grammar with 2 and 5 histories (Catalan's numbers), aⁿbⁿcⁿ by rules with context; and refusals with their place.

## 4. Dependencies

`kinesis::language` (and through it `kinesis::terms` and `kinesis::kernel`); physis's `dynamics` (`flows.hpp`: the Runge–Kutta step) and `quantum` (`lindblad.hpp`); calx's `linalg`, `random` (streams keyed by seed, path and source) and `integrators`; morphe's `formal` (`grammar.hpp`, `earley.hpp`). The checks use calx's `verification`.

## 5. References

Arnold (1989), *Mathematical Methods of Classical Mechanics*, 2nd ed., §12–13 and §20; Kloeden and Platen (1992), *Numerical Solution of Stochastic Differential Equations*, Eq. 4.9.15, §10.2 and §10.3; Chambers, Mallows and Stuck (1976), symmetric stable draws; Mjolsness and Yosiphon (2006), a process semantics for dynamical grammars; the idea of 2026-11-12, §3.1 and §6. In the checks: Jaynes and Cummings (1963).
