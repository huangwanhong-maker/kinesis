# A Catalogue of Models

Status: **draft** (2026-11-13). Module `kinesis::catalogue`.

## 1. Concepts

Many fields have written down laws that recur: a pair of equations that oscillates, a map that doubles its period, a neuron that fires and is reset, a process pulled toward a level, a grammar of nested words. A model of the catalogue is such a form of law, written as a declaration in the language of dynamical systems, with the values of its parameters in common use and with where it comes from.

A model serves two purposes. It is a starting point for a study: it is read, run and analysed as any declaration is. And it is a hypothesis for a record: with its parameters freed (induction's `as_template`) it is a template, and a record is matched against the catalogue by fitting each template whose shape the record could have and ranking them (induction's `match`).

A record's components have no names a catalogue knows. A model's states are therefore named anew after the record's, in each order that is tried: all orders for three components or fewer, the record's own order above that. Theory: [docs/theory.md](../docs/theory.md).

*Caveat.* Nothing in an entry names a place or a measurement. The values given are those in common use, which are not always those of the first paper, and an entry says so where it matters. Measured values belong to the domain that owns the object a model describes.

## 2. Content

- `models.hpp`: `Entry{id, title, field, about, source, text}`; `entries()`, the thirty models; `find(id)`; `system(entry)` and `system(id)`, the declaration an entry holds; `states_of(s)`; `plain(s)`, whether a model is a flow or a map given by first-order equations alone; `discrete(s)`; `renamed(model, names)`, a plain model with its states named anew, a parameter that shares a name with a new state moved aside first; `candidates(names, maps)`, the templates a record of those components could be a path of.
- `turbulence.hpp`: `ShellModel{shells, nu, forcing, k0, ratio, a, b, c, flux}`; `wavenumber(m, n)`; `shell_model(m)`, the Sabra model written out as a declaration, two real components a shell, with what is observed of it (each shell's energy, the whole energy, the second quadratic quantity, the power put in and dissipated, the flux of energy through the shells named); `start_of(m)`, a start with every shell moved, since a path from zero stays on the forced shell.

The models, by form:

| Form | Models |
|---|---|
| a flow of twenty-eight components | `sabra`, a shell model of turbulence |
| flows of three components | `lorenz`, `rossler`, `chen`, `thomas`, `sir`, `hindmarsh_rose` |
| flows of two | `van_der_pol`, `fitzhugh_nagumo`, `lotka_volterra`, `brusselator`, `selkov`, `duffing`, `stuart_landau` |
| a flow of one | `logistic_growth` |
| maps | `logistic_map`, `henon`, `standard_map` |
| flows with events | `integrate_and_fire`, `izhikevich` |
| systems given by an action | `pendulum`, `kepler` |
| systems with noise | `ornstein_uhlenbeck`, `double_well`, `geometric_brownian`, `square_root_diffusion` |
| open quantum systems | `driven_atom`, `jaynes_cummings` |
| rules | `anbn`, `balanced` |

*Caveat.* Only plain models are offered as candidates for a record: a model with events, noise, inputs or an action is not named anew. A ranking is among the candidates tried, and says nothing of a law that is not in the catalogue.

## 3. Checks

60: every entry declares one system of its own name, of the form it is listed under, and is written out and read back the same. Then what is known of each. The exponents of Lorenz's, Chen's and Thomas's systems sum to their constant divergences, Hénon's to ln b and the standard map's to zero; Rössler's two equilibria, the Brusselator's with the real part (b − 1 − a²)/2 of its eigenvalues, Sel'kov's and Duffing's three are at their closed forms; the epidemic keeps s + i + r and its final size solves ln(s/s₀) = −(β/γ)(1 − s); Van der Pol's cycle at μ = 0.1 has the amplitude 2 and the period 2π(1 + μ²/16); Lotka and Volterra's quantity is kept; Hopf's normal form ends on the circle of radius √μ; Verhulst's growth follows its closed form; the logistic map at r = 4 has the exponent ln 2; the integrate-and-fire neuron fires with the period τ ln((RI + E − v_r)/(RI + E − θ)); the pendulum's small swings have the period 2π√(l/g), and Kepler's circular orbit closes after 2π with its energy and angular momentum kept; Ornstein and Uhlenbeck's stationary variance is σ²/(2θ), within its standard error; the driven atom's stationary excited population is Ω²/(γ² + 2Ω²); balanced brackets of 2n letters number as the Catalan numbers, each with one history. The shell model: its nonlinear term gives the energy and the second quadratic quantity no rate of change at five states, to rounding, and the energy is kept along a path; in its steady state the power put in is the power dissipated within three of their errors, the flux of energy through shells 3 and 6 is the power put in less what the shells up to each dissipate, the energy of the shells falls as a power of the wavenumber steeper than 2/3 and less than 0.9 over eight octaves, the flatness grows from about 2 at shell 4 to half as much again or more at shell 9, and the largest exponent is positive. And the templates: thirty-six candidates for a record of three components, fourteen for two, four for a map of two; a model whose parameters share the record's names has the same field when named anew; and a record of Rössler's system, under other names and in another order, is matched to Rössler's among all thirty-six, with its three numbers recovered.

## 4. Dependencies

`kinesis::language` (and through it `kinesis::terms` and `kinesis::kernel`). The checks use `kinesis::systems`, `kinesis::analysis` (its exponents and its averages) and `kinesis::induction`, physis's `quantum` and calx's `verification`.

## 5. References

Each entry names its source: L'vov, Podivilov, Pomyalov, Procaccia and Vandembroucq (1998) and the review of Biferale (2003) for the shell model; Lorenz (1963); Rössler (1976); Chen and Ueta (1999); Thomas (1999); Kermack and McKendrick (1927); Hindmarsh and Rose (1984); van der Pol (1926); FitzHugh (1961) and Nagumo, Arimoto and Yoshizawa (1962); Lotka (1925) and Volterra (1926); Prigogine and Lefever (1968); Sel'kov (1968); Duffing (1918); Landau (1944) and Stuart (1960); Verhulst (1838); May (1976); Hénon (1976); Chirikov (1979); Lapicque (1907); Izhikevich (2003); Lagrange (1788); Newton (1687); Uhlenbeck and Ornstein (1930); Kramers (1940); Samuelson (1965); Cox, Ingersoll and Ross (1985); Breuer and Petruccione (2002); Jaynes and Cummings (1963); Chomsky (1956); Chomsky and Schützenberger (1963). In the checks: Strogatz (2015), *Nonlinear Dynamics and Chaos*, 2nd ed., for the series of Van der Pol's period.
