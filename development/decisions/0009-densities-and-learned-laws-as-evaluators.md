# 0009. Densities Moved by a Generator, a Value Equation Whose Policy Drives Them, and Learned Laws as Evaluators

- **Status:** Accepted, provisionally (2026-11-21).
- **Date:** 2026-11-21
- **Decided by:** Claude, at the infrastructure's Stage 7 of evolutionary allocation (its ADR 0139), after its survey's Parts VII and VIII (*Evolutionary heterogeneous resource allocation*): declared density laws give the exact effect of a sequence of interventions under their assumptions; learned laws give rankings that must be held to baselines and withheld outside the record's support.

## 1. Context

`balance` holds balance laws of a density on a graph whose coefficients are constants, linear in what is learned. A density whose drift and diffusion vary from cell to cell (Bouchaud and Mézard's wealth, a policy's drift), in types between which it switches, is moved instead by a Markov generator, and a value equation solved by the upwind scheme gives the drift. Laws learned from a record (`induction`) had no discipline for being used to evaluate interventions: no baselines, no support, no reading across an ensemble.

## 2. Decision

A new module, `densities`:

- `generator.hpp`: `Grid`; `Generator` and `generator(grid, drift, diffusion, switching)`, the upwind chain (Kushner and Dupuis 2001), rates out of a state removed where they would leave (reflection), rows summing to zero; `Flow` and `flow(G, dt)`, the implicit step (I − dt Aᵀ) factorised once by calx's band LU (its decision 0020), keeping mass and positivity; `stationary(G)` by inverse iteration with a small shift; readings (`total`, `mean`, `below`, `type_mass`, `gini`); interventions as maps on densities, `carry` (split between the two nodes about the destination: mass and mean kept within the grid, what is carried beyond reported) and `transfer`;
- `value.hpp`: `Stock` (a stock with an inflow in each regime, a return and a control drawn, valued by CRRA utility and discounted; economics reads it as saving) and `solve`, Achdou et al.'s implicit upwind scheme with the state constraints at both ends; its generator is the one the density's forward equation transposes, exactly.

`induction` gains `evaluators.hpp`:

- `held_to_baselines`: a law's rolling one-step forecasts against tyche's random walk, historical mean and fitted autoregression by the test of equal accuracy (tyche's decision 0031), with a verdict that says when the law is beaten and when it is not shown better;
- `Support`, `support_of`, `outside`: the box of a record's ranges, a point's variables beyond it named;
- `Member`, `Run`, `run`, `Ranking` and `rank`: an ensemble of laws evaluating two sequences of interventions, the claim that the first is better holding, failing, arguable or undecided across the laws, undecided when either leaves the support or the order is immaterial under every law, with the difference by law;
- `bagged`: sequentially thresholded least squares on resamples of the rows, each term's inclusion probability (Fasel et al. 2022).

## 3. Alternatives Considered

- **Extending `balance`.** Its laws are linear in coefficients and learned by least squares; a generator with variable coefficients and a value equation are another object, run differently.
- **Reference values from the authors' codes.** Achdou et al.'s are MATLAB and were not run; the checks use the closed form of the permanent-income control and the propositions on the mass at the limit, and say so.
- **A support finer than a box.** A point inside the box but far from the record is not flagged, and the header says so; a nearest-neighbour or convex-hull support is a later refinement.

## 4. Consequences

`check_kinesis_densities` (23 checks): rows summing to zero; mass kept over 300 implicit steps; the mean moving exactly as the step says under a linear drift, the rate reflection removed counted (to 10⁻¹⁴); Bouchaud and Mézard's inverse gamma reached at first order in the step; the order of a transfer and a flow δ(1 − q^N) exactly; translations commuting and transfers below a line not; the Gini coefficient; the permanent-income control within 1%, nearer as the grid is refined; with two regimes the lower holding at the limit with a mass there (by extrapolation in √h), the higher growing, the regimes' masses the switching chain's, the stationary density solving Aᵀg = 0. `check_kinesis_evaluators` (17 checks): baselines' verdicts (a true law not beaten, a fitted autoregression not shown better than itself, a wrong law beaten by every baseline); the four standings over laws x' = a x + u with the difference a³(1 − a) exactly; the order immaterial under the integrator and not under a law whose effect depends on its state; inputs and states beyond the support flagged at their steps; bagged fits keeping the true terms. The economics domain does not use these; the allocation domain's `evaluators` does.
