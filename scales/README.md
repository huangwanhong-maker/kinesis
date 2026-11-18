# Systems Across Scales

Status: **draft** (2026-11-18). Module `kinesis::scales` ([decision 0007](../development/decisions/0007-systems-across-scales.md)).

## 1. Concepts

A system across scales lives on a **hierarchy of units**: levels, the finest first, each unit of a level below the top with one parent in the level above, so that the levels are a chain of partitions (cells within districts within a city). A level may carry what its units are measured by, a graph of neighbours with conductances, and matrices of pairwise costs.

A **quantity** has a value for every unit of its level. It is a state (advanced by a law), a driver (given over time) or fixed. Its kind says how it is carried up: an extensive quantity is summed, an intensive one averaged with a declared weight, an ordinal one read through thresholds on a latent value.

What crosses a level, or depends on more than one unit, is a **derived quantity**, made by one of a few operators, so that no law names another level by accident: a value of the level's own quantities; a sum, mean, weighted mean, minimum, maximum or reduction (by a measure the caller registers) of the units beneath; a parent's value read by each unit beneath it; a share or a logit among siblings (allocation of a parent's total); z-scores among siblings; the graph Laplacian and the mean of neighbours; and the inflow of the production-constrained spatial interaction model. A model in which coarse quantities are derived up and fine ones allocated down agrees across its levels by construction, and a ratio defined on summed quantities is the ratio of sums at every level.

A **law** advances one state: as a flow over its own step, as a map once a period, or as an equilibrium solved at every instant (the quasi-steady state of a fast variable, Tikhonov 1952; Fenichel 1979). Its expression is a term of the language. Laws with different steps advance side by side and exchange their states at the start of each of their steps (the explicit Jacobi scheme of co-simulation); with one step the whole is one system under Runge–Kutta's classical step.

Every value a model is given says **how it is known**: observed, estimated (with a standard error), declared (with a range) or modelled. A model **lists what it needs and has not been given**, and refuses to run with each item named, its level and the units that lack it: a quantity no record holds and no estimator gives is declared, never filled in. Declared ranges and estimated errors are **drawn on a Latin hypercube** and the outputs summarised by quantiles; the declarations are **ranked by Morris's screening**. Unknowns are **fitted** by Levenberg–Marquardt on runs from the initial values, to observations at any level; an unknown the record does not identify is reported so, and the fit is judged on what it was held from against persistence. Whether a linear fine law **aggregates perfectly** under a partition is decided by Iwasa, Andreasen and Levin's condition. Theory: [docs/theory.md](../docs/theory.md) §9.

## 2. Content

- `hierarchy.hpp`: `Level{name, units, parent, measures, edges, costs}`; `Hierarchy` with `children`, `parent`, `ancestor`, `groups` (siblings), `group`, `neighbours`, `measure`, `costs`, refusing a parent that does not exist, a coarse unit with nothing beneath it and costs of the wrong size.
- `model.hpp`: `Quantity`, `Derived`, `Law`, `Parameter`, `Model` (with its `step` and its registered `reducers`); `Kind`, `Role`, `Op`, `Form` and their names; `compile(hierarchy, model)`, every name resolved, every expression compiled for its level, the derived quantities ordered by dependence, refusing with the way to mend it a law that names another level, derived quantities in a circle, a state without a law, a step that is not a whole multiple, a reducer not given.
- `run.hpp`: `Source`, `Value` (with `centre` and `uncertain`), `Series`, `Inputs`; `missing(compiled, inputs, t0, t1)` and `sentence`; `Options`, `Run`, `Stepper`; `simulate(compiled, inputs, options)`.
- `ensemble.hpp`: `Factor`, `factors(inputs)`, `at_quantiles`; `ensemble(compiled, inputs, options, {draws, seed, low, high, outputs})`, quantile bands; `screen(compiled, inputs, options, output, unit)`, Morris's μ*, μ and σ for each declaration.
- `fit.hpp`: `Observation`, `FitOptions{held_after, most}`, `Fitted` (values, standard errors, correlation, condition, what is not identified, predictions, skill against persistence on what was held aside, the inputs with the unknowns set as estimated); `fit(compiled, inputs, options, observations)`.
- `aggregation.hpp`: `perfect_aggregation(M, A)`, the coarse matrix M̃ and the relative residual of A M = M̃ A.

*Caveat.* The multirate scheme is of the first order in a slow law's step wherever a fast state varies within it; laws that must agree more closely share a step. An equilibrium is found by damped iteration of x = g, which converges when g is a contraction about the solution; a run counts the instants at which it did not, and an ensemble leaves those draws out.

*Caveat.* Standard errors come from the Jacobian at the optimum (a local, linearised statement), and the Jacobian is by forward differences of whole runs. A fit that stops at a bound, or whose correlation matrix is singular, reports its unknowns as not identified rather than giving numbers.

*Caveat.* A declared range is drawn as uniform, a standard error as normal, and the values of one quantity at one quantile across its units (comonotone): assumptions, said where they are used.

## 3. Checks

66: the hierarchy's children, groups and ancestors, and its refusals; sums across three levels, a ratio index on sums equal to the ratio of sums and to the ratios weighted by their denominators, and not to the mean of ratios; exponential growth in every block with the city's total Σ x₀ e^{rt}; a logistic total allocated by fixed shares, each block w_i/Σw of it and the sum equal to the total to rounding; a map of period two steps; diffusion between two cells, the difference falling as e^{−2Dt} and the content kept; two laws at two steps against the exact exponential of their system, the error halving with the slow step, and Runge–Kutta's error with one step; a fast variable at rest giving e^{−t/2} exactly, and as a fast flow the error falling in proportion to ε; the inflow against the formula by hand, its total conserved, and β = 0 sharing by attractiveness; Harris and Wilson's centres reaching D_j = κW_j with ΣW = ΣO/κ; logit shares and z-scores among siblings; a registered weighted median of ordinal levels read through steps; what is missing named with its level and count, an equilibrium needing no initial value; an ensemble over a declared rate, its median and quantiles e^{0.2}, e^{0.11}, e^{0.29}; a quantity's units drawn at one quantile; screening ranking the wider declaration first; a logistic's r and K recovered exactly, skill against persistence on what was held aside, noisy records within four standard errors, two unknowns that enter only as a product reported as not identified; a lumpable law aggregating perfectly into its coarse matrix and a perturbed one not; and the compiler's refusals.

## 4. Dependencies

`kinesis::language` (an expression read as a term), `kinesis::terms` (terms compiled to programs), `kinesis::kernel` (refusals); calx's `linalg`; tyche's `uncertainty` (Latin hypercubes, Morris), `descriptive` (quantiles) and `probability` (the normal quantile). The checks use calx's `random`.

## 5. References

Wilson (1971), Environment and Planning 3, 1–32; Harris and Wilson (1978), Environment and Planning A 10, 371–388; Tikhonov (1952), Mat. Sb. 31, 575–586; Fenichel (1979), J. Differ. Equ. 31, 53–98; Kübler and Schiehlen (2000), Math. Comput. Model. Dyn. Syst. 6, 93–113; Iwasa, Andreasen and Levin (1987), Ecological Modelling 37, 287–302; Kemeny and Snell (1960), Finite Markov Chains, §6.3; McKay, Beckman and Conover (1979), Technometrics 21, 239–245; Morris (1991), Technometrics 33, 161–174; Levenberg (1944), Q. Appl. Math. 2, 164–168; Marquardt (1963), J. SIAM 11, 431–441. The infrastructure's survey *Indexes of cities across scales* (2026) collects what the module rests on and what it does not yet do.
