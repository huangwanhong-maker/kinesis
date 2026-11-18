# 0007. Systems Across Scales: Quantities on a Hierarchy of Units, Carried Between Levels by Operators, Advanced Each at Its Own Pace

- **Status:** Accepted, provisionally
- **Date:** 2026-11-18
- **Decided by:** the user ("make urban science studio can build multi-scale dynamical system model(s) and visualize them on top on them. Also build estimators for some indexes if can … if cannot well estimate, them should take up when apply model") and Claude; the infrastructure's ADR 0128

## 1. Context

A city's indexes live at many scales at once: a block's floor space, a district's vacancy, a city's total. Models of them are multi-scale in two ways the literature treats together (the infrastructure's survey *Indexes of cities across scales*, Part I): levels of units (blocks within districts within a city), each with laws of its own, coupled by aggregation upward and allocation downward under control totals; and paces (fast flows equilibrating within slow change of structure, as in Harris and Wilson's retail centres, or laws stepped yearly beside laws stepped monthly). kinesis held systems of systems (decision 0003's `coupling`: flows merged, parts co-simulated) and balance laws of one density on a graph (decision 0005), but no family of components over a hierarchy of units, no operators between levels, and no record of how each input is known. The language reads no indexed families (the theory note, §11).

## 2. Decision

1. **A module `scales`**, on the language (an expression read as a term), the terms (compiled programs), calx's linear algebra and tyche's uncertainty, descriptive statistics and normal quantile.
2. **A hierarchy** is a chain of partitions, the finest level first; levels may carry measures, a graph of neighbours and pairwise costs.
3. **Quantities on levels**, states, drivers or fixed, extensive, intensive or ordinal. **Every crossing of a level is a derived quantity by a named operator** (value; sum, mean, weighted mean, minimum, maximum, reduction by a registered measure, over the units beneath; parent; share, logit and z-scores among siblings; Laplacian and neighbours on the graph; the production-constrained inflow). A law that names a quantity of another level is refused, with the operator that would bring it.
4. **Laws** are flows over their own steps, maps over their own periods, or equilibria solved at every instant (a fast variable at its quasi-steady state). Paces are exchanged at the start of each step (Jacobi); one shared step gives Runge–Kutta on the whole.
5. **Provenance on every input** (observed, estimated with a standard error, declared with a range, modelled). **What a model needs and has not been given is refused**, item by item with its level and count, never filled in.
6. **Ensembles over what is declared or estimated** (a Latin hypercube; one factor per quantity, comonotone across its units), with **declarations ranked** by Morris's screening.
7. **Unknowns fitted** to observations at any level by Levenberg–Marquardt on runs; an unknown the record does not identify (a relative standard error above a half, a correlation above 0.95 with another, a singular normal matrix) is reported so; the fit is judged on what it was held from against persistence.
8. **Perfect aggregation** of linear laws decided by A M = M̃ A.
9. **Declarations as structures, not grammar.** A system across scales is a C++ structure (and JSON in the infrastructure's API); its laws and operators' arguments are expressions of the language. A clause of the grammar for indexed families is left for when the form has settled.

## 3. Alternatives Considered

- **Families of components in the grammar** (`state x[i in block]`): the language's natural extension, but it changes the grammar, the reader, the writer and `form_of` together (CLAUDE.md §5.2), for a form not yet settled. The structure comes first; the clause can read into it later.
- **Extending `balance`**: balance laws are one density on a graph, with fitted coefficients; a hierarchy of levels with operators between them and laws of any form is a different object, which may hold balance laws as one of its laws (the Laplacian operator does).
- **Measures (Gini, HHI, entropy) written here**: they are nomos's, which kinesis may not depend on (invariant 1). The module owns the interface (a reducer of values and weights), and a caller registers nomos's measures.
- **Gauss–Seidel or extrapolated multirate schemes** (MRI-GARK): of higher order, and the survey ranks them; the Jacobi scheme is what `coupling` already states and checks, and is sufficient for laws stepped yearly beside fast equilibria. Higher order is left for a case that needs it.
- **Profile likelihood for identifiability**: exact where the linearisation is not, at the cost of a fit per point of each profile; the linearised gate is the first form, and says that it is one.

## 4. Consequences

- The infrastructure's API composes the module (`api/cpp/science/scales.hpp`) with nomos's measures as reducers and tyche's estimators, for the Urban Studio's plug-in *Cities across scales*.
- 66 checks; kinesis has 562.
- Not built: indexed families in the grammar; Gauss–Seidel and higher-order multirate schemes; the profile likelihood; history matching and ABC; a Kalman filter on a tree; reconciliation across levels (MinT); discrete events of development (an allocation with capacities); the Bayesian inference of Harris–Wilson.
