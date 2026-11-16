# 0005. Balance Laws on Graphs, Run and Learned from a Few Snapshots

- **Status:** Accepted, provisionally
- **Date:** 2026-11-16
- **Decided by:** the user ("Please utilize dynamical system library to learn population dynamics with various models, like fluid dynamics, and simulate, visualization, in the urban studio") and Claude; the infrastructure's ADR 0124

## 1. Context

A city's population recorded cell by cell (GHS-POP's 1 km grid, every fifth year since 1975) is a field of thousands of components known at a dozen times. The library could not learn its laws. Its declarations are for small systems: no families of components, no coupling of neighbours, no sum over a graph (decision 0004 deferred indexed families). Its fits read derivatives from a record, which a dozen snapshots cannot give, and its dynamic modes fit an n × n operator, which needs more snapshots than components. The theory note listed "a discretised fluid" as not built.

## 2. Decision

1. **A module `balance`**: balance laws of a density on a graph of cells, in finite volumes. A law is a sum of terms with coefficients: growth (linear, logistic, Gompertz's), a constant, a source in a covariate; diffusion, crowding (the porous-medium flux) and drift up a covariate's potential, each conservative by construction. A square grid with a mask is the first graph, so that a coast is a boundary of no flux.
2. **Runs** by the classical fourth-order Runge–Kutta step, its step bounded by Gershgorin's estimate of the stiffness; a density made negative is set to zero and the quantity added is reported.
3. **Learning from a few snapshots.** The integral form by the trapezoid (linear least squares, or thresholded with induction's stlsq on scaled columns); then shooting, the law run from the first snapshot and fitted to the rest by Levenberg–Marquardt. A diffusion or crowding found negative is set aside (it is ill-posed), and said to be.
4. **Judgement on snapshots kept aside**, against carrying the last one forward: the skill, and the skill of the pattern once each forecast's total is set to the record's. The two separate how much from where.
5. **Baselines beside the laws**: a rate of its own for each cell, and the record's dynamic modes by the method of snapshots (exact DMD on the snapshots' span).
6. **The name.** The README had planned `fields` for statistical fields that run, which is another thing; this module is `balance`, and `fields` stays planned.

## 3. Alternatives Considered

- **Indexed families in the grammar** (`state u[cells]` with a neighbour sum). The general answer, and a change to the grammar, the reader, the writer, `form_of` and every tool. A balance law on a graph is a structured class that needs none of it, and its laws remain readable as sums of named terms. The grammar's families stay wanted.
- **physis's FieldKinetics and GraphTransport for the runs.** They hold reaction–diffusion with any kinetics and upwind transport. Learning needs the law linear in its coefficients and its fluxes consistent with its regression, which the centred drift gives and the upwind scheme does not (its flux changes form with the sign of the coefficient). The operators are a few lines on a graph and are written here; physis remains where the laws of physics are written.
- **Derivatives of the record.** Five years apart they are not derivatives.
- **Plain DMD.** It needs more snapshots than cells.

## 4. Consequences

- 31 checks of the module; kinesis's 464 in all.
- The infrastructure's API composes the module (`science/population.hpp`) and the Urban Studio's plug-in *Population dynamics* works with it.
- **What the first records show:** at 1 km, Melbourne's and Lagos's fields hold no diffusion: one is found negative and set aside, or too small to move anything. A city's density is pinned by its buildings, and a law that would flatten its centre is not in its record. The pattern of change is placed best by the record's modes in Melbourne, and by a drift along the settled landscape in Lagos, slightly better than uniform growth.
- **Not built:** graphs other than grids (zones of any shape, with their faces' widths) are supported by the operators but not yet built from data; upwind drift for runs far from regression; an implicit step for stiff laws; spatially varying coefficients; standard errors by a spatial bootstrap.
