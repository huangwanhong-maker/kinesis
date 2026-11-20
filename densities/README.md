# Densities Moved by a Generator

Status: **draft** (2026-11-21). Module `kinesis::densities` ([decision 0009](../development/decisions/0009-densities-and-learned-laws-as-evaluators.md)).

## 1. Concepts

A density here is the mass of a population over a line of positions, cut into cells, in types between which it switches: units over their holdings, in regimes of their inflow. It is moved by a diffusion whose drift and diffusion vary from position to position, dX = μ_k(X) dt + √(2D_k(X)) dB in type k, and by switching from type k to l at rate λ_kl. On the grid this is a continuous-time Markov chain whose rates out of a cell go up by max(μ, 0)/h + D/h², down by max(−μ, 0)/h + D/h², and across to another type by λ_kl: the upwind scheme, locally consistent with the diffusion (Kushner and Dupuis 2001). Its generator has non-negative rates off the diagonal and rows summing to zero, and the density's forward equation is dg/dt = Aᵀg, which keeps the mass exactly. Rates that would leave the grid are removed, a reflection; the mean moves by Σ g μ exactly where no rate was removed, the diffusion's two rates cancelling in it.

Time is stepped implicitly, (I − dt Aᵀ)g′ = g, an M-matrix that keeps mass and positivity whatever the step; the stationary density solves Aᵀg = 0 by inverse iteration. Both rest on calx's band LU, the matrix having bandwidth equal to the number of types.

An intervention is a map on densities: each state's mass carried to another position in its type, split between the two nodes about its destination, so that mass and mean are kept wherever the destination lies within the grid. Composed with flows in both orders, two interventions give the exact effect of order under the declared law; what the law implies, not that anything behaves so.

A value equation gives a drift. A stock with an inflow in each regime, a return r and a control c drawn from it, valued by CRRA utility discounted at ρ, has the value of Hamilton, Jacobi and Bellman, solved by Achdou, Han, Lasry, Lions and Moll's implicit upwind scheme (2022), with the state constraints at both ends; its policy's drift moves the density, and the density's forward equation is the value equation's generator transposed, exactly. Economics reads it as a household's saving (Huggett 1993); the domain that uses it names it so.

## 2. Content

- `generator.hpp`: `Grid{lower, upper, cells}`; `Generator{grid, types, A}`, `generator(grid, drift, diffusion, switching)`, `add_rate`, `row_defect`; `Flow{dt, F}` with `step` and `run`, `flow(G, dt)`; `stationary(G)`; `total`, `mean`, `below`, `type_mass`, `gini`; `Moved{density, beyond}`, `carry(G, g, destination)`, `transfer(G, g, amount, line)`.
- `value.hpp`: `Stock{rho, r, gamma, inflow, switching, grid}`; `Solution{value, control, drift, G, iterations, change, converged}`; `solve(stock, delta, most, tolerance)`.

*Caveat.* The scheme is accurate to first order in the step. Reflection at the grid's far end removes a diffusion rate of order D/h, which a heavy tail feels; the grid must reach where the density is negligible. The reference values of Achdou et al.'s own codes (MATLAB) were not generated; the checks use the closed form and their propositions.

## 3. Checks

`check_kinesis_densities` (23): a generator's rows summing to zero; mass kept over 300 implicit steps, and non-negative; under a linear drift the mean moving exactly as the implicit step says, the rate reflection removed counted (to 10⁻¹⁴); Bouchaud and Mézard's stationary inverse gamma law (α = 4, β = 3) reached at first order in the step; a transfer and a flow in both orders, their means differing by δ(1 − q^N) exactly; translations commuting, transfers below a line not; carrying keeping mass and adding exactly the transfer to the mean; the Gini coefficient of two equal masses; with one regime the control within 1% of the permanent-income line κ(a + z/r), nearer as the grid is refined; with two regimes and r < ρ the lower regime holding at the limit with a mass there (by Richardson's extrapolation in √h) and the higher growing there, the regimes' masses the switching chain's exactly, and the stationary density solving Aᵀg = 0.

## 4. Dependencies

`kinesis::kernel` (refusals); calx's `linalg` (band matrices and their LU). The checks use calx's `verification`.
