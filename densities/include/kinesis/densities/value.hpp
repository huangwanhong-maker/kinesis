// A value equation whose policy is the drift of a density (decision 0009).
//
// A stock a ≥ a_min receives an inflow z_k in its regime k, returns r on itself, and is drawn on at the rate c, the
// control; the regime switches at rates λ_kl; what is drawn is valued by u(c) = c^{1−γ}/(1−γ) (log c when γ = 1),
// discounted at ρ. (Economics reads it as a household's saving out of an income that changes, Huggett 1993; the domain
// that uses it names it so.) The value solves the Hamilton–Jacobi–Bellman equation
//   ρ v_k(a) = max_c { u(c) + v_k'(a)(z_k + r a − c) } + Σ_l λ_kl (v_l(a) − v_k(a)),
// whose maximiser is c = (u')⁻¹(v_k'(a)) = v_k'(a)^{−1/γ}, and the stock then drifts at s_k(a) = z_k + r a − c. Solved by
// Achdou et al.'s implicit upwind scheme (2022, Rev. Econ. Stud. 89, 45–86, and its numerical appendix): the derivative
// is taken forward where the drift it implies is positive, backward where negative, and otherwise the stock is held
// (c = z_k + r a); at the grid's ends the state constraints impose v'(a_min) = u'(z_k + r a_min) and
// v'(a_max) = u'(z_k + r a_max), so that the stock never leaves. Each iteration solves
//   ((1/Δ + ρ) I − A) v_{n+1} = u(c_n) + v_n/Δ,
// A the generator of the chain whose drift is the policy (generator.hpp, with no diffusion), until v stops changing.
// The stationary density of the units solves Aᵀ g = 0: the forward equation is the transpose of the value equation's
// generator, exactly.
//
// Known results it is held to (the checks): with one regime and the natural limit a_min = −z/r, the control is
// c = κ (a + z/r) with κ = (ρ − r(1 − γ))/γ (the permanent-income solution of the Euler equation ċ/c = (r − ρ)/γ); with
// two regimes and r < ρ the stock in the lower regime reaches a_min in finite time and stays, so its stationary density has
// a mass at a_min, while in the higher regime it grows there (Achdou et al. 2022, their Propositions 1 and 2, after Huggett 1993). The
// reference values of the authors' own codes were not generated here (they are MATLAB); the checks use the closed form
// and the propositions.
#pragma once

#include <kinesis/densities/generator.hpp>

#include <calx/linalg/band.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace kinesis::densities {

struct Stock {
  double rho = 0.05, r = 0.03, gamma = 2;
  std::vector<double> inflow;                   // z_k, in each regime
  std::vector<std::vector<double>> switching;   // λ_kl (one regime: empty)
  Grid grid;                                    // the stock a, from a_min
};

struct Solution {
  std::vector<std::vector<double>> value, control, drift;  // [regime][node]: v, c and s = z + r a − c
  Generator G;                                             // the chain moved by the policy's drift
  std::size_t iterations = 0;
  double change = 0;
  bool converged = false;
};

namespace detail {
inline double utility(double c, double gamma) { return gamma == 1 ? std::log(c) : std::pow(c, 1 - gamma) / (1 - gamma); }
inline double marginal(double c, double gamma) { return std::pow(c, -gamma); }
inline double inverse_marginal(double p, double gamma) { return std::pow(std::max(p, 1e-300), -1 / gamma); }
}  // namespace detail

inline Solution solve(const Stock& m, double delta = 1000, std::size_t most = 500, double tolerance = 1e-9) {
  require(m.grid);
  const std::size_t K = m.inflow.size(), n = m.grid.cells;
  if (K == 0) throw Refusal("an inflow for each regime");
  if (K > 1 && m.switching.size() != K) throw Refusal("the rates of switching between the regimes");
  if (!(m.rho > 0) || !(m.gamma > 0)) throw Refusal("a positive discount rate and curvature");
  const double h = m.grid.step();
  for (std::size_t k = 0; k < K; k++)
    if (!(m.inflow[k] + m.r * m.grid.lower > 0)) throw Refusal("the inflow must cover the return at the lowest stock, or the control there is not positive");
  Solution S;
  S.value.assign(K, std::vector<double>(n));
  for (std::size_t k = 0; k < K; k++)
    for (std::size_t i = 0; i < n; i++) S.value[k][i] = detail::utility(m.inflow[k] + m.r * m.grid.at(i), m.gamma) / m.rho;
  std::vector<std::vector<double>> drift(K, std::vector<double>(n)), zero(K, std::vector<double>(n, 0.0));
  auto policy = [&](const std::vector<std::vector<double>>& v) {
    for (std::size_t k = 0; k < K; k++)
      for (std::size_t i = 0; i < n; i++) {
        const double a = m.grid.at(i), hold = m.inflow[k] + m.r * a;
        const double vf = i + 1 < n ? (v[k][i + 1] - v[k][i]) / h : detail::marginal(hold, m.gamma);
        const double vb = i > 0 ? (v[k][i] - v[k][i - 1]) / h : detail::marginal(hold, m.gamma);
        const double sf = hold - detail::inverse_marginal(vf, m.gamma), sb = hold - detail::inverse_marginal(vb, m.gamma);
        double c, s;
        if (sf > 0 && i + 1 < n) c = detail::inverse_marginal(vf, m.gamma), s = sf;
        else if (sb < 0 && i > 0) c = detail::inverse_marginal(vb, m.gamma), s = sb;
        else c = hold, s = 0;
        S.control[k][i] = c, drift[k][i] = s;
      }
  };
  S.control.assign(K, std::vector<double>(n));
  for (S.iterations = 0; S.iterations < most; S.iterations++) {
    policy(S.value);
    const Generator G = generator(m.grid, drift, zero, m.switching);
    Band M = G.A;
    for (double& x : M.a) x = -x;
    for (std::size_t s = 0; s < G.states(); s++) M(s, s) += 1 / delta + m.rho;
    std::vector<double> rhs(G.states());
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t k = 0; k < K; k++) rhs[G.index(i, k)] = detail::utility(S.control[k][i], m.gamma) + S.value[k][i] / delta;
    const auto x = calx::linalg::band_solve(M, rhs);
    S.change = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t k = 0; k < K; k++) {
        S.change = std::max(S.change, std::fabs(x[G.index(i, k)] - S.value[k][i]));
        S.value[k][i] = x[G.index(i, k)];
      }
    if (S.change < tolerance) {
      S.converged = true;
      S.iterations++;
      break;
    }
  }
  policy(S.value);
  S.drift = drift;
  S.G = generator(m.grid, drift, zero, m.switching);
  return S;
}

}  // namespace kinesis::densities
