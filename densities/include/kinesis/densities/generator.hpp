// A density on a grid of cells moved by a Markov generator: drift and diffusion that vary from cell to cell, in types
// between which it switches (decision 0009).
//
// A state is a cell i of a uniform grid x_0 < … < x_{n−1} (step h) and a type k of K, numbered cell-major, s = iK + k.
// A diffusion dX = μ_k(X) dt + √(2 D_k(X)) dB in type k, switching from type k to l at rate λ_kl, is approximated by the
// continuous-time Markov chain whose rates out of (i, k) are
//   to (i + 1, k):  max(μ, 0)/h + D/h²        to (i − 1, k):  max(−μ, 0)/h + D/h²        to (i, l):  λ_kl,
// the rates that would leave the grid removed (reflection). This is the upwind scheme: its generator A has non-negative
// rates off the diagonal and rows summing to zero, the density's forward (Kolmogorov) equation is dg/dt = Aᵀ g, and its
// mass Σ g is kept exactly, since every column of Aᵀ sums to zero (Kushner and Dupuis 2001, Numerical Methods for
// Stochastic Control Problems in Continuous Time, ch. 5, the locally consistent chain; Achdou, Han, Lasry, Lions and
// Moll 2022, "Income and wealth distribution in macroeconomics: a continuous-time approach", Rev. Econ. Stud. 89, 45–86,
// whose forward equation is the transpose of their value equation's generator). The mean moves by Σ g μ exactly where no
// rate was removed, since the diffusion's two rates cancel in it; the scheme is accurate to first order in h.
//
// A step of time dt is implicit, (I − dt Aᵀ) g' = g, an M-matrix: it keeps mass and positivity for any dt. A
// stationary density solves Aᵀ g = 0 with Σ g = 1, found by inverse iteration with a small shift. Both use calx's band
// LU, the matrix having bandwidth K.
//
// An intervention is a map on densities: each state's mass carried to another place in the same type (a transfer of
// wealth to those below a line, say), split between the two nodes about its destination in proportion to nearness, so
// that mass and mean are kept exactly wherever the destination lies within the grid (what is carried beyond is put on
// the last node and reported). Two interventions with a flow between them, composed in both orders, give the exact effect
// of order under the declared law: what the law implies, not that anything behaves so.
#pragma once

#include <kinesis/kernel/answer.hpp>

#include <calx/linalg/band.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace kinesis::densities {

using calx::linalg::Band;
using kernel::Refusal;

struct Grid {
  double lower = 0, upper = 1;
  std::size_t cells = 101;
  double step() const { return (upper - lower) / static_cast<double>(cells - 1); }
  double at(std::size_t i) const { return lower + step() * static_cast<double>(i); }
};

inline void require(const Grid& g) {
  if (g.cells < 3) throw Refusal("a grid needs at least three cells");
  if (!(g.upper > g.lower)) throw Refusal("a grid's upper end must lie above its lower");
}

struct Generator {
  Grid grid;
  std::size_t types = 1;
  Band A;  // row s: the rates of leaving s for each state, the diagonal minus their sum
  std::size_t states() const { return grid.cells * types; }
  std::size_t index(std::size_t i, std::size_t k) const { return i * types + k; }
};

// The rates of one state, added into A (the diagonal kept as minus their sum).
inline void add_rate(Band& A, std::size_t from, std::size_t to, double rate) {
  if (!(rate >= 0) || !std::isfinite(rate)) throw Refusal("a rate must be finite and non-negative");
  if (rate == 0 || from == to) return;
  A(from, to) += rate;
  A(from, from) -= rate;
}

// drift[k][i] and diffusion[k][i] (D ≥ 0) at each node of each type; switching[k][l] the rate from type k to type l
// (the diagonal ignored; empty with one type).
inline Generator generator(const Grid& g, const std::vector<std::vector<double>>& drift, const std::vector<std::vector<double>>& diffusion,
                           const std::vector<std::vector<double>>& switching = {}) {
  require(g);
  const std::size_t K = drift.size(), n = g.cells;
  if (K == 0 || diffusion.size() != K) throw Refusal("a drift and a diffusion for each type");
  if (K > 1 && switching.size() != K) throw Refusal("the rates of switching between the types, a row for each");
  Generator G;
  G.grid = g, G.types = K;
  G.A = Band(n * K, K, K);
  const double h = g.step();
  for (std::size_t k = 0; k < K; k++) {
    if (drift[k].size() != n || diffusion[k].size() != n) throw Refusal("a drift and a diffusion at every node of type " + std::to_string(k));
    for (std::size_t i = 0; i < n; i++) {
      const double mu = drift[k][i], D = diffusion[k][i];
      if (!(D >= 0) || !std::isfinite(D) || !std::isfinite(mu)) throw Refusal("a diffusion must be finite and non-negative, a drift finite");
      const std::size_t s = G.index(i, k);
      if (i + 1 < n) add_rate(G.A, s, G.index(i + 1, k), std::max(mu, 0.0) / h + D / (h * h));
      if (i > 0) add_rate(G.A, s, G.index(i - 1, k), std::max(-mu, 0.0) / h + D / (h * h));
      if (K > 1) {
        if (switching[k].size() != K) throw Refusal("a rate of switching to each type");
        for (std::size_t l = 0; l < K; l++)
          if (l != k) add_rate(G.A, s, G.index(i, l), switching[k][l]);
      }
    }
  }
  return G;
}

// The largest |row sum| of A: zero up to rounding for a generator.
inline double row_defect(const Generator& G) {
  double worst = 0;
  const std::size_t N = G.states();
  for (std::size_t s = 0; s < N; s++) {
    double sum = 0;
    for (std::size_t t = (s >= G.A.kl ? s - G.A.kl : 0); t <= std::min(N - 1, s + G.A.ku); t++) sum += G.A(s, t);
    worst = std::max(worst, std::fabs(sum));
  }
  return worst;
}

// The implicit flow of time dt: (I − dt Aᵀ) factorised once.
struct Flow {
  double dt = 0;
  calx::linalg::BandLu F;
  std::vector<double> step(const std::vector<double>& g) const { return calx::linalg::band_solve(F, g); }
  std::vector<double> run(std::vector<double> g, std::size_t steps) const {
    for (std::size_t k = 0; k < steps; k++) g = step(g);
    return g;
  }
};

inline Flow flow(const Generator& G, double dt) {
  if (!(dt > 0)) throw Refusal("a step of time must be positive");
  Band M = calx::linalg::transpose(G.A);
  for (double& v : M.a) v *= -dt;
  for (std::size_t s = 0; s < G.states(); s++) M(s, s) += 1;
  return {dt, calx::linalg::band_lu(std::move(M))};
}

// The stationary density (mass per state, summing to one): inverse iteration on Aᵀ with a small shift.
inline std::vector<double> stationary(const Generator& G, double tolerance = 1e-14, std::size_t most = 200) {
  const std::size_t N = G.states();
  double scale = 0;
  for (double v : G.A.a) scale = std::max(scale, std::fabs(v));
  if (!(scale > 0)) throw Refusal("a generator with no rate has no single stationary density");
  Band M = calx::linalg::transpose(G.A);
  for (double& v : M.a) v = -v;
  const double shift = 1e-10 * scale;
  for (std::size_t s = 0; s < N; s++) M(s, s) += shift;
  const auto F = calx::linalg::band_lu(std::move(M));
  std::vector<double> g(N, 1.0 / static_cast<double>(N));
  for (std::size_t it = 0; it < most; it++) {
    auto y = calx::linalg::band_solve(F, g);
    double sum = 0;
    for (double& v : y) v = std::max(v, 0.0), sum += v;
    if (!(sum > 0)) throw Refusal("inverse iteration lost the density");
    double change = 0;
    for (std::size_t s = 0; s < N; s++) {
      y[s] /= sum;
      change += std::fabs(y[s] - g[s]);
    }
    g.swap(y);
    if (change < tolerance) break;
  }
  return g;
}

// ---- readings of a density ----------------------------------------------------------------------------------------------

inline double total(const std::vector<double>& g) {
  double s = 0;
  for (double v : g) s += v;
  return s;
}
inline double mean(const Generator& G, const std::vector<double>& g) {
  double s = 0, m = 0;
  for (std::size_t i = 0; i < G.grid.cells; i++)
    for (std::size_t k = 0; k < G.types; k++) s += G.grid.at(i) * g[G.index(i, k)], m += g[G.index(i, k)];
  return m > 0 ? s / m : 0;
}
// The mass at nodes strictly below `line`.
inline double below(const Generator& G, const std::vector<double>& g, double line) {
  double m = 0;
  for (std::size_t i = 0; i < G.grid.cells && G.grid.at(i) < line; i++)
    for (std::size_t k = 0; k < G.types; k++) m += g[G.index(i, k)];
  return m;
}
inline double type_mass(const Generator& G, const std::vector<double>& g, std::size_t k) {
  double m = 0;
  for (std::size_t i = 0; i < G.grid.cells; i++) m += g[G.index(i, k)];
  return m;
}
// The Gini coefficient of the positions held with these masses (all types together; positions must be non-negative):
// G = 1 − Σ_i p_i (S_{i−1} + S_i)/S_n, S_i the cumulative Σ_{j ≤ i} p_j x_j, over the nodes in increasing order (the
// trapezoid under the Lorenz curve).
inline double gini(const Generator& G, const std::vector<double>& g) {
  double S = 0, prev = 0, area = 0, m = 0;
  std::vector<double> p(G.grid.cells, 0);
  for (std::size_t i = 0; i < G.grid.cells; i++)
    for (std::size_t k = 0; k < G.types; k++) p[i] += g[G.index(i, k)];
  for (double v : p) m += v;
  double top = 0;
  for (std::size_t i = 0; i < G.grid.cells; i++) {
    if (G.grid.at(i) < 0 && p[i] > 0) throw Refusal("a Gini coefficient of negative positions");
    top += p[i] * G.grid.at(i);
  }
  if (!(top > 0) || !(m > 0)) return 0;
  for (std::size_t i = 0; i < G.grid.cells; i++) {
    S = prev + p[i] * G.grid.at(i);
    area += p[i] / m * (prev + S) / top;
    prev = S;
  }
  return 1 - area;
}

// ---- interventions ------------------------------------------------------------------------------------------------------

struct Moved {
  std::vector<double> density;
  double beyond = 0;  // mass whose destination lay outside the grid, put on its nearest end
};

// Each state's mass carried to destination(x, k), split linearly between the two nodes about it.
inline Moved carry(const Generator& G, const std::vector<double>& g, const std::function<double(double, std::size_t)>& destination) {
  if (g.size() != G.states()) throw Refusal("a density with a mass for every state");
  Moved out;
  out.density.assign(g.size(), 0);
  const double h = G.grid.step();
  const std::size_t n = G.grid.cells;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k < G.types; k++) {
      const double m = g[G.index(i, k)];
      if (m == 0) continue;
      double y = destination(G.grid.at(i), k);
      if (!std::isfinite(y)) throw Refusal("a destination must be finite");
      if (y < G.grid.lower) out.beyond += m, y = G.grid.lower;
      if (y > G.grid.upper) out.beyond += m, y = G.grid.upper;
      const double u = (y - G.grid.lower) / h;
      std::size_t j = static_cast<std::size_t>(std::floor(u));
      if (j >= n - 1) j = n - 2;
      const double f = u - static_cast<double>(j);
      out.density[G.index(j, k)] += m * (1 - f);
      out.density[G.index(j + 1, k)] += m * f;
    }
  return out;
}

// A transfer of `amount` to every unit of mass at positions below `line`.
inline Moved transfer(const Generator& G, const std::vector<double>& g, double amount, double line) {
  return carry(G, g, [amount, line](double x, std::size_t) { return x < line ? x + amount : x; });
}

}  // namespace kinesis::densities
