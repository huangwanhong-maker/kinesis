// Landscapes: a function on the state space whose valleys are where a system stays.
//
// The word covers several functions, which agree only in special cases (Zhou, Aliyu, Aurell and Huang 2012,
// J. R. Soc. Interface 9, 3539–3553). Each one here says which it is.
//
// The potential of a gradient system. For x' = −∇U the landscape is U itself. A flow of one dimension is always
// one, U(x) = −∫ f dx; in more dimensions a flow is a gradient exactly when its Jacobian is symmetric everywhere
// on a simply connected domain, which is tested at points, and U is then the line integral of −f from a base
// point. The minima of U are the stable equilibria and the saddles between them the barriers.
//
// The stationary landscape of a stochastic system, −ln p_ss(x). For one dimension it is exact: an Itô equation
// dX = a dt + b dW has the stationary density p ∝ b⁻² exp(∫ 2a/b² dx) when the integral is normalisable
// (Gardiner 2009, Stochastic Methods, 4th ed., §5.3.3). For additive noise b = σ this is p ∝ exp(−2U/σ²): the
// stationary landscape is the potential times 2/σ². In more dimensions it is estimated here from long paths
// (four by default), as minus the logarithm of a histogram, with a ceiling where no sample fell; it is not a
// Lyapunov function of the deterministic part in general, and is exactly twice the Freidlin–Wentzell
// quasi-potential over σ² only in the limit of small noise within a basin (Zhou et al. 2012).
//
// The time to cross a barrier. For one dimension with additive noise, Kramers's formula for the mean time from a
// well at x_a over a barrier at x_b (Kramers 1940; Gardiner 2009, §5.5.3, in the limit ΔU ≫ σ²/2):
//     τ = 2π / √(U''(x_a) |U''(x_b)|) · exp(2 (U(x_b) − U(x_a)) / σ²).
#pragma once

#include <kinesis/systems/dynamics.hpp>
#include <kinesis/systems/stochastic.hpp>

#include <calx/random/stream.hpp>

#include <cmath>
#include <numbers>
#include <vector>

namespace kinesis::analysis {

using systems::Dynamics;
using systems::State;
using systems::Stochastic;

struct Potential {
  std::vector<double> x, u;
  std::vector<double> minima, maxima;  // the places of its wells and of its barriers, within the grid
};

namespace detail {
// The extrema of a sampled function, each refined by the parabola through its three samples.
inline void extrema(Potential& p) {
  for (std::size_t i = 1; i + 1 < p.x.size(); i++) {
    const double a = p.u[i - 1], b = p.u[i], c = p.u[i + 1];
    if (!((b < a && b <= c) || (b > a && b >= c))) continue;
    const double h = p.x[i + 1] - p.x[i], curvature = a - 2 * b + c;
    const double place = curvature == 0 ? p.x[i] : p.x[i] - 0.5 * h * (c - a) / curvature;
    (b < a ? p.minima : p.maxima).push_back(place);
  }
}
}  // namespace detail

// U(x) = −∫ f dx of a flow of one dimension on [lower, upper], zero at the lower end; by Simpson's rule.
inline Potential potential(const Dynamics& d, double lower, double upper, std::size_t cells = 400) {
  if (d.dimension() != 1 || d.discrete) throw kernel::Refusal("the potential −∫ f dx is of a flow of one dimension; " + d.name + " has " + std::to_string(d.dimension()));
  Potential p;
  const double h = (upper - lower) / static_cast<double>(cells);
  State x(1), f(1);
  const auto field = [&](double at) {
    x[0] = at;
    d.field(0, x, f);
    return f[0];
  };
  double u = 0, previous = field(lower);
  p.x.push_back(lower), p.u.push_back(0);
  for (std::size_t i = 1; i <= cells; i++) {
    const double a = lower + h * static_cast<double>(i - 1), mid = field(a + 0.5 * h), end = field(a + h);
    u -= h / 6 * (previous + 4 * mid + end);
    previous = end;
    p.x.push_back(a + h), p.u.push_back(u);
  }
  detail::extrema(p);
  return p;
}

// Whether a flow is a gradient: the largest asymmetry of its Jacobian at `points` points of a box, against the
// largest entry there. Zero (to rounding) at every point is the condition on a simply connected domain.
inline double asymmetry(const Dynamics& d, const std::vector<double>& lower, const std::vector<double>& upper, std::size_t points = 200, std::uint64_t seed = 1) {
  const std::size_t n = d.dimension();
  calx::random::Stream r(seed, 0x67726164u, 0, 0);  // 'grad'
  calx::linalg::Matrix J(n);
  double worst = 0;
  for (std::size_t k = 0; k < points; k++) {
    State x(n);
    for (std::size_t i = 0; i < n; i++) x[i] = r.uniform(lower[i], upper[i]);
    d.jacobian(0, x, J);
    double scale = 1e-300, gap = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) scale = std::max(scale, std::fabs(J(i, j))), gap = std::max(gap, std::fabs(J(i, j) - J(j, i)));
    worst = std::max(worst, gap / scale);
  }
  return worst;
}

// U(x) − U(base) = −∫ f · dl along the straight line from `base`, for a gradient flow (Simpson's rule).
inline double potential_between(const Dynamics& d, const State& base, const State& x, std::size_t cells = 200) {
  const std::size_t n = d.dimension();
  State at(n), f(n);
  const auto along = [&](double s) {
    for (std::size_t i = 0; i < n; i++) at[i] = base[i] + s * (x[i] - base[i]);
    d.field(0, at, f);
    double dot = 0;
    for (std::size_t i = 0; i < n; i++) dot += f[i] * (x[i] - base[i]);
    return dot;
  };
  double u = 0;
  const double h = 1.0 / static_cast<double>(cells);
  for (std::size_t i = 0; i < cells; i++) {
    const double a = h * static_cast<double>(i);
    u -= h / 6 * (along(a) + 4 * along(a + 0.5 * h) + along(a + h));
  }
  return u;
}

// −ln p_ss of a stochastic system of one dimension on a grid, exact up to the quadrature and to a constant: zero
// at its least value.
inline Potential stationary_landscape(const Stochastic& d, double lower, double upper, std::size_t cells = 400) {
  if (d.dimension() != 1 || d.discrete || d.sources.size() != 1) throw kernel::Refusal("the exact stationary landscape is of a stochastic differential of one dimension with one source");
  Potential p;
  const double h = (upper - lower) / static_cast<double>(cells);
  State x(1), a(1);
  std::vector<double> B;
  const auto ratio = [&](double at) {  // 2a / b²
    x[0] = at;
    d.drift(0, x, a);
    d.diffusion(0, x, B);
    return 2 * a[0] / (B[0] * B[0]);
  };
  const auto log_b2 = [&](double at) {
    x[0] = at;
    d.diffusion(0, x, B);
    return std::log(B[0] * B[0]);
  };
  double integral = 0, previous = ratio(lower);
  p.x.push_back(lower), p.u.push_back(log_b2(lower));
  for (std::size_t i = 1; i <= cells; i++) {
    const double s = lower + h * static_cast<double>(i - 1), mid = ratio(s + 0.5 * h), end = ratio(s + h);
    integral += h / 6 * (previous + 4 * mid + end);
    previous = end;
    p.x.push_back(s + h), p.u.push_back(log_b2(s + h) - integral);
  }
  const double least = *std::min_element(p.u.begin(), p.u.end());
  for (double& v : p.u) v -= least;
  detail::extrema(p);
  return p;
}

// Kramers's mean time over a barrier, from the well at `well` to the barrier at `barrier`, for dx = −U' dt + σ dW.
inline double kramers_time(const Dynamics& gradient, double well, double barrier, double sigma) {
  if (gradient.dimension() != 1) throw kernel::Refusal("Kramers's formula is for one dimension");
  calx::linalg::Matrix J(1);
  gradient.jacobian(0, {well}, J);
  const double curvature_well = -J(0, 0);       // U'' = −f'
  gradient.jacobian(0, {barrier}, J);
  const double curvature_barrier = -J(0, 0);
  if (!(curvature_well > 0) || !(curvature_barrier < 0)) throw kernel::Refusal("Kramers's formula needs a well (U'' > 0) and a barrier (U'' < 0) at the places given");
  const double height = potential_between(gradient, {well}, {barrier});
  return 2 * std::numbers::pi / std::sqrt(curvature_well * -curvature_barrier) * std::exp(2 * height / (sigma * sigma));
}

// The stationary landscape −ln p of two components of a stochastic system, from a long path.
struct DensityLandscape {
  std::vector<double> xs, ys;      // the centres of the cells
  std::vector<double> u;           // u[iy * nx + ix], zero at its least; `ceiling` where no sample fell
  std::vector<double> count;
  double ceiling = 0;
  std::size_t samples = 0;
};

struct DensityOptions {
  std::size_t component_x = 0, component_y = 1;
  double x0 = -2, x1 = 2, y0 = -2, y1 = 2;
  std::size_t nx = 60, ny = 60;
  double time = 2000, step = 0.01, transient = 50;
  std::uint64_t seed = 1;
  std::uint32_t paths = 4;
};

inline DensityLandscape density_landscape(const Stochastic& d, const State& start, const DensityOptions& o = {}) {
  DensityLandscape L;
  for (std::size_t i = 0; i < o.nx; i++) L.xs.push_back(o.x0 + (o.x1 - o.x0) * (static_cast<double>(i) + 0.5) / static_cast<double>(o.nx));
  for (std::size_t i = 0; i < o.ny; i++) L.ys.push_back(o.y0 + (o.y1 - o.y0) * (static_cast<double>(i) + 0.5) / static_cast<double>(o.ny));
  L.count.assign(o.nx * o.ny, 0.0);
  const bool one = d.dimension() == 1;
  for (std::uint32_t path = 0; path < o.paths; path++) {
    systems::StochasticOptions so;
    so.seed = o.seed, so.path = path, so.step = o.step, so.t1 = o.transient + o.time / static_cast<double>(o.paths);
    const auto tr = systems::simulate(d, start, so);
    for (std::size_t k = 0; k < tr.x.size(); k++) {
      if (tr.t[k] < o.transient) continue;
      const double x = tr.x[k][o.component_x], y = one ? 0.5 * (o.y0 + o.y1) : tr.x[k][o.component_y];
      const double fx = (x - o.x0) / (o.x1 - o.x0), fy = (y - o.y0) / (o.y1 - o.y0);
      if (fx < 0 || fx >= 1 || fy < 0 || fy >= 1) continue;
      L.count[static_cast<std::size_t>(fy * static_cast<double>(o.ny)) * o.nx + static_cast<std::size_t>(fx * static_cast<double>(o.nx))] += 1;
      L.samples++;
    }
  }
  double most = 0;
  for (double c : L.count) most = std::max(most, c);
  if (most == 0) throw kernel::Refusal("no sample of the path fell in the box of the landscape");
  L.ceiling = std::log(most);  // where one sample would stand, and beyond
  L.u.resize(L.count.size());
  for (std::size_t i = 0; i < L.count.size(); i++) L.u[i] = L.count[i] > 0 ? std::log(most / L.count[i]) : L.ceiling;
  return L;
}

}  // namespace kinesis::analysis
