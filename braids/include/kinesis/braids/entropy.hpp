// The topological entropy of a braid, between a bound from below and an estimate.
//
// A braid on n strands is a mapping class of the disc with n punctures, and its topological entropy is the least
// entropy of the maps in the class: the logarithm of the dilatation when the class is pseudo-Anosov, and zero when
// it is of finite order (Fathi, Laudenbach and Poénaru 1979, Travaux de Thurston sur les surfaces; Boyland 1994,
// "Topological methods in surface dynamics", Topology Appl. 58, 223–298). A flow that moves n points along a
// braid has at least that entropy, which is why braids of trajectories measure mixing (Boyland, Aref and Stremler
// 2000, J. Fluid Mech. 403, 277–304; Thiffeault 2005, Phys. Rev. Lett. 94, 084502).
//
// The estimate. The braid is applied again and again to the curve system of the generating loops, in Dynnikov's
// coordinates (morphe/algebra/braid.hpp), and the growth of their norm per application is followed until two
// successive rates agree to a tolerance (Moussafir 2006; Thiffeault 2010, §III). The integers are exact; the
// estimate is a limit, and says whether it was reached. A rate that settles below the tolerance is reported as
// zero growth, which is what a braid of finite order or a reducible one with no growing part gives.
//
// The bound. For |t| = 1 the spectral radius of the Burau matrix bounds the entropy from below,
// h ≥ log sup_{|t|=1} ρ(ψ_t(β)) (Fried 1986; Kolev 1989, C. R. Acad. Sci. Paris 309, 835–838; Band and Boyland
// 2007, Algebr. Geom. Topol. 7, 1345–1378). It is taken over a grid of t, so it is itself a bound from below on
// the supremum. For three strands it is attained at t = −1.
//
// The exponent of a finite time. For a braid read over a time T, (1/T) log(|β·ℓ| / |ℓ|) with ℓ the generating
// loops (Budišić and Thiffeault 2015, Chaos 25, 087407): here with the norm of Dynnikov's coordinates in the place
// of their count of intersections, which has the same exponential rate and not the same value at short times.
#pragma once

#include <kinesis/kernel/answer.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>
#include <morphe/algebra/braid.hpp>

#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

namespace kinesis::braids {

struct Entropy {
  double estimate = 0;       // per application of the braid
  double lower = 0;          // the Burau bound
  std::size_t applications = 0;
  bool settled = false;      // two successive rates agreed to the tolerance
};

// log sup ρ(ψ_t(β)) over `grid` values of t on the unit circle.
inline double burau_bound(const morphe::algebra::Braid& b, std::size_t grid = 180) {
  const std::size_t n = b.strands;
  if (n < 3) return 0;
  const auto psi = morphe::algebra::burau(b);
  double best = 1;
  for (std::size_t g = 1; g < grid; g++) {
    const std::complex<double> t = std::polar(1.0, 2 * std::numbers::pi * static_cast<double>(g) / static_cast<double>(grid));
    // A complex matrix M as the real matrix [[Re M, −Im M], [Im M, Re M]], whose eigenvalues are M's and their conjugates.
    calx::linalg::Matrix real(2 * n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) {
        const std::complex<double> z = psi[i][j].at(t);
        real(i, j) = z.real(), real(i, j + n) = -z.imag();
        real(i + n, j) = z.imag(), real(i + n, j + n) = z.real();
      }
    best = std::max(best, calx::linalg::spectral_radius(real));
  }
  return std::log(best);
}

inline Entropy entropy(const morphe::algebra::Braid& b, double tolerance = 1e-10, std::size_t most = 400) {
  Entropy e;
  e.lower = burau_bound(b);
  if (b.strands < 3 || b.word.empty()) {
    e.settled = true;
    return e;
  }
  morphe::algebra::Dynnikov u = morphe::algebra::generating_loops(b.strands);
  std::vector<double> logs{std::log(u.norm().to_double())};
  double rate = 0, last = -1;
  for (std::size_t k = 1; k <= most; k++) {
    u = morphe::algebra::act(b, u);
    const double size = u.norm().to_double();
    if (!std::isfinite(size) || size <= 0) break;
    logs.push_back(std::log(size));
    rate = logs[k] - logs[k - 1], e.applications = k;
    if (k > 3 && std::fabs(rate - last) < tolerance) {
      e.settled = true;
      break;
    }
    last = rate;
  }
  if (!e.settled && e.applications >= 8) {
    // Not settled: the mean rate over the second half of the applications. A norm that grows as a power of k
    // (a reducible braid with twists, no growing part) has a mean rate below 3 log(k) / k, and is given as zero.
    const std::size_t k = e.applications, half = k / 2;
    rate = (logs[k] - logs[half]) / static_cast<double>(k - half);
    if (rate < 3 * std::log(static_cast<double>(k)) / static_cast<double>(k)) rate = 0;
  }
  e.estimate = std::max(0.0, rate);
  return e;
}

// The braiding exponent of a braid read over a time T.
inline double braiding_exponent(const morphe::algebra::Braid& b, double duration) {
  if (b.strands < 3 || !(duration > 0)) return 0;
  const morphe::algebra::Dynnikov u = morphe::algebra::generating_loops(b.strands);
  return std::log(morphe::algebra::act(b, u).norm().to_double() / u.norm().to_double()) / duration;
}

}  // namespace kinesis::braids
