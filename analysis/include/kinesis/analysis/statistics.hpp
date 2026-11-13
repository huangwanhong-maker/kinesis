// What a system does on average: the means of its observations along a run, with their errors, and the law by
// which they scale.
//
// A system that neither settles nor repeats is described by its statistics: the mean energy at each scale of a
// turbulent flow, the mean abundance of a species that fluctuates. A declaration names what is observed
// (`difference { observe … }`), and averages() follows those observations along one long path.
//
// The mean and its error. The path is sampled at equal times after a transient, the samples are cut into blocks,
// and the standard error of the mean is that of the block means (Flyvbjerg and Petersen 1989, J. Chem. Phys. 91,
// 461–466): the standard deviation of the block means over the square root of their number. *Caveat.* It is right
// when a block is long against the time the observation remembers; where it is not, the error is understated,
// and no test of that is made here.
//
// The second moment. With each mean ⟨o⟩ comes ⟨o²⟩/⟨o⟩². It is 1 for an observation that does not fluctuate.
// For o = |u|², the squared size of a complex amplitude, it is the flatness ⟨|u|⁴⟩/⟨|u|²⟩² of u: 2 when u is a
// complex Gaussian variable, and more when u is intermittent, quiet for long and then large.
//
// Scaling. Means that follow a power law of a scale, ⟨o_n⟩ ∝ k_n^{−ζ}, lie on a line in logarithms. scaling()
// fits that line by least squares and gives the exponent with the standard error of the slope from the residuals
// (Press et al. 2007, Numerical Recipes, 3rd ed., §15.2). *Caveat.* A line fitted over a few octaves measures an
// exponent over those octaves: it does not show that a law holds beyond them, and the error is of the fit, not
// of the means.
#pragma once

#include <kinesis/systems/dynamics.hpp>
#include <kinesis/systems/simulate.hpp>

#include <cmath>
#include <span>
#include <string>
#include <vector>

namespace kinesis::analysis {

using systems::Dynamics;
using systems::State;

// The mean of a series and the standard error of that mean from `blocks` blocks of equal length; samples beyond
// the last whole block are left out of the error and kept in the mean.
struct BlockMean {
  double mean = 0, error = 0;
  std::size_t samples = 0, blocks = 0;
};
inline BlockMean block_mean(std::span<const double> series, std::size_t blocks = 10) {
  BlockMean out;
  out.samples = series.size();
  if (series.empty()) return out;
  for (double v : series) out.mean += v / static_cast<double>(series.size());
  const std::size_t length = series.size() / std::max<std::size_t>(blocks, 1);
  if (blocks < 2 || length == 0) return out;
  std::vector<double> means(blocks, 0.0);
  for (std::size_t b = 0; b < blocks; b++)
    for (std::size_t k = 0; k < length; k++) means[b] += series[b * length + k] / static_cast<double>(length);
  double centre = 0, spread = 0;
  for (double m : means) centre += m / static_cast<double>(blocks);
  for (double m : means) spread += (m - centre) * (m - centre) / static_cast<double>(blocks - 1);
  out.error = std::sqrt(spread / static_cast<double>(blocks)), out.blocks = blocks;
  return out;
}

struct Statistic {
  std::string name;
  double mean = 0, error = 0;
  double second = 0;   // ⟨o²⟩ / ⟨o⟩², or 0 where the mean is zero
};

struct StatisticsOptions {
  double step = 0.01;        // a flow's step
  double transient = 100;    // time (iterations) before the average begins
  double time = 1000;        // time (iterations) averaged over
  std::size_t every = 10;    // steps between samples
  std::size_t blocks = 10;
};

struct Statistics {
  std::vector<Statistic> of;   // one for each observation the system declares, in their order
  double time = 0;
  std::size_t samples = 0, blocks = 0;
  State final_state;
};

// The averages of a system's observations along the path from x0.
inline Statistics averages(const Dynamics& d, State x, const StatisticsOptions& o = {}) {
  if (x.size() != d.dimension()) throw kernel::Refusal(d.name + " has " + std::to_string(d.dimension()) + " components; " + std::to_string(x.size()) + " were given");
  if (d.observables.empty()) throw kernel::Refusal(d.name + " declares nothing to observe: name what is averaged in `difference { observe … }`");
  if (!d.discrete && !(o.step > 0)) throw kernel::Refusal("the step of a flow is positive");
  systems::Options run;
  run.step = o.step;
  double t = 0;
  if (d.discrete) {
    for (double k = 0; k < o.transient; k += 1) systems::iterate(d, t, x);
  } else {
    for (double s = 0; s < o.transient - 0.5 * o.step; s += o.step) systems::advance(d, t, x, o.step, run);
  }
  const std::size_t every = std::max<std::size_t>(o.every, 1);
  const std::size_t steps = static_cast<std::size_t>(d.discrete ? std::max(0.0, o.time) : std::ceil(o.time / o.step - 1e-9));
  const std::size_t m = d.observables.size();
  std::vector<std::vector<double>> series(m);
  for (auto& s : series) s.reserve(steps / every + 1);
  std::vector<double> squares(m, 0.0);
  for (std::size_t k = 1; k <= steps; k++) {
    if (d.discrete) systems::iterate(d, t, x);
    else systems::advance(d, t, x, o.step, run);
    if (k % every) continue;
    for (std::size_t i = 0; i < m; i++) {
      const double v = d.observe(i, t, x);
      if (!std::isfinite(v)) throw kernel::Refusal(d.name + " left the finite numbers at t = " + terms::number_text(t));
      series[i].push_back(v), squares[i] += v * v;
    }
  }
  Statistics out;
  out.time = d.discrete ? static_cast<double>(steps) : static_cast<double>(steps) * o.step;
  out.samples = series.front().size(), out.final_state = x;
  if (out.samples == 0) throw kernel::Refusal("no sample was taken: the time averaged over is shorter than the time between samples");
  for (std::size_t i = 0; i < m; i++) {
    const BlockMean b = block_mean(series[i], o.blocks);
    out.blocks = b.blocks;
    const double second = squares[i] / static_cast<double>(out.samples);
    out.of.push_back({d.observables[i], b.mean, b.error, b.mean != 0 ? second / (b.mean * b.mean) : 0.0});
  }
  return out;
}

// The power law y ∝ x^exponent through points with positive coordinates, by least squares on their logarithms.
struct Scaling {
  double exponent = 0, error = 0;   // the slope, and its standard error from the residuals
  double prefactor = 0;             // y at x = 1
  std::size_t points = 0;
};
inline Scaling scaling(std::span<const double> x, std::span<const double> y) {
  if (x.size() != y.size() || x.size() < 3) throw kernel::Refusal("a power law is fitted through three points or more, as many abscissas as ordinates");
  const double n = static_cast<double>(x.size());
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  std::vector<double> lx(x.size()), ly(x.size());
  for (std::size_t i = 0; i < x.size(); i++) {
    if (!(x[i] > 0) || !(y[i] > 0)) throw kernel::Refusal("a power law is fitted through positive values; one is not");
    lx[i] = std::log(x[i]), ly[i] = std::log(y[i]);
    sx += lx[i], sy += ly[i], sxx += lx[i] * lx[i], sxy += lx[i] * ly[i];
  }
  const double spread = sxx - sx * sx / n;
  if (!(spread > 0)) throw kernel::Refusal("a power law needs abscissas that differ");
  Scaling s;
  s.points = x.size();
  s.exponent = (sxy - sx * sy / n) / spread;
  const double intercept = (sy - s.exponent * sx) / n;
  s.prefactor = std::exp(intercept);
  double residual = 0;
  for (std::size_t i = 0; i < x.size(); i++) {
    const double r = ly[i] - intercept - s.exponent * lx[i];
    residual += r * r;
  }
  s.error = std::sqrt(residual / (n - 2) / spread);
  return s;
}

}  // namespace kinesis::analysis
