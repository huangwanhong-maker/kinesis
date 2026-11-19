// A drifting restoring rate against a constant one, by likelihood, with the test's error rates from the models themselves.
//
// Boettiger and Hastings (2012, J. R. Soc. Interface 9, 2527–2539; read in full by the survey, Part I §4) set the question
// of an early warning as a choice between two models of the record, a stationary Ornstein–Uhlenbeck process and one whose
// stability is being lost, compared by their likelihood ratio, and they measure the test's reliability by simulating the
// ratio under each fitted model (a parametric bootstrap), so that the warning comes with its false-alarm and detection
// rates and its ROC curve. A record sampled at any times has exact OU transitions,
//
//     x_{i+1} | x_i ~ N(θ + (x_i − θ) e^{−κΔ_i},  σ² (1 − e^{−2κΔ_i}) / (2κ)),   Δ_i = t_{i+1} − t_i,
//
// with the random walk's variance σ²Δ_i as κ → 0. The alternatives here keep θ and σ and let κ move in time:
//
//   linear   κ(t) = κ₀ + κ₁ (t − t₀), as before a transcritical or pitchfork bifurcation under linear forcing;
//   fold     κ(t)² = κ₀² + c (t − t₀), as before a fold (κ = 2√(μ_c − μ) with μ linear in time),
//
// κ taken constant over each interval at its value at the interval's start, which is exact in the limit of slow change and
// otherwise an approximation, said so. Their own alternative, the linearised saddle-node, also moves the noise and the mean
// with the parameter; these keep them fixed, so that a fall of κ is not confused with a change of noise.
//
// The time at which the rate would reach zero is reported only when the alternative fits better than the null at the
// stated level and the bootstrap of the alternative gives a bounded interval: tipping times extrapolated from records are
// unidentifiable more often than not (Ben-Yami et al. 2024, read in full by the survey, Part I §4: under reasonable
// alternatives the estimated times of the Atlantic overturning's collapse range over centuries, with up to a quarter
// infinite).
#pragma once

#include <kinesis/transitions/space.hpp>

#include <calx/optimisation/nelder_mead.hpp>
#include <calx/random/stream.hpp>
#include <tyche/descriptive/summary.hpp>

#include <numbers>
#include <span>

namespace kinesis::transitions {

enum class Drift { linear, fold };

inline const char* name(Drift d) { return d == Drift::linear ? "a restoring rate falling linearly in time" : "a restoring rate whose square falls linearly in time (a fold)"; }

struct OuFit {
  double theta = 0, kappa0 = 0, sigma = 0, slope = 0;  // slope: κ₁ (linear) or c (fold); 0 under the null
  double log_likelihood = 0;
  bool converged = false;
};

struct LikelihoodTest {
  Drift drift = Drift::linear;
  OuFit null, alternative;
  double deviance = 0;              // 2(ℓ₁ − ℓ₀)
  double p = 1;                     // the share of the null's replicates with a deviance at least as large
  double false_alarm = 0, detection = 0, auc = 0.5;
  std::vector<double> null_deviances, alternative_deviances;
  bool tipping_reported = false;
  double tipping = std::numeric_limits<double>::infinity(), tipping_lower = 0, tipping_upper = std::numeric_limits<double>::infinity();
  std::string tipping_sentence;
  std::vector<std::string> assumptions;
};

namespace detail {

inline double kappa_at(Drift d, const OuFit& f, double t, double t0) {
  if (d == Drift::linear) return f.kappa0 + f.slope * (t - t0);
  const double k2 = f.kappa0 * f.kappa0 + f.slope * (t - t0);
  return k2 > 0 ? std::sqrt(k2) : 0.0;
}

inline void transition(double x, double kappa, double theta, double sigma, double dt, double& mean, double& variance) {
  if (std::abs(kappa * dt) < 1e-8) {
    mean = x;
    variance = sigma * sigma * dt;
    return;
  }
  const double e = std::exp(-kappa * dt);
  mean = theta + (x - theta) * e;
  variance = sigma * sigma * (1 - e * e) / (2 * kappa);
}

inline double log_likelihood(const std::vector<double>& t, const std::vector<double>& x, Drift d, const OuFit& f, bool drifting) {
  double ll = 0;
  for (std::size_t i = 0; i + 1 < x.size(); i++) {
    const double k = drifting ? kappa_at(d, f, t[i], t[0]) : f.kappa0;
    double m, v;
    transition(x[i], k, f.theta, f.sigma, t[i + 1] - t[i], m, v);
    if (!(v > 0) || !std::isfinite(v)) return -1e300;
    ll += -0.5 * (std::log(2 * std::numbers::pi * v) + (x[i + 1] - m) * (x[i + 1] - m) / v);
  }
  return ll;
}

inline OuFit fit(const std::vector<double>& t, const std::vector<double>& x, Drift d, bool drifting, const OuFit* start) {
  const double duration = t.back() - t.front();
  OuFit init;
  if (start) init = *start;
  else {
    init.theta = tyche::descriptive::mean(x);
    double sxx = 0, sxy = 0;
    for (std::size_t i = 0; i + 1 < x.size(); i++) sxx += (x[i] - init.theta) * (x[i] - init.theta), sxy += (x[i + 1] - init.theta) * (x[i] - init.theta);
    const double phi = std::clamp(sxx > 0 ? sxy / sxx : 0.5, 0.01, 0.999);
    const double dt = duration / static_cast<double>(x.size() - 1);
    init.kappa0 = -std::log(phi) / dt;
    init.sigma = std::sqrt(2 * init.kappa0 * std::max(1e-300, tyche::descriptive::variance(x, 1)));
  }
  // Parameters: θ, ln κ₀, ln σ, and the slope scaled by the record's duration.
  auto unpack = [&](std::span<const double> p) {
    OuFit f;
    f.theta = p[0], f.kappa0 = std::exp(p[1]), f.sigma = std::exp(p[2]);
    f.slope = drifting ? p[3] / duration * (d == Drift::fold ? f.kappa0 : 1.0) : 0.0;
    return f;
  };
  auto objective = [&](std::span<const double> p) {
    const OuFit f = unpack(p);
    return -log_likelihood(t, x, d, f, drifting);
  };
  std::vector<double> p0 = {init.theta, std::log(std::max(1e-12, init.kappa0)), std::log(std::max(1e-300, init.sigma))};
  if (drifting) p0.push_back(start && start->slope != 0 ? start->slope * duration / (d == Drift::fold ? init.kappa0 : 1.0) : 0.0);
  calx::optimisation::NelderMeadOptions o;
  o.max_evaluations = 2000;
  o.initial_step = drifting ? std::vector<double>{0.1 * std::max(1e-6, init.sigma / std::sqrt(2 * init.kappa0)), 0.2, 0.2, 0.2 * init.kappa0}
                            : std::vector<double>{0.1 * std::max(1e-6, init.sigma / std::sqrt(2 * init.kappa0)), 0.2, 0.2};
  const auto r = calx::optimisation::nelder_mead(objective, p0, o);
  OuFit out = unpack(r.x);
  out.log_likelihood = -r.value;
  out.converged = r.converged;
  return out;
}

// A record of the OU process with its rate at the times given, exact for a rate constant over each interval.
inline std::vector<double> simulate(const std::vector<double>& t, Drift d, const OuFit& f, bool drifting, calx::random::Stream& rng) {
  std::vector<double> x(t.size());
  const double k0 = drifting ? kappa_at(d, f, t[0], t[0]) : f.kappa0;
  x[0] = f.theta + (k0 > 0 ? f.sigma / std::sqrt(2 * k0) : f.sigma) * rng.normal();
  for (std::size_t i = 0; i + 1 < t.size(); i++) {
    const double k = drifting ? kappa_at(d, f, t[i], t[0]) : f.kappa0;
    double m, v;
    transition(x[i], k, f.theta, f.sigma, t[i + 1] - t[i], m, v);
    x[i + 1] = m + std::sqrt(std::max(0.0, v)) * rng.normal();
  }
  return x;
}

inline double crossing(Drift d, const OuFit& f, double t0) {
  if (d == Drift::linear) return f.slope < 0 ? t0 - f.kappa0 / f.slope : std::numeric_limits<double>::infinity();
  return f.slope < 0 ? t0 - f.kappa0 * f.kappa0 / f.slope : std::numeric_limits<double>::infinity();
}

}  // namespace detail

struct LikelihoodOptions {
  Drift drift = Drift::linear;
  std::size_t replicates = 100;
  double alpha = 0.05;
  std::uint64_t seed = 1;
};

inline LikelihoodTest likelihood_test(const std::vector<double>& times, const std::vector<double>& x, const LikelihoodOptions& o = {}) {
  if (times.size() != x.size()) throw Refusal("a likelihood test needs one value for each time");
  if (x.size() < 40) throw Refusal("a likelihood test of a drifting rate needs a record of at least 40 states");
  for (std::size_t i = 1; i < times.size(); i++)
    if (!(times[i] > times[i - 1])) throw Refusal("the times of a record must increase");
  LikelihoodTest out;
  out.drift = o.drift;
  out.null = detail::fit(times, x, o.drift, false, nullptr);
  OuFit start = out.null;
  start.slope = 0;
  out.alternative = detail::fit(times, x, o.drift, true, &start);
  if (out.alternative.log_likelihood < out.null.log_likelihood) {
    out.alternative = out.null;
    out.alternative.slope = 0;
  }
  out.deviance = 2 * (out.alternative.log_likelihood - out.null.log_likelihood);

  calx::random::Stream rng(o.seed, 0x6c696b65u, 0, 0);  // 'like'
  std::vector<double> crossings;
  auto deviance_of = [&](const std::vector<double>& y, double* t_c) {
    const OuFit f0 = detail::fit(times, y, o.drift, false, nullptr);
    OuFit s = f0;
    s.slope = 0;
    OuFit f1 = detail::fit(times, y, o.drift, true, &s);
    if (f1.log_likelihood < f0.log_likelihood) f1 = f0, f1.slope = 0;
    if (t_c) *t_c = detail::crossing(o.drift, f1, times.front());
    return 2 * (f1.log_likelihood - f0.log_likelihood);
  };
  for (std::size_t r = 0; r < o.replicates; r++) {
    out.null_deviances.push_back(deviance_of(detail::simulate(times, o.drift, out.null, false, rng), nullptr));
    double t_c = 0;
    out.alternative_deviances.push_back(deviance_of(detail::simulate(times, o.drift, out.alternative, true, rng), &t_c));
    crossings.push_back(t_c);
  }
  std::size_t above = 0;
  for (double v : out.null_deviances) above += v >= out.deviance;
  out.p = static_cast<double>(above + 1) / static_cast<double>(o.replicates + 1);
  std::vector<double> sorted = out.null_deviances;
  std::sort(sorted.begin(), sorted.end());
  const double critical = tyche::descriptive::quantile_sorted(sorted, 1 - o.alpha);
  std::size_t fa = 0, hit = 0;
  for (double v : out.null_deviances) fa += v > critical;
  for (double v : out.alternative_deviances) hit += v > critical;
  out.false_alarm = static_cast<double>(fa) / static_cast<double>(o.replicates);
  out.detection = static_cast<double>(hit) / static_cast<double>(o.replicates);
  double s = 0;
  for (double a : out.alternative_deviances)
    for (double b : out.null_deviances) s += a > b ? 1.0 : (a == b ? 0.5 : 0.0);
  out.auc = s / (static_cast<double>(o.replicates) * static_cast<double>(o.replicates));

  out.tipping = detail::crossing(o.drift, out.alternative, times.front());
  std::size_t unbounded = 0;
  std::vector<double> finite;
  for (double c : crossings) {
    if (std::isfinite(c)) finite.push_back(c);
    else unbounded++;
  }
  const bool significant = out.p < o.alpha;
  const bool bounded = static_cast<double>(unbounded) <= 0.025 * static_cast<double>(crossings.size()) && finite.size() >= 10;
  if (significant && bounded && std::isfinite(out.tipping)) {
    std::sort(finite.begin(), finite.end());
    out.tipping_lower = tyche::descriptive::quantile_sorted(finite, 0.025);
    out.tipping_upper = tyche::descriptive::quantile_sorted(finite, 0.975);
    out.tipping_reported = true;
    out.tipping_sentence = "under " + std::string(name(o.drift)) + " the rate would reach zero at " + std::to_string(out.tipping) +
                           " (95% of the bootstrap between " + std::to_string(out.tipping_lower) + " and " + std::to_string(out.tipping_upper) +
                           "), an extrapolation that holds only while the forcing goes on as it has";
  } else if (!significant) {
    out.tipping_sentence = "no time of tipping is given: the drifting rate does not fit better than a constant one at the stated level";
  } else {
    out.tipping_sentence = "no time of tipping is given: in " + std::to_string(unbounded) + " of " + std::to_string(crossings.size()) +
                           " replicates of the fitted alternative the rate never reaches zero, so the interval is unbounded";
  }
  out.assumptions = {"an Ornstein–Uhlenbeck process with constant mean and noise", std::string("the alternative: ") + name(o.drift),
                     "the rate taken constant over each sampling interval"};
  return out;
}

}  // namespace kinesis::transitions
