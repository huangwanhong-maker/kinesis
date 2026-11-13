// Signs of criticality, read from records: slowing down before a transition, correlations without a scale, and
// avalanches.
//
// These work on a record alone. Each is an indicator with known ways of misleading, which are stated; none is a
// certificate of a transition or of a critical state.
//
// Slowing down (Wissel 1984; Held and Kleinen 2004, Geophys. Res. Lett. 31, L23207; Scheffer et al. 2009, Nature
// 461, 53–59; the methods as Dakos et al. 2012, PLoS ONE 7, e41010). Near a bifurcation where an equilibrium
// loses stability, the rate λ at which a perturbation decays goes to zero. Fluctuations about the equilibrium are
// then an Ornstein–Uhlenbeck process to first order, whose variance σ²/(2λ) and autocorrelation at the sampling
// lag e^{−λΔt} both rise. early_warnings() removes a slow trend (a Gaussian kernel smoother), takes the variance,
// the autocorrelation at lag 1 and the skewness of the residuals in a rolling window, and measures each one's
// trend by Kendall's τ against time. For the variance and the autocorrelation a significance is given: the
// share of surrogate series with a τ at least as large, the surrogates being red noise fitted to the residuals
// (tyche's), the null of a stationary linear process with the series' own correlation. The skewness is given
// with its τ alone.
// *Caveat.* A rising indicator is not a transition: trends arise by chance in correlated series (which the
// surrogates measure), a window chosen after looking inflates them (Boettiger and Hastings 2012, J. R. Soc.
// Interface 9, 2527–2539), and transitions driven by noise or by the rate of change of a parameter need give no
// warning at all (Ditlevsen and Johnsen 2010, Geophys. Res. Lett. 37, L19703; Ashwin et al. 2012, Phil. Trans.
// R. Soc. A 370, 1166–1184).
//
// Correlations without a scale. Detrended fluctuation analysis (Peng et al. 1994, Phys. Rev. E 49, 1685–1689;
// Hardstone et al. 2012, Front. Physiol. 3, 450): the profile of the series is cut into windows of n samples, a
// line is removed in each, and the root-mean-square residual F(n) grows as n^α: α = ½ for a series without
// memory, 1 for 1/f noise, 3/2 for a random walk.
//
// Avalanches (Beggs and Plenz 2003, J. Neurosci. 23, 11167–11177). An avalanche is a run of consecutive bins in
// which the activity is above a threshold; its size is the activity summed over the run. At a critical branching
// process sizes follow a power law of exponent 3/2. The exponent of a power law above x_min is fitted by maximum
// likelihood (Clauset, Shalizi and Newman 2009, SIAM Rev. 51, 661–703, Eqs. 3.1 and 3.7), with its standard
// error (α − 1)/√n and the Kolmogorov–Smirnov distance of the fit.
// *Caveat.* A power law over a decade or two does not show criticality: thresholds and mechanisms that are not
// critical give the same exponents (Touboul and Destexhe 2017, Phys. Rev. E 95, 012413).
//
// The branching ratio m is the mean number of units one active unit activates in the next step; m = 1 is
// critical. The slope r_1 of the regression of A_{t+1} on A_t estimates it when every unit is observed, and is
// biased when only a share is. The slopes r_k at lags k decay as r_k = b mᵏ, with only b depending on the share
// observed, so that m is read from the decay (Wilting and Priesemann 2018, Nat. Commun. 9, 2325). Here the fit is
// a least-squares line through log r_k, on the lags at which r_k is positive.
#pragma once

#include <kinesis/kernel/answer.hpp>

#include <calx/random/stream.hpp>
#include <tyche/descriptive/rank.hpp>
#include <tyche/descriptive/summary.hpp>
#include <tyche/time_series/autocorrelation.hpp>
#include <tyche/time_series/spectral.hpp>

#include <algorithm>
#include <cmath>
#include <span>
#include <string>
#include <vector>

namespace kinesis::analysis {

// ---- slowing down ----------------------------------------------------------------------------------------

struct WarningOptions {
  double window = 0.5;          // the rolling window, as a share of the series
  double bandwidth = 0.1;       // the smoother's standard deviation, as a share of the series; 0 removes no trend
  std::size_t surrogates = 200; // 0 gives no significance
  std::uint64_t seed = 1;
};

struct Warnings {
  std::vector<double> at;                                  // the index at which each window ends
  std::vector<double> variance, autocorrelation, skewness; // in each window
  double tau_variance = 0, tau_autocorrelation = 0, tau_skewness = 0;
  double p_variance = 1, p_autocorrelation = 1;            // the share of surrogates with τ at least as large
  std::vector<double> trend, residual;
};

// The series smoothed by a Gaussian kernel, the slow trend that early_warnings() takes away; `sd` is the kernel's
// standard deviation in samples.
inline std::vector<double> smooth(std::span<const double> x, double sd) {
  const std::size_t n = x.size();
  std::vector<double> out(n);
  const long reach = static_cast<long>(std::ceil(4 * sd));
  for (std::size_t i = 0; i < n; i++) {
    double sum = 0, weight = 0;
    for (long k = -reach; k <= reach; k++) {
      const long j = static_cast<long>(i) + k;
      if (j < 0 || j >= static_cast<long>(n)) continue;
      const double w = std::exp(-0.5 * static_cast<double>(k * k) / (sd * sd));
      sum += w * x[static_cast<std::size_t>(j)], weight += w;
    }
    out[i] = sum / weight;
  }
  return out;
}

namespace detail {
inline double lag_one(std::span<const double> x) {
  const auto r = tyche::time_series::acf(x, 1);
  return r.size() > 1 ? r[1] : 0;
}
struct Rolling {
  std::vector<double> at, variance, autocorrelation, skewness;
};
inline Rolling rolling(std::span<const double> residual, std::size_t window) {
  Rolling r;
  for (std::size_t end = window; end <= residual.size(); end++) {
    const std::span<const double> w = residual.subspan(end - window, window);
    r.at.push_back(static_cast<double>(end - 1));
    r.variance.push_back(tyche::descriptive::variance(w));
    r.autocorrelation.push_back(lag_one(w));
    r.skewness.push_back(tyche::descriptive::skewness(w));
  }
  return r;
}
inline double trend_of(const std::vector<double>& at, const std::vector<double>& y) { return tyche::descriptive::kendall_tau_b(at, y); }
}  // namespace detail

inline Warnings early_warnings(std::span<const double> x, const WarningOptions& o = {}) {
  const std::size_t n = x.size();
  const std::size_t window = static_cast<std::size_t>(std::round(o.window * static_cast<double>(n)));
  if (n < 40 || window < 10 || window >= n) throw kernel::Refusal("early warnings need a series of at least 40 samples and a window within it");
  Warnings w;
  w.trend = o.bandwidth > 0 ? smooth(x, o.bandwidth * static_cast<double>(n)) : std::vector<double>(n, tyche::descriptive::mean(x));
  w.residual.resize(n);
  for (std::size_t i = 0; i < n; i++) w.residual[i] = x[i] - w.trend[i];
  const auto r = detail::rolling(w.residual, window);
  w.at = r.at, w.variance = r.variance, w.autocorrelation = r.autocorrelation, w.skewness = r.skewness;
  w.tau_variance = detail::trend_of(w.at, w.variance);
  w.tau_autocorrelation = detail::trend_of(w.at, w.autocorrelation);
  w.tau_skewness = detail::trend_of(w.at, w.skewness);
  if (o.surrogates > 0) {
    const auto red = tyche::time_series::red_noise(w.residual);
    calx::random::Stream rng(o.seed, 0x7761726eu, 0, 0);  // 'warn'
    std::size_t above_variance = 0, above_autocorrelation = 0;
    for (std::size_t k = 0; k < o.surrogates; k++) {
      const std::vector<double> s = tyche::time_series::red_noise_surrogate(red, n, rng);
      const auto rs = detail::rolling(s, window);
      above_variance += detail::trend_of(rs.at, rs.variance) >= w.tau_variance;
      above_autocorrelation += detail::trend_of(rs.at, rs.autocorrelation) >= w.tau_autocorrelation;
    }
    // With the series itself counted among the draws, so that a significance is never zero.
    w.p_variance = static_cast<double>(above_variance + 1) / static_cast<double>(o.surrogates + 1);
    w.p_autocorrelation = static_cast<double>(above_autocorrelation + 1) / static_cast<double>(o.surrogates + 1);
  }
  return w;
}

// ---- detrended fluctuation analysis ----------------------------------------------------------------------

struct Fluctuation {
  std::vector<double> scales, fluctuation;  // n and F(n)
  double exponent = 0;                      // the slope of log F against log n
};

inline Fluctuation dfa(std::span<const double> x, std::size_t smallest = 8, std::size_t count = 12) {
  const std::size_t N = x.size();
  if (N < 16 * smallest) throw kernel::Refusal("detrended fluctuation analysis needs a longer series");
  const double mean = tyche::descriptive::mean(x);
  std::vector<double> profile(N);
  double run = 0;
  for (std::size_t i = 0; i < N; i++) run += x[i] - mean, profile[i] = run;
  Fluctuation f;
  const double largest = static_cast<double>(N) / 8;
  for (std::size_t k = 0; k < count; k++) {
    const std::size_t n = static_cast<std::size_t>(std::round(static_cast<double>(smallest) * std::pow(largest / static_cast<double>(smallest), static_cast<double>(k) / static_cast<double>(count - 1))));
    if (!f.scales.empty() && static_cast<double>(n) == f.scales.back()) continue;
    double total = 0;
    std::size_t windows = 0;
    for (std::size_t start = 0; start + n <= N; start += n, windows++) {
      // The least-squares line through the profile in the window, and the mean square about it.
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      for (std::size_t i = 0; i < n; i++) {
        const double t = static_cast<double>(i), y = profile[start + i];
        sx += t, sy += y, sxx += t * t, sxy += t * y;
      }
      const double m = static_cast<double>(n), slope = (m * sxy - sx * sy) / (m * sxx - sx * sx), intercept = (sy - slope * sx) / m;
      for (std::size_t i = 0; i < n; i++) {
        const double e = profile[start + i] - (intercept + slope * static_cast<double>(i));
        total += e * e;
      }
    }
    f.scales.push_back(static_cast<double>(n));
    f.fluctuation.push_back(std::sqrt(total / static_cast<double>(windows * n)));
  }
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (std::size_t k = 0; k < f.scales.size(); k++) {
    const double a = std::log(f.scales[k]), b = std::log(f.fluctuation[k]);
    sx += a, sy += b, sxx += a * a, sxy += a * b;
  }
  const double m = static_cast<double>(f.scales.size());
  f.exponent = (m * sxy - sx * sy) / (m * sxx - sx * sx);
  return f;
}

// ---- avalanches and power laws ---------------------------------------------------------------------------

struct Avalanches {
  std::vector<double> sizes, durations;
};

inline Avalanches avalanches(std::span<const double> activity, double threshold = 0) {
  Avalanches a;
  double size = 0, duration = 0;
  for (double v : activity) {
    if (v > threshold) {
      size += v, duration += 1;
    } else if (duration > 0) {
      a.sizes.push_back(size), a.durations.push_back(duration);
      size = 0, duration = 0;
    }
  }
  if (duration > 0) a.sizes.push_back(size), a.durations.push_back(duration);
  return a;
}

struct PowerLaw {
  double exponent = 0, error = 0;   // α and (α − 1)/√n
  double xmin = 0;
  std::size_t count = 0;            // the values at or above xmin
  double distance = 0;              // Kolmogorov and Smirnov's, between the values and the fitted law
};

// The power law p(x) ∝ x^{−α} above xmin, by maximum likelihood. `discrete` uses the approximation for integer
// values, α = 1 + n / Σ ln(x_i / (xmin − ½)), good for xmin of about 6 and more.
inline PowerLaw power_law(std::span<const double> values, double xmin, bool discrete = false) {
  std::vector<double> x;
  for (double v : values)
    if (v >= xmin) x.push_back(v);
  if (x.size() < 10) throw kernel::Refusal("a power law is fitted to at least 10 values at or above xmin; there are " + std::to_string(x.size()));
  std::sort(x.begin(), x.end());
  const double base = discrete ? xmin - 0.5 : xmin;
  double logs = 0;
  for (double v : x) logs += std::log(v / base);
  PowerLaw p;
  p.xmin = xmin, p.count = x.size();
  p.exponent = 1 + static_cast<double>(x.size()) / logs;
  p.error = (p.exponent - 1) / std::sqrt(static_cast<double>(x.size()));
  for (std::size_t i = 0; i < x.size(); i++) {
    const double model = 1 - std::pow(x[i] / base, 1 - p.exponent);
    p.distance = std::max({p.distance, std::fabs(model - static_cast<double>(i) / static_cast<double>(x.size())), std::fabs(model - static_cast<double>(i + 1) / static_cast<double>(x.size()))});
  }
  return p;
}

// ---- the branching ratio ---------------------------------------------------------------------------------

struct Branching {
  double conventional = 0;          // r_1
  double ratio = 0;                 // m, from the decay of r_k
  double amplitude = 0;             // b
  std::vector<double> slopes;       // r_1 … r_K
  std::size_t used = 0;             // the lags in the fit
};

inline Branching branching_ratio(std::span<const double> activity, std::size_t lags = 40) {
  const std::size_t n = activity.size();
  if (n < 10 * lags) throw kernel::Refusal("the branching ratio needs a series at least ten times its longest lag");
  Branching b;
  for (std::size_t k = 1; k <= lags; k++) {
    // The slope of the least-squares line of A_{t+k} on A_t.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const std::size_t m = n - k;
    for (std::size_t t = 0; t < m; t++) sx += activity[t], sy += activity[t + k], sxx += activity[t] * activity[t], sxy += activity[t] * activity[t + k];
    const double md = static_cast<double>(m);
    b.slopes.push_back((md * sxy - sx * sy) / (md * sxx - sx * sx));
  }
  b.conventional = b.slopes.front();
  double sk = 0, sl = 0, skk = 0, skl = 0;
  for (std::size_t k = 1; k <= lags; k++) {
    const double r = b.slopes[k - 1];
    if (!(r > 0)) break;  // the decay is read as far as it stays above zero
    const double kd = static_cast<double>(k), l = std::log(r);
    sk += kd, sl += l, skk += kd * kd, skl += kd * l, b.used++;
  }
  if (b.used < 3) throw kernel::Refusal("the autocorrelation of the activity does not decay over three lags or more: no branching ratio is read");
  const double u = static_cast<double>(b.used), slope = (u * skl - sk * sl) / (u * skk - sk * sk);
  b.ratio = std::exp(slope), b.amplitude = std::exp((sl - slope * sk) / u);
  return b;
}

}  // namespace kinesis::analysis
