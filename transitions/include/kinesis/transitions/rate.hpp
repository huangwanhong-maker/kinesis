// The restoring rate κ, measured directly from a window of one quantity.
//
// Near a stable equilibrium a quantity is, to first order, an Ornstein–Uhlenbeck process dx = −κ x dt + σ dW, whose
// variance σ²/(2κ) and lag-1 autocorrelation e^{−κΔt} are the classical indicators. Both rise as κ falls, but both also
// rise when the noise grows or reddens with κ unchanged (Morr and Boers 2024, Phys. Rev. X 14, 021037, read in full
// by the survey, Part I §3). Estimating κ itself separates the two. Four estimators share one answer:
//
//   ar1        the least-squares lag-1 coefficient φ̂ of the centred window, κ̂ = −ln φ̂ / Δt; the white-noise case.
//   increment  the regression of increments on the state, Δx_t = λ x_t + e_t, with e_t an AR(1) process estimated with it by
//              Cochrane and Orcutt's iteration (1949, J. Am. Stat. Assoc. 44, 32–61): the generalised least squares Boers
//              (2021, Nat. Clim. Change 11, 680–688) fits; κ̂ = −ln(1 + λ̂)/Δt.
//   red        the autocorrelation of an OU process driven by OU noise of rate θ, fitted over lags 1…L:
//                  r(τ) = (θ e^{−κτ} − κ e^{−θτ}) / (θ − κ),   with r(τ) = (1 + κτ) e^{−κτ} at θ = κ
//              (the spectrum σ²/((κ² + ω²)(θ² + ω²)) split into partial fractions; the theory note, §6.1). r is symmetric in
//              κ and θ, so the record identifies the pair of rates and not which is the system's: the answer gives both,
//              and takes the slower as the restoring rate only under the stated assumption that the noise is the faster
//              (Morr and Boers's estimator from the autocorrelation; their estimator from the spectrum carries the same
//              information).
//   drift      the local slope b of the conditional mean of increments, E[Δx | x], about the window's mode (a Kramers–
//              Moyal estimate of the drift's derivative, Friedrich et al. 2011, Phys. Rep. 506, 87–162): for OU
//              sampled at Δt the slope is e^{−κΔt} − 1 exactly, so κ̂ = −ln(1 + b)/Δt. It reads the restoring force where
//              the state stays, which matters when the drift is not linear across the window.
//
// Every answer carries the diffusion σ² that makes the window's variance σ²/(2κ) when the noise is white: comparing how
// variance and κ̂ move across windows tells slowing down (κ falls, variance rises in proportion) from a change of noise
// (variance rises, κ does not fall).
#pragma once

#include <kinesis/transitions/space.hpp>

#include <calx/optimisation/nelder_mead.hpp>
#include <tyche/descriptive/summary.hpp>
#include <tyche/time_series/autocorrelation.hpp>

#include <numbers>
#include <span>

namespace kinesis::transitions {

enum class RateEstimator { ar1, increment, red, drift };

inline const char* name(RateEstimator e) {
  switch (e) {
    case RateEstimator::ar1: return "restoring rate (AR(1))";
    case RateEstimator::increment: return "restoring rate (increments, AR(1) noise)";
    case RateEstimator::red: return "restoring rate (red noise)";
    default: return "restoring rate (drift slope)";
  }
}

struct Rate {
  RateEstimator estimator = RateEstimator::ar1;
  double kappa = 0;                       // κ̂; ≤ 0 when no restoring force is measured
  double lower = 0, upper = 0;            // an approximate 95% interval
  double phi = 0;                         // the lag-1 coefficient the estimate rests on
  double diffusion = 0;                   // σ̂² with variance σ²/(2κ) under white noise
  double noise_rate = 0;                  // red: the other rate of the pair (θ̂); 0 otherwise
  bool red_noise = false;                 // red: the red-noise autocorrelation fits markedly better than white noise's e^{−κτ}
  double fit_red = 0, fit_white = 0;      // red: the weighted residuals of the two fits
  bool restoring = true;                  // false when φ̂ ≥ 1: no restoring force is measured
  std::vector<std::string> assumptions;
};

namespace detail {

inline std::vector<double> centred(std::span<const double> x) {
  const double m = tyche::descriptive::mean(x);
  std::vector<double> out(x.begin(), x.end());
  for (double& v : out) v -= m;
  return out;
}

inline double kappa_of(double phi, double dt) { return phi > 0 && phi < 1 ? -std::log(phi) / dt : (phi >= 1 ? 0.0 : std::numeric_limits<double>::infinity()); }

inline void interval_from_phi(Rate& r, double phi, double se, double dt) {
  const double hi_phi = std::min(phi + 1.96 * se, 1 - 1e-12), lo_phi = std::max(phi - 1.96 * se, 1e-12);
  r.lower = kappa_of(hi_phi, dt);
  r.upper = kappa_of(lo_phi, dt);
}

inline double diffusion_of(double kappa, double variance) { return kappa > 0 ? 2 * kappa * variance : 0; }

}  // namespace detail

inline Rate ar1_rate(std::span<const double> x, double dt) {
  const std::size_t n = x.size();
  if (n < 10) throw Refusal("a restoring rate needs a window of at least 10 states");
  const auto c = detail::centred(x);
  double sxx = 0, sxy = 0;
  for (std::size_t t = 0; t + 1 < n; t++) sxx += c[t] * c[t], sxy += c[t + 1] * c[t];
  if (!(sxx > 0)) throw Refusal("a window that does not vary has no restoring rate");
  Rate r;
  r.estimator = RateEstimator::ar1;
  r.phi = sxy / sxx;
  r.restoring = r.phi < 1;
  r.kappa = detail::kappa_of(r.phi, dt);
  const double se = std::sqrt(std::max(1e-12, 1 - r.phi * r.phi) / static_cast<double>(n - 1));
  detail::interval_from_phi(r, r.phi, se, dt);
  r.diffusion = detail::diffusion_of(r.kappa, tyche::descriptive::variance(c, 1));
  r.assumptions = {"a linear restoring force with white, stationary noise within the window"};
  return r;
}

inline Rate increment_rate(std::span<const double> x, double dt, int iterations = 20) {
  const std::size_t n = x.size();
  if (n < 12) throw Refusal("a restoring rate from increments needs a window of at least 12 states");
  const auto c = detail::centred(x);
  // Δx_t = x_{t+1} − x_t against x_t, t = 0 … n − 2.
  std::vector<double> dx(n - 1), xs(n - 1);
  for (std::size_t t = 0; t + 1 < n; t++) dx[t] = c[t + 1] - c[t], xs[t] = c[t];
  auto fit = [](const std::vector<double>& y, const std::vector<double>& z) {
    double szz = 0, szy = 0;
    for (std::size_t i = 0; i < y.size(); i++) szz += z[i] * z[i], szy += z[i] * y[i];
    return szz > 0 ? szy / szz : 0.0;
  };
  double lambda = fit(dx, xs), rho = 0, se = 0;
  for (int k = 0; k < iterations; k++) {
    std::vector<double> e(n - 1);
    for (std::size_t t = 0; t + 1 < n; t++) e[t] = dx[t] - lambda * xs[t];
    double see = 0, sel = 0;
    for (std::size_t t = 1; t + 1 < n; t++) see += e[t - 1] * e[t - 1], sel += e[t] * e[t - 1];
    const double rho_new = see > 0 ? std::clamp(sel / see, -0.99, 0.99) : 0.0;
    std::vector<double> y(n - 2), z(n - 2);
    for (std::size_t t = 1; t + 1 < n; t++) y[t - 1] = dx[t] - rho_new * dx[t - 1], z[t - 1] = xs[t] - rho_new * xs[t - 1];
    const double lambda_new = fit(y, z);
    double szz = 0, ssr = 0;
    for (std::size_t i = 0; i < y.size(); i++) {
      szz += z[i] * z[i];
      const double u = y[i] - lambda_new * z[i];
      ssr += u * u;
    }
    se = szz > 0 ? std::sqrt(ssr / static_cast<double>(y.size() - 1) / szz) : 0;
    const bool settled = std::abs(lambda_new - lambda) < 1e-10 && std::abs(rho_new - rho) < 1e-10;
    lambda = lambda_new, rho = rho_new;
    if (settled) break;
  }
  Rate r;
  r.estimator = RateEstimator::increment;
  r.phi = 1 + lambda;
  r.restoring = r.phi < 1;
  r.kappa = detail::kappa_of(r.phi, dt);
  detail::interval_from_phi(r, r.phi, se, dt);
  r.noise_rate = 0;
  r.diffusion = detail::diffusion_of(r.kappa, tyche::descriptive::variance(c, 1));
  r.assumptions = {"a linear restoring force; the noise an AR(1) process of lag-1 coefficient " + std::to_string(rho)};
  return r;
}

// The autocorrelation at lag τ of an OU process of rate κ driven by OU noise of rate θ.
inline double red_autocorrelation(double kappa, double theta, double tau) {
  if (std::abs(theta - kappa) <= 1e-9 * std::max(kappa, theta)) return (1 + kappa * tau) * std::exp(-kappa * tau);
  return (theta * std::exp(-kappa * tau) - kappa * std::exp(-theta * tau)) / (theta - kappa);
}

inline Rate red_rate(std::span<const double> x, double dt, std::size_t lags = 0, bool noise_faster = true) {
  const std::size_t n = x.size();
  if (n < 30) throw Refusal("a restoring rate under red noise needs a window of at least 30 states");
  // The lags fitted span the record's correlation: up to where its autocorrelation first falls below 0.1, at least 4
  // and at most a quarter of the window; the starting rates are scaled by the lag τ_e at which it falls below 1/e.
  const auto probe = tyche::time_series::acf(x, n / 4);
  std::size_t L = lags, crossing = probe.size() - 1;
  for (std::size_t k = 1; k < probe.size(); k++)
    if (probe[k] < std::exp(-1.0)) {
      crossing = k;
      break;
    }
  if (L == 0) {
    L = n / 4;
    for (std::size_t k = 1; k < probe.size(); k++)
      if (probe[k] < 0.1) {
        L = k;
        break;
      }
    L = std::clamp<std::size_t>(L, 4, n / 4);
  }
  const auto r_hat = tyche::time_series::acf(x, L);
  const double k0 = 1 / (static_cast<double>(crossing) * dt);
  // Each lag weighted by the inverse of Bartlett's variance of its estimate, (1 + 2 Σ_{j<k} r_j²)/n (Bartlett 1946, J. R.
  // Stat. Soc. Suppl. 8, 27–41), so that the long lags, estimated worst, count least.
  std::vector<double> weight(L + 1, 1);
  double cumulative = 0;
  for (std::size_t k = 1; k <= L; k++) {
    weight[k] = 1 / (1 + 2 * cumulative);
    cumulative += r_hat[k] * r_hat[k];
  }
  auto loss = [&](std::span<const double> p) {
    const double a = std::exp(p[0]), b = std::exp(p[1]);
    double s = 0;
    for (std::size_t k = 1; k <= L; k++) {
      const double e = r_hat[k] - red_autocorrelation(a, b, static_cast<double>(k) * dt);
      s += weight[k] * e * e;
    }
    // Rates far above the sampling rate are not resolved: hold them below 20/Δt.
    for (double q : {a, b})
      if (q > 20 / dt) s += (q * dt - 20) * (q * dt - 20);
    return s;
  };
  calx::optimisation::NelderMeadOptions o;
  o.initial_step = {0.5, 0.5};
  calx::optimisation::NelderMeadResult best;
  best.value = std::numeric_limits<double>::infinity();
  for (const auto& start : {std::pair{0.5, 3.0}, std::pair{1.0, 10.0}, std::pair{0.3, 1.5}, std::pair{0.8, 40.0}}) {
    const auto r = calx::optimisation::nelder_mead(loss, {std::log(start.first * k0), std::log(std::min(start.second * k0, 15 / dt))}, o);
    if (r.value < best.value) best = r;
  }
  const double a = std::exp(best.x[0]), b = std::exp(best.x[1]);
  // White noise's autocorrelation e^{−κτ}, fitted the same way: the red fit is taken to show red noise when it leaves less than
  // three quarters of the white fit's weighted residual and its noise rate is resolved by the sampling (θΔt < 3). Both bounds are
  // declared, not derived: provisional, and stated wherever the judgement is used.
  auto white = [&](std::span<const double> p) {
    const double k = std::exp(p[0]);
    double s = 0;
    for (std::size_t j = 1; j <= L; j++) {
      const double e = r_hat[j] - std::exp(-k * static_cast<double>(j) * dt);
      s += weight[j] * e * e;
    }
    return s;
  };
  calx::optimisation::NelderMeadOptions one;
  one.initial_step = {0.5};
  const auto w_fit = calx::optimisation::nelder_mead(white, {std::log(k0)}, one);
  Rate r;
  r.red_noise = best.value < 0.75 * w_fit.value && std::max(a, b) * dt < 3;
  r.fit_red = best.value, r.fit_white = w_fit.value;
  r.estimator = RateEstimator::red;
  const double slow = std::min(a, b), fast = std::max(a, b);
  r.kappa = noise_faster ? slow : fast;
  r.noise_rate = noise_faster ? fast : slow;
  r.phi = std::exp(-r.kappa * dt);
  r.restoring = true;
  // The interval: the ar1 interval's relative width, carried over (an approximation, said so).
  const Rate w = ar1_rate(x, dt);
  const double rel = w.kappa > 0 && std::isfinite(w.upper) ? (w.upper - w.lower) / (2 * w.kappa) : 1.0;
  r.lower = r.kappa * std::max(0.0, 1 - rel);
  r.upper = r.kappa * (1 + rel);
  r.diffusion = detail::diffusion_of(r.kappa, tyche::descriptive::variance(x, 1));
  r.assumptions = {"a linear restoring force driven by Ornstein–Uhlenbeck noise",
                   "the record identifies the pair of rates {" + std::to_string(slow) + ", " + std::to_string(fast) + "}; " +
                       (noise_faster ? "the noise is taken as the faster" : "the noise is taken as the slower"),
                   "the interval is the AR(1) interval's relative width, an approximation"};
  return r;
}

inline Rate drift_rate(std::span<const double> x, double dt) {
  const std::size_t n = x.size();
  if (n < 30) throw Refusal("a drift slope needs a window of at least 30 states");
  const auto c = detail::centred(x);
  const double sd = std::sqrt(tyche::descriptive::variance(c, 1));
  if (!(sd > 0)) throw Refusal("a window that does not vary has no drift");
  // The window's mode: the peak of a Gaussian kernel density on a grid (Silverman's bandwidth).
  const double h = 1.06 * sd * std::pow(static_cast<double>(n), -0.2);
  double lo = *std::min_element(c.begin(), c.end()), hi = *std::max_element(c.begin(), c.end()), mode = 0, top = -1;
  for (int g = 0; g <= 200; g++) {
    const double z = lo + (hi - lo) * g / 200.0;
    double f = 0;
    for (double v : c) f += std::exp(-0.5 * (z - v) * (z - v) / (h * h));
    if (f > top) top = f, mode = z;
  }
  // Increments against the state, weighted by a Gaussian of one standard deviation about the mode.
  double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
  for (std::size_t t = 0; t + 1 < n; t++) {
    const double w = std::exp(-0.5 * (c[t] - mode) * (c[t] - mode) / (sd * sd));
    const double y = c[t + 1] - c[t];
    sw += w, swx += w * c[t], swy += w * y, swxx += w * c[t] * c[t], swxy += w * c[t] * y;
  }
  const double denom = sw * swxx - swx * swx;
  if (!(denom > 0)) throw Refusal("the drift's slope is not determined: the window's states do not spread about its mode");
  const double b = (sw * swxy - swx * swy) / denom;
  double ssr = 0;
  const double a0 = (swy - b * swx) / sw;
  for (std::size_t t = 0; t + 1 < n; t++) {
    const double w = std::exp(-0.5 * (c[t] - mode) * (c[t] - mode) / (sd * sd));
    const double e = c[t + 1] - c[t] - a0 - b * c[t];
    ssr += w * e * e;
  }
  const double effective = sw * sw / [&] {
    double s2 = 0;
    for (std::size_t t = 0; t + 1 < n; t++) {
      const double w = std::exp(-0.5 * (c[t] - mode) * (c[t] - mode) / (sd * sd));
      s2 += w * w;
    }
    return s2;
  }();
  const double se = std::sqrt(ssr / sw / std::max(1.0, effective - 2) * sw / (swxx - swx * swx / sw));
  Rate r;
  r.estimator = RateEstimator::drift;
  r.phi = 1 + b;
  r.restoring = r.phi < 1;
  r.kappa = detail::kappa_of(r.phi, dt);
  detail::interval_from_phi(r, r.phi, se, dt);
  r.diffusion = detail::diffusion_of(r.kappa, tyche::descriptive::variance(c, 1));
  r.assumptions = {"the drift linear about the window's mode, read with a Gaussian weight of one standard deviation"};
  return r;
}

inline Rate restoring_rate(std::span<const double> x, double dt, RateEstimator e) {
  switch (e) {
    case RateEstimator::ar1: return ar1_rate(x, dt);
    case RateEstimator::increment: return increment_rate(x, dt);
    case RateEstimator::red: return red_rate(x, dt);
    default: return drift_rate(x, dt);
  }
}

}  // namespace kinesis::transitions
