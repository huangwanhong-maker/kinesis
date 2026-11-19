// The indicators, each with what it estimates, the transitions it can see and those it is blind to.
//
// Every indicator here is a functional of a window of residual coordinates (the survey, Part II §11.2):
//   of one coordinate   variance; lag-1 autocorrelation; skewness; the restoring rate by each estimator of rate.hpp; the
//                       modulus of the dominant eigenvalue of an AR(1) or AR(2);
//   of many             the largest eigenvalue of the shrunk covariance; its growth relative to the first window's,
//                       λ_max(Ĉ₀⁻¹Ĉ), invariant under linear changes of coordinates; the modulus of the propagator's leading
//                       eigenvalue; the mean absolute correlation; the share of the variance on the first principal axis;
//   on a neighbourhood  the spatial variance, skewness and Moran's I of the residual field at each instant, averaged over the
//                       window (Kéfi et al. 2014, PLoS ONE 9, e92097; Dakos et al. 2010, Theor. Ecol. 3, 163–174), and each
//                       unit's variance relative to its own first window, averaged over units (baseline referencing, the
//                       survey Part II §2).
//
// What each is blind to is stated from the structure, before any data (Part II §1–2): every slowing-down indicator misses
// transitions without a local loss of stability; a covariance indicator misses a critical mode the noise does not excite;
// a spatial variance misses a critical mode uniform across the units, which shows in their mean and not in their spread.
#pragma once

#include <kinesis/transitions/kinds.hpp>
#include <kinesis/transitions/rate.hpp>
#include <kinesis/transitions/trend.hpp>

namespace kinesis::transitions {

inline const std::vector<std::string>& blind_to_all() {
  static const std::vector<std::string> b = {
      "rate-induced tipping with the frozen system's stability unchanged",
      "noise-induced tipping at fixed parameters",
      "global bifurcations (saddle-node on a cycle, homoclinic) and crises",
      "non-smooth transitions (a border collision, a threshold)"};
  return b;
}

inline std::vector<std::string> with(std::vector<std::string> a, const std::vector<std::string>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

namespace detail {

inline std::vector<double> column_rows(const Rectangular& r, std::size_t c, std::size_t from, std::size_t to) {
  std::vector<double> out(to - from);
  for (std::size_t t = from; t < to; t++) out[t - from] = r(t, c);
  return out;
}

inline const std::vector<std::string>& slowing_assumes() {
  static const std::vector<std::string> a = {"the parameter moves slowly against the system's own relaxation",
                                             "the noise is stationary in size and colour"};
  return a;
}

}  // namespace detail

// ---- one coordinate ---------------------------------------------------------------------------------------------

inline std::vector<Indicator> univariate_indicators(std::size_t c = 0, const std::vector<RateEstimator>& rates = {RateEstimator::ar1, RateEstimator::red}) {
  std::vector<Indicator> out;
  const std::vector<std::string> stationary = {"fold", "transcritical", "pitchfork"};
  out.push_back({"variance", "the variance σ²/(2κ) of the fluctuations", Direction::rising,
                 [c](const Rectangular& r, std::size_t a, std::size_t b, double) { return tyche::descriptive::variance(detail::column_rows(r, c, a, b), 1); },
                 with(stationary, {"Hopf"}), with(blind_to_all(), {"slowing down masked by noise that shrinks"}),
                 with(detail::slowing_assumes(), {"misled by noise that grows: a rise in variance alone is not slowing down"})});
  out.back().positive = true;
  out.push_back({"lag-1 autocorrelation", "the lag-1 coefficient e^{−κΔt}", Direction::rising,
                 [c](const Rectangular& r, std::size_t a, std::size_t b, double) {
                   const auto w = detail::column_rows(r, c, a, b);
                   return tyche::time_series::acf(w, 1)[1];
                 },
                 stationary, with(blind_to_all(), {"period-doubling (it tends to −1)", "a Hopf bifurcation sampled coarsely (it may fall)"}),
                 with(detail::slowing_assumes(), {"misled by noise that reddens"})});
  for (RateEstimator e : rates) {
    out.push_back({name(e), "the restoring rate κ, directly", Direction::falling,
                   [c, e](const Rectangular& r, std::size_t a, std::size_t b, double dt) {
                     const auto rate = restoring_rate(detail::column_rows(r, c, a, b), dt, e);
                     return rate.kappa;
                   },
                   stationary, with(blind_to_all(), {"period-doubling", "an oscillatory loss of stability"}),
                   e == RateEstimator::red ? std::vector<std::string>{"the parameter moves slowly against the system's own relaxation",
                                                                      "the noise is an Ornstein–Uhlenbeck process, the faster of the pair of rates"}
                                           : detail::slowing_assumes()});
  }
  out.push_back({"dominant eigenvalue modulus", "|μ| of an AR(1) or AR(2), approaching 1 at any local bifurcation of codimension one",
                 Direction::rising,
                 [c](const Rectangular& r, std::size_t a, std::size_t b, double) { return std::abs(ar_dominant(detail::column_rows(r, c, a, b)).mu); },
                 with(stationary, {"period-doubling", "Hopf", "Neimark–Sacker"}), blind_to_all(), detail::slowing_assumes()});
  out.push_back({"skewness", "the asymmetry of the fluctuations, which grows near a fold and with flickering", Direction::either,
                 [c](const Rectangular& r, std::size_t a, std::size_t b, double) { return tyche::descriptive::skewness(detail::column_rows(r, c, a, b), false); },
                 {"fold (asymmetric potential)", "flickering between states"}, with(blind_to_all(), {"symmetric bifurcations (pitchfork)"}),
                 {"misled by skewed noise"}});
  return out;
}

// ---- many coordinates -------------------------------------------------------------------------------------------

inline std::vector<Indicator> multivariate_indicators() {
  std::vector<Indicator> out;
  const std::vector<std::string> stationary = {"fold", "transcritical", "pitchfork"};
  const std::vector<std::string> unexcited = {"a critical mode the noise does not excite (wᵀB = 0), which shows no rise in variance"};
  auto window = [](const Rectangular& r, std::size_t a, std::size_t b) { return rows(r, a, b); };
  out.push_back({"largest covariance eigenvalue", "λ_max of the shrunk covariance, ≈ wᵀBBᵀw / (2|λ₁|) along the critical direction",
                 Direction::rising,
                 [window](const Rectangular& r, std::size_t a, std::size_t b, double) { return covariance_spectrum(covariance(window(r, a, b), true).c).largest; },
                 with(stationary, {"Hopf"}), with(blind_to_all(), unexcited),
                 with(detail::slowing_assumes(), {"depends on the units of the coordinates"})});
  out.back().positive = true;
  out.push_back({"relative variance growth", "λ_max(Ĉ₀⁻¹Ĉ) against the first window, invariant under linear changes of coordinates",
                 Direction::rising,
                 [window](const Rectangular& r, std::size_t a, std::size_t b, double) {
                   const std::size_t w = b - a;
                   const Matrix C0 = covariance(window(r, 0, w), true).c;
                   return relative_growth(C0, covariance(window(r, a, b), true).c);
                 },
                 with(stationary, {"Hopf"}), with(blind_to_all(), unexcited), detail::slowing_assumes()});
  out.back().positive = true;
  out.push_back({"leading propagator modulus", "|μ₁| of the lag-1 propagator Â = Ĉ(1)Ĉ(0)⁻¹, approaching 1 at any local bifurcation",
                 Direction::rising,
                 [window](const Rectangular& r, std::size_t a, std::size_t b, double) { return std::abs(propagator(window(r, a, b), 1e-4).spectrum.front()); },
                 with(stationary, {"period-doubling", "Hopf", "Neimark–Sacker"}), with(blind_to_all(), unexcited),
                 with(detail::slowing_assumes(), {"needs a window several times longer than the dimension"})});
  out.push_back({"mean absolute correlation", "the coordinates' mean |correlation|, which grows as one mode dominates", Direction::rising,
                 [window](const Rectangular& r, std::size_t a, std::size_t b, double) { return mean_absolute_correlation(covariance(window(r, a, b), false).c); },
                 stationary, with(blind_to_all(), with(unexcited, {"a critical mode confined to one coordinate"})), detail::slowing_assumes()});
  out.push_back({"first principal share", "λ_max / tr Ĉ, the share of the variance along the leading axis", Direction::rising,
                 [window](const Rectangular& r, std::size_t a, std::size_t b, double) { return covariance_spectrum(covariance(window(r, a, b), true).c).share; },
                 stationary, with(blind_to_all(), unexcited), with(detail::slowing_assumes(), {"depends on the units of the coordinates"})});
  out.back().positive = true;
  return out;
}

// ---- a neighbourhood of units -----------------------------------------------------------------------------------

namespace detail {

struct Weights {
  std::vector<std::vector<std::pair<std::size_t, double>>> rows;
  double total = 0;
};

inline Weights weights_of(const StateSpace& s) {
  Weights w;
  w.rows.assign(s.units, {});
  for (const auto& e : s.neighbourhood()) {
    w.rows[e.a].push_back({e.b, e.conductance});
    w.rows[e.b].push_back({e.a, e.conductance});
    w.total += 2 * e.conductance;
  }
  return w;
}

// Moran's I of a field z over units: (n / W) Σ w_ij z̃_i z̃_j / Σ z̃_i².
inline double moran(const std::vector<double>& z, const Weights& w) {
  const std::size_t n = z.size();
  const double m = std::accumulate(z.begin(), z.end(), 0.0) / static_cast<double>(n);
  double num = 0, den = 0;
  for (std::size_t i = 0; i < n; i++) {
    den += (z[i] - m) * (z[i] - m);
    for (const auto& [j, v] : w.rows[i]) num += v * (z[i] - m) * (z[j] - m);
  }
  return den > 0 && w.total > 0 ? static_cast<double>(n) / w.total * num / den : 0.0;
}

}  // namespace detail

inline std::vector<Indicator> spatial_indicators(const StateSpace& s) {
  std::vector<Indicator> out;
  const std::size_t per = s.coordinates_per_unit(), n = s.units;
  const auto w = detail::weights_of(s);
  const std::vector<std::string> sees = {"fold", "transcritical", "pitchfork", "with a critical mode that is not uniform across the units"};
  const std::vector<std::string> uniform = {"a critical mode uniform across the units, which shows in their mean and not in their spread"};
  // The field of coordinate k over units at time t.
  auto field = [per, n](const Rectangular& r, std::size_t t, std::size_t k) {
    std::vector<double> z(n);
    for (std::size_t u = 0; u < n; u++) z[u] = r(t, u * per + k);
    return z;
  };
  out.push_back({"spatial variance", "the variance of the residual field across units, averaged over the window", Direction::rising,
                 [field, per](const Rectangular& r, std::size_t a, std::size_t b, double) {
                   double s2 = 0;
                   for (std::size_t t = a; t < b; t++)
                     for (std::size_t k = 0; k < per; k++) s2 += tyche::descriptive::variance(field(r, t, k), 1);
                   return s2 / static_cast<double>((b - a) * per);
                 },
                 sees, with(blind_to_all(), uniform), detail::slowing_assumes()});
  out.back().positive = true;
  out.push_back({"spatial skewness", "the asymmetry of the residual field across units", Direction::either,
                 [field, per](const Rectangular& r, std::size_t a, std::size_t b, double) {
                   double s3 = 0;
                   for (std::size_t t = a; t < b; t++)
                     for (std::size_t k = 0; k < per; k++) s3 += tyche::descriptive::skewness(field(r, t, k), false);
                   return s3 / static_cast<double>((b - a) * per);
                 },
                 {"fold (asymmetric potential) in a spatial field", "patches flipping state"}, with(blind_to_all(), uniform), {"misled by skewed noise"}});
  if (w.total > 0) {
    out.push_back({"Moran's I", "the correlation of neighbouring units' fluctuations, which grows as the critical mode's reach does",
                   Direction::rising,
                   [field, per, w](const Rectangular& r, std::size_t a, std::size_t b, double) {
                     double m = 0;
                     for (std::size_t t = a; t < b; t++)
                       for (std::size_t k = 0; k < per; k++) m += detail::moran(field(r, t, k), w);
                     return m / static_cast<double>((b - a) * per);
                   },
                   sees, with(blind_to_all(), {"a critical mode that alternates between neighbours (it lowers Moran's I)"}), detail::slowing_assumes()});
  }
  out.push_back({"relative unit variance", "each coordinate's variance relative to its own first window, averaged", Direction::rising,
                 [](const Rectangular& r, std::size_t a, std::size_t b, double) {
                   const std::size_t w0 = b - a;
                   double s = 0;
                   std::size_t count = 0;
                   for (std::size_t c = 0; c < r.cols; c++) {
                     const double base = tyche::descriptive::variance(detail::column_rows(r, c, 0, w0), 1);
                     if (base > 0) s += tyche::descriptive::variance(detail::column_rows(r, c, a, b), 1) / base, count++;
                   }
                   return count ? s / static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
                 },
                 {"fold", "transcritical", "pitchfork", "Hopf"}, with(blind_to_all(), {"a critical mode the noise does not excite"}),
                 detail::slowing_assumes()});
  out.back().positive = true;
  return out;
}

}  // namespace kinesis::transitions
