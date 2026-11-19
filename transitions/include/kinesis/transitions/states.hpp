// How many states a record visits, and how often it switches between them: the one family that speaks to noise-induced
// transitions (the survey, Part I §3).
//
// A record that dwells in one well of a landscape and is pushed over a barrier by noise gives no slowing down: its
// restoring rate in the well is unchanged until it leaves (Ditlevsen and Johnsen 2010, Geophys. Res. Lett. 37, L19703).
// What can be read is the landscape. The density p of the window's values is estimated by a Gaussian kernel (Silverman's
// bandwidth), and −ln p, which is 2U/σ² for a gradient system with additive noise, is fitted by polynomials of even
// degree 2 to 8 by weighted least squares (the weights p, so that the sparse tails count less); the degree is chosen by
// Schwarz's criterion among those whose leading coefficient is positive, and the wells counted as the fitted
// polynomial's minima within the data (after Livina, Kwasniok and Lenton 2010, Clim. Past 6, 77–82, whose method fits the
// potential's degree to the record). With two wells the record's switches across the barrier are counted with a
// hysteresis of a quarter of a standard deviation, and Kramers's mean time to leave the deeper well, in the fitted
// landscape Φ = −ln p and with σ² from the increments, is τ = 2π / ((σ²/2) √(Φ''(a) |Φ''(b)|)) e^{ΔΦ}.
//
// It reports proximity to switching, never the approach to a bifurcation: a record may switch with no parameter moving.
// The bimodality coefficient (Pfister et al. 2013, Front. Psychol. 4, 700) is given beside it, above 5/9 for a
// bimodal sample.
#pragma once

#include <kinesis/transitions/space.hpp>

#include <tyche/descriptive/summary.hpp>

#include <numbers>
#include <span>

namespace kinesis::transitions {

struct Wells {
  std::size_t wells = 1;
  std::size_t degree = 2;                // the fitted potential's degree
  std::vector<double> minima, maxima;    // in the window's values
  std::vector<double> grid, landscape;   // −ln p on a grid, and the fitted polynomial on it
  std::vector<double> fitted;
  std::size_t switches = 0;              // crossings of the barrier, with hysteresis
  double bimodality = 0;
  double escape_time = std::numeric_limits<double>::infinity();  // Kramers's time to leave the deeper well
  std::string sentence;
};

namespace detail {

// Weighted least squares of y on 1, z, …, z^k; the coefficients, lowest first, and the weighted residual sum.
inline std::pair<std::vector<double>, double> weighted_polynomial(const std::vector<double>& z, const std::vector<double>& y,
                                                                 const std::vector<double>& w, std::size_t k) {
  const std::size_t m = k + 1;
  Matrix A(m), B(m);
  for (std::size_t i = 0; i < z.size(); i++) {
    std::vector<double> p(m, 1);
    for (std::size_t j = 1; j < m; j++) p[j] = p[j - 1] * z[i];
    for (std::size_t a = 0; a < m; a++) {
      for (std::size_t b = 0; b < m; b++) A(a, b) += w[i] * p[a] * p[b];
      B(a, 0) += w[i] * p[a] * y[i];
    }
  }
  for (std::size_t a = 0; a < m; a++) A(a, a) += 1e-12 * (A(a, a) + 1);
  const Matrix X = calx::linalg::solve(A, B);
  std::vector<double> c(m);
  for (std::size_t a = 0; a < m; a++) c[a] = X(a, 0);
  double rss = 0;
  for (std::size_t i = 0; i < z.size(); i++) {
    double v = 0, p = 1;
    for (std::size_t a = 0; a < m; a++) v += c[a] * p, p *= z[i];
    rss += w[i] * (y[i] - v) * (y[i] - v);
  }
  return {c, rss};
}

inline double polynomial(const std::vector<double>& c, double z) {
  double v = 0;
  for (std::size_t a = c.size(); a-- > 0;) v = v * z + c[a];
  return v;
}

}  // namespace detail

inline Wells wells(std::span<const double> x, double dt) {
  const std::size_t n = x.size();
  if (n < 50) throw Refusal("counting wells needs a window of at least 50 states");
  const double mean = tyche::descriptive::mean(x), sd = std::sqrt(tyche::descriptive::variance(x, 1));
  if (!(sd > 0)) throw Refusal("a window that does not vary visits one state");
  std::vector<double> z(n);
  for (std::size_t i = 0; i < n; i++) z[i] = (x[i] - mean) / sd;
  const double h = 1.06 * std::pow(static_cast<double>(n), -0.2);
  const double lo = *std::min_element(z.begin(), z.end()), hi = *std::max_element(z.begin(), z.end());
  Wells out;
  const std::size_t G = 120;
  std::vector<double> g(G), y, w, gz;
  for (std::size_t i = 0; i < G; i++) {
    g[i] = lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(G - 1);
    double f = 0;
    for (double v : z) f += std::exp(-0.5 * (g[i] - v) * (g[i] - v) / (h * h));
    f /= static_cast<double>(n) * h * std::sqrt(2 * std::numbers::pi);
    if (f > 1e-12) gz.push_back(g[i]), y.push_back(-std::log(f)), w.push_back(f);
  }
  double best = std::numeric_limits<double>::infinity();
  std::vector<double> coefficients;
  const double N = static_cast<double>(n);
  double wsum = 0;
  for (double v : w) wsum += v;
  for (std::size_t k : {2u, 4u, 6u, 8u}) {
    auto [c, rss] = detail::weighted_polynomial(gz, y, w, k);
    if (!(c.back() > 0)) continue;
    const double bic = N * std::log(std::max(1e-300, rss / wsum)) + static_cast<double>(k + 1) * std::log(N);
    if (bic < best - 1e-9) best = bic, coefficients = c, out.degree = k;
  }
  if (coefficients.empty()) coefficients = detail::weighted_polynomial(gz, y, w, 2).first, out.degree = 2;
  for (std::size_t i = 0; i < gz.size(); i++) {
    out.grid.push_back(mean + sd * gz[i]);
    out.landscape.push_back(y[i]);
    out.fitted.push_back(detail::polynomial(coefficients, gz[i]));
  }
  // Minima and maxima of the fitted landscape strictly inside the data.
  std::vector<double> zmin, zmax;
  for (std::size_t i = 1; i + 1 < out.fitted.size(); i++) {
    if (out.fitted[i] < out.fitted[i - 1] && out.fitted[i] <= out.fitted[i + 1]) zmin.push_back(gz[i]);
    if (out.fitted[i] > out.fitted[i - 1] && out.fitted[i] >= out.fitted[i + 1]) zmax.push_back(gz[i]);
  }
  out.wells = std::max<std::size_t>(1, zmin.size());
  for (double v : zmin) out.minima.push_back(mean + sd * v);
  for (double v : zmax) out.maxima.push_back(mean + sd * v);

  const double s = tyche::descriptive::skewness(z, true), kurt = tyche::descriptive::kurtosis(z, true, true);
  out.bimodality = (s * s + 1) / (kurt + 3 * (N - 1) * (N - 1) / ((N - 2) * (N - 3)));

  if (zmin.size() >= 2 && !zmax.empty()) {
    const double barrier = zmax.front();
    int side = z[0] < barrier ? -1 : 1;
    for (double v : z) {
      if (side < 0 && v > barrier + 0.25) side = 1, out.switches++;
      else if (side > 0 && v < barrier - 0.25) side = -1, out.switches++;
    }
    // Kramers's time from the deeper well over the barrier, in the fitted landscape (in standardised units).
    auto second = [&](double at) {
      const double e = 1e-3;
      return (detail::polynomial(coefficients, at + e) - 2 * detail::polynomial(coefficients, at) + detail::polynomial(coefficients, at - e)) / (e * e);
    };
    const double a = detail::polynomial(coefficients, zmin.front()) <= detail::polynomial(coefficients, zmin.back()) ? zmin.front() : zmin.back();
    const double dphi = detail::polynomial(coefficients, barrier) - detail::polynomial(coefficients, a);
    double vi = 0;
    for (std::size_t i = 0; i + 1 < n; i++) vi += (z[i + 1] - z[i]) * (z[i + 1] - z[i]);
    const double sigma2 = vi / static_cast<double>(n - 1) / dt;
    const double ca = second(a), cb = std::abs(second(barrier));
    if (ca > 0 && cb > 0 && sigma2 > 0) out.escape_time = 2 * std::numbers::pi / ((sigma2 / 2) * std::sqrt(ca * cb)) * std::exp(dphi);
  }
  if (out.wells >= 2)
    out.sentence = "the window's values fall in " + std::to_string(out.wells) + " wells, crossed " + std::to_string(out.switches) +
                   " times: the record is near switching between states, which says nothing of whether a parameter is moving";
  else
    out.sentence = "the window's values fall in one well";
  return out;
}

}  // namespace kinesis::transitions
