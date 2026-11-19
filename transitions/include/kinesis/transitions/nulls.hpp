// Records that have everything of the observed one except what a test is about: the nulls against which a trend in an
// indicator is judged, and the declared alternative against which the test's power is measured.
//
// Each null answers a different question (the survey, Part I §4):
//   red (fitted)  a stationary linear process fitted to the record's residuals: AR(1) for one coordinate (Dakos et al. 2012,
//                 PLoS ONE 7, e41010), VAR(1) for several when the record is long enough against the dimension, otherwise
//                 independent AR(1) coordinates with the innovations' correlation kept. It is the parametric null, and the
//                 one from which the operating characteristics are computed.
//   phase         the same power spectrum with random phases (Theiler et al. 1992, Physica D 58, 77–94): any stationary
//                 linear Gaussian process with the record's correlation; for several coordinates the same random phase at
//                 each frequency, keeping the cross-spectrum (Prichard and Theiler 1994, Phys. Rev. Lett. 73, 951–954).
//   iaaft         the same spectrum and the same distribution of values, by iterated amplitude adjustment (Schreiber and
//                 Schmitz 1996, Phys. Rev. Lett. 77, 635–638): a static nonlinear transformation of such a process.
//   block         the stationary bootstrap (Politis and Romano 1994, J. Am. Stat. Assoc. 89, 1303–1313): blocks of random
//                 geometric length, keeping dependence within a block and assuming none of the linear forms above.
//
// The declared alternative is the fitted null with its leading mode slowed: its rate falls linearly over the record to a
// declared share of its initial value (a tenth by default) while the noise is unchanged, the slowing down a bifurcation
// would bring. A test that seldom rejects under it at the record's length cannot report that there is no signal.
#pragma once

#include <kinesis/transitions/linear.hpp>

#include <calx/random/stream.hpp>
#include <calx/transforms/fft.hpp>
#include <tyche/time_series/spectral.hpp>

namespace kinesis::transitions {

enum class Null { red, phase, iaaft, block };

inline const char* name(Null n) {
  switch (n) {
    case Null::red: return "fitted red noise";
    case Null::phase: return "random phases";
    case Null::iaaft: return "iterated amplitude-adjusted phases";
    default: return "stationary block bootstrap";
  }
}

inline const char* description(Null n) {
  switch (n) {
    case Null::red: return "a stationary linear process fitted to the residuals (AR(1), or VAR(1) for several coordinates)";
    case Null::phase: return "any stationary linear Gaussian process with the residuals' spectrum";
    case Null::iaaft: return "a static transformation of a stationary linear Gaussian process with the residuals' spectrum and values";
    default: return "a stationary process whose dependence is shorter than the bootstrap's blocks";
  }
}

// ---- the fitted linear null and its alternative -----------------------------------------------------------------

struct LinearNull {
  Matrix a;          // the propagator (diagonal when independent)
  Matrix noise_root; // a square root of the innovations' covariance
  bool full = true;  // VAR(1) when true, independent AR(1) coordinates otherwise
  std::string described;
};

inline LinearNull fit_null(const Rectangular& residual) {
  const std::size_t T = residual.rows, d = residual.cols;
  LinearNull out;
  Matrix Q;
  if (d == 1 || T >= 10 * d + 10) {
    const Propagator P = propagator(residual, d == 1 ? 0.0 : 1e-4);
    out.a = P.a;
    Q = P.noise;
    out.full = true;
    out.described = d == 1 ? "AR(1) fitted to the residuals" : "VAR(1) fitted to the residuals";
  } else {
    out.a = Matrix(d);
    std::vector<std::vector<double>> innovations(d);
    for (std::size_t j = 0; j < d; j++) {
      const auto x = column(residual, j);
      const auto r = tyche::time_series::red_noise(x);
      const double phi = std::clamp(r.rho, -0.999, 0.999);
      out.a(j, j) = phi;
      for (std::size_t t = 0; t + 1 < T; t++) innovations[j].push_back(x[t + 1] - phi * x[t]);
    }
    const Covariance c = covariance(Rectangular::from_columns(innovations), true);
    Q = c.c;
    out.full = false;
    out.described = "independent AR(1) coordinates fitted to the residuals, their innovations' correlation kept (the record is short against its dimension)";
  }
  // Ensure stability of the fitted null.
  double rho = 0;
  for (const auto& z : calx::linalg::eigenvalues(out.a)) rho = std::max(rho, std::abs(z));
  if (rho >= 0.999) {
    for (double& v : out.a.a) v *= 0.999 / rho;
    out.described += " (scaled to be stable)";
  }
  out.noise_root = detail::symmetric_function(Q, [](double v) { return std::sqrt(std::max(0.0, v)); });
  return out;
}

namespace detail {

inline std::vector<double> draw_noise(const LinearNull& m, calx::random::Stream& rng) {
  const std::size_t d = m.a.n;
  std::vector<double> z(d), e(d, 0);
  for (auto& v : z) v = rng.normal();
  for (std::size_t i = 0; i < d; i++)
    for (std::size_t j = 0; j < d; j++) e[i] += m.noise_root(i, j) * z[j];
  return e;
}

// The propagator with its leading mode's modulus set so that its rate is `share` of the rate it has in A.
inline Matrix slowed(const Matrix& A, double share) {
  const auto spectrum = calx::linalg::eigenvalues(A);
  std::size_t lead = 0;
  for (std::size_t k = 1; k < spectrum.size(); k++)
    if (std::abs(spectrum[k]) > std::abs(spectrum[lead])) lead = k;
  const Complex mu = spectrum[lead];
  const double modulus = std::abs(mu);
  if (!(modulus > 0)) return A;
  const double target = std::pow(modulus, share);  // e^{−share·κΔt}
  const Complex mu_new = std::polar(target, std::arg(mu));
  if (A.n == 1) {
    Matrix out(1);
    out(0, 0) = mu_new.real();
    return out;
  }
  const auto v = eigenvector(A, mu, false);
  auto w = eigenvector(A, mu, true);
  for (auto& z : w) z = std::conj(z);
  Complex wv = 0;
  for (std::size_t i = 0; i < A.n; i++) wv += std::conj(w[i]) * v[i];
  Matrix out = A;
  const bool real = std::abs(mu.imag()) <= 1e-12 * std::max(1.0, modulus);
  const Complex delta = mu_new - mu;
  for (std::size_t i = 0; i < A.n; i++)
    for (std::size_t j = 0; j < A.n; j++) {
      const Complex term = delta * v[i] * std::conj(w[j]) / wv;
      out(i, j) += real ? term.real() : 2 * term.real();
    }
  return out;
}

}  // namespace detail

// A record of the fitted null, or of its alternative when `share` < 1 (the leading rate falling linearly to `share` of
// its value over the record), started from the stationary state of its first propagator.
inline Rectangular simulate_null(const LinearNull& m, std::size_t T, calx::random::Stream& rng, double share = 1) {
  const std::size_t d = m.a.n;
  Rectangular out(T, d);
  const Matrix end = share < 1 ? detail::slowed(m.a, share) : m.a;
  std::vector<double> x(d, 0);
  for (int burn = 0; burn < 300; burn++) {
    const auto e = detail::draw_noise(m, rng);
    std::vector<double> next(d, 0);
    for (std::size_t i = 0; i < d; i++) {
      for (std::size_t j = 0; j < d; j++) next[i] += m.a(i, j) * x[j];
      next[i] += e[i];
    }
    x = next;
  }
  for (std::size_t t = 0; t < T; t++) {
    for (std::size_t i = 0; i < d; i++) out(t, i) = x[i];
    const double s = T > 1 ? static_cast<double>(t) / static_cast<double>(T - 1) : 0;
    const auto e = detail::draw_noise(m, rng);
    std::vector<double> next(d, 0);
    for (std::size_t i = 0; i < d; i++) {
      for (std::size_t j = 0; j < d; j++) next[i] += ((1 - s) * m.a(i, j) + s * end(i, j)) * x[j];
      next[i] += e[i];
    }
    x = next;
  }
  return out;
}

// ---- surrogates of an observed record ---------------------------------------------------------------------------

namespace detail {

inline std::vector<double> iaaft(const std::vector<double>& x, calx::random::Stream& rng, int iterations = 100) {
  const std::size_t n = x.size();
  std::vector<double> sorted = x;
  std::sort(sorted.begin(), sorted.end());
  std::vector<calx::transforms::Complex> X(x.begin(), x.end());
  calx::transforms::fft(X);
  std::vector<double> amplitude(n);
  for (std::size_t k = 0; k < n; k++) amplitude[k] = std::abs(X[k]);
  std::vector<double> r = x;
  for (std::size_t i = n; i > 1; i--) std::swap(r[i - 1], r[rng.below(i)]);
  std::vector<std::size_t> order(n), previous;
  for (int it = 0; it < iterations; it++) {
    std::vector<calx::transforms::Complex> R(r.begin(), r.end());
    calx::transforms::fft(R);
    for (std::size_t k = 0; k < n; k++) {
      const double a = std::abs(R[k]);
      R[k] = a > 0 ? R[k] * (amplitude[k] / a) : calx::transforms::Complex(amplitude[k], 0);
    }
    calx::transforms::fft(R, true);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return R[a].real() < R[b].real(); });
    for (std::size_t k = 0; k < n; k++) r[order[k]] = sorted[k];
    if (order == previous) break;
    previous = order;
  }
  return r;
}

inline std::vector<std::size_t> stationary_bootstrap_indices(std::size_t n, double mean_block, calx::random::Stream& rng) {
  std::vector<std::size_t> idx(n);
  std::size_t at = rng.below(n);
  for (std::size_t t = 0; t < n; t++) {
    if (t > 0) at = rng.uniform() < 1 / mean_block ? rng.below(n) : (at + 1) % n;
    idx[t] = at;
  }
  return idx;
}

}  // namespace detail

// A surrogate of a record of residual coordinates under a null. Phase surrogates of several coordinates share their random
// phases; IAAFT is applied coordinate by coordinate (keeping each one's spectrum and values, not the cross-spectrum exactly).
inline Rectangular surrogate(const Rectangular& residual, Null null, calx::random::Stream& rng, const LinearNull* fitted = nullptr) {
  const std::size_t T = residual.rows, d = residual.cols;
  switch (null) {
    case Null::red: {
      if (!fitted) throw Refusal("a fitted red-noise surrogate needs its fitted null");
      return simulate_null(*fitted, T, rng);
    }
    case Null::phase: {
      std::vector<double> phases(T, 0);
      for (std::size_t j = 1; j + j < T; j++) phases[j] = rng.uniform(0, 2 * std::numbers::pi);
      Rectangular out(T, d);
      for (std::size_t c = 0; c < d; c++) {
        const auto x = column(residual, c);
        std::vector<calx::transforms::Complex> X(x.begin(), x.end());
        calx::transforms::fft(X);
        for (std::size_t j = 1; j + j < T; j++) {
          X[j] *= std::polar(1.0, phases[j]);
          X[T - j] = std::conj(X[j]);
        }
        calx::transforms::fft(X, true);
        for (std::size_t t = 0; t < T; t++) out(t, c) = X[t].real();
      }
      return out;
    }
    case Null::iaaft: {
      Rectangular out(T, d);
      for (std::size_t c = 0; c < d; c++) {
        const auto s = detail::iaaft(column(residual, c), rng);
        for (std::size_t t = 0; t < T; t++) out(t, c) = s[t];
      }
      return out;
    }
    default: {
      const double block = std::max(5.0, std::cbrt(static_cast<double>(T)) * 2);
      const auto idx = detail::stationary_bootstrap_indices(T, block, rng);
      Rectangular out(T, d);
      for (std::size_t t = 0; t < T; t++)
        for (std::size_t c = 0; c < d; c++) out(t, c) = residual(idx[t], c);
      return out;
    }
  }
}

}  // namespace kinesis::transitions
