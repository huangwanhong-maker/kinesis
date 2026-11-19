// Which kind of transition a slowing record approaches: where its dominant eigenvalue goes.
//
// At a local bifurcation of codimension one an eigenvalue of the linearisation crosses the stability boundary, and where
// it crosses tells the kind (the survey, Part I §2, after Kuehn 2013, J. Nonlinear Sci. 23, 457–510): sampled at a step
// Δt, the dominant eigenvalue μ of the discrete propagator approaches +1 at a fold, a transcritical or a pitchfork
// bifurcation (a stationary loss of stability), −1 at a period-doubling of a map, and the unit circle away from the real
// axis at a Hopf or Neimark–Sacker bifurcation (an oscillatory one). The dominant eigenvalue of a record ("DEV",
// Grziwotz et al. 2023, Sci. Adv. 9, eabq4558) is read here from a linear autoregression in each window: for one
// coordinate an AR(p) of order 1 or 2 chosen by Akaike's criterion (an AR(2) with complex roots is a spectrum with a peak
// away from zero frequency, the distinction Bury, Bauch and Anand 2020, J. R. Soc. Interface 17, 20200482, draw between
// fold and Hopf spectra), for several the lag-1 propagator.
//
// Lag-1 autocorrelation sees only the first case: it tends to −1 before a period-doubling and may fall before a Hopf
// sampled coarsely (ω₀Δt > π/2), so a test of "autocorrelation rising" misses both, while the modulus |μ| rises in all three.
//
// The shape of the rate's fall is a weaker hint, valid only when the parameter moves linearly in time: at a fold the rate
// is κ = 2√(μ_c − μ), so κ² falls linearly; at a transcritical or pitchfork bifurcation κ = |μ − μ_c| falls linearly.
#pragma once

#include <kinesis/transitions/linear.hpp>

#include <tyche/descriptive/rank.hpp>

#include <numbers>
#include <span>

namespace kinesis::transitions {

struct Dominant {
  Complex mu;
  std::size_t order = 1;   // the autoregression's order (1 or 2) or the coordinates of a propagator
};

// The dominant eigenvalue of an AR(1) or AR(2) fitted to a centred window by least squares, the order chosen by AIC.
inline Dominant ar_dominant(std::span<const double> x) {
  const std::size_t n = x.size();
  if (n < 12) throw Refusal("a dominant eigenvalue needs a window of at least 12 states");
  const double m = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(n);
  std::vector<double> c(n);
  for (std::size_t t = 0; t < n; t++) c[t] = x[t] - m;
  // AR(1) on t = 2 … n − 1 so that both orders use the same observations.
  double s11 = 0, s1y = 0, s22 = 0, s12 = 0, s2y = 0, syy = 0;
  for (std::size_t t = 2; t < n; t++) {
    const double y = c[t], a = c[t - 1], b = c[t - 2];
    s11 += a * a, s1y += a * y, s22 += b * b, s12 += a * b, s2y += b * y, syy += y * y;
  }
  const double N = static_cast<double>(n - 2);
  const double phi = s11 > 0 ? s1y / s11 : 0;
  const double rss1 = std::max(1e-300, syy - phi * s1y);
  const double det = s11 * s22 - s12 * s12;
  Dominant d;
  if (det > 0) {
    const double a1 = (s1y * s22 - s2y * s12) / det, a2 = (s2y * s11 - s1y * s12) / det;
    const double rss2 = std::max(1e-300, syy - a1 * s1y - a2 * s2y);
    const double aic1 = N * std::log(rss1 / N) + 2, aic2 = N * std::log(rss2 / N) + 4;
    if (aic2 < aic1) {
      // Roots of z² − a₁ z − a₂ = 0, the companion's eigenvalues.
      const Complex disc = std::sqrt(Complex(a1 * a1 + 4 * a2, 0));
      const Complex z1 = (a1 + disc) / 2.0, z2 = (a1 - disc) / 2.0;
      d.mu = std::abs(z1) >= std::abs(z2) ? z1 : z2;
      if (d.mu.imag() < 0) d.mu = std::conj(d.mu);
      d.order = 2;
      return d;
    }
  }
  d.mu = Complex(phi, 0);
  d.order = 1;
  return d;
}

enum class Approach { stationary, flip, oscillatory, undetermined };

inline const char* name(Approach a) {
  switch (a) {
    case Approach::stationary: return "towards +1: a stationary loss of stability (fold, transcritical or pitchfork)";
    case Approach::flip: return "towards −1: a period-doubling (flip)";
    case Approach::oscillatory: return "towards the unit circle off the real axis: an oscillatory loss of stability (Hopf or Neimark–Sacker)";
    default: return "undetermined";
  }
}

struct KindReading {
  Approach approach = Approach::undetermined;
  std::vector<double> at;
  std::vector<Complex> path;        // the dominant eigenvalue in each window
  double tau_modulus = 0;           // Kendall's τ of |μ| against time
  double change_modulus = 0;        // |μ| in the last window less |μ| in the first
  double frequency = 0;             // the dominant mode's frequency at the end, arg μ / Δt
  std::string shape;                // the rate's shape, a weak hint
  std::string sentence;
};

// The kind a path of dominant eigenvalues points to, read over its last third.
inline KindReading classify(const std::vector<double>& at, const std::vector<Complex>& path, double dt) {
  KindReading k;
  k.at = at;
  k.path = path;
  if (path.size() < 6) {
    k.sentence = "too few windows to read where the dominant eigenvalue goes";
    return k;
  }
  std::vector<double> modulus(path.size());
  for (std::size_t i = 0; i < path.size(); i++) modulus[i] = std::abs(path[i]);
  k.tau_modulus = tyche::descriptive::kendall_tau_b(at, modulus);
  k.change_modulus = modulus.back() - modulus.front();
  const std::size_t from = path.size() - std::max<std::size_t>(3, path.size() / 3);
  std::size_t negative = 0, complex_count = 0;
  double argument = 0;
  for (std::size_t i = from; i < path.size(); i++) {
    const double a = std::abs(std::arg(path[i]));
    argument += a;
    if (a > 0.9 * std::numbers::pi) negative++;
    else if (a > std::numbers::pi / 12) complex_count++;
  }
  const double count = static_cast<double>(path.size() - from);
  argument /= count;
  k.frequency = argument / dt;
  if (static_cast<double>(negative) >= 0.6 * count) k.approach = Approach::flip;
  else if (static_cast<double>(complex_count) >= 0.6 * count) k.approach = Approach::oscillatory;
  else k.approach = Approach::stationary;

  if (k.approach == Approach::stationary) {
    // The rate's shape: R² of κ and of κ² against time, over windows with a positive rate.
    std::vector<double> t, kappa, kappa2;
    for (std::size_t i = 0; i < path.size(); i++)
      if (modulus[i] > 0 && modulus[i] < 1) {
        t.push_back(at[i]);
        const double r = -std::log(modulus[i]) / dt;
        kappa.push_back(r);
        kappa2.push_back(r * r);
      }
    auto r2 = [&](const std::vector<double>& y) {
      const std::size_t n = y.size();
      if (n < 4) return 0.0;
      const double mt = std::accumulate(t.begin(), t.end(), 0.0) / n, my = std::accumulate(y.begin(), y.end(), 0.0) / n;
      double stt = 0, sty = 0, syy = 0;
      for (std::size_t i = 0; i < n; i++) stt += (t[i] - mt) * (t[i] - mt), sty += (t[i] - mt) * (y[i] - my), syy += (y[i] - my) * (y[i] - my);
      return stt > 0 && syy > 0 ? sty * sty / (stt * syy) : 0.0;
    };
    const double a = r2(kappa), b = r2(kappa2);
    if (std::abs(a - b) > 0.02 && std::max(a, b) > 0.3)
      k.shape = b > a ? "the rate's square falls more nearly linearly than the rate, as at a fold (weak: it assumes a parameter moving linearly)"
                      : "the rate falls more nearly linearly than its square, as at a transcritical or pitchfork bifurcation (weak: it assumes a parameter moving linearly)";
  }
  k.sentence = std::string("the dominant eigenvalue goes ") + name(k.approach) +
               (k.change_modulus > 0 ? ", its modulus rising from " : ", its modulus not rising, from ") + std::to_string(modulus.front()).substr(0, 5) +
               " to " + std::to_string(modulus.back()).substr(0, 5);
  return k;
}

}  // namespace kinesis::transitions
