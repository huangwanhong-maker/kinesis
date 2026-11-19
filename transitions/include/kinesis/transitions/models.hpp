// Records with known answers: the positive cases a detector must see and the negative cases on which it must stay silent
// (the survey, Part I §9.4, Part II §11.3).
//
// Positive: the normal forms of the fold (ẋ = μ − x²), the transcritical (ẋ = μx − x²) and the pitchfork (ẋ = μx − x³)
// bifurcations and of the Hopf bifurcation (ż = (μ + iω)z − |z|²z) driven slowly towards their critical points; a linear
// map approaching a period-doubling (x_{n+1} = −(1 − δ)x_n); an Ornstein–Uhlenbeck process whose rate falls; a network
// Ornstein–Uhlenbeck process dx = −(Γ + κL)x dt + σ dW whose rates Γ fall.
// Negative: a stationary Ornstein–Uhlenbeck process; rate-induced tipping in ẋ = (x + rt)² − 1, where in the moving
// coordinate y = x + rt the stable state y = −√(1 − r) has the constant rate 2√(1 − r) for constant r < 1 and none for
// r > 1, so the system tips at r > 1 with nothing slowing (the survey, Part I §9.3); noise-induced switching in the double
// well ẋ = x − x³ at fixed parameters; and an OU process of constant rate driven by noise whose correlation grows while its
// variance is held, which raises variance and autocorrelation with no slowing at all (Morr and Boers 2024's spurious case).
//
// Stochastic differentials are stepped by Euler and Maruyama with substeps; linear ones exactly.
#pragma once

#include <kinesis/transitions/linear.hpp>

#include <calx/random/stream.hpp>

#include <functional>

namespace kinesis::transitions::models {

// Times 0, Δt, …, (T − 1)Δt.
inline std::vector<double> times(std::size_t T, double dt) {
  std::vector<double> t(T);
  for (std::size_t i = 0; i < T; i++) t[i] = static_cast<double>(i) * dt;
  return t;
}

// An OU process dx = −κ(t) x dt + σ dW at the times given, exact for κ constant over each interval.
inline std::vector<double> ou(const std::vector<double>& t, const std::function<double(double)>& kappa, double sigma, calx::random::Stream& rng) {
  std::vector<double> x(t.size());
  const double k0 = kappa(t[0]);
  x[0] = sigma / std::sqrt(2 * k0) * rng.normal();
  for (std::size_t i = 0; i + 1 < t.size(); i++) {
    const double k = kappa(t[i]), h = t[i + 1] - t[i], e = std::exp(-k * h);
    x[i + 1] = x[i] * e + sigma * std::sqrt((1 - e * e) / (2 * k)) * rng.normal();
  }
  return x;
}

// A one-dimensional SDE dx = f(x, t) dt + σ dW from x0, by Euler and Maruyama with `sub` substeps a sample.
inline std::vector<double> sde(std::size_t T, double dt, double x0, const std::function<double(double, double)>& f, double sigma,
                               calx::random::Stream& rng, int sub = 20) {
  std::vector<double> x(T);
  double v = x0, t = 0;
  const double h = dt / sub, s = sigma * std::sqrt(h);
  for (std::size_t i = 0; i < T; i++) {
    x[i] = v;
    for (int k = 0; k < sub; k++) {
      v += f(v, t) * h + s * rng.normal();
      t += h;
    }
  }
  return x;
}

// μ moving linearly from μ₀ to μ₁ over T samples.
inline std::function<double(double)> ramp(double mu0, double mu1, std::size_t T, double dt) {
  const double duration = static_cast<double>(T - 1) * dt;
  return [=](double t) { return mu0 + (mu1 - mu0) * std::min(1.0, t / duration); };
}

inline std::vector<double> fold(std::size_t T, double dt, double mu0, double mu1, double sigma, calx::random::Stream& rng) {
  const auto mu = ramp(mu0, mu1, T, dt);
  return sde(T, dt, std::sqrt(mu0), [mu](double x, double t) { return mu(t) - x * x; }, sigma, rng);
}

inline std::vector<double> transcritical(std::size_t T, double dt, double mu0, double mu1, double sigma, calx::random::Stream& rng) {
  const auto mu = ramp(mu0, mu1, T, dt);
  return sde(T, dt, 0, [mu](double x, double t) { return mu(t) * x - x * x; }, sigma, rng);
}

inline std::vector<double> pitchfork(std::size_t T, double dt, double mu0, double mu1, double sigma, calx::random::Stream& rng) {
  const auto mu = ramp(mu0, mu1, T, dt);
  return sde(T, dt, 0, [mu](double x, double t) { return mu(t) * x - x * x * x; }, sigma, rng);
}

inline std::vector<double> double_well(std::size_t T, double dt, double sigma, calx::random::Stream& rng) {
  return sde(T, dt, -1, [](double x, double) { return x - x * x * x; }, sigma, rng);
}

// ẋ = (x + rt)² − 1 from the stable state x = −1 at t = 0 (r < 1 tracks y = −√(1 − r) after a transient).
inline std::vector<double> rate_tipping(std::size_t T, double dt, double r, double sigma, calx::random::Stream& rng) {
  const double y0 = r < 1 ? -std::sqrt(1 - r) : -1.0;
  return sde(T, dt, y0, [r](double x, double t) {
    const double y = x + r * t;
    return std::max(-50.0, std::min(50.0, y * y - 1));
  }, sigma, rng);
}

// The Hopf normal form ż = (μ + iω) z − |z|² z + σ(dW₁ + i dW₂), as two coordinates (x, y) per sample.
inline std::vector<std::vector<double>> hopf(std::size_t T, double dt, double mu0, double mu1, double omega, double sigma,
                                            calx::random::Stream& rng, int sub = 20) {
  const auto mu = ramp(mu0, mu1, T, dt);
  std::vector<std::vector<double>> out(T, std::vector<double>(2));
  double x = 0, y = 0, t = 0;
  const double h = dt / sub, s = sigma * std::sqrt(h);
  for (std::size_t i = 0; i < T; i++) {
    out[i] = {x, y};
    for (int k = 0; k < sub; k++) {
      const double m = mu(t), r2 = x * x + y * y;
      const double fx = m * x - omega * y - r2 * x, fy = omega * x + m * y - r2 * y;
      x += fx * h + s * rng.normal();
      y += fy * h + s * rng.normal();
      t += h;
    }
  }
  return out;
}

// x_{n+1} = −(1 − δ_n) x_n + σ ε_n, δ moving linearly from δ₀ to δ₁.
inline std::vector<double> flip(std::size_t T, double delta0, double delta1, double sigma, calx::random::Stream& rng) {
  std::vector<double> x(T);
  double v = 0;
  for (std::size_t i = 0; i < T; i++) {
    x[i] = v;
    const double d = delta0 + (delta1 - delta0) * static_cast<double>(i) / static_cast<double>(T - 1);
    v = -(1 - d) * v + sigma * rng.normal();
  }
  return x;
}

// dx = −κ x dt + η dt, dη = −θ(t) η dt + √(2θ(t)) s dW: the noise η keeps its variance s² while its correlation time
// 1/θ grows from 1/θ₀ to 1/θ₁ over the record; κ does not change.
inline std::vector<double> reddening(std::size_t T, double dt, double kappa, double theta0, double theta1, double s,
                                     calx::random::Stream& rng, int sub = 20) {
  std::vector<double> x(T);
  double v = 0, eta = s * rng.normal(), t = 0;
  const double duration = static_cast<double>(T - 1) * dt, h = dt / sub;
  for (std::size_t i = 0; i < T; i++) {
    x[i] = v;
    for (int k = 0; k < sub; k++) {
      const double th = theta0 + (theta1 - theta0) * std::min(1.0, t / duration);
      v += (-kappa * v + eta) * h;
      eta += -th * eta * h + std::sqrt(2 * th * h) * s * rng.normal();
      t += h;
    }
  }
  return x;
}

// The network OU process dx = −(Γ(t) + κL) x dt + σ dW on n units, L the graph Laplacian of the edges, Γ(t) diagonal,
// stepped exactly at Δt (M = Γ + κL symmetric: x ↦ e^{−MΔt}x plus noise of covariance V diag(σ²(1 − e^{−2λΔt})/(2λ)) Vᵀ).
inline std::vector<std::vector<double>> network_ou(std::size_t T, double dt, std::size_t n,
                                                  const std::vector<std::pair<std::size_t, std::size_t>>& edges, double coupling,
                                                  const std::function<std::vector<double>(double)>& gamma, double sigma,
                                                  calx::random::Stream& rng) {
  Matrix L(n);
  for (const auto& [a, b] : edges) {
    L(a, a) += 1, L(b, b) += 1, L(a, b) -= 1, L(b, a) -= 1;
  }
  std::vector<std::vector<double>> out(T, std::vector<double>(n, 0));
  std::vector<double> x(n, 0);
  for (std::size_t i = 0; i < T + 200; i++) {
    const double t = i < 200 ? 0.0 : static_cast<double>(i - 200) * dt;
    if (i >= 200) out[i - 200] = x;
    Matrix M(n);
    const auto g = gamma(t);
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t b = 0; b < n; b++) M(a, b) = coupling * L(a, b) + (a == b ? g[a] : 0);
    const auto e = calx::linalg::symmetric_eigen(M);
    std::vector<double> y(n, 0), z(n);
    for (std::size_t k = 0; k < n; k++) {
      double proj = 0;
      for (std::size_t a = 0; a < n; a++) proj += e.vectors(a, k) * x[a];
      const double lam = e.values[k], decay = std::exp(-lam * dt);
      const double var = lam > 1e-12 ? sigma * sigma * (1 - decay * decay) / (2 * lam) : sigma * sigma * dt;
      z[k] = proj * decay + std::sqrt(var) * rng.normal();
    }
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t k = 0; k < n; k++) y[a] += e.vectors(a, k) * z[k];
    x = y;
  }
  return out;
}

// A VAR(1) record x_{t+1} = A(s) x_t + B ε with s = t/(T − 1), started after a burn-in at A(0).
inline std::vector<std::vector<double>> var1(std::size_t T, const std::function<Matrix(double)>& A, const Matrix& B, calx::random::Stream& rng) {
  const std::size_t n = B.n;
  std::vector<std::vector<double>> out(T);
  std::vector<double> x(n, 0);
  auto step = [&](const Matrix& a) {
    std::vector<double> e(n), next(n, 0);
    for (auto& v : e) v = rng.normal();
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) next[i] += a(i, j) * x[j] + B(i, j) * e[j];
    x = next;
  };
  const Matrix a0 = A(0);
  for (int i = 0; i < 300; i++) step(a0);
  for (std::size_t t = 0; t < T; t++) {
    out[t] = x;
    step(A(T > 1 ? static_cast<double>(t) / static_cast<double>(T - 1) : 0));
  }
  return out;
}

}  // namespace kinesis::transitions::models
