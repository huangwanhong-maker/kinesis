// The linear core of every detector on a state of many coordinates.
//
// Near an equilibrium a stochastic system is, to first order, x_{t+1} = A x_t + ε_t with Cov ε = Q, sampled at a step Δt
// from dx = J x dt + B dW with A = e^{JΔt}. Its stationary covariance solves C = A C Aᵀ + Q (the discrete Lyapunov
// equation), and every indicator of slowing down estimates one of three objects (the survey, Part II §1): the spectrum of
// the lag-1 propagator Â = Ĉ(1) Ĉ(0)⁻¹, the leading eigenpair of the covariance Ĉ, or a projection on an estimate of the
// critical direction. As an eigenvalue λ₁ of J approaches zero, μ₁ = e^{λ₁Δt} approaches one, and near a simple fold
// C ≈ (wᵀBBᵀw / 2|λ₁|) v vᵀ, where v and w are J's right and left eigenvectors of λ₁: w says whether the noise excites
// the critical mode, v where it shows. A record whose noise does not excite the mode (wᵀB = 0) shows no rise in variance
// whatever λ₁ does; a record that observes only coordinates orthogonal to v shows none either.
//
// Covariances are shrunk towards a multiple of the identity by Ledoit and Wolf's rule (2004, J. Multivariate Anal. 88,
// 365–411, in the normalised Frobenius norm), so that a window short against the
// dimension still gives a well-conditioned estimate; the propagator is regularised by a ridge of the same scale.
#pragma once

#include <kinesis/transitions/space.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>

#include <complex>
#include <limits>
#include <numeric>

namespace kinesis::transitions {

using Complex = std::complex<double>;

// ---- covariance -------------------------------------------------------------------------------------------------

struct Covariance {
  Matrix c;            // the estimate (shrunk when asked)
  double shrinkage = 0;  // the weight on the target m I
  std::vector<double> mean;
};

// The covariance of the rows of x (each row an observation), shrunk by Ledoit and Wolf's rule when `shrink`.
inline Covariance covariance(const Rectangular& x, bool shrink = true) {
  const std::size_t n = x.rows, p = x.cols;
  if (n < 2) throw Refusal("a covariance needs two observations or more");
  Covariance out;
  out.mean.assign(p, 0);
  for (std::size_t t = 0; t < n; t++)
    for (std::size_t j = 0; j < p; j++) out.mean[j] += x(t, j) / static_cast<double>(n);
  Matrix S(p);
  for (std::size_t t = 0; t < n; t++)
    for (std::size_t i = 0; i < p; i++)
      for (std::size_t j = 0; j < p; j++) S(i, j) += (x(t, i) - out.mean[i]) * (x(t, j) - out.mean[j]) / static_cast<double>(n);
  if (!shrink || p == 1) {
    for (double& v : S.a) v *= static_cast<double>(n) / static_cast<double>(n - 1);
    out.c = S;
    return out;
  }
  const double dp = static_cast<double>(p);
  double m = 0;
  for (std::size_t i = 0; i < p; i++) m += S(i, i) / dp;
  double d2 = 0;
  for (std::size_t i = 0; i < p; i++)
    for (std::size_t j = 0; j < p; j++) {
      const double v = S(i, j) - (i == j ? m : 0);
      d2 += v * v / dp;
    }
  double b2 = 0;
  for (std::size_t t = 0; t < n; t++) {
    double s = 0;
    for (std::size_t i = 0; i < p; i++)
      for (std::size_t j = 0; j < p; j++) {
        const double v = (x(t, i) - out.mean[i]) * (x(t, j) - out.mean[j]) - S(i, j);
        s += v * v / dp;
      }
    b2 += s / (static_cast<double>(n) * static_cast<double>(n));
  }
  b2 = std::min(b2, d2);
  const double w = d2 > 0 ? b2 / d2 : 1;
  out.shrinkage = w;
  out.c = Matrix(p);
  for (std::size_t i = 0; i < p; i++)
    for (std::size_t j = 0; j < p; j++) out.c(i, j) = (1 - w) * S(i, j) + (i == j ? w * m : 0);
  return out;
}

// ---- the lag-1 propagator ---------------------------------------------------------------------------------------

struct Propagator {
  Matrix a;                    // Â, the least-squares map x_t ↦ x_{t+1} of the centred window
  Matrix noise;                // the covariance of its residuals, Q̂
  std::vector<Complex> spectrum;  // Â's eigenvalues, largest modulus first
  double ridge = 0;
};

// Â = C₁ (C₀ + δ I)⁻¹ with C₀ = Σ x_t x_tᵀ and C₁ = Σ x_{t+1} x_tᵀ over the centred window; δ = ridge · tr C₀ / p.
inline Propagator propagator(const Rectangular& x, double ridge = 1e-6) {
  const std::size_t n = x.rows, p = x.cols;
  if (n < p + 3) throw Refusal("a propagator of " + std::to_string(p) + " coordinates needs a window of at least " + std::to_string(p + 3) + " states");
  std::vector<double> mean(p, 0);
  for (std::size_t t = 0; t < n; t++)
    for (std::size_t j = 0; j < p; j++) mean[j] += x(t, j) / static_cast<double>(n);
  Matrix C0(p), C1(p);
  for (std::size_t t = 0; t + 1 < n; t++)
    for (std::size_t i = 0; i < p; i++)
      for (std::size_t j = 0; j < p; j++) {
        C0(i, j) += (x(t, i) - mean[i]) * (x(t, j) - mean[j]);
        C1(i, j) += (x(t + 1, i) - mean[i]) * (x(t, j) - mean[j]);
      }
  double trace = 0;
  for (std::size_t i = 0; i < p; i++) trace += C0(i, i);
  Propagator out;
  out.ridge = ridge * trace / static_cast<double>(p);
  for (std::size_t i = 0; i < p; i++) C0(i, i) += out.ridge;
  // Â = C₁ C₀⁻¹: solve C₀ Âᵀ = C₁ᵀ (C₀ symmetric).
  Matrix C1t(p);
  for (std::size_t i = 0; i < p; i++)
    for (std::size_t j = 0; j < p; j++) C1t(i, j) = C1(j, i);
  const Matrix At = calx::linalg::solve(C0, C1t);
  out.a = Matrix(p);
  for (std::size_t i = 0; i < p; i++)
    for (std::size_t j = 0; j < p; j++) out.a(i, j) = At(j, i);
  out.noise = Matrix(p);
  for (std::size_t t = 0; t + 1 < n; t++) {
    std::vector<double> e(p);
    for (std::size_t i = 0; i < p; i++) {
      double s = x(t + 1, i) - mean[i];
      for (std::size_t j = 0; j < p; j++) s -= out.a(i, j) * (x(t, j) - mean[j]);
      e[i] = s;
    }
    for (std::size_t i = 0; i < p; i++)
      for (std::size_t j = 0; j < p; j++) out.noise(i, j) += e[i] * e[j] / static_cast<double>(n - 1 - p);
  }
  out.spectrum = calx::linalg::eigenvalues(out.a);
  std::sort(out.spectrum.begin(), out.spectrum.end(), [](Complex a, Complex b) {
    if (std::abs(a) != std::abs(b)) return std::abs(a) > std::abs(b);
    return a.imag() > b.imag();
  });
  return out;
}

// The rate and frequency of a discrete eigenvalue μ at a step Δt: μ = e^{(−κ + iω)Δt}.
struct Mode {
  Complex mu;
  double modulus = 0, argument = 0;
  double rate = 0;        // κ = −ln|μ| / Δt, the rate at which the mode decays (zero at a bifurcation)
  double frequency = 0;   // ω = arg μ / Δt
};

inline Mode mode(Complex mu, double dt) {
  Mode m;
  m.mu = mu;
  m.modulus = std::abs(mu);
  m.argument = std::arg(mu);
  m.rate = m.modulus > 0 ? -std::log(m.modulus) / dt : std::numeric_limits<double>::infinity();
  m.frequency = m.argument / dt;
  return m;
}

// ---- eigenvectors -----------------------------------------------------------------------------------------------

namespace detail {

// Solves M z = b for a complex matrix by elimination with partial pivoting; returns false when M is singular to
// working precision.
inline bool complex_solve(std::vector<Complex> M, std::vector<Complex>& b, std::size_t n) {
  for (std::size_t c = 0; c < n; c++) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < n; r++)
      if (std::abs(M[r * n + c]) > std::abs(M[piv * n + c])) piv = r;
    if (std::abs(M[piv * n + c]) == 0) return false;
    if (piv != c) {
      for (std::size_t j = 0; j < n; j++) std::swap(M[c * n + j], M[piv * n + j]);
      std::swap(b[c], b[piv]);
    }
    for (std::size_t r = c + 1; r < n; r++) {
      const Complex f = M[r * n + c] / M[c * n + c];
      if (f == Complex(0)) continue;
      for (std::size_t j = c; j < n; j++) M[r * n + j] -= f * M[c * n + j];
      b[r] -= f * b[c];
    }
  }
  for (std::size_t r = n; r-- > 0;) {
    Complex s = b[r];
    for (std::size_t k = r + 1; k < n; k++) s -= M[r * n + k] * b[k];
    b[r] = s / M[r * n + r];
  }
  return true;
}

}  // namespace detail

// The right eigenvector of A for its eigenvalue μ (or of Aᵀ, giving the conjugate of the left eigenvector, when
// `transpose`), by inverse iteration, normalised to unit length with its largest entry real and positive.
inline std::vector<Complex> eigenvector(const Matrix& A, Complex mu, bool transpose = false) {
  const std::size_t n = A.n;
  const double scale = std::max(1e-300, calx::linalg::norm1(A));
  const Complex shift = mu + Complex(1e-10 * scale, 1e-10 * scale);
  std::vector<Complex> M(n * n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) M[i * n + j] = (transpose ? A(j, i) : A(i, j)) - (i == j ? shift : Complex(0));
  std::vector<Complex> v(n, Complex(1.0 / std::sqrt(static_cast<double>(n))));
  for (std::size_t i = 0; i < n; i++) v[i] += Complex(0.01 * static_cast<double>(i % 7), 0);
  for (int iteration = 0; iteration < 8; iteration++) {
    if (!detail::complex_solve(M, v, n)) break;
    double norm = 0;
    for (const auto& z : v) norm += std::norm(z);
    norm = std::sqrt(norm);
    for (auto& z : v) z /= norm;
  }
  std::size_t big = 0;
  for (std::size_t i = 1; i < n; i++)
    if (std::abs(v[i]) > std::abs(v[big])) big = i;
  const Complex phase = std::abs(v[big]) > 0 ? std::conj(v[big]) / std::abs(v[big]) : Complex(1);
  for (auto& z : v) z *= phase;
  return v;
}

// The leading mode of a propagator with its right eigenvector (where it shows) and left eigenvector (what excites it).
struct Leading {
  Mode mode;
  std::vector<Complex> right, left;
  double excitation = 0;  // wᴴ Q̂ w / ‖w‖², the noise's variance along the left eigenvector
};

inline Leading leading(const Propagator& P, double dt) {
  if (P.spectrum.empty()) throw Refusal("a propagator of no coordinates has no leading mode");
  Leading L;
  L.mode = mode(P.spectrum.front(), dt);
  L.right = eigenvector(P.a, L.mode.mu, false);
  auto w = eigenvector(P.a, L.mode.mu, true);
  for (auto& z : w) z = std::conj(z);
  L.left = w;
  const std::size_t n = P.a.n;
  Complex s = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) s += std::conj(w[i]) * P.noise(i, j) * w[j];
  L.excitation = s.real();
  return L;
}

// ---- covariance spectra -----------------------------------------------------------------------------------------

struct CovarianceSpectrum {
  double largest = 0, trace = 0, share = 0;  // λ_max, tr C, λ_max / tr C
  std::vector<double> direction;             // the leading eigenvector
};

inline CovarianceSpectrum covariance_spectrum(const Matrix& C) {
  const auto e = calx::linalg::symmetric_eigen(C);
  CovarianceSpectrum s;
  const std::size_t n = C.n;
  s.largest = e.values.back();
  for (double v : e.values) s.trace += v;
  s.share = s.trace > 0 ? s.largest / s.trace : 0;
  s.direction.resize(n);
  for (std::size_t i = 0; i < n; i++) s.direction[i] = e.vectors(i, n - 1);
  return s;
}

// The largest eigenvalue of C₀⁻¹ C, the factor by which the variance along the most changed direction has grown since a
// baseline C₀, invariant under any linear change of coordinates (and so under the units each coordinate is measured in).
inline double relative_growth(const Matrix& C0, const Matrix& C) {
  const Matrix H = detail::symmetric_function(C0, [](double v) {
    if (!(v > 0)) throw Refusal("the baseline covariance is singular; shrink it or observe more states");
    return 1 / std::sqrt(v);
  });
  const Matrix M = H * C * H;
  return calx::linalg::symmetric_eigen(M).values.back();
}

// The mean absolute correlation between distinct coordinates.
inline double mean_absolute_correlation(const Matrix& C) {
  const std::size_t n = C.n;
  if (n < 2) return 0;
  double s = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      const double d = std::sqrt(C(i, i) * C(j, j));
      s += d > 0 ? std::abs(C(i, j)) / d : 0;
    }
  return s / (static_cast<double>(n) * static_cast<double>(n - 1) / 2);
}

// ---- Lyapunov equations, for exact answers ----------------------------------------------------------------------

// The C with C = A C Aᵀ + Q (A stable), by its Kronecker form (I − A ⊗ A) vec C = vec Q.
inline Matrix discrete_lyapunov(const Matrix& A, const Matrix& Q) {
  const std::size_t n = A.n, N = n * n;
  Matrix K(N), R(N);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      for (std::size_t k = 0; k < n; k++)
        for (std::size_t l = 0; l < n; l++) K(i * n + j, k * n + l) = (i == k && j == l ? 1.0 : 0.0) - A(i, k) * A(j, l);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) R(i * n + j, 0) = Q(i, j);
  const Matrix X = calx::linalg::solve(K, R);
  Matrix C(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) C(i, j) = X(i * n + j, 0);
  return C;
}

// The C with J C + C Jᵀ + B Bᵀ = 0 (J stable), the stationary covariance of dx = J x dt + B dW.
inline Matrix continuous_lyapunov(const Matrix& J, const Matrix& BBt) {
  const std::size_t n = J.n, N = n * n;
  Matrix K(N), R(N);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      for (std::size_t k = 0; k < n; k++) {
        K(i * n + j, k * n + j) += J(i, k);
        K(i * n + j, i * n + k) += J(j, k);
      }
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) R(i * n + j, 0) = -BBt(i, j);
  const Matrix X = calx::linalg::solve(K, R);
  Matrix C(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) C(i, j) = X(i * n + j, 0);
  return C;
}

}  // namespace kinesis::transitions
