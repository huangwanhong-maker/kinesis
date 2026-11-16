// The modes of a field recorded at a few equal steps: dynamic mode decomposition by the method of snapshots (decision 0005).
//
// A field of n cells recorded m times, with n far larger than m, is too short a record for induction's dmd, which fits
// an n × n matrix. Exact DMD (Schmid 2010, J. Fluid Mech. 656, 5–28; Tu, Rowley, Luchtenburg, Brunton and Kutz 2014,
// J. Comput. Dyn. 1, 391–421, Algorithm 2) fits the operator only on the span of the snapshots. With X = [u_0 … u_{m−2}]
// and Y = [u_1 … u_{m−1}]:
//
//   the snapshots' Gram matrix XᵀX = V Σ² Vᵀ (Sirovich's method of snapshots, 1987, Q. Appl. Math. 45, 561–571), its
//   r largest values kept; U = X V Σ⁻¹ the proper orthogonal modes; Ã = Uᵀ Y V Σ⁻¹ = Σ⁻¹ Vᵀ (XᵀY) V Σ⁻¹, r × r;
//   Ã w_k = λ_k w_k; the modes φ_k = Y V Σ⁻¹ w_k / λ_k; the amplitudes b from Φ b = u_0 by least squares;
//   u(t) ≈ Re Σ_k b_k λ_k^{(t − t_0)/h} φ_k.
//
// Each mode grows or decays at the rate Re(ln λ)/h and turns at Im(ln λ)/h. For a field that a linear law moves, the
// modes are the law's own (Tu et al., Theorem 3); for one that a nonlinear law moves they are the best linear account of
// the record, and a forecast from them is a linear extrapolation in the span of what was seen. *Caveat.* With few
// snapshots every mode is fitted to few numbers, and the record's error biases the eigenvalues toward zero (Dawson,
// Hemati, Williams and Rowley 2016, Exp. Fluids 57, 42); the residual says how well the modes carry each snapshot to
// the next.
#pragma once

#include <kinesis/balance/learn.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/linalg/eigenvalues.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

namespace kinesis::balance {

using cd = std::complex<double>;

struct Mode {
  cd eigenvalue, rate;     // λ, and ln(λ) / h
  cd amplitude;            // b
  std::vector<cd> shape;   // φ over the cells
  double size = 0;         // |b| ‖φ‖: how much of the first snapshot the mode carries
};

struct Modes {
  double t0 = 0, step = 0;
  std::vector<double> singular_values;  // of the snapshots X, largest first
  std::size_t rank = 0;
  std::vector<Mode> modes;              // by decreasing size
  double residual = 0;                  // ‖Y − Ũ Ã Σ Vᵀ‖ / ‖Y‖: what the fitted operator does not carry, relative
  std::vector<double> at(double t) const {
    if (modes.empty()) return {};
    std::vector<double> u(modes[0].shape.size(), 0.0);
    const double k = (t - t0) / step;
    for (const Mode& m : modes) {
      if (std::abs(m.eigenvalue) == 0 && k != 0) continue;  // a mode that vanishes after one step
      const cd f = std::abs(m.eigenvalue) == 0 ? m.amplitude : m.amplitude * std::exp(k * std::log(m.eigenvalue));
      for (std::size_t i = 0; i < u.size(); i++) u[i] += (f * m.shape[i]).real();
    }
    return u;
  }
};

namespace detail {
// x solving A x = b, A small and complex, by Gaussian elimination with partial pivoting; A singular is refused.
inline std::vector<cd> csolve(std::vector<std::vector<cd>> A, std::vector<cd> b) {
  const std::size_t n = b.size();
  for (std::size_t c = 0; c < n; c++) {
    std::size_t p = c;
    for (std::size_t i = c + 1; i < n; i++)
      if (std::abs(A[i][c]) > std::abs(A[p][c])) p = i;
    if (std::abs(A[p][c]) == 0) throw Refusal("a singular system in the modes");
    std::swap(A[p], A[c]), std::swap(b[p], b[c]);
    for (std::size_t i = c + 1; i < n; i++) {
      const cd f = A[i][c] / A[c][c];
      for (std::size_t j = c; j < n; j++) A[i][j] -= f * A[c][j];
      b[i] -= f * b[c];
    }
  }
  std::vector<cd> x(n);
  for (std::size_t i = n; i-- > 0;) {
    cd s = b[i];
    for (std::size_t j = i + 1; j < n; j++) s -= A[i][j] * x[j];
    x[i] = s / A[i][i];
  }
  return x;
}

// An eigenvector of the small real matrix A for its eigenvalue λ, by inverse iteration with a slight shift.
inline std::vector<cd> eigenvector(const calx::linalg::Matrix& A, cd lambda) {
  const std::size_t n = A.n;
  double scale = 0;
  for (double v : A.a) scale = std::max(scale, std::abs(v));
  const cd shift = lambda + cd(1e-10 * std::max(scale, 1e-300), 1e-10 * std::max(scale, 1e-300));
  std::vector<std::vector<cd>> M(n, std::vector<cd>(n));
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) M[i][j] = A(i, j) - (i == j ? shift : cd(0));
  std::vector<cd> x(n, cd(1.0, 0.3));
  for (int it = 0; it < 6; it++) {
    x = csolve(M, x);
    double norm = 0;
    for (const cd& v : x) norm += std::norm(v);
    norm = std::sqrt(norm);
    for (cd& v : x) v /= norm;
  }
  return x;
}
}  // namespace detail

// The modes of a record taken at equal steps; `rank` 0 keeps every singular value above `tolerance` times the largest.
// The Gram matrix squares the singular values, so that a value below about √ε of the largest is rounding, not a mode:
// the tolerance is 10⁻⁷ by default.
inline Modes modes(const Record& r, std::size_t rank = 0, double tolerance = 1e-7) {
  if (r.t.size() < 3) throw Refusal("modes need at least three times");
  if (r.t.size() != r.u.size()) throw Refusal("a record needs one field per time");
  const double h = r.t[1] - r.t[0];
  for (std::size_t k = 1; k < r.t.size(); k++)
    if (std::abs((r.t[k] - r.t[k - 1]) - h) > 1e-9 * std::max(1.0, std::abs(h))) throw Refusal("modes need the record's times equally spaced");
  const std::size_t n = r.u[0].size(), m = r.t.size() - 1;  // m pairs
  for (const auto& u : r.u)
    if (u.size() != n) throw Refusal("every field of the record must have the same cells");
  // Gram matrices XᵀX and XᵀY, m × m.
  calx::linalg::Matrix G(m), H(m);
  for (std::size_t a = 0; a < m; a++)
    for (std::size_t b = 0; b < m; b++) {
      double g = 0, y = 0;
      for (std::size_t i = 0; i < n; i++) g += r.u[a][i] * r.u[b][i], y += r.u[a][i] * r.u[b + 1][i];
      G(a, b) = g, H(a, b) = y;
    }
  const calx::linalg::SymmetricEigen se = calx::linalg::symmetric_eigen(G);
  std::vector<std::size_t> order(m);
  for (std::size_t k = 0; k < m; k++) order[k] = m - 1 - k;  // largest first
  Modes out;
  out.t0 = r.t[0], out.step = h;
  const double top = std::sqrt(std::max(se.values[order[0]], 0.0));
  for (std::size_t k : order) out.singular_values.push_back(std::sqrt(std::max(se.values[k], 0.0)));
  std::size_t keep = 0;
  for (double s : out.singular_values)
    if (s > tolerance * top) keep++;
  if (rank > 0) keep = std::min(keep, rank);
  if (keep == 0) throw Refusal("the record is zero");
  out.rank = keep;
  // Ã = Σ⁻¹ Vᵀ H V Σ⁻¹.
  calx::linalg::Matrix At(keep);
  for (std::size_t p = 0; p < keep; p++)
    for (std::size_t q = 0; q < keep; q++) {
      double s = 0;
      for (std::size_t a = 0; a < m; a++)
        for (std::size_t b = 0; b < m; b++) s += se.vectors(a, order[p]) * H(a, b) * se.vectors(b, order[q]);
      At(p, q) = s / (out.singular_values[p] * out.singular_values[q]);
    }
  const std::vector<cd> lambdas = calx::linalg::eigenvalues(At);
  // Y V Σ⁻¹, n × keep.
  std::vector<std::vector<double>> YV(keep, std::vector<double>(n, 0.0));
  for (std::size_t q = 0; q < keep; q++)
    for (std::size_t b = 0; b < m; b++) {
      const double c = se.vectors(b, order[q]) / out.singular_values[q];
      if (c == 0) continue;
      for (std::size_t i = 0; i < n; i++) YV[q][i] += c * r.u[b + 1][i];
    }
  for (const cd& l : lambdas) {
    Mode md;
    md.eigenvalue = l;
    md.rate = std::abs(l) > 0 ? std::log(l) / h : cd(-std::numeric_limits<double>::infinity(), 0);
    const std::vector<cd> w = detail::eigenvector(At, l);
    md.shape.assign(n, cd(0));
    for (std::size_t q = 0; q < keep; q++)
      for (std::size_t i = 0; i < n; i++) md.shape[i] += YV[q][i] * w[q];
    if (std::abs(l) > 0)
      for (cd& v : md.shape) v /= l;
    out.modes.push_back(std::move(md));
  }
  // Amplitudes: Φ b = u_0 in least squares, by the normal equations (keep × keep).
  const std::size_t K = out.modes.size();
  std::vector<std::vector<cd>> N(K, std::vector<cd>(K, cd(0)));
  std::vector<cd> rhs(K, cd(0));
  for (std::size_t a = 0; a < K; a++) {
    for (std::size_t b = 0; b < K; b++)
      for (std::size_t i = 0; i < n; i++) N[a][b] += std::conj(out.modes[a].shape[i]) * out.modes[b].shape[i];
    for (std::size_t i = 0; i < n; i++) rhs[a] += std::conj(out.modes[a].shape[i]) * r.u[0][i];
  }
  const std::vector<cd> bs = detail::csolve(N, rhs);
  for (std::size_t a = 0; a < K; a++) {
    out.modes[a].amplitude = bs[a];
    double norm = 0;
    for (const cd& v : out.modes[a].shape) norm += std::norm(v);
    out.modes[a].size = std::abs(bs[a]) * std::sqrt(norm);
  }
  std::sort(out.modes.begin(), out.modes.end(), [](const Mode& a, const Mode& b) { return a.size > b.size; });
  // How well the fitted operator carries each snapshot to the next: ‖Y − Ỹ‖ / ‖Y‖, Ỹ = U Ã Uᵀ X.
  double num = 0, den = 0;
  for (std::size_t b = 0; b < m; b++) {
    // Uᵀ x_b = Σ Vᵀ e_b (coordinates of x_b on U).
    std::vector<double> z(keep);
    for (std::size_t p = 0; p < keep; p++) z[p] = out.singular_values[p] * se.vectors(b, order[p]);
    std::vector<double> az(keep, 0.0);
    for (std::size_t p = 0; p < keep; p++)
      for (std::size_t q = 0; q < keep; q++) az[p] += At(p, q) * z[q];
    // U az = X V Σ⁻¹ az.
    std::vector<double> coef(m, 0.0);
    for (std::size_t p = 0; p < keep; p++)
      for (std::size_t a = 0; a < m; a++) coef[a] += se.vectors(a, order[p]) / out.singular_values[p] * az[p];
    for (std::size_t i = 0; i < n; i++) {
      double v = 0;
      for (std::size_t a = 0; a < m; a++) v += coef[a] * r.u[a][i];
      num += (r.u[b + 1][i] - v) * (r.u[b + 1][i] - v);
      den += r.u[b + 1][i] * r.u[b + 1][i];
    }
  }
  out.residual = den > 0 ? std::sqrt(num / den) : 0;
  return out;
}

}  // namespace kinesis::balance
