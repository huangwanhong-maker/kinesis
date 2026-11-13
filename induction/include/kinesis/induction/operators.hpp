// Linear operators fitted to a record: dynamic mode decomposition, and its extension on a dictionary.
//
// A record sampled at equal steps gives pairs (x_k, x_{k+1}). Dynamic mode decomposition finds the matrix A that
// best carries each sample to the next, A = argmin Σ |x_{k+1} − A x_k|² (Schmid 2010, J. Fluid Mech. 656, 5–28;
// the exact form of Tu et al. 2014, J. Comput. Dyn. 1, 391–421). For a linear system sampled at a step h it is
// the propagator e^{hL}, and its eigenvalues λ give the rates log(λ)/h of L exactly.
//
// On a dictionary of functions ψ_1 … ψ_M of the state the same least squares gives a matrix K that carries the
// values of the dictionary one step: a finite section of Koopman's operator, the lift of the system to its
// observables (extended DMD: Williams, Kevrekidis and Rowley 2015, J. Nonlinear Sci. 25, 1307–1346). When the
// span of the dictionary is invariant under the operator, the eigenvalues of K are eigenvalues of the operator.
// *Caveat.* When it is not, K has eigenvalues that belong to no eigenfunction (spectral pollution; Colbrook and
// Townsend 2024, Commun. Pure Appl. Math. 77, 221–283, give residuals that tell them apart, which are not
// computed here). A record with noise biases the eigenvalues of plain DMD towards zero.
#pragma once

#include <kinesis/induction/sparse.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>
#include <calx/linalg/qr.hpp>

#include <complex>
#include <vector>

namespace kinesis::induction {

struct Linear {
  calx::linalg::Matrix matrix;                    // the values at one step from those before: y ≈ matrix · x
  std::vector<std::complex<double>> eigenvalues;  // by decreasing modulus
  std::vector<std::complex<double>> rates;        // log(λ) / h: growth rates and angular frequencies
  double residual = 0;                            // the root-mean-square of what one step leaves unexplained
};

namespace detail {
// The matrix K with rows[k + 1] ≈ K rows[k], by least squares, row by row.
inline Linear propagator(const Rectangular& values, double step) {
  const std::size_t n = values.rows, m = values.cols;
  Rectangular now(n - 1, m);
  for (std::size_t k = 0; k + 1 < n; k++)
    for (std::size_t j = 0; j < m; j++) now(k, j) = values(k, j);
  const calx::linalg::Qr qr = calx::linalg::qr(now);
  Linear out;
  out.matrix = calx::linalg::Matrix(m);
  double sum = 0;
  std::vector<double> next(n - 1);
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t k = 0; k + 1 < n; k++) next[k] = values(k + 1, i);
    const std::vector<double> row = calx::linalg::least_squares(qr, next);
    for (std::size_t j = 0; j < m; j++) out.matrix(i, j) = row[j];
    for (std::size_t k = 0; k + 1 < n; k++) {
      double model = 0;
      for (std::size_t j = 0; j < m; j++) model += row[j] * values(k, j);
      sum += (model - next[k]) * (model - next[k]);
    }
  }
  out.residual = std::sqrt(sum / static_cast<double>((n - 1) * m));
  out.eigenvalues = calx::linalg::eigenvalues(out.matrix);
  std::sort(out.eigenvalues.begin(), out.eigenvalues.end(), [](const std::complex<double>& a, const std::complex<double>& b) { return std::abs(a) > std::abs(b); });
  for (const auto& l : out.eigenvalues) out.rates.push_back(std::log(l) / step);
  return out;
}
}  // namespace detail

// Dynamic mode decomposition of a record sampled at equal steps.
inline Linear dmd(const Trajectory& record) {
  if (record.x.size() < record.x[0].size() + 2) throw kernel::Refusal("dynamic mode decomposition needs more samples than components");
  Rectangular values(record.x.size(), record.x[0].size());
  for (std::size_t k = 0; k < values.rows; k++)
    for (std::size_t j = 0; j < values.cols; j++) values(k, j) = record.x[k][j];
  return detail::propagator(values, record.t[1] - record.t[0]);
}

// Extended dynamic mode decomposition on a dictionary of terms of the record's names.
inline Linear edmd(const Trajectory& record, const std::vector<TermPtr>& dictionary) {
  if (record.x.size() < dictionary.size() + 2) throw kernel::Refusal("the record has fewer samples than the dictionary has functions");
  return detail::propagator(evaluate(dictionary, record.names, record.x), record.t[1] - record.t[0]);
}

}  // namespace kinesis::induction
