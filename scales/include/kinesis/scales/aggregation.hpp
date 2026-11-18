// When a coarse law is exactly what a fine law does to the sums of its parts (decision 0007).
//
// A linear law x' = M x on n fine quantities and an aggregation y = A x onto k coarse ones (A a k × n matrix: rows of ones on
// each group for sums, or weights) aggregate perfectly when the coarse quantities obey a closed law of their own, y' = M̃ y,
// for every initial x. That holds if and only if A M = M̃ A for some M̃, that is, when the row space of A is invariant under
// Mᵀ (Iwasa, Andreasen and Levin 1987, Ecological Modelling 37, 287–302, "Aggregation in model ecosystems. I. Perfect
// aggregation"). For a partition summed, the condition reads: for every two groups I, J, the total Σ_{i∈I} M_ij that a
// unit j of J sends into I is the same for every j in J (the column form of Kemeny and Snell's lumpability, 1960, Finite
// Markov Chains, §6.3). The least-squares M̃ = A M Aᵀ (A Aᵀ)⁻¹ is reported with the relative residual
// ‖A M − M̃ A‖_F / ‖A M‖_F: zero exactly when the aggregation is perfect, and otherwise how far a coarse law of the sums
// misdescribes what the fine law does to them.
#pragma once

#include <cmath>

#include "calx/linalg/dense.hpp"
#include "calx/linalg/rectangular.hpp"
#include "kinesis/kernel/answer.hpp"

namespace kinesis::scales {

struct Aggregation {
  calx::linalg::Rectangular coarse;  // M̃, k × k
  double residual = 0;               // ‖A M − M̃ A‖_F / ‖A M‖_F
  bool perfect = false;
};

inline Aggregation perfect_aggregation(const calx::linalg::Rectangular& M, const calx::linalg::Rectangular& A, double tolerance = 1e-10) {
  if (M.rows != M.cols) throw kernel::Refusal("a linear law's matrix is square");
  if (A.cols != M.rows || A.rows == 0) throw kernel::Refusal("the aggregation must have one column per fine quantity");
  const std::size_t n = M.rows, k = A.rows;
  const calx::linalg::Rectangular AM = A * M;
  calx::linalg::Matrix G(k);
  for (std::size_t a = 0; a < k; a++)
    for (std::size_t b = 0; b < k; b++) {
      double s = 0;
      for (std::size_t j = 0; j < n; j++) s += A(a, j) * A(b, j);
      G(a, b) = s;
    }
  const calx::linalg::Matrix Gi = calx::linalg::inverse(G);
  calx::linalg::Rectangular AMAt(k, k);
  for (std::size_t a = 0; a < k; a++)
    for (std::size_t b = 0; b < k; b++) {
      double s = 0;
      for (std::size_t j = 0; j < n; j++) s += AM(a, j) * A(b, j);
      AMAt(a, b) = s;
    }
  Aggregation out;
  out.coarse = calx::linalg::Rectangular(k, k);
  for (std::size_t a = 0; a < k; a++)
    for (std::size_t b = 0; b < k; b++) {
      double s = 0;
      for (std::size_t c = 0; c < k; c++) s += AMAt(a, c) * Gi(c, b);
      out.coarse(a, b) = s;
    }
  const calx::linalg::Rectangular MA = out.coarse * A;
  double num = 0, den = 0;
  for (std::size_t i = 0; i < AM.a.size(); i++) num += (AM.a[i] - MA.a[i]) * (AM.a[i] - MA.a[i]), den += AM.a[i] * AM.a[i];
  out.residual = den > 0 ? std::sqrt(num / den) : std::sqrt(num);
  out.perfect = out.residual <= tolerance;
  return out;
}

}  // namespace kinesis::scales
