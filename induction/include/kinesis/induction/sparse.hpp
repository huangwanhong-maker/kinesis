// Libraries of terms and sparse regression on them: the linear algebra of going from a record to a law.
//
// A record of a system is a table of states at its times (systems::Trajectory). Its derivatives in time are
// estimated by central differences of the second order (one-sided at the two ends; the steps need not be equal),
// which amplify noise by the inverse of the step; the weak form below avoids them, and asks for equal steps.
//
// A library is a list of terms of the record's names, θ_1 … θ_K: by default the monomials up to a degree. A law
// x' = Σ ξ_k θ_k(x) that holds on the record is a solution of the linear system Θ ξ = x', Θ the library evaluated
// on the record. Among the solutions a sparse one is sought by sequentially thresholded least squares (Brunton,
// Proctor and Kutz 2016, Proc. Natl Acad. Sci. 113, 3932–3937, their algorithm in the supporting text): solve by
// least squares, set to zero the coefficients below a threshold in magnitude, solve again on what is left, until
// the support does not change. It converges in at most K steps to a local minimiser of the least-squares error
// with a penalty on the number of terms (Zhang and Schaeffer 2019, Multiscale Model. Simul. 17, 948–972).
// *Caveat.* The answer depends on the library holding the true terms and on the threshold; with a library that
// does not, a sparse law is still returned, and fits. A record does not say which.
//
// The weak form (Messenger and Bortz 2021, Multiscale Model. Simul. 19, 1474–1497; Schaeffer and McCalla 2017).
// For a test function φ that vanishes at the ends of a window, ∫ φ x' dt = −∫ φ' x dt, so the law tested against
// φ is −∫ φ' x dt = Σ ξ_k ∫ φ θ_k(x) dt: an equation between integrals of the record, with no derivative of it.
// Here φ(t) = (1 − s²)^p on windows that overlap, with the trapezoid rule; the same sparse regression follows.
#pragma once

#include <kinesis/systems/simulate.hpp>
#include <kinesis/terms/program.hpp>

#include <calx/linalg/qr.hpp>
#include <calx/linalg/rectangular.hpp>

#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace kinesis::induction {

using calx::linalg::Rectangular;
using systems::State;
using systems::Trajectory;
using terms::TermPtr;

// The derivative in time of each component of a record at each of its times, by differences of second order.
inline std::vector<State> derivatives(const Trajectory& r) {
  const std::size_t n = r.x.size();
  if (n < 3) throw kernel::Refusal("a derivative needs at least three samples");
  std::vector<State> d(n, State(r.x[0].size()));
  for (std::size_t i = 0; i < r.x[0].size(); i++) {
    for (std::size_t k = 1; k + 1 < n; k++) {
      // The derivative of the parabola through three samples, at the middle one, for steps that may differ.
      const double h0 = r.t[k] - r.t[k - 1], h1 = r.t[k + 1] - r.t[k];
      d[k][i] = (-h1 / (h0 * (h0 + h1))) * r.x[k - 1][i] + ((h1 - h0) / (h0 * h1)) * r.x[k][i] + (h0 / (h1 * (h0 + h1))) * r.x[k + 1][i];
    }
    const double h0 = r.t[1] - r.t[0], h1 = r.t[2] - r.t[1];
    d[0][i] = -(2 * h0 + h1) / (h0 * (h0 + h1)) * r.x[0][i] + (h0 + h1) / (h0 * h1) * r.x[1][i] - h0 / (h1 * (h0 + h1)) * r.x[2][i];
    const double g0 = r.t[n - 2] - r.t[n - 3], g1 = r.t[n - 1] - r.t[n - 2];
    d[n - 1][i] = g1 / (g0 * (g0 + g1)) * r.x[n - 3][i] - (g0 + g1) / (g0 * g1) * r.x[n - 2][i] + (g0 + 2 * g1) / (g1 * (g0 + g1)) * r.x[n - 1][i];
  }
  return d;
}

// The monomials of the variables up to a degree, the constant first, by degree and then in the order given.
inline std::vector<TermPtr> monomials(const std::vector<std::string>& variables, unsigned degree, bool constant = true) {
  std::vector<TermPtr> out;
  std::vector<unsigned> power(variables.size(), 0);
  std::function<void(std::size_t, unsigned, unsigned)> make = [&](std::size_t at, unsigned left, unsigned total) {
    if (at == variables.size()) {
      if (left != 0) return;
      if (total == 0 && !constant) return;
      std::vector<TermPtr> factors;
      for (std::size_t i = 0; i < variables.size(); i++)
        if (power[i]) factors.push_back(terms::power(terms::symbol(variables[i]), terms::number(power[i])));
      out.push_back(factors.empty() ? terms::number(1) : terms::product(std::move(factors)));
      return;
    }
    for (unsigned p = left + 1; p-- > 0;) {
      power[at] = p;
      make(at + 1, left - p, total);
    }
    power[at] = 0;
  };
  for (unsigned d = 0; d <= degree; d++) make(0, d, d);
  return out;
}

// The library a declaration names: poly(x, y, 3) is the monomials of x and y up to the third degree; any other
// term is a library of that one term; a sum is the library of its terms.
inline std::vector<TermPtr> library_of(const TermPtr& spec) {
  if (spec->kind == terms::Kind::call && spec->name == "poly" && spec->args.size() >= 2 && terms::is_integer(spec->args.back())) {
    std::vector<std::string> variables;
    for (std::size_t i = 0; i + 1 < spec->args.size(); i++) {
      if (spec->args[i]->kind != terms::Kind::symbol) throw kernel::Refusal("poly(…, degree) lists names and then a degree");
      variables.push_back(spec->args[i]->name);
    }
    return monomials(variables, static_cast<unsigned>(spec->args.back()->value));
  }
  if (spec->kind == terms::Kind::sum) return spec->args;
  return {spec};
}

// The terms evaluated on a record: one row for each time, one column for each term. `names` are the record's.
inline Rectangular evaluate(const std::vector<TermPtr>& library, const std::vector<std::string>& names, const std::vector<State>& rows, const std::vector<double>* times = nullptr, const std::string& time = "t") {
  terms::Frame frame;
  for (const auto& n : names) frame.add(n);
  const int slot_t = frame.add(time);
  std::vector<terms::Program> programs;
  for (const auto& t : library) programs.push_back(terms::compile(t, frame.resolver()));
  Rectangular theta(rows.size(), library.size());
  std::vector<double> values(frame.size(), 0.0);
  for (std::size_t k = 0; k < rows.size(); k++) {
    std::copy(rows[k].begin(), rows[k].end(), values.begin());
    values[static_cast<std::size_t>(slot_t)] = times ? (*times)[k] : 0.0;
    for (std::size_t j = 0; j < library.size(); j++) theta(k, j) = programs[j](values);
  }
  return theta;
}

struct SparseOptions {
  double threshold = 0.05;     // coefficients below it in magnitude are set to zero
  std::size_t most = 20;       // rounds of thresholding
  std::vector<char> keep;      // keep[k] != 0: the coefficient is never thresholded (empty: none is kept)
};

struct SparseFit {
  std::vector<double> coefficients;
  std::size_t terms = 0;       // not zero
  double residual = 0;         // the root-mean-square of Θ ξ − y
  std::size_t rounds = 0;
};

// A sparse solution of Θ ξ ≈ y by sequentially thresholded least squares.
inline SparseFit stlsq(const Rectangular& theta, std::span<const double> y, const SparseOptions& o = {}) {
  const std::size_t K = theta.cols, N = theta.rows;
  std::vector<char> active(K, 1);
  SparseFit fit;
  fit.coefficients.assign(K, 0.0);
  for (fit.rounds = 0; fit.rounds < o.most; fit.rounds++) {
    std::vector<std::size_t> columns;
    for (std::size_t k = 0; k < K; k++)
      if (active[k]) columns.push_back(k);
    std::fill(fit.coefficients.begin(), fit.coefficients.end(), 0.0);
    if (columns.empty()) break;
    Rectangular sub(N, columns.size());
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < columns.size(); j++) sub(i, j) = theta(i, columns[j]);
    const std::vector<double> xi = calx::linalg::least_squares(sub, y);
    bool changed = false;
    for (std::size_t j = 0; j < columns.size(); j++) {
      const std::size_t k = columns[j];
      const bool kept = !o.keep.empty() && o.keep[k];
      if (!kept && std::fabs(xi[j]) < o.threshold) active[k] = 0, changed = true;
      else fit.coefficients[k] = xi[j];
    }
    if (!changed) break;
  }
  double sum = 0;
  for (std::size_t i = 0; i < N; i++) {
    double model = 0;
    for (std::size_t k = 0; k < K; k++) model += theta(i, k) * fit.coefficients[k];
    sum += (model - y[i]) * (model - y[i]);
  }
  fit.residual = std::sqrt(sum / static_cast<double>(std::max<std::size_t>(N, 1)));
  for (double c : fit.coefficients) fit.terms += c != 0;
  return fit;
}

// The weak form of x_i' = Σ ξ_k θ_k on a record sampled at equal steps: for each window, the row of ∫ φ θ_k dt
// and the value −∫ φ' x_i dt. `half` is the half-width of a window in samples; windows advance by `stride`.
struct Weak {
  Rectangular matrix;                 // windows × terms
  std::vector<std::vector<double>> targets;  // for each component, a value for each window
};
inline Weak weak_form(const Trajectory& r, const Rectangular& theta, std::size_t half = 20, std::size_t stride = 5, unsigned power = 4) {
  const std::size_t n = r.x.size(), dim = r.x[0].size();
  if (n < 2 * half + 1) throw kernel::Refusal("the record is shorter than a window of the weak form");
  const double h = r.t[1] - r.t[0], width = static_cast<double>(half) * h;
  std::vector<double> phi(2 * half + 1), dphi(2 * half + 1);
  for (std::size_t j = 0; j <= 2 * half; j++) {
    const double s = (static_cast<double>(j) - static_cast<double>(half)) / static_cast<double>(half);
    phi[j] = std::pow(1 - s * s, static_cast<double>(power));
    dphi[j] = -2 * static_cast<double>(power) * s * std::pow(1 - s * s, static_cast<double>(power) - 1) / width;
  }
  std::vector<std::size_t> starts;
  for (std::size_t start = 0; start + 2 * half < n; start += stride) starts.push_back(start);
  Weak w;
  w.matrix = Rectangular(starts.size(), theta.cols);
  w.targets.assign(dim, std::vector<double>(starts.size(), 0.0));
  for (std::size_t row = 0; row < starts.size(); row++) {
    for (std::size_t j = 0; j <= 2 * half; j++) {
      const double weight = (j == 0 || j == 2 * half ? 0.5 : 1.0) * h;
      const std::size_t k = starts[row] + j;
      for (std::size_t c = 0; c < theta.cols; c++) w.matrix(row, c) += weight * phi[j] * theta(k, c);
      for (std::size_t i = 0; i < dim; i++) w.targets[i][row] -= weight * dphi[j] * r.x[k][i];
    }
  }
  return w;
}

}  // namespace kinesis::induction
