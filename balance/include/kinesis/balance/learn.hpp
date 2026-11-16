// Learning a balance law on a graph from a record of the field at a few times (decision 0005).
//
// A record is the field at times t_0 < t_1 < … < t_m, which need not be equally spaced and may be few (a census every
// five or ten years). Since every law of law.hpp is linear in its coefficients, two estimates are made.
//
// THE INTEGRAL FORM. Between two times the law integrates to u(t_{k+1}) − u(t_k) = Σ_j θ_j ∫ T_j[u] dt, and the
// integral is taken by the trapezoidal rule, (t_{k+1} − t_k)(T_j[u(t_k)] + T_j[u(t_{k+1})])/2. Every cell and every
// interval gives one equation, linear in θ, solved by weighted least squares (QR, calx). It needs no derivative of the
// record, which with a few samples cannot be had, and its error is the trapezoid's, of the second order in the step:
// small when the field changes little between the times, biased when it changes much. It is the integral form of sparse
// regression for coarse samples (Schaeffer and McCalla 2017, Phys. Rev. E 96, 023302), on a graph.
//
// SHOOTING. The law is run from the first time and compared with the record at every later time; the coefficients that
// make the run closest are found by the Levenberg–Marquardt iteration (Levenberg 1944; Marquardt 1963, J. SIAM 11,
// 431–441), from the integral estimate, the Jacobian by differences of runs. Its residual is the run's distance from
// the record, which is what a forecast is judged by; it has no bias of the step. It finds a local minimum.
//
// WEIGHTS. A census count's error grows with the count, so each equation is weighed by 1/√(ū + c), ū the mean of its
// cell's values over the interval (or its value) and c a floor, the field's mean unless stated: Poisson's variance
// with room for the empty cells. `weights = "none"` weighs all alike.
//
// SPARSITY. With a threshold, the integral form's coefficients are found by sequentially thresholded least squares
// (induction's stlsq) on columns scaled to unit root-mean-square, so that the threshold is a share of the change
// explained, and the terms it keeps are then refitted.
//
// WELL-POSEDNESS. A diffusion with a negative coefficient is the backward heat equation, which is ill-posed (Hadamard 1923,
// Lectures on Cauchy's Problem, §I.3): its run sharpens every small wiggle without bound, and nothing it forecasts can be
// trusted. So by default a diffusion or a crowding that the least squares find negative is set aside, said so, and the
// other terms are fitted again without it (an active set, as in non-negative least squares; Lawson and Hanson 1974,
// Solving Least Squares Problems, ch. 23), until none is negative.
//
// STANDARD ERRORS are those of least squares, s²(XᵀX)⁻¹ with the residual variance s², of the last linearisation for
// shooting. *Caveat.* They treat the equations as independent; the errors of neighbouring cells and of one cell in
// successive intervals are not, so they understate the uncertainty, by a factor that a spatial bootstrap would measure.
// A coefficient two standard errors from zero is not thereby shown to be non-zero.
#pragma once

#include <kinesis/balance/law.hpp>
#include <kinesis/induction/sparse.hpp>

#include <calx/linalg/qr.hpp>
#include <calx/linalg/rectangular.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace kinesis::balance {

using calx::linalg::Rectangular;

struct Record {
  std::vector<double> t;
  std::vector<std::vector<double>> u;  // the field at each time
};

inline void require_record(const Graph& g, const Record& r, std::size_t least = 2) {
  if (r.t.size() != r.u.size()) throw Refusal("a record needs one field per time");
  if (r.t.size() < least) throw Refusal("the record has " + std::to_string(r.t.size()) + " times; this needs at least " + std::to_string(least));
  for (std::size_t k = 0; k < r.t.size(); k++) {
    require_size(g, r.u[k], "a field of the record");
    if (k > 0 && !(r.t[k] > r.t[k - 1])) throw Refusal("the record's times must increase");
  }
}

inline Record first(const Record& r, std::size_t count) {
  Record out;
  for (std::size_t k = 0; k < count && k < r.t.size(); k++) out.t.push_back(r.t[k]), out.u.push_back(r.u[k]);
  return out;
}

struct LearnOptions {
  std::string weights = "poisson";  // or "none"
  double floor = -1;                // c in 1/√(ū + c); below zero: the record's mean
  double threshold = 0;             // sparsity on scaled columns; 0: none
  std::vector<char> keep;           // terms never thresholded
  std::size_t most = 60;            // Levenberg–Marquardt iterations
  bool well_posed = true;           // a diffusion or crowding found negative is dropped, and the rest refitted
};

struct Estimate {
  Law law;
  std::vector<double> standard_errors;
  double residual = 0;    // root-mean-square of the weighted residuals
  double explained = 0;   // 1 − (residual sum of squares) / (total sum of squares about the mean), weighted
  std::size_t rows = 0;   // equations
  std::size_t terms = 0;  // coefficients not zero
  std::size_t iterations = 0;
  bool converged = true;
  std::string method;
  std::vector<std::string> dropped;  // terms set aside, each with why
};

namespace detail {

inline double floor_of(const Record& r, const LearnOptions& o) {
  if (o.floor >= 0) return o.floor;
  double s = 0;
  std::size_t n = 0;
  for (const auto& u : r.u)
    for (double v : u) s += v, n++;
  return n ? std::max(s / static_cast<double>(n), 1e-12) : 1.0;
}

inline double weight(double mean, double c, const LearnOptions& o) {
  return o.weights == "none" ? 1.0 : 1.0 / std::sqrt(std::max(mean, 0.0) + c);
}

// Weighted least squares with standard errors, on the columns that are active.
inline void solve(const Rectangular& X, const std::vector<double>& y, const std::vector<char>& active, Estimate& e) {
  std::vector<std::size_t> cols;
  for (std::size_t j = 0; j < X.cols; j++)
    if (active[j]) cols.push_back(j);
  e.law.coefficients.assign(X.cols, 0.0);
  e.standard_errors.assign(X.cols, 0.0);
  double tss = 0, mean = 0;
  for (double v : y) mean += v;
  mean /= static_cast<double>(y.size());
  for (double v : y) tss += (v - mean) * (v - mean);
  if (cols.empty()) {
    double rss = 0;
    for (double v : y) rss += v * v;
    e.residual = std::sqrt(rss / static_cast<double>(y.size()));
    e.explained = tss > 0 ? 1 - rss / tss : 0;
    e.terms = 0;
    return;
  }
  Rectangular A(X.rows, cols.size());
  for (std::size_t i = 0; i < X.rows; i++)
    for (std::size_t c = 0; c < cols.size(); c++) A(i, c) = X(i, cols[c]);
  const calx::linalg::Qr f = calx::linalg::qr(A);
  const std::vector<double> b = calx::linalg::least_squares(f, y);
  double rss = 0;
  for (std::size_t i = 0; i < A.rows; i++) {
    double m = 0;
    for (std::size_t c = 0; c < cols.size(); c++) m += A(i, c) * b[c];
    rss += (y[i] - m) * (y[i] - m);
  }
  const double dof = std::max(1.0, static_cast<double>(A.rows) - static_cast<double>(f.rank));
  const calx::linalg::Matrix C = calx::linalg::normal_inverse(f);
  for (std::size_t c = 0; c < cols.size(); c++) {
    e.law.coefficients[cols[c]] = b[c];
    e.standard_errors[cols[c]] = std::sqrt(std::max(0.0, rss / dof * C(c, c)));
  }
  e.residual = std::sqrt(rss / static_cast<double>(A.rows));
  e.explained = tss > 0 ? 1 - rss / tss : 0;
  e.terms = cols.size();
}

}  // namespace detail

// The integral form: weighted least squares (or thresholded, with o.threshold) on increments between the record's times.
inline Estimate integral(const Graph& g, const std::vector<Term>& terms, const Record& r, const Covariates& cov = {},
                         const LearnOptions& o = {}) {
  require_record(g, r);
  if (terms.empty()) throw Refusal("a law needs at least one term");
  { Law probe{terms, std::vector<double>(terms.size(), 0.0)}; require_law(g, probe, cov); }
  const std::size_t n = g.size(), m = r.t.size() - 1, p = terms.size();
  const double c = detail::floor_of(r, o);
  std::vector<std::vector<std::vector<double>>> T(r.t.size(), std::vector<std::vector<double>>(p));
  for (std::size_t k = 0; k < r.t.size(); k++)
    for (std::size_t j = 0; j < p; j++) T[k][j] = evaluate(g, terms[j], r.u[k], cov);
  Rectangular X(n * m, p);
  std::vector<double> y(n * m);
  for (std::size_t k = 0; k < m; k++) {
    const double dt = r.t[k + 1] - r.t[k];
    for (std::size_t i = 0; i < n; i++) {
      const std::size_t row = k * n + i;
      const double w = detail::weight(0.5 * (r.u[k][i] + r.u[k + 1][i]), c, o);
      y[row] = w * (r.u[k + 1][i] - r.u[k][i]);
      for (std::size_t j = 0; j < p; j++) X(row, j) = w * dt * 0.5 * (T[k][j][i] + T[k + 1][j][i]);
    }
  }
  Estimate e;
  e.law.terms = terms;
  e.rows = n * m;
  std::vector<char> active(p, 1);
  e.method = "least squares on the integral form";
  if (o.threshold > 0) {
    // Scale each column and the target to unit root-mean-square, threshold there, and refit what is kept.
    Rectangular S(X.rows, p);
    std::vector<double> scale(p, 1.0), ys(y.size());
    double yr = 0;
    for (double v : y) yr += v * v;
    yr = std::sqrt(yr / static_cast<double>(y.size()));
    if (yr == 0) yr = 1;
    for (std::size_t j = 0; j < p; j++) {
      double s = 0;
      for (std::size_t i = 0; i < X.rows; i++) s += X(i, j) * X(i, j);
      scale[j] = std::sqrt(s / static_cast<double>(X.rows));
      if (scale[j] == 0) scale[j] = 1, active[j] = 0;
      for (std::size_t i = 0; i < X.rows; i++) S(i, j) = X(i, j) / scale[j];
    }
    for (std::size_t i = 0; i < y.size(); i++) ys[i] = y[i] / yr;
    induction::SparseOptions so;
    so.threshold = o.threshold;
    so.keep = o.keep;
    const induction::SparseFit sf = induction::stlsq(S, ys, so);
    for (std::size_t j = 0; j < p; j++) active[j] = active[j] && sf.coefficients[j] != 0;
    e.method = "thresholded least squares on the integral form (threshold " + std::to_string(o.threshold) + " of the change on scaled terms), refitted";
  }
  detail::solve(X, y, active, e);
  while (o.well_posed) {
    bool again = false;
    for (std::size_t j = 0; j < p; j++)
      if (active[j] && (terms[j].kind == Kind::diffusion || terms[j].kind == Kind::crowding) && e.law.coefficients[j] < 0) {
        active[j] = 0, again = true;
        e.dropped.push_back(std::string(name(terms[j].kind)) + ": found negative (" + std::to_string(e.law.coefficients[j]) + "), an ill-posed backward diffusion, and set aside");
      }
    if (!again) break;
    detail::solve(X, y, active, e);
  }
  return e;
}

// The weighted residuals of a run of the law from the record's first time, at every later time.
namespace detail {
inline bool run_residuals(const Graph& g, const Law& law, const Record& r, const Covariates& cov, const std::vector<std::vector<double>>& w,
                          std::vector<double>& out) {
  try {
    const std::vector<double> times(r.t.begin() + 1, r.t.end());
    const Run run = simulate(g, law, r.u[0], r.t[0], times, cov);
    out.clear();
    for (std::size_t k = 0; k < times.size(); k++)
      for (std::size_t i = 0; i < g.size(); i++) out.push_back(w[k][i] * (run.u[k][i] - r.u[k + 1][i]));
    return true;
  } catch (const Refusal&) {
    return false;
  }
}
}  // namespace detail

// Shooting: the coefficients whose run from the record's first time comes closest to it, from `start`.
inline Estimate shooting(const Graph& g, const Estimate& start, const Record& r, const Covariates& cov = {}, const LearnOptions& o = {}) {
  require_record(g, r);
  require_law(g, start.law, cov);
  const double c = detail::floor_of(r, o);
  std::vector<std::vector<double>> w(r.t.size() - 1, std::vector<double>(g.size()));
  for (std::size_t k = 1; k < r.t.size(); k++)
    for (std::size_t i = 0; i < g.size(); i++) w[k - 1][i] = detail::weight(r.u[k][i], c, o);
  // The coefficients that move are those not zero in the start (a sparse start stays sparse).
  std::vector<std::size_t> free;
  for (std::size_t j = 0; j < start.law.terms.size(); j++)
    if (start.law.coefficients[j] != 0) free.push_back(j);
  Estimate e = start;
  e.method = "shooting: runs from the first time fitted to the record by Levenberg–Marquardt, from the integral form";
  std::vector<double> res, trial;
  if (!detail::run_residuals(g, e.law, r, cov, w, res)) throw Refusal("the law of the integral form blows up when run from the record's first time");
  auto sse = [](const std::vector<double>& v) { double s = 0; for (double x : v) s += x * x; return s; };
  double best = sse(res);
  const std::size_t p = free.size(), N = res.size();
  Rectangular J(N, p);
  double lambda = -1;
  e.iterations = 0;
  e.converged = false;
  auto jacobian = [&]() -> bool {
    for (std::size_t q = 0; q < p; q++) {
      const std::size_t j = free[q];
      const double th = e.law.coefficients[j];
      const double se = j < start.standard_errors.size() ? start.standard_errors[j] : 0;
      double h = std::max(1e-6 * std::abs(th), 1e-3 * se);
      if (!(h > 0)) h = 1e-8;
      Law l = e.law;
      l.coefficients[j] = th + h;
      if (!detail::run_residuals(g, l, r, cov, w, trial)) {
        l.coefficients[j] = th - h, h = -h;
        if (!detail::run_residuals(g, l, r, cov, w, trial)) return false;
      }
      for (std::size_t i = 0; i < N; i++) J(i, q) = (trial[i] - res[i]) / h;
    }
    return true;
  };
  if (p > 0 && jacobian()) {
    for (; e.iterations < o.most; e.iterations++) {
      std::vector<double> diag(p, 0.0);
      for (std::size_t q = 0; q < p; q++)
        for (std::size_t i = 0; i < N; i++) diag[q] += J(i, q) * J(i, q);
      if (lambda < 0) lambda = 1e-3;
      bool accepted = false;
      for (int tries = 0; tries < 12 && !accepted; tries++) {
        // Solve [J; √(λ D)] δ = [−r; 0] by least squares.
        Rectangular A(N + p, p);
        std::vector<double> b(N + p, 0.0);
        for (std::size_t i = 0; i < N; i++) {
          for (std::size_t q = 0; q < p; q++) A(i, q) = J(i, q);
          b[i] = -res[i];
        }
        for (std::size_t q = 0; q < p; q++) A(N + q, q) = std::sqrt(lambda * std::max(diag[q], 1e-300));
        const std::vector<double> d = calx::linalg::least_squares(A, b);
        Law l = e.law;
        for (std::size_t q = 0; q < p; q++) {
          l.coefficients[free[q]] += d[q];
          const Kind kq = l.terms[free[q]].kind;
          if (o.well_posed && (kq == Kind::diffusion || kq == Kind::crowding) && l.coefficients[free[q]] < 0) l.coefficients[free[q]] = 0;  // projected
        }
        if (detail::run_residuals(g, l, r, cov, w, trial) && sse(trial) < best) {
          const double gain = (best - sse(trial)) / std::max(best, 1e-300);
          double step = 0, size = 0;
          for (std::size_t q = 0; q < p; q++) step += d[q] * d[q], size += l.coefficients[free[q]] * l.coefficients[free[q]];
          e.law = l, res = trial, best = sse(trial), lambda = std::max(lambda / 3, 1e-12), accepted = true;
          if (gain < 1e-10 || step <= 1e-20 * std::max(size, 1e-300)) e.converged = true;
        } else {
          lambda *= 4;
        }
      }
      if (!accepted) { e.converged = true; break; }  // no step lowers the misfit: a minimum, as far as the iteration can tell
      if (e.converged) break;
      if (!jacobian()) break;
    }
  } else if (p == 0) {
    e.converged = true;
  }
  // Standard errors from the last linearisation; the explained share is of the states, about their mean.
  if (p > 0 && jacobian()) {
    const calx::linalg::Qr f = calx::linalg::qr(J);
    const calx::linalg::Matrix C = calx::linalg::normal_inverse(f);
    const double dof = std::max(1.0, static_cast<double>(N) - static_cast<double>(f.rank));
    e.standard_errors.assign(e.law.terms.size(), 0.0);
    for (std::size_t q = 0; q < p; q++) e.standard_errors[free[q]] = std::sqrt(std::max(0.0, best / dof * C(q, q)));
  }
  double mean = 0, tss = 0;
  std::size_t count = 0;
  for (std::size_t k = 1; k < r.t.size(); k++)
    for (std::size_t i = 0; i < g.size(); i++) mean += w[k - 1][i] * r.u[k][i], count++;
  mean /= static_cast<double>(count);
  for (std::size_t k = 1; k < r.t.size(); k++)
    for (std::size_t i = 0; i < g.size(); i++) tss += (w[k - 1][i] * r.u[k][i] - mean) * (w[k - 1][i] * r.u[k][i] - mean);
  for (std::size_t q = 0; q < p; q++)
    if (e.law.coefficients[free[q]] == 0) e.dropped.push_back(std::string(name(e.law.terms[free[q]].kind)) + ": held at zero by shooting, which would have made it negative");
  e.rows = N;
  e.residual = std::sqrt(best / static_cast<double>(N));
  e.explained = tss > 0 ? 1 - best / tss : 0;
  e.terms = 0;
  for (double c : e.law.coefficients) e.terms += c != 0;
  return e;
}

// How a forecast fared at times kept aside, against carrying the last field forward unchanged.
struct Score {
  std::vector<double> t;                       // the times kept aside
  std::vector<std::vector<double>> predicted;  // the forecast at each
  double rmse = 0;                             // root-mean-square error of the density, all cells and times
  double persistence_rmse = 0;                 // the same for the last field carried forward
  double skill = 0;                            // 1 − SSE / SSE of persistence: above zero, better than no change
  double total_error = 0;                      // the forecast content against the record's, relative, at the last time
  double change_correlation = 0;               // between predicted and recorded change since the last time used
  double pattern_skill = 0;                    // the skill of the forecast rescaled, at each time, to the record's content:
                                               // where the quantity went, apart from how much of it there was
};

inline Score score(const Graph& g, const std::vector<std::vector<double>>& predicted, const Record& r, std::size_t from) {
  require_record(g, r);
  if (from + 1 >= r.t.size()) throw Refusal("no time is left after the last one used, to judge a forecast by");
  if (predicted.size() != r.t.size() - from - 1) throw Refusal("the forecast is not one field per time kept aside");
  Score s;
  double sse = 0, pse = 0, n = 0, sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0, shape = 0, pshape = 0;
  for (std::size_t k = from + 1, q = 0; k < r.t.size(); k++, q++) {
    s.t.push_back(r.t[k]);
    s.predicted.push_back(predicted[q]);
    require_size(g, predicted[q], "a forecast field");
    const double cp = content(g, predicted[q]), cr = content(g, r.u[k]), cl = content(g, r.u[from]);
    const double fp = cp > 0 ? cr / cp : 1, fl = cl > 0 ? cr / cl : 1;  // each rescaled to the record's content then
    for (std::size_t i = 0; i < g.size(); i++) {
      const double e = predicted[q][i] - r.u[k][i], pe = r.u[from][i] - r.u[k][i];
      sse += e * e, pse += pe * pe, n++;
      const double es = fp * predicted[q][i] - r.u[k][i], ps = fl * r.u[from][i] - r.u[k][i];
      shape += es * es, pshape += ps * ps;
      const double a = predicted[q][i] - r.u[from][i], b = r.u[k][i] - r.u[from][i];
      sa += a, sb += b, saa += a * a, sbb += b * b, sab += a * b;
    }
  }
  s.rmse = std::sqrt(sse / n);
  s.persistence_rmse = std::sqrt(pse / n);
  s.skill = pse > 0 ? 1 - sse / pse : 0;
  s.pattern_skill = pshape > 0 ? 1 - shape / pshape : 0;
  const double va = saa - sa * sa / n, vb = sbb - sb * sb / n;
  s.change_correlation = va > 0 && vb > 0 ? (sab - sa * sb / n) / std::sqrt(va * vb) : 0;
  const double rec = content(g, r.u.back()), pre = content(g, predicted.back());
  s.total_error = rec != 0 ? (pre - rec) / rec : 0;
  return s;
}

// A law's forecast from the record's time `from` to every later time.
inline std::vector<std::vector<double>> forecast(const Graph& g, const Law& law, const Record& r, std::size_t from, const Covariates& cov = {},
                                                 double* clipped = nullptr) {
  const std::vector<double> times(r.t.begin() + static_cast<long>(from) + 1, r.t.end());
  const Run run = simulate(g, law, r.u[from], r.t[from], times, cov);
  if (clipped) *clipped = run.clipped;
  return run.u;
}

// The baseline of a rate of its own for every cell: u_i(t) = (u_i(t_a) + c) e^{r_i (t − t_a)} − c, r_i from the first and
// the last times used, c a floor of one unit of density so that an empty cell has a rate. It has a parameter per cell,
// and no cell affects another.
struct OwnRates {
  std::vector<double> rate;
  double t0 = 0, floor = 1;
  std::vector<double> u0;
  std::vector<double> at(double t) const {
    std::vector<double> u(u0.size());
    for (std::size_t i = 0; i < u.size(); i++) u[i] = std::max(0.0, (u0[i] + floor) * std::exp(rate[i] * (t - t0)) - floor);
    return u;
  }
};

inline OwnRates own_rates(const Record& r, double floor = 1) {
  if (r.t.size() < 2) throw Refusal("a rate needs two times");
  OwnRates o;
  o.t0 = r.t.back(), o.floor = floor, o.u0 = r.u.back();
  const double span = r.t.back() - r.t.front();
  for (std::size_t i = 0; i < r.u.back().size(); i++)
    o.rate.push_back(std::log((std::max(r.u.back()[i], 0.0) + floor) / (std::max(r.u.front()[i], 0.0) + floor)) / span);
  return o;
}

}  // namespace kinesis::balance
