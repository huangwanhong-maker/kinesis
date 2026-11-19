// What a detector knows of a record and of the space its states live in.
//
// A record is a sequence of states at increasing times. A state need not be a number: it may be a field over units, a
// composition (a population by age, an economy by sector), a distribution given by a histogram, or a covariance matrix
// (the coupling of channels of a recording). A detector measures fluctuations, and fluctuations are measured in linear
// coordinates with an inner product; so the space declares a chart, a map from a native state to coordinates in which
// Euclidean geometry is the space's own (the infrastructure's survey *Critical transitions in structured systems*,
// Part II §11):
//
//   identity     a vector or a field: each part of each unit is a coordinate.
//   composition  each unit's parts are a composition, read in isometric log-ratio coordinates (Egozcue et al. 2003,
//                Math. Geol. 35, 279–300), so that the Euclidean distance is Aitchison's; with `with_total`, the logarithm
//                of the unit's total is a further coordinate, so that what changes the size and what changes the shape are
//                kept apart.
//   histogram    each unit's parts are counts in ordered bins with declared edges, read as the quantile function at
//                declared probabilities: on the real line the Wasserstein distance W₂ between two distributions is the L²
//                distance between their quantile functions, which the coordinates approximate on the grid.
//   covariance   each unit's parts are an m × m symmetric positive definite matrix C, read as the upper triangle of
//                log(R^{-1/2} C R^{-1/2}), the off-diagonal entries multiplied by √2, for a reference R (the log-Euclidean
//                mean of the record's first states). At R the chart is an isometry of the affine-invariant metric
//                (Pennec, Fillard and Ayache 2006, Int. J. Comput. Vis. 66, 41–66): the Euclidean norm of the coordinates
//                of C is the affine-invariant distance from R to C. Away from R it is a first approximation.
//
// The units may carry a neighbourhood (a weighted graph, for spatial statistics) and a hierarchy (kinesis::scales'
// hierarchy, whose finest level is the units, for indicators at every level). A quantity is extensive (summed when
// units are grouped) or intensive (averaged).
#pragma once

#include <kinesis/kernel/answer.hpp>
#include <kinesis/scales/hierarchy.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/linalg/rectangular.hpp>
#include <tyche/compositions/balances.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace kinesis::transitions {

using calx::linalg::Matrix;
using calx::linalg::Rectangular;
using kernel::Refusal;

// ---- records ----------------------------------------------------------------------------------------------------

struct Record {
  std::vector<double> times;                 // strictly increasing
  std::vector<std::vector<double>> states;   // one native state a time, all of one length
};

inline void check(const Record& r) {
  if (r.times.size() != r.states.size()) throw Refusal("a record needs one state for each time");
  if (r.times.size() < 2) throw Refusal("a record needs two states or more");
  for (std::size_t t = 1; t < r.times.size(); t++)
    if (!(r.times[t] > r.times[t - 1])) throw Refusal("the times of a record must increase");
  for (const auto& s : r.states)
    if (s.size() != r.states[0].size()) throw Refusal("the states of a record must have one length");
}

// The step of a regularly sampled record, or nothing when the steps differ by more than a relative 10⁻⁶.
inline std::optional<double> regular_step(const std::vector<double>& times) {
  if (times.size() < 2) return std::nullopt;
  const double dt = (times.back() - times.front()) / static_cast<double>(times.size() - 1);
  for (std::size_t t = 1; t < times.size(); t++)
    if (std::abs(times[t] - times[t - 1] - dt) > 1e-6 * std::abs(dt)) return std::nullopt;
  return dt;
}

inline double require_regular(const Record& r, const std::string& tool) {
  const auto dt = regular_step(r.times);
  if (!dt) throw Refusal(tool + " needs a regularly sampled record; this one's steps differ");
  return *dt;
}

// A record of one scalar quantity.
inline Record scalar_record(const std::vector<double>& times, const std::vector<double>& x) {
  Record r;
  r.times = times;
  for (double v : x) r.states.push_back({v});
  return r;
}

// ---- the state space --------------------------------------------------------------------------------------------

enum class Chart { identity, composition, histogram, covariance };
enum class Kind { intensive, extensive };

inline const char* name(Chart c) {
  switch (c) {
    case Chart::identity: return "identity";
    case Chart::composition: return "composition";
    case Chart::histogram: return "histogram";
    default: return "covariance";
  }
}

struct StateSpace {
  std::size_t units = 1;                     // the native state is units × parts, unit by unit
  std::size_t parts = 1;
  Chart chart = Chart::identity;
  Kind kind = Kind::intensive;
  bool with_total = false;                   // composition, histogram: the log of each unit's total is a coordinate too
  std::vector<double> bin_edges;             // histogram: parts + 1 increasing edges
  std::vector<double> probabilities;         // histogram: the quantile grid, in (0, 1)
  std::size_t reference = 0;                 // covariance: states averaged for R (0: a tenth of the record, at least 2)
  std::vector<std::string> unit_names, part_names;
  std::vector<scales::Level::Edge> edges;    // the neighbourhood of the units, when there is no hierarchy
  std::optional<scales::Hierarchy> hierarchy;  // its finest level is the units; its edges are their neighbourhood

  std::size_t native_size() const { return units * parts; }
  std::size_t matrix_order() const {
    const auto m = static_cast<std::size_t>(std::llround(std::sqrt(static_cast<double>(parts))));
    return m;
  }
  // Coordinates a unit has in the chart.
  std::size_t coordinates_per_unit() const {
    switch (chart) {
      case Chart::identity: return parts;
      case Chart::composition: return parts - 1 + (with_total ? 1 : 0);
      case Chart::histogram: return probabilities.size() + (with_total ? 1 : 0);
      default: {
        const std::size_t m = matrix_order();
        return m * (m + 1) / 2;
      }
    }
  }
  std::size_t dimension() const { return units * coordinates_per_unit(); }

  const std::vector<scales::Level::Edge>& neighbourhood() const { return hierarchy ? hierarchy->level(0).edges : edges; }
  std::string unit_name(std::size_t u) const { return u < unit_names.size() ? unit_names[u] : "unit " + std::to_string(u + 1); }
};

inline void check(const StateSpace& s) {
  if (s.units == 0 || s.parts == 0) throw Refusal("a state space needs at least one unit and one part");
  if (!s.unit_names.empty() && s.unit_names.size() != s.units) throw Refusal("a state space names every unit or none");
  if (!s.part_names.empty() && s.part_names.size() != s.parts) throw Refusal("a state space names every part or none");
  if (s.chart == Chart::composition && s.parts < 2) throw Refusal("a composition needs two parts or more");
  if (s.chart == Chart::histogram) {
    if (s.bin_edges.size() != s.parts + 1) throw Refusal("a histogram needs parts + 1 bin edges");
    for (std::size_t k = 1; k < s.bin_edges.size(); k++)
      if (!(s.bin_edges[k] > s.bin_edges[k - 1])) throw Refusal("the bin edges of a histogram must increase");
    if (s.probabilities.empty()) throw Refusal("a histogram is read at declared probabilities; none were given");
    for (double p : s.probabilities)
      if (!(p > 0 && p < 1)) throw Refusal("the probabilities of a histogram's quantiles lie in (0, 1)");
  }
  if (s.chart == Chart::covariance) {
    const std::size_t m = s.matrix_order();
    if (m * m != s.parts || m < 1) throw Refusal("a covariance state needs m × m parts a unit");
  }
  if (s.hierarchy && s.hierarchy->size(0) != s.units)
    throw Refusal("the finest level of the hierarchy must be the state's units");
  for (const auto& e : s.neighbourhood())
    if (e.a >= s.units || e.b >= s.units || e.a == e.b || !(e.conductance >= 0))
      throw Refusal("an edge of the neighbourhood joins units that are not two of the state's, or has a negative weight");
}

// A record in the chart's coordinates: one row a time, one column a coordinate, with the unit of each coordinate.
struct Charted {
  Rectangular x;                         // times × coordinates
  std::vector<std::size_t> unit_of;      // the unit of each coordinate
  std::vector<std::string> names;        // each coordinate's name
  std::vector<double> times;
};

namespace detail {

// The symmetric function f(S) = V f(Λ) Vᵀ of a symmetric matrix.
template <class F>
Matrix symmetric_function(const Matrix& S, F&& f) {
  const auto e = calx::linalg::symmetric_eigen(S);
  const std::size_t n = S.n;
  Matrix out(n);
  for (std::size_t k = 0; k < n; k++) {
    const double fk = f(e.values[k]);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) out(i, j) += e.vectors(i, k) * fk * e.vectors(j, k);
  }
  return out;
}

inline Matrix unit_matrix(const std::vector<double>& state, std::size_t unit, std::size_t m) {
  Matrix C(m);
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = 0; j < m; j++) C(i, j) = 0.5 * (state[unit * m * m + i * m + j] + state[unit * m * m + j * m + i]);
  return C;
}

inline void require_positive_definite(const Matrix& C, double t, std::size_t unit) {
  const auto e = calx::linalg::symmetric_eigen(C);
  if (!(e.values.front() > 0))
    throw Refusal("the covariance of unit " + std::to_string(unit + 1) + " at time " + std::to_string(t) + " is not positive definite");
}

// The quantile function of a histogram at p, linear within each bin (a density constant in each bin).
inline double histogram_quantile(const std::vector<double>& counts, std::size_t offset, std::size_t bins,
                                 const std::vector<double>& edges, double total, double p) {
  const double target = p * total;
  double cumulative = 0;
  for (std::size_t k = 0; k < bins; k++) {
    const double c = counts[offset + k];
    if (cumulative + c >= target && c > 0) {
      const double share = (target - cumulative) / c;
      return edges[k] + share * (edges[k + 1] - edges[k]);
    }
    cumulative += c;
  }
  return edges[bins];
}

}  // namespace detail

// The record read in the space's chart.
inline Charted chart(const StateSpace& s, const Record& r) {
  check(s);
  check(r);
  if (r.states[0].size() != s.native_size())
    throw Refusal("the record's states have " + std::to_string(r.states[0].size()) + " values; the space declares " +
                  std::to_string(s.units) + " units × " + std::to_string(s.parts) + " parts");
  const std::size_t T = r.times.size(), per = s.coordinates_per_unit(), d = s.dimension();
  Charted out;
  out.x = Rectangular(T, d);
  out.times = r.times;
  out.unit_of.resize(d);
  out.names.resize(d);
  for (std::size_t u = 0; u < s.units; u++)
    for (std::size_t k = 0; k < per; k++) {
      out.unit_of[u * per + k] = u;
      std::string coordinate;
      switch (s.chart) {
        case Chart::identity: coordinate = k < s.part_names.size() ? s.part_names[k] : "part " + std::to_string(k + 1); break;
        case Chart::composition:
          coordinate = (s.with_total && k == 0) ? "log total" : "balance " + std::to_string(k + (s.with_total ? 0 : 1));
          break;
        case Chart::histogram:
          coordinate = (s.with_total && k == 0) ? "log total" : "quantile " + std::to_string(k + (s.with_total ? 0 : 1));
          break;
        default: coordinate = "log-covariance " + std::to_string(k + 1);
      }
      out.names[u * per + k] = (s.units > 1 ? s.unit_name(u) + ": " : std::string()) + coordinate;
    }

  switch (s.chart) {
    case Chart::identity:
      for (std::size_t t = 0; t < T; t++)
        for (std::size_t c = 0; c < d; c++) out.x(t, c) = r.states[t][c];
      break;
    case Chart::composition: {
      const auto basis = tyche::compositions::pivot_basis(s.parts);
      for (std::size_t t = 0; t < T; t++)
        for (std::size_t u = 0; u < s.units; u++) {
          std::vector<double> parts(r.states[t].begin() + static_cast<std::ptrdiff_t>(u * s.parts),
                                    r.states[t].begin() + static_cast<std::ptrdiff_t>((u + 1) * s.parts));
          double total = 0;
          for (double v : parts) {
            if (!(v > 0))
              throw Refusal("a part of " + s.unit_name(u) + " at time " + std::to_string(r.times[t]) +
                            " is not positive; a composition is read in log-ratios, so a zero must be replaced by a declared value first");
            total += v;
          }
          const auto z = tyche::compositions::ilr(parts, basis);
          std::size_t c = u * per;
          if (s.with_total) out.x(t, c++) = std::log(total);
          for (double v : z) out.x(t, c++) = v;
        }
      break;
    }
    case Chart::histogram:
      for (std::size_t t = 0; t < T; t++)
        for (std::size_t u = 0; u < s.units; u++) {
          double total = 0;
          for (std::size_t k = 0; k < s.parts; k++) {
            const double v = r.states[t][u * s.parts + k];
            if (!(v >= 0)) throw Refusal("a count of " + s.unit_name(u) + " at time " + std::to_string(r.times[t]) + " is negative");
            total += v;
          }
          if (!(total > 0)) throw Refusal("the histogram of " + s.unit_name(u) + " at time " + std::to_string(r.times[t]) + " is empty");
          std::size_t c = u * per;
          if (s.with_total) out.x(t, c++) = std::log(total);
          for (double p : s.probabilities)
            out.x(t, c++) = detail::histogram_quantile(r.states[t], u * s.parts, s.parts, s.bin_edges, total, p);
        }
      break;
    case Chart::covariance: {
      const std::size_t m = s.matrix_order();
      const std::size_t n_ref = s.reference > 0 ? std::min(s.reference, T) : std::max<std::size_t>(2, T / 10);
      for (std::size_t u = 0; u < s.units; u++) {
        // R = exp(mean log C) over the reference states, the log-Euclidean mean (Arsigny et al. 2007).
        Matrix mean_log(m);
        for (std::size_t t = 0; t < n_ref; t++) {
          const Matrix C = detail::unit_matrix(r.states[t], u, m);
          detail::require_positive_definite(C, r.times[t], u);
          const Matrix L = detail::symmetric_function(C, [](double v) { return std::log(v); });
          for (std::size_t i = 0; i < m * m; i++) mean_log.a[i] += L.a[i] / static_cast<double>(n_ref);
        }
        const Matrix R = detail::symmetric_function(mean_log, [](double v) { return std::exp(v); });
        const Matrix R_half_inv = detail::symmetric_function(R, [](double v) { return 1 / std::sqrt(v); });
        for (std::size_t t = 0; t < T; t++) {
          const Matrix C = detail::unit_matrix(r.states[t], u, m);
          detail::require_positive_definite(C, r.times[t], u);
          Matrix W = R_half_inv * C * R_half_inv;
          for (std::size_t i = 0; i < m; i++)
            for (std::size_t j = i + 1; j < m; j++) W(i, j) = W(j, i) = 0.5 * (W(i, j) + W(j, i));
          const Matrix L = detail::symmetric_function(W, [](double v) { return std::log(v); });
          std::size_t c = u * per;
          for (std::size_t i = 0; i < m; i++)
            for (std::size_t j = i; j < m; j++) out.x(t, c++) = i == j ? L(i, i) : std::sqrt(2.0) * L(i, j);
        }
      }
      break;
    }
  }
  return out;
}

// The columns of a charted record that belong to some units, or one column.
inline Rectangular columns(const Rectangular& x, const std::vector<std::size_t>& which) {
  Rectangular out(x.rows, which.size());
  for (std::size_t t = 0; t < x.rows; t++)
    for (std::size_t k = 0; k < which.size(); k++) out(t, k) = x(t, which[k]);
  return out;
}

inline std::vector<double> column(const Rectangular& x, std::size_t c) {
  std::vector<double> out(x.rows);
  for (std::size_t t = 0; t < x.rows; t++) out[t] = x(t, c);
  return out;
}

inline Rectangular rows(const Rectangular& x, std::size_t from, std::size_t to) {
  Rectangular out(to - from, x.cols);
  for (std::size_t t = from; t < to; t++)
    for (std::size_t c = 0; c < x.cols; c++) out(t - from, c) = x(t, c);
  return out;
}

}  // namespace kinesis::transitions
