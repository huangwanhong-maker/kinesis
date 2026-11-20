// Laws learned from a record, used to evaluate sequences of interventions, with what the record cannot support said
// (decision 0009).
//
// A learned law is first held to the forecasts that need no law (tyche's random walk, historical mean and fitted
// first-order autoregression), one step ahead from rolling origins, by Diebold and Mariano's test with Harvey, Leybourne
// and Newbold's correction: a law that a baseline beats is reported so, and one not shown better than a baseline is not
// called better (the survey of evolutionary allocation, Part VIII §6, recommendation A.3).
//
// A record supports what it varied. Its support here is the range of each state and input it holds; a sequence that
// drives a state or an input outside it is flagged at the step and the variable, and its value is not taken as the law's
// word (Lucas 1976: a law fitted under one policy need not hold under another; tyche's estimates withheld outside their
// area of applicability are the spatial form of the same discipline). The ranges are a box, the least a support can be:
// a point inside the box but far from every record is not flagged, and says so.
//
// An ensemble of laws (posterior draws, bootstrap refits, rival forms, and the linear benchmark among them) evaluates two
// sequences of interventions. The order of two interventions matters under a law only through its state-dependent terms
// (under a linear, time-invariant law the effects add, and order enters only through timing), and those are the terms a
// record identifies worst. So the claim "the first sequence is better" is given a standing across the ensemble: it holds
// when every law says so, fails when every law says the opposite, is arguable when a weighted majority says so, and is
// undecided otherwise or when either sequence leaves the record's support; beside it the difference by law and its range.
//
// Bagged sparse laws (Fasel, Kutz, Brunton and Brunton 2022, "Ensemble-SINDy", Proc. R. Soc. A 478, 20210904):
// sequentially thresholded least squares refitted on resamples of the rows, each term's inclusion probability the share of
// refits that keep it, so that an order effect resting on a term kept one time in three is seen to rest on little.
#pragma once

#include <kinesis/induction/sparse.hpp>

#include <calx/random/stream.hpp>
#include <tyche/time_series/forecasts.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace kinesis::induction {

using kernel::Refusal;

// ---- a law held to the forecasts that need none -------------------------------------------------------------------------

struct Baseline {
  std::string name;
  double mse = 0, statistic = 0, p = 1;
  std::string verdict;  // "the law forecasts better", "the baseline forecasts better", "no difference shown"
};

struct Held {
  double mse = 0;
  std::size_t forecasts = 0;
  std::vector<Baseline> against;
  bool beaten = false, better_than_all = false;
  std::string verdict;
};

inline Held held_to_baselines(std::span<const double> y, std::size_t start, const tyche::time_series::Forecaster& law, double level = 0.05) {
  namespace ts = tyche::time_series;
  const auto L = ts::rolling(y, start, law);
  Held h;
  h.mse = L.mse(), h.forecasts = L.truth.size();
  const std::pair<const char*, ts::Forecaster> bases[] = {{"random walk", ts::random_walk}, {"historical mean", ts::historical_mean}, {"first-order autoregression", ts::ar1}};
  std::string unbeaten;
  h.better_than_all = true;
  for (const auto& [name, f] : bases) {
    const auto B = ts::rolling(y, start, f);
    const auto e = ts::diebold_mariano(L.truth, L.forecast, B.forecast);
    Baseline b;
    b.name = name, b.mse = B.mse(), b.statistic = e.statistic, b.p = e.p;
    if (e.p < level && e.statistic < 0) b.verdict = "the law forecasts better";
    else if (e.p < level && e.statistic > 0) b.verdict = "the baseline forecasts better", h.beaten = true;
    else b.verdict = "no difference shown";
    if (b.verdict != "the law forecasts better") {
      h.better_than_all = false;
      if (b.verdict == "no difference shown") unbeaten += (unbeaten.empty() ? "" : ", ") + std::string(name);
    }
    h.against.push_back(b);
  }
  if (h.beaten) {
    std::string by;
    for (const auto& b : h.against)
      if (b.verdict == "the baseline forecasts better") by += (by.empty() ? "" : ", ") + b.name;
    h.verdict = "beaten by the " + by;
  } else if (h.better_than_all) {
    h.verdict = "better than every baseline";
  } else {
    h.verdict = "not shown better than the " + unbeaten;
  }
  return h;
}

// ---- the record's support -------------------------------------------------------------------------------------------------

struct Support {
  std::vector<std::string> names;
  std::vector<double> lower, upper;
};

inline Support support_of(const std::vector<std::vector<double>>& rows, std::vector<std::string> names) {
  if (rows.empty()) throw Refusal("a support needs at least one record");
  const std::size_t d = rows[0].size();
  if (names.size() != d) throw Refusal("a name for each variable of the record");
  Support S;
  S.names = std::move(names);
  S.lower = rows[0], S.upper = rows[0];
  for (const auto& r : rows) {
    if (r.size() != d) throw Refusal("every record has the same variables");
    for (std::size_t j = 0; j < d; j++) S.lower[j] = std::min(S.lower[j], r[j]), S.upper[j] = std::max(S.upper[j], r[j]);
  }
  return S;
}

// What of a point lies outside the support's box (with a margin of each range's share): one message per variable.
inline std::vector<std::string> outside(const Support& S, const std::vector<double>& x, double margin = 0) {
  if (x.size() != S.names.size()) throw Refusal("a point with every variable of the support");
  std::vector<std::string> out;
  for (std::size_t j = 0; j < x.size(); j++) {
    const double pad = margin * (S.upper[j] - S.lower[j]);
    if (x[j] < S.lower[j] - pad || x[j] > S.upper[j] + pad)
      out.push_back(S.names[j] + " = " + std::to_string(x[j]) + " beyond [" + std::to_string(S.lower[j]) + ", " + std::to_string(S.upper[j]) + "]");
  }
  return out;
}

// ---- an ensemble of laws evaluating sequences ------------------------------------------------------------------------------

// A law as a map: the next state from the state and the input applied in that step.
using Step = std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)>;

struct Member {
  std::string name;
  double weight = 1;
  Step step;
};

struct Run {
  std::vector<std::vector<double>> states;  // x_0 … x_T
  std::vector<std::string> flags;           // where the run left the record's support
};

inline Run run(const Member& m, const std::vector<double>& x0, const std::vector<std::vector<double>>& inputs, const Support* states = nullptr,
               const Support* in = nullptr, double margin = 0) {
  Run r;
  r.states.push_back(x0);
  for (std::size_t t = 0; t < inputs.size(); t++) {
    if (in)
      for (const auto& f : outside(*in, inputs[t], margin)) r.flags.push_back("step " + std::to_string(t) + ": input " + f);
    r.states.push_back(m.step(r.states.back(), inputs[t]));
    if (states)
      for (const auto& f : outside(*states, r.states.back(), margin)) r.flags.push_back("step " + std::to_string(t + 1) + ": state " + f);
  }
  return r;
}

struct Ranking {
  std::string claim = "the first sequence is better";
  std::string standing;  // holds, fails, arguable, undecided
  std::string reason;
  std::vector<double> difference;  // by law: value(first) − value(second)
  double mean = 0, lower = 0, upper = 0, share_first = 0;
  bool immaterial = false;         // every law gives the two the same value: the order does not matter under any
  std::vector<std::string> flags;
};

inline Ranking rank(const std::vector<Member>& laws, const std::vector<double>& x0, const std::vector<std::vector<double>>& first,
                    const std::vector<std::vector<double>>& second, const std::function<double(const Run&)>& value, const Support* states = nullptr,
                    const Support* inputs = nullptr, double tolerance = 1e-12, double margin = 0) {
  if (laws.empty()) throw Refusal("an ensemble needs at least one law");
  Ranking R;
  double W = 0, for_first = 0, against = 0;
  R.lower = std::numeric_limits<double>::infinity(), R.upper = -std::numeric_limits<double>::infinity();
  R.immaterial = true;
  for (const auto& m : laws) {
    if (!(m.weight > 0)) throw Refusal("every law in an ensemble has a positive weight");
    const Run a = run(m, x0, first, states, inputs, margin), b = run(m, x0, second, states, inputs, margin);
    for (const auto& f : a.flags) R.flags.push_back(m.name + ", the first: " + f);
    for (const auto& f : b.flags) R.flags.push_back(m.name + ", the second: " + f);
    const double d = value(a) - value(b);
    R.difference.push_back(d);
    W += m.weight;
    R.mean += m.weight * d;
    R.lower = std::min(R.lower, d), R.upper = std::max(R.upper, d);
    if (d > tolerance) for_first += m.weight;
    if (d < -tolerance) against += m.weight;
    if (std::fabs(d) > tolerance) R.immaterial = false;
  }
  R.mean /= W;
  R.share_first = for_first / W;
  if (!R.flags.empty()) {
    R.standing = "undecided";
    R.reason = "out of the record's support (" + R.flags.front() + (R.flags.size() > 1 ? ", and " + std::to_string(R.flags.size() - 1) + " more" : std::string()) + ")";
  } else if (R.immaterial) {
    R.standing = "undecided";
    R.reason = "the order is immaterial under every law";
  } else if (for_first == W) {
    R.standing = "holds", R.reason = "every law says so";
  } else if (against == W) {
    R.standing = "fails", R.reason = "every law says the opposite";
  } else if (for_first > W / 2) {
    R.standing = "arguable", R.reason = "a weighted majority of the laws says so, not every one";
  } else {
    R.standing = "undecided", R.reason = "the laws disagree";
  }
  return R;
}

// ---- bagged sparse laws -------------------------------------------------------------------------------------------------

struct Bagged {
  std::vector<std::vector<double>> coefficients;  // one refit per resample
  std::vector<double> inclusion;                  // per term, the share of refits that keep it
  std::vector<double> median;                     // per term, the median coefficient over the refits
};

inline Bagged bagged(const Rectangular& theta, std::span<const double> y, std::size_t resamples, std::uint64_t seed, const SparseOptions& o = {}) {
  const std::size_t N = theta.rows, K = theta.cols;
  if (y.size() != N || N == 0) throw Refusal("a value for each row of the library");
  if (resamples == 0) throw Refusal("at least one resample");
  Bagged B;
  B.inclusion.assign(K, 0);
  for (std::size_t b = 0; b < resamples; b++) {
    calx::random::Stream r(seed, 0x62616767u /* "bagg" */, static_cast<std::uint32_t>(b), 0);
    Rectangular t(N, K);
    std::vector<double> v(N);
    for (std::size_t i = 0; i < N; i++) {
      const std::size_t row = static_cast<std::size_t>(r.below(N));
      for (std::size_t k = 0; k < K; k++) t(i, k) = theta(row, k);
      v[i] = y[row];
    }
    const auto fit = stlsq(t, v, o);
    B.coefficients.push_back(fit.coefficients);
    for (std::size_t k = 0; k < K; k++) B.inclusion[k] += fit.coefficients[k] != 0 ? 1.0 / static_cast<double>(resamples) : 0.0;
  }
  B.median.assign(K, 0);
  for (std::size_t k = 0; k < K; k++) {
    std::vector<double> c;
    for (const auto& f : B.coefficients) c.push_back(f[k]);
    std::nth_element(c.begin(), c.begin() + static_cast<std::ptrdiff_t>(c.size() / 2), c.end());
    B.median[k] = c[c.size() / 2];
  }
  return B;
}

}  // namespace kinesis::induction
