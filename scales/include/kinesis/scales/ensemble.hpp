// What the uncertainty of a system's inputs does to what it shows, and which input it answers to most (decision 0007).
//
// Factors. Every input that states an uncertainty is a factor: a value declared with a range, or estimated with a
// standard error. The values of one quantity are one factor, drawn together across its units at one quantile (comonotone:
// a range declared for a quantity states doubt about its level, the same doubt in every unit), and so are the values of a
// driver's series and of a parameter. A declared range is read as uniform over it, a standard error as a normal law about
// the estimate; both are assumptions, and are said to be.
//
// Ensembles. The factors are drawn on a Latin hypercube (McKay, Beckman and Conover 1979, Technometrics 21, 239–245; tyche's
// uncertainty), the model is run once per draw, and every value asked for is summarised at each time and unit by its
// quantiles over the draws (Hyndman and Fan's type 7). A draw the model refuses (an equilibrium not reached, a value that
// cannot be) is counted and left out, never filled in.
//
// Screening. Which factor moves an outcome most is answered by Morris's elementary effects (Morris 1991, Technometrics 33,
// 161–174), μ* the mean absolute effect of moving the factor across half its range: the order in which declarations are
// worth making better.
#pragma once

#include <cmath>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <vector>

#include "kinesis/scales/run.hpp"
#include "tyche/descriptive/summary.hpp"
#include "tyche/probability/distributions.hpp"
#include "tyche/uncertainty/uncertainty.hpp"

namespace kinesis::scales {

struct Factor {
  enum What { value, series, parameter } what = value;
  std::string name;
  Source source = Source::declared;
  std::string law;  // "uniform over the declared range", "normal about the estimate", or both
};

inline const char* name(Factor::What w) {
  switch (w) {
    case Factor::value: return "value";
    case Factor::series: return "series";
    default: return "parameter";
  }
}

namespace detail {

inline double drawn(const Value& v, double u) {
  u = std::clamp(u, 1e-12, 1 - 1e-12);
  if (std::isfinite(v.lower) && std::isfinite(v.upper) && v.upper > v.lower) return v.lower + u * (v.upper - v.lower);
  if (v.sigma > 0 && std::isfinite(v.value)) return v.value + v.sigma * tyche::probability::Normal::standard_quantile(u);
  return centre(v);
}

inline void describe(Factor& f, const Value& v, bool& ranged, bool& normal) {
  if (std::isfinite(v.lower) && std::isfinite(v.upper) && v.upper > v.lower) ranged = true;
  else if (v.sigma > 0) normal = true;
  if (v.source == Source::declared) f.source = Source::declared;
}

inline std::string law_words(bool ranged, bool normal) {
  if (ranged && normal) return "uniform over a declared range, or normal about an estimate";
  return ranged ? "uniform over the declared range" : "normal about the estimate";
}

}  // namespace detail

// Every input that states an uncertainty, one factor per quantity, series or parameter.
inline std::vector<Factor> factors(const Inputs& in) {
  std::vector<Factor> out;
  for (const auto& [n, vs] : in.values) {
    bool ranged = false, normal = false;
    Factor f{Factor::value, n, Source::estimated, ""};
    for (const Value& v : vs)
      if (uncertain(v)) detail::describe(f, v, ranged, normal);
    if (ranged || normal) f.law = detail::law_words(ranged, normal), out.push_back(f);
  }
  for (const auto& [n, s] : in.series) {
    bool ranged = false, normal = false;
    Factor f{Factor::series, n, Source::estimated, ""};
    for (const auto& row : s.at)
      for (const Value& v : row)
        if (uncertain(v)) detail::describe(f, v, ranged, normal);
    if (ranged || normal) f.law = detail::law_words(ranged, normal), out.push_back(f);
  }
  for (const auto& [n, v] : in.parameters) {
    if (!uncertain(v)) continue;
    bool ranged = false, normal = false;
    Factor f{Factor::parameter, n, Source::estimated, ""};
    detail::describe(f, v, ranged, normal);
    f.law = detail::law_words(ranged, normal);
    out.push_back(f);
  }
  return out;
}

// The inputs with each factor set at its quantile u_k: every value fixed, so that a run uses exactly it.
inline Inputs at_quantiles(const Inputs& in, const std::vector<Factor>& fs, std::span<const double> u) {
  Inputs out = in;
  for (std::size_t k = 0; k < fs.size(); k++) {
    const Factor& f = fs[k];
    const auto fix = [&](Value& v) {
      if (!uncertain(v)) return;
      v.value = detail::drawn(v, u[k]);
      v.lower = v.upper = missing_value;
      v.sigma = 0;
    };
    if (f.what == Factor::value)
      for (Value& v : out.values[f.name]) fix(v);
    else if (f.what == Factor::series)
      for (auto& row : out.series[f.name].at)
        for (Value& v : row) fix(v);
    else
      fix(out.parameters[f.name]);
  }
  return out;
}

struct Band {
  std::vector<std::vector<double>> low, middle, high;  // [time][unit]
};

struct EnsembleOptions {
  std::size_t draws = 64;
  std::uint64_t seed = 1;
  double low = 0.05, high = 0.95;
  std::vector<std::string> outputs;  // the quantities summarised; empty: every one the run keeps
};

struct Ensemble {
  std::vector<double> times;
  std::vector<Factor> factors;
  std::map<std::string, Band> bands;
  std::map<std::string, std::string> level;
  std::size_t draws = 0, refused = 0;
  std::vector<std::string> notes;
};

// The model run over draws of its uncertain inputs, each output summarised by its quantiles.
inline Ensemble ensemble(const Compiled& c, const Inputs& in, Options o, const EnsembleOptions& e) {
  Ensemble out;
  out.factors = factors(in);
  if (!e.outputs.empty()) o.record = e.outputs;
  if (out.factors.empty()) {
    const Run r = simulate(c, in, o);
    out.times = r.times;
    out.level = r.level;
    out.draws = 1;
    for (const auto& [n, v] : r.values) out.bands[n] = Band{v, v, v};
    out.notes.push_back("no input states an uncertainty, so one run stands for the ensemble");
    return out;
  }
  if (e.draws < 2) throw Refusal("an ensemble needs at least two draws");
  const auto design = tyche::uncertainty::latin_hypercube(e.draws, out.factors.size(), e.seed);
  std::map<std::string, std::vector<std::vector<std::vector<double>>>> kept;  // name → draw → time → unit
  for (std::size_t d = 0; d < e.draws; d++) {
    std::vector<double> u(out.factors.size());
    for (std::size_t k = 0; k < u.size(); k++) u[k] = design(d, k);
    try {
      const Run r = simulate(c, at_quantiles(in, out.factors, u), o);
      if (r.counts.unconverged) {
        out.refused++;
        continue;
      }
      if (out.times.empty()) out.times = r.times, out.level = r.level;
      for (const auto& [n, v] : r.values) kept[n].push_back(v);
    } catch (const Refusal&) {
      out.refused++;
    }
  }
  out.draws = e.draws - out.refused;
  if (out.draws < 2) throw Refusal("fewer than two draws of the ensemble ran; the declared ranges may hold values the model cannot take");
  for (const auto& [n, draws] : kept) {
    Band b;
    const std::size_t T = draws.front().size(), U = draws.front().front().size();
    b.low.assign(T, std::vector<double>(U));
    b.middle = b.high = b.low;
    std::vector<double> xs(draws.size());
    for (std::size_t t = 0; t < T; t++)
      for (std::size_t uu = 0; uu < U; uu++) {
        for (std::size_t d = 0; d < draws.size(); d++) xs[d] = draws[d][t][uu];
        b.low[t][uu] = tyche::descriptive::quantile(xs, e.low);
        b.middle[t][uu] = tyche::descriptive::quantile(xs, 0.5);
        b.high[t][uu] = tyche::descriptive::quantile(xs, e.high);
      }
    out.bands[n] = std::move(b);
  }
  out.notes.push_back("each factor is drawn at one quantile across all its units; ranges are read as uniform, standard errors as normal");
  if (out.refused) out.notes.push_back(std::to_string(out.refused) + " draws were refused by the model and left out");
  return out;
}

struct Screening {
  std::vector<Factor> factors;
  std::vector<double> mu_star, mu, sigma;  // per factor, in the outcome's units, for a move across half its range
  std::string outcome;
  std::size_t evaluations = 0;
};

// Which uncertain input moves an outcome most: the outcome is the sum over the units of `output` at the run's last time
// (or its value at one unit when `unit` names one), screened by Morris's elementary effects.
inline Screening screen(const Compiled& c, const Inputs& in, Options o, const std::string& output, const std::string& unit = "",
                        std::size_t trajectories = 10, std::uint64_t seed = 1) {
  Screening s;
  s.factors = factors(in);
  if (s.factors.empty()) throw Refusal("no input states an uncertainty, so there is nothing to screen");
  const auto lv = c.level_of.find(output);
  if (lv == c.level_of.end()) throw Refusal("the model has no quantity " + output);
  o.record = {output};
  o.report = {o.t1};
  const std::size_t which = unit.empty() ? 0 : c.hierarchy.unit(lv->second, unit);
  s.outcome = unit.empty() ? "the sum of " + output + " over its units at the last time" : output + " at " + unit + " at the last time";
  const tyche::uncertainty::Model f = [&](std::span<const double> u) {
    const Run r = simulate(c, at_quantiles(in, s.factors, u), o);
    const auto& v = r.values.at(output).back();
    if (!unit.empty()) return v[which];
    double total = 0;
    for (double x : v) total += x;
    return total;
  };
  std::vector<tyche::uncertainty::Quantile> inputs(s.factors.size(), [](double u) { return u; });
  const auto m = tyche::uncertainty::morris(f, inputs, trajectories, 4, seed);
  s.mu_star = m.mu_star;
  s.mu = m.mu;
  s.sigma = m.sigma;
  s.evaluations = m.evaluations;
  return s;
}

}  // namespace kinesis::scales
