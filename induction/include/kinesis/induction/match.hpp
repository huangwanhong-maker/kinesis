// Matching models to a record: known models made templates, stochastic templates fitted, and candidates ranked
// by how well each continues a part of the record it was not fitted to.
//
// A model as a template. A model of a catalogue is a declaration with parameters. as_template() frees some or all
// of them: each freed parameter becomes an unknown number, and the model's structure is kept. Matching a record
// against a catalogue is then fitting each template (fit.hpp) and comparing.
//
// Stochastic templates. For dX_i = a_i dt + b_i dW_i with unknowns, the increments of the record give the drift
// by least squares on ΔX/h (the first Kramers–Moyal coefficient; Risken 1989, The Fokker–Planck Equation, 2nd
// ed., §4.1) and the size of the noise from what the drift leaves, E[(ΔX − aΔt)²] = b² h: an unknown s in
// b = s·g(x) is √⟨(ΔX − aΔt)² / (g² h)⟩. Both are of the first order in the step h of the record: a record
// sampled coarsely against the system's own rates biases them (Friedrich, Peinke, Sahimi and Tabar 2011, Phys.
// Rep. 506, 87–162, §3). kramers_moyal() gives the same two coefficients without a template, in bins.
//
// Ranking. Each candidate is fitted to the first part of the record and run on the rest: the error of one step
// from every held-out sample, and the horizon over which a run from the first held-out sample stays near the
// record. Candidates are ordered by Schwarz's criterion on the one-step errors, N ln(mean square error) +
// k ln N, k the coefficients fitted (Schwarz 1978, Ann. Statist. 6, 461–464), so that a model with more freedom
// must earn it. *Caveat.* On a chaotic record the horizon is bounded by the largest Lyapunov exponent whatever
// the model, and one step says little about the attractor; a ranking is among the candidates given, and says
// nothing of a model that was not proposed.
#pragma once

#include <kinesis/induction/fit.hpp>
#include <kinesis/induction/surrogate.hpp>
#include <kinesis/systems/dynamics.hpp>
#include <kinesis/systems/simulate.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace kinesis::induction {

// The model with the parameters named (or all, when none is named) made unknown numbers.
inline System as_template(const System& model, const std::vector<std::string>& free = {}) {
  System t = model;
  std::vector<std::string> freed;
  t.params.clear();
  for (const auto& p : model.params) {
    const bool chosen = p.name.arguments.empty() && (free.empty() || std::find(free.begin(), free.end(), p.name.name) != free.end());
    if (chosen) freed.push_back(p.name.name), t.unknowns.push_back({language::Unknown::Kind::number, p.name.name, nullptr});
    else t.params.push_back(p);
  }
  const auto open = [&](const TermPtr& x) {
    return terms::rewrite(x, [&](const TermPtr& s) -> TermPtr {
      return s->kind == terms::Kind::symbol && s->primes == 0 && std::find(freed.begin(), freed.end(), s->name) != freed.end() ? terms::hole(s->name) : nullptr;
    });
  };
  for (auto& e : t.change.equations) e.value = open(e.value);
  for (auto& a : t.change.actions) a.value = open(a.value);
  for (auto& p : t.params) p.value = open(p.value);
  t.constraint.conserved.clear();
  return t;
}

// ---- stochastic records ----------------------------------------------------------------------------------

struct Moments {
  std::vector<double> x, drift, diffusion, count;  // in each bin: its centre, ⟨Δx⟩/h, ⟨Δx²⟩/h, the samples
};

// The first two Kramers–Moyal coefficients of a scalar series sampled at the step h, in `bins` bins of [lower, upper].
inline Moments kramers_moyal(std::span<const double> series, double h, std::size_t bins, double lower, double upper) {
  Moments m;
  m.x.resize(bins), m.drift.assign(bins, 0.0), m.diffusion.assign(bins, 0.0), m.count.assign(bins, 0.0);
  for (std::size_t b = 0; b < bins; b++) m.x[b] = lower + (upper - lower) * (static_cast<double>(b) + 0.5) / static_cast<double>(bins);
  for (std::size_t k = 0; k + 1 < series.size(); k++) {
    const double f = (series[k] - lower) / (upper - lower);
    if (f < 0 || f >= 1) continue;
    const std::size_t b = static_cast<std::size_t>(f * static_cast<double>(bins));
    const double dx = series[k + 1] - series[k];
    m.drift[b] += dx, m.diffusion[b] += dx * dx, m.count[b] += 1;
  }
  for (std::size_t b = 0; b < bins; b++)
    if (m.count[b] > 0) m.drift[b] /= m.count[b] * h, m.diffusion[b] /= m.count[b] * h;
  return m;
}

// The unknowns of a stochastic differential with diagonal noise, from a record of it.
inline Fit fit_stochastic(const System& tmpl, const Trajectory& record) {
  if (systems::form_of(tmpl) != systems::Form::stochastic || systems::discrete_time(tmpl)) throw kernel::Refusal(tmpl.name + " is not a stochastic differential");
  const systems::Scope sc = systems::scope(tmpl);
  const detail::Working w = detail::open(tmpl, sc);
  const std::size_t U = w.coefficients.size(), n = record.x.size();
  if (n < 10) throw kernel::Refusal("a fit needs a record of at least ten samples");
  const double h = record.t[1] - record.t[0];
  terms::Frame frame;
  std::vector<std::size_t> columns;
  for (const auto& st : sc.states) {
    const auto it = std::find(record.names.begin(), record.names.end(), st);
    if (it == record.names.end()) throw kernel::Refusal("the record has no column '" + st + "', a state of " + tmpl.name);
    frame.add(st), columns.push_back(static_cast<std::size_t>(it - record.names.begin()));
  }
  const std::size_t first_param = frame.size();
  for (const auto& p : sc.params) frame.add(p);
  const int slot_t = frame.add(sc.time);
  const std::size_t first_unknown = frame.size();
  for (const auto& c : w.coefficients) frame.add(c.symbol);
  std::vector<double> values(frame.size(), 0.0);
  std::copy(sc.param_values.begin(), sc.param_values.end(), values.begin() + static_cast<std::ptrdiff_t>(first_param));
  const auto load = [&](std::size_t k) {
    for (std::size_t c = 0; c < columns.size(); c++) values[c] = record.x[k][columns[c]];
    values[static_cast<std::size_t>(slot_t)] = record.t[k];
  };
  // The drift and the one noise coefficient of each equation.
  struct Law {
    std::size_t column;
    TermPtr drift, noise;
  };
  std::vector<Law> laws;
  const std::string dt = "d" + sc.time;
  for (const auto& e : tmpl.change.equations) {
    if (e.name.size() < 2 || e.name[0] != 'd') throw kernel::Refusal("a stochastic differential is written dX = a * " + dt + " + b * dW");
    const std::string target = e.name.substr(1);
    const auto v = sc.vectors.find(target);
    const std::size_t count = v == sc.vectors.end() ? 1 : v->second.size();
    for (std::size_t k = 0; k < count; k++) {
      const std::string state = v == sc.vectors.end() ? target : v->second[k];
      const auto it = std::find(sc.states.begin(), sc.states.end(), state);
      if (it == sc.states.end()) throw kernel::Refusal("'" + state + "' is not a state of " + tmpl.name);
      const TermPtr rhs = w.opened(sc.resolve(v == sc.vectors.end() ? e.value : systems::component(sc, e.value, k, "d")));
      Law law{static_cast<std::size_t>(it - sc.states.begin()), terms::derivative(rhs, dt), terms::number(0)};
      std::size_t driven = 0;
      for (const auto& src : sc.sources) {
        const TermPtr b = terms::derivative(rhs, "d" + src.name);
        if (!terms::is_number(b, 0)) law.noise = b, driven++;
      }
      if (driven > 1) throw kernel::Refusal("the noise of " + state + " comes from several sources: only diagonal noise is fitted");
      laws.push_back(std::move(law));
    }
  }
  Fit out;
  std::vector<double> solution(U, 0.0);
  // The drift: affine in its unknowns, by least squares on the increments over the step.
  std::vector<char> in_drift(U, 0);
  for (const auto& law : laws)
    for (std::size_t u = 0; u < U; u++)
      if (terms::depends_on(law.drift, w.coefficients[u].symbol, 0)) in_drift[u] = 1;
  std::vector<std::size_t> drift_unknowns;
  for (std::size_t u = 0; u < U; u++)
    if (in_drift[u]) drift_unknowns.push_back(u);
  std::vector<terms::Program> rest;
  std::vector<std::vector<terms::Program>> coefficient(laws.size());
  for (std::size_t e = 0; e < laws.size(); e++) {
    TermPtr r0 = laws[e].drift;
    for (std::size_t u : drift_unknowns) {
      const TermPtr c = terms::derivative(laws[e].drift, w.coefficients[u].symbol);
      for (std::size_t v : drift_unknowns)
        if (terms::depends_on(c, w.coefficients[v].symbol, 0)) throw kernel::Refusal("the drift of " + tmpl.name + " is not linear in its unknowns");
      coefficient[e].push_back(terms::compile(c, frame.resolver()));
      r0 = terms::substitute(r0, w.coefficients[u].symbol, 0, terms::number(0));
    }
    rest.push_back(terms::compile(r0, frame.resolver()));
  }
  const std::size_t rows = n - 1;
  if (!drift_unknowns.empty()) {
    Rectangular theta(rows * laws.size(), drift_unknowns.size());
    std::vector<double> y(rows * laws.size());
    for (std::size_t k = 0; k < rows; k++) {
      load(k);
      for (std::size_t e = 0; e < laws.size(); e++) {
        for (std::size_t j = 0; j < drift_unknowns.size(); j++) theta(e * rows + k, j) = coefficient[e][j](values);
        y[e * rows + k] = (record.x[k + 1][columns[laws[e].column]] - record.x[k][columns[laws[e].column]]) / h - rest[e](values);
      }
    }
    const std::vector<double> xi = calx::linalg::least_squares(theta, y);
    for (std::size_t j = 0; j < drift_unknowns.size(); j++) solution[drift_unknowns[j]] = xi[j];
  }
  // The noise: b = s g with one unknown s, from what the drift leaves of the increments.
  std::copy(solution.begin(), solution.end(), values.begin() + static_cast<std::ptrdiff_t>(first_unknown));
  double mismatch = 0, size = 0;
  for (std::size_t e = 0; e < laws.size(); e++) {
    std::size_t unknown = U;
    for (std::size_t u = 0; u < U; u++)
      if (terms::depends_on(laws[e].noise, w.coefficients[u].symbol, 0)) {
        if (unknown != U || in_drift[u]) throw kernel::Refusal("the noise of a state has one unknown of its own, which the drift does not hold");
        unknown = u;
      }
    const terms::Program drift = terms::compile(laws[e].drift, frame.resolver());
    terms::Program shape;
    if (unknown != U) {
      const TermPtr g = terms::derivative(laws[e].noise, w.coefficients[unknown].symbol);
      if (terms::depends_on(g, w.coefficients[unknown].symbol, 0)) throw kernel::Refusal("the noise is not proportional to its unknown");
      shape = terms::compile(g, frame.resolver());
    }
    double sum = 0;
    for (std::size_t k = 0; k < rows; k++) {
      load(k);
      const double dx = record.x[k + 1][columns[laws[e].column]] - record.x[k][columns[laws[e].column]], left = dx - drift(values) * h;
      mismatch += left * left, size += dx * dx;
      if (unknown != U) {
        const double g = shape(values);
        sum += left * left / (g * g * h);
      }
    }
    if (unknown != U) solution[unknown] = std::sqrt(sum / static_cast<double>(rows));
  }
  out.residual = std::sqrt(mismatch / std::max(size, 1e-300));  // near 1: the increments of a diffusion are mostly noise
  out.samples = rows * laws.size();
  out.method = "Kramers-Moyal coefficients: least squares for the drift, the quadratic variation for the noise";
  out.system = detail::filled(tmpl, w, solution, out);
  return out;
}

// ---- ranking ---------------------------------------------------------------------------------------------

struct Candidate {
  std::string name;
  bool fitted = false;
  std::string refused;           // why it could not be fitted or run
  Fit fit;
  double one_step = std::numeric_limits<double>::infinity();   // rms error of one step on the held-out part, against the rms increment
  std::size_t horizon = 0;       // held-out samples a free run stays near the record
  double criterion = std::numeric_limits<double>::infinity();  // Schwarz's, on the one-step errors
};

struct MatchOptions {
  double train = 0.7;            // the share of the record fitted to
  double tolerance = 0.3;        // of the horizon, against the spread of the record
  FitOptions fit;
};

inline std::vector<Candidate> match(const std::vector<System>& templates, const Trajectory& record, const MatchOptions& o = {}) {
  const std::size_t n = record.x.size(), cut = static_cast<std::size_t>(o.train * static_cast<double>(n));
  if (cut < 5 || n - cut < 5) throw kernel::Refusal("the record is too short to keep a part aside");
  Trajectory head, tail;
  head.names = tail.names = record.names;
  head.t.assign(record.t.begin(), record.t.begin() + static_cast<std::ptrdiff_t>(cut)), head.x.assign(record.x.begin(), record.x.begin() + static_cast<std::ptrdiff_t>(cut));
  tail.t.assign(record.t.begin() + static_cast<std::ptrdiff_t>(cut), record.t.end()), tail.x.assign(record.x.begin() + static_cast<std::ptrdiff_t>(cut), record.x.end());
  const double h = record.t[1] - record.t[0];
  std::vector<Candidate> out;
  for (const System& tmpl : templates) {
    Candidate c;
    c.name = tmpl.name;
    try {
      c.fit = fit(tmpl, head, o.fit);
      const System& model = c.fit.system;
      const systems::Dynamics d = systems::dynamics(systems::form_of(model) == systems::Form::lagrangian ? systems::varied(model, model.name) : model);
      // The record's columns in the order of the model's components; a model of higher order needs its
      // derivatives, which the record may not hold.
      std::vector<std::size_t> column;
      for (const auto& st : d.states) {
        const auto it = std::find(record.names.begin(), record.names.end(), st);
        if (it == record.names.end()) throw kernel::Refusal("the record has no column '" + st + "' to start the model from");
        column.push_back(static_cast<std::size_t>(it - record.names.begin()));
      }
      const auto state_at = [&](const State& row) {
        State x(column.size());
        for (std::size_t i = 0; i < column.size(); i++) x[i] = row[column[i]];
        return x;
      };
      systems::Options so;
      so.step = d.discrete ? 1 : h / std::max(1.0, std::ceil(h / 0.01)), so.every = 100000000;
      double error = 0, size = 0;
      for (std::size_t k = 0; k + 1 < tail.x.size(); k++) {
        so.t0 = d.discrete ? 0 : tail.t[k], so.t1 = d.discrete ? 1 : tail.t[k + 1];
        const State next = systems::simulate(d, state_at(tail.x[k]), so).x.back(), truth = state_at(tail.x[k + 1]), before = state_at(tail.x[k]);
        for (std::size_t i = 0; i < next.size(); i++) error += (next[i] - truth[i]) * (next[i] - truth[i]), size += (truth[i] - before[i]) * (truth[i] - before[i]);
      }
      const double N = static_cast<double>((tail.x.size() - 1) * column.size());
      c.one_step = std::sqrt(error / std::max(size, 1e-300));
      c.criterion = N * std::log(std::max(error / N, 1e-300)) + static_cast<double>(c.fit.parameters) * std::log(N);
      // A free run from the first held-out sample.
      std::vector<State> run, truth;
      State x = state_at(tail.x[0]);
      for (std::size_t k = 0; k + 1 < tail.x.size(); k++) {
        so.t0 = d.discrete ? 0 : tail.t[k], so.t1 = d.discrete ? 1 : tail.t[k + 1];
        try {
          x = systems::simulate(d, x, so).x.back();
        } catch (const kernel::Refusal&) {
          break;
        }
        run.push_back(x), truth.push_back(state_at(tail.x[k + 1]));
      }
      if (!run.empty()) {
        std::vector<State> all;
        for (const auto& row : tail.x) all.push_back(state_at(row));
        truth.resize(run.size());
        c.horizon = induction::horizon(run, truth, o.tolerance);
      }
      c.fitted = true;
    } catch (const std::exception& e) {
      c.refused = e.what();
    }
    out.push_back(std::move(c));
  }
  std::stable_sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b) { return a.criterion < b.criterion; });
  return out;
}

}  // namespace kinesis::induction
