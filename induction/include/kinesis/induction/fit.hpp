// Filling the unknowns of a template from a record: the one mechanism by which a model is matched, discovered or
// learned (the idea of 2026-11-12, §5.4).
//
// A template is a declaration with unknowns: numbers (?a), and sparse combinations of a library of terms (?f,
// declared `unknown f : sparse poly(x, y, 3)`). fit() returns the same declaration with the unknowns replaced by
// what the record gives, so that the result is a system of the language, and runs, and is analysed, as any other.
//
// The equations a record must satisfy. For a flow, each equation q⁽ᵏ⁾ = rhs is held at every sample, with the
// derivatives of the record on the left (differences of second order, taken k times; a record that holds q' as a
// column is believed). For a map, the left side is the next sample. For a system given by an action, the
// equations are those of its stationary paths (systems/lagrangian.hpp), so that an unknown inside an action is
// fitted through the Euler–Lagrange equations it implies.
//
// Linear unknowns. Each unknown is made a symbol, and the right side is differentiated in it. When no derivative
// holds an unknown, every equation is affine in them, rhs = r₀ + Σ c_u u, and the unknowns solve a linear
// least-squares problem over all equations and samples, each equation weighed by the size of its left side. The
// coefficients of sparse unknowns are thresholded (sparse.hpp); numbers never are. With `weak`, equations of the
// first order are tested against windows instead, and the regression takes no derivative of the record; the
// differences of the record are still taken to weigh the equations and to report the residual.
//
// Unknowns that are not linear are found by the Levenberg–Marquardt iteration (Marquardt 1963, J. Soc. Indust.
// Appl. Math. 11, 431–441) on the same residuals, with the Jacobian derived from the terms, from the values 1.
// It finds a local minimum, and says how many iterations it took; another start may find another.
//
// What a fit does not say. The residual measures how well the equations hold on the record, not whether the
// template is right: a template that cannot express the law still returns its best coefficients. Unknown
// functions (`unknown V : function`) are not filled here.
#pragma once

#include <kinesis/induction/sparse.hpp>
#include <kinesis/language/declaration.hpp>
#include <kinesis/systems/lagrangian.hpp>
#include <kinesis/systems/model.hpp>

#include <calx/linalg/dense.hpp>

#include <map>
#include <string>
#include <vector>

namespace kinesis::induction {

using language::System;

struct FitOptions {
  double threshold = 0.05;      // of the coefficients of sparse unknowns
  bool weak = false;            // test equations of the first order against windows
  std::size_t half = 20, stride = 5;
  std::size_t most = 100;       // iterations, when the unknowns are not linear
};

struct Fit {
  System system;                           // the template with its unknowns filled
  std::map<std::string, double> numbers;   // the numbers found
  std::map<std::string, TermPtr> terms;    // the sparse unknowns, as terms
  std::size_t parameters = 0;              // coefficients found that are not zero
  double residual = 0;                     // the root-mean-square mismatch of the equations, against that of their left sides
  std::size_t samples = 0;                 // rows of the regression
  bool linear = true;
  std::size_t iterations = 0;
  std::string method;
};

namespace detail {

// An unknown coefficient of a template: a number, or the k-th coefficient of a sparse unknown.
struct Coefficient {
  std::string symbol;   // its name as a symbol of the working terms
  std::string unknown;  // the unknown it belongs to
  int term = -1;        // its place in the library, or −1 for a number
};

struct Working {
  std::vector<Coefficient> coefficients;
  std::map<std::string, std::vector<TermPtr>> libraries;
  // The term with each unknown replaced by symbols.
  TermPtr opened(const TermPtr& t) const {
    return terms::rewrite(t, [&](const TermPtr& x) -> TermPtr {
      if (x->kind != terms::Kind::hole) return nullptr;
      const auto lib = libraries.find(x->name);
      if (lib == libraries.end()) return terms::symbol("$" + x->name);
      std::vector<TermPtr> parts;
      for (std::size_t k = 0; k < lib->second.size(); k++) parts.push_back(terms::symbol("$" + x->name + "." + std::to_string(k)) * lib->second[k]);
      return terms::sum(std::move(parts));
    });
  }
};

inline Working open(const System& tmpl, const systems::Scope& sc) {
  Working w;
  std::set<std::string> declared;
  for (const auto& u : tmpl.unknowns) {
    declared.insert(u.name);
    if (u.kind == language::Unknown::Kind::function) throw kernel::Refusal("the unknown function ?" + u.name + " is not filled here: only numbers and sparse combinations of a library are");
    if (u.kind == language::Unknown::Kind::number) {
      w.coefficients.push_back({"$" + u.name, u.name, -1});
    } else {
      auto& lib = w.libraries[u.name];
      for (const auto& t : library_of(sc.resolve(u.library))) lib.push_back(t);
      for (std::size_t k = 0; k < lib.size(); k++) w.coefficients.push_back({"$" + u.name + "." + std::to_string(k), u.name, static_cast<int>(k)});
    }
  }
  // An unknown that is used and not declared is a number.
  std::set<std::string> used;
  for (const auto& e : tmpl.change.equations) terms::holes(e.value, used);
  for (const auto& a : tmpl.change.actions) terms::holes(a.value, used);
  for (const auto& n : used)
    if (!declared.count(n)) w.coefficients.push_back({"$" + n, n, -1});
  return w;
}

// The template with its unknowns given their values.
inline System filled(const System& tmpl, const Working& w, const std::vector<double>& values, Fit& fit) {
  std::map<std::string, TermPtr> by_name;
  for (const auto& [name, lib] : w.libraries) {
    std::vector<TermPtr> parts;
    for (std::size_t k = 0; k < w.coefficients.size(); k++)
      if (w.coefficients[k].unknown == name && values[k] != 0) parts.push_back(terms::number(values[k]) * lib[static_cast<std::size_t>(w.coefficients[k].term)]);
    by_name[name] = terms::sum(std::move(parts));
    fit.terms[name] = by_name[name];
  }
  for (std::size_t k = 0; k < w.coefficients.size(); k++) {
    fit.parameters += values[k] != 0;
    if (w.coefficients[k].term < 0) by_name[w.coefficients[k].unknown] = terms::number(values[k]), fit.numbers[w.coefficients[k].unknown] = values[k];
  }
  const auto fill = [&](const TermPtr& t) {
    return terms::rewrite(t, [&](const TermPtr& x) -> TermPtr {
      if (x->kind != terms::Kind::hole) return nullptr;
      const auto it = by_name.find(x->name);
      return it == by_name.end() ? nullptr : it->second;
    });
  };
  System out = tmpl;
  out.unknowns.clear();
  for (auto& e : out.change.equations) e.value = fill(e.value);
  for (auto& a : out.change.actions) a.value = fill(a.value);
  for (auto& c : out.constraint.conserved) c = fill(c);
  for (auto& ob : out.difference.observations) ob.value = fill(ob.value);
  return out;
}

inline double rms(const std::vector<double>& v) {
  double s = 0;
  for (double x : v) s += x * x;
  return std::sqrt(s / static_cast<double>(std::max<std::size_t>(v.size(), 1)));
}

}  // namespace detail

namespace detail {
inline Fit fit_equations(const System& equations, const System& tmpl, const Working& w, const Trajectory& record, const FitOptions& o);
}

inline Fit fit(const System& tmpl, const Trajectory& record, const FitOptions& o = {}) {
  const systems::Form form = systems::form_of(tmpl);
  const detail::Working w = detail::open(tmpl, systems::scope(tmpl));
  if (form == systems::Form::lagrangian) {
    // The unknowns of an action are made symbols, the action is varied with them as constants, and they are
    // fitted through the equations of its stationary paths; they are returned in the action.
    System opened = tmpl;
    opened.unknowns.clear();
    for (auto& a : opened.change.actions) a.value = w.opened(a.value);
    opened.constraint.conserved.clear();
    return detail::fit_equations(systems::varied(opened, tmpl.name), tmpl, w, record, o);
  }
  if (form != systems::Form::flow && form != systems::Form::map) throw kernel::Refusal("a template of equations or of an action is fitted here; " + tmpl.name + " is a " + systems::name(form) + " system");
  return detail::fit_equations(tmpl, tmpl, w, record, o);
}

inline Fit detail::fit_equations(const System& equations, const System& tmpl, const Working& w, const Trajectory& record, const FitOptions& o) {
  const systems::Scope sc = systems::scope(equations);
  const std::size_t U = w.coefficients.size(), n = record.x.size();
  if (n < 5) throw kernel::Refusal("a fit needs a record of at least five samples");

  // The columns of the record, and its derivatives as the equations need them.
  std::map<std::pair<std::string, unsigned>, std::vector<double>> column;
  const auto record_column = [&](const std::string& name) -> const std::vector<double>* {
    const auto it = std::find(record.names.begin(), record.names.end(), name);
    if (it == record.names.end()) return nullptr;
    auto& c = column[{name, 0u}];
    if (c.empty()) c = record.column(static_cast<std::size_t>(it - record.names.begin()));
    return &c;
  };
  struct Law {
    std::string state;
    unsigned order;
    TermPtr rhs;
  };
  std::vector<Law> laws;
  for (const auto& e : equations.change.equations) {
    const auto v = sc.vectors.find(e.name);
    const std::size_t count = v == sc.vectors.end() ? 1 : v->second.size();
    for (std::size_t k = 0; k < count; k++) laws.push_back({v == sc.vectors.end() ? e.name : v->second[k], e.primes, w.opened(sc.resolve(v == sc.vectors.end() ? e.value : systems::component(sc, e.value, k)))});
  }
  for (const auto& st : sc.states)
    if (!record_column(st)) throw kernel::Refusal("the record has no column '" + st + "', a state of " + tmpl.name);
  const bool discrete = sc.discrete;
  std::size_t rows = discrete ? n - 1 : n;
  if (!discrete) {
    for (const auto& law : laws)
      for (unsigned k = 1; k <= law.order; k++) {
        auto& c = column[{law.state, k}];
        if (!c.empty()) continue;
        const std::string primed = law.state + std::string(k, '\'');
        const auto it = std::find(record.names.begin(), record.names.end(), primed);
        if (it != record.names.end()) {
          c = record.column(static_cast<std::size_t>(it - record.names.begin()));
          continue;
        }
        Trajectory one;
        one.t = record.t;
        for (double v : column[{law.state, k - 1}]) one.x.push_back({v});
        for (const auto& d : derivatives(one)) c.push_back(d[0]);
      }
  }

  // The frame of the working terms: the columns, the parameters, the time, the unknowns.
  terms::Frame frame;
  std::vector<const std::vector<double>*> columns;
  for (const auto& [key, values] : column) frame.add(key.first, key.second), columns.push_back(&values);
  const std::size_t first_param = frame.size();
  for (const auto& p : sc.params) frame.add(p);
  const int slot_t = frame.add(sc.time);
  const std::size_t first_unknown = frame.size();
  for (const auto& c : w.coefficients) frame.add(c.symbol);
  std::vector<double> values(frame.size(), 0.0);
  std::copy(sc.param_values.begin(), sc.param_values.end(), values.begin() + static_cast<std::ptrdiff_t>(first_param));
  const auto load = [&](std::size_t k) {
    for (std::size_t c = 0; c < columns.size(); c++) values[c] = (*columns[c])[k];
    values[static_cast<std::size_t>(slot_t)] = record.t[k];
  };
  const auto compiled = [&](const TermPtr& t, const std::string& what) {
    try {
      return terms::compile(t, frame.resolver());
    } catch (const kernel::Refusal& e) {
      throw kernel::Refusal(what + ": " + e.what());
    }
  };

  // Each law: its left side on the record, and the derivatives of its right side in the unknowns.
  const std::size_t E = laws.size();
  std::vector<std::vector<double>> left(E);
  std::vector<std::vector<TermPtr>> gradient(E, std::vector<TermPtr>(U));
  Fit out;
  for (std::size_t e = 0; e < E; e++) {
    if (discrete) {
      const auto& c = column[{laws[e].state, 0u}];
      left[e].assign(c.begin() + 1, c.end());
    } else {
      left[e] = column[{laws[e].state, laws[e].order}];
    }
    for (std::size_t u = 0; u < U; u++) {
      gradient[e][u] = terms::derivative(laws[e].rhs, w.coefficients[u].symbol);
      for (const auto& c : w.coefficients)
        if (terms::depends_on(gradient[e][u], c.symbol, 0)) out.linear = false;
    }
  }
  std::vector<double> weight(E);
  for (std::size_t e = 0; e < E; e++) weight[e] = 1 / std::max(detail::rms(left[e]), 1e-300);
  std::vector<double> solution(U, 0.0);
  out.samples = rows * E;

  if (U == 0) {
    out.method = "nothing to fit";
  } else if (out.linear) {
    // rhs = r₀ + Σ c_u u: the columns c_u and the remainder r₀, on the record.
    std::vector<terms::Program> rest;
    std::vector<std::vector<terms::Program>> coefficient(E);
    for (std::size_t e = 0; e < E; e++) {
      TermPtr r0 = laws[e].rhs;
      for (const auto& c : w.coefficients) r0 = terms::substitute(r0, c.symbol, 0, terms::number(0));
      rest.push_back(compiled(r0, "the equation of " + laws[e].state));
      for (std::size_t u = 0; u < U; u++) coefficient[e].push_back(compiled(gradient[e][u], "the equation of " + laws[e].state));
    }
    const bool weak = o.weak && !discrete && std::all_of(laws.begin(), laws.end(), [](const Law& l) { return l.order == 1; });
    Rectangular theta;
    std::vector<double> y;
    if (weak) {
      // For each law, the columns c_u and r₀ against the windows; the left side is −∫ φ' x dt.
      for (std::size_t e = 0; e < E; e++) {
        Rectangular wide(n, U + 1);
        for (std::size_t k = 0; k < n; k++) {
          load(k);
          for (std::size_t u = 0; u < U; u++) wide(k, u) = coefficient[e][u](values);
          wide(k, U) = rest[e](values);
        }
        Trajectory one;
        one.t = record.t;
        for (double v : column[{laws[e].state, 0u}]) one.x.push_back({v});
        const Weak tested = weak_form(one, wide, o.half, o.stride);
        const std::size_t windows = tested.matrix.rows, base = theta.rows;
        Rectangular grown(base + windows, U);
        std::copy(theta.a.begin(), theta.a.end(), grown.a.begin());
        // Windows are weighed by the size of the left side of this law across them.
        const double scale = 1 / std::max(detail::rms(tested.targets[0]), 1e-300);
        for (std::size_t r = 0; r < windows; r++) {
          for (std::size_t u = 0; u < U; u++) grown(base + r, u) = scale * tested.matrix(r, u);
          y.push_back(scale * (tested.targets[0][r] - tested.matrix(r, U)));
        }
        theta = std::move(grown);
      }
      out.samples = theta.rows;
    } else {
      theta = Rectangular(rows * E, U);
      y.resize(rows * E);
      for (std::size_t k = 0; k < rows; k++) {
        load(k);
        for (std::size_t e = 0; e < E; e++) {
          for (std::size_t u = 0; u < U; u++) theta(e * rows + k, u) = weight[e] * coefficient[e][u](values);
          y[e * rows + k] = weight[e] * (left[e][k] - rest[e](values));
        }
      }
    }
    SparseOptions so;
    so.threshold = o.threshold;
    so.keep.assign(U, 0);
    bool sparse = false;
    for (std::size_t u = 0; u < U; u++) so.keep[u] = w.coefficients[u].term < 0, sparse = sparse || w.coefficients[u].term >= 0;
    const SparseFit sf = stlsq(theta, y, so);
    solution = sf.coefficients;
    out.iterations = sf.rounds;
    out.method = std::string(weak ? "the weak form, " : "") + (sparse ? "sequentially thresholded least squares" : "least squares");
  } else {
    // Levenberg–Marquardt on the residuals rhs − left, with the Jacobian from the terms.
    std::vector<terms::Program> right;
    std::vector<std::vector<terms::Program>> jacobian(E);
    for (std::size_t e = 0; e < E; e++) {
      right.push_back(compiled(laws[e].rhs, "the equation of " + laws[e].state));
      for (std::size_t u = 0; u < U; u++) jacobian[e].push_back(compiled(gradient[e][u], "the equation of " + laws[e].state));
    }
    std::fill(solution.begin(), solution.end(), 1.0);
    const auto cost = [&](const std::vector<double>& p) {
      std::copy(p.begin(), p.end(), values.begin() + static_cast<std::ptrdiff_t>(first_unknown));
      double c = 0;
      for (std::size_t k = 0; k < rows; k++) {
        load(k);
        for (std::size_t e = 0; e < E; e++) {
          const double r = weight[e] * (right[e](values) - left[e][k]);
          c += r * r;
        }
      }
      return c;
    };
    double mu = 1e-3, current = cost(solution);
    for (out.iterations = 0; out.iterations < o.most; out.iterations++) {
      calx::linalg::Matrix A(U), b(U);
      std::copy(solution.begin(), solution.end(), values.begin() + static_cast<std::ptrdiff_t>(first_unknown));
      std::vector<double> row(U);
      for (std::size_t k = 0; k < rows; k++) {
        load(k);
        for (std::size_t e = 0; e < E; e++) {
          const double r = weight[e] * (right[e](values) - left[e][k]);
          for (std::size_t u = 0; u < U; u++) row[u] = weight[e] * jacobian[e][u](values);
          for (std::size_t u = 0; u < U; u++) {
            b(u, 0) -= row[u] * r;
            for (std::size_t v = 0; v < U; v++) A(u, v) += row[u] * row[v];
          }
        }
      }
      bool improved = false;
      for (int attempt = 0; attempt < 30 && !improved; attempt++) {
        calx::linalg::Matrix damped = A;
        for (std::size_t u = 0; u < U; u++) damped(u, u) += mu * (A(u, u) > 0 ? A(u, u) : 1.0);
        std::vector<double> trial = solution;
        try {
          const calx::linalg::Matrix step = calx::linalg::solve(damped, b);
          for (std::size_t u = 0; u < U; u++) trial[u] += step(u, 0);
        } catch (const std::exception&) {
          mu *= 10;
          continue;
        }
        const double c = cost(trial);
        if (std::isfinite(c) && c < current) {
          const bool settled = current - c < 1e-14 * (1 + current);
          solution = trial, current = c, mu = std::max(mu / 10, 1e-12), improved = true;
          if (settled) out.iterations = o.most;
        } else {
          mu *= 10;
        }
      }
      if (!improved) break;
    }
    out.method = "Levenberg-Marquardt";
  }

  // The residual of the equations, with the values found.
  {
    std::copy(solution.begin(), solution.end(), values.begin() + static_cast<std::ptrdiff_t>(first_unknown));
    double mismatch = 0, size = 0;
    for (std::size_t e = 0; e < E; e++) {
      const terms::Program right = compiled(laws[e].rhs, "the equation of " + laws[e].state);
      for (std::size_t k = 0; k < rows; k++) {
        load(k);
        const double r = weight[e] * (right(values) - left[e][k]);
        mismatch += r * r, size += weight[e] * weight[e] * left[e][k] * left[e][k];
      }
    }
    out.residual = std::sqrt(mismatch / std::max(size, 1e-300));
  }
  out.system = detail::filled(tmpl, w, solution, out);
  return out;
}

// The template that asks for a sparse law of every state: x_i' = ?f_i, each a combination of the monomials of
// the record's names up to `degree`. Fitting it is the discovery of a law with no law proposed.
inline System discovery_template(const Trajectory& record, unsigned degree, const std::string& name = "discovered", bool discrete = false) {
  System s;
  s.name = name;
  language::States st;
  st.space.name = "real";
  for (const auto& n : record.names) st.names.push_back({n, {}});
  s.states.push_back(st);
  std::vector<TermPtr> arguments;
  for (const auto& n : record.names) arguments.push_back(terms::symbol(n));
  arguments.push_back(terms::number(degree));
  for (std::size_t i = 0; i < record.names.size(); i++) {
    const std::string f = "f" + std::to_string(i + 1);
    s.unknowns.push_back({language::Unknown::Kind::sparse, f, terms::call("poly", arguments)});
    s.change.equations.push_back({record.names[i], 1, {}, terms::hole(f), {}});
  }
  if (discrete) s.structure.timed = true, s.structure.time = "n", s.structure.time_domain.kind = language::Domain::Kind::natural, s.structure.time_domain.name = "natural";
  return s;
}

}  // namespace kinesis::induction
