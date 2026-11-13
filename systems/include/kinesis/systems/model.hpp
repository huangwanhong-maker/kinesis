// What a declaration means before it runs: its form, and the names its terms may use.
//
// The form of a system is read from what its `change` holds and how its `structure` reads it (form_of):
//
//   rules                                  a rewriting system: its histories are derivations
//   a Hamiltonian or jumps                 an open quantum system
//   an action, histories stationary        a Lagrangian system: its paths make the action stationary
//   an action, histories weight            a statistical field: its configurations are weighted
//   equations with a source of noise       a stochastic system
//   equations, time in the naturals        a map
//   equations                              a flow
//   parts                                  a system composed of others
//   = derivation source                    a system derived from another
//
// A scope resolves the names of a declaration once: the time, the scalar components of the states (a state in
// real^n is n components X1 … Xn), the parameters with their values, the functions (parameters with arguments,
// which are inlined where they are applied), the inputs, and the sources of noise. A name that a term uses and
// the scope does not hold is an error of the script, reported with the place of the item.
#pragma once

#include <kinesis/language/declaration.hpp>
#include <kinesis/terms/program.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace kinesis::systems {

using language::Script;
using language::ScriptError;
using language::System;
using terms::TermPtr;

using State = std::vector<double>;

enum class Form { empty, flow, map, stochastic, rules, lagrangian, field, quantum, composite, derived };

inline const char* name(Form f) {
  switch (f) {
    case Form::flow: return "flow";
    case Form::map: return "map";
    case Form::stochastic: return "stochastic";
    case Form::rules: return "rules";
    case Form::lagrangian: return "lagrangian";
    case Form::field: return "field";
    case Form::quantum: return "quantum";
    case Form::composite: return "composite";
    case Form::derived: return "derived";
    default: return "empty";
  }
}

inline bool discrete_time(const System& s) {
  return s.structure.timed && (s.structure.time_domain.kind == language::Domain::Kind::natural || s.structure.time_domain.kind == language::Domain::Kind::integer);
}

inline Form form_of(const System& s) {
  if (s.derived()) return Form::derived;
  if (!s.parts.empty()) return Form::composite;
  if (!s.change.rules.empty()) return Form::rules;
  if (!s.change.hamiltonians.empty() || !s.change.jumps.empty()) return Form::quantum;
  if (!s.change.actions.empty()) return s.structure.stationary.empty() ? Form::field : Form::lagrangian;
  if (!s.change.equations.empty()) {
    for (const auto& n : s.noises)
      if (!n.law) return Form::stochastic;
    return discrete_time(s) ? Form::map : Form::flow;
  }
  return Form::empty;
}

// A source of noise: a process that drives the system (a Wiener process for a stochastic differential; a draw at
// each step for a map), with its kind and the numbers of its law.
struct Source {
  std::string name;
  std::string kind;            // wiener, normal, uniform, levy, cauchy, poisson, exponential
  std::vector<double> numbers;
};

// The law a state starts with: `noise X(0) ~ normal(0, 1)`. `data` is a law known only by samples.
struct InitialLaw {
  std::string state;
  std::string kind;
  std::vector<double> numbers;
};

struct Scope {
  const System* system = nullptr;
  std::string time = "t";
  bool discrete = false;
  std::vector<std::string> states;                             // the scalar components, in the order declared
  std::map<std::string, std::vector<std::string>> vectors;    // a vector's name → its components (states and sources)
  std::vector<std::string> params;
  std::vector<double> param_values;
  std::map<std::string, const language::Binding*> functions;
  std::vector<std::string> inputs;
  std::vector<double> input_values;
  std::vector<Source> sources;
  std::vector<InitialLaw> initial;

  bool is_state(const std::string& n) const { return std::find(states.begin(), states.end(), n) != states.end(); }
  int param_index(const std::string& n) const {
    const auto it = std::find(params.begin(), params.end(), n);
    return it == params.end() ? -1 : static_cast<int>(it - params.begin());
  }
  bool declares(const std::string& n) const {
    if (n == time || is_state(n) || param_index(n) >= 0 || vectors.count(n)) return true;
    if (std::find(inputs.begin(), inputs.end(), n) != inputs.end()) return true;
    for (const auto& s : sources)
      if (s.name == n) return true;
    return false;
  }

  // The term with the functions of the scope inlined and the constants pi and e given their values.
  TermPtr resolve(const TermPtr& t, int depth = 0) const {
    if (depth > 32) throw kernel::Refusal("the functions of '" + system->name + "' apply one another without end");
    return terms::rewrite(t, [&](const TermPtr& x) -> TermPtr {
      if (x->kind == terms::Kind::call) {
        const auto f = functions.find(x->name);
        if (f != functions.end() && f->second->name.arguments.size() == x->args.size()) {
          TermPtr body = f->second->value;
          // Formals are replaced together: first by placeholders no script can write, then by the arguments.
          for (std::size_t i = 0; i < x->args.size(); i++) body = terms::substitute(body, f->second->name.arguments[i], 0, terms::symbol("#" + std::to_string(i)));
          for (std::size_t i = 0; i < x->args.size(); i++) body = terms::substitute(body, "#" + std::to_string(i), 0, resolve(x->args[i], depth + 1));
          return resolve(body, depth + 1);
        }
      }
      if (x->kind == terms::Kind::symbol && x->primes == 0 && !declares(x->name)) {
        if (x->name == "pi") return terms::number(std::numbers::pi);
        if (x->name == "e") return terms::number(std::numbers::e);
      }
      return nullptr;
    });
  }

  // The value of a term of numbers and parameters, or an error of the script at `at`.
  double constant(const TermPtr& t, language::Place at, const std::string& what) const {
    try {
      return terms::evaluate(resolve(t), [&](const std::string& n, unsigned primes) -> double {
        const int i = param_index(n);
        if (i < 0 || primes) throw kernel::Refusal("'" + n + std::string(primes, '\'') + "' is not a parameter");
        return param_values[static_cast<std::size_t>(i)];
      });
    } catch (const ScriptError&) {
      throw;
    } catch (const kernel::Refusal& e) {
      throw ScriptError(at, what + ": " + e.what());
    }
  }
};

namespace detail {
inline std::vector<double> numbers_of(const Scope& sc, const TermPtr& law, language::Place at) {
  std::vector<double> out;
  if (law->kind == terms::Kind::call)
    for (const auto& a : law->args) out.push_back(sc.constant(a, at, "a number of the law " + law->name));
  return out;
}
inline std::size_t count_of(const Scope& sc, const std::string& text, language::Place at) {
  if (text.empty()) return 1;
  const double v = sc.constant(language::TermPtr(std::isdigit(static_cast<unsigned char>(text[0])) ? terms::number(std::strtod(text.c_str(), nullptr)) : terms::symbol(text)), at, "the dimension '" + text + "'");
  if (v < 1 || std::floor(v) != v || v > 4096) throw ScriptError(at, "the dimension '" + text + "' is " + terms::number_text(v) + ", not a count");
  return static_cast<std::size_t>(v);
}
}  // namespace detail

inline Scope scope(const System& s) {
  Scope sc;
  sc.system = &s;
  if (s.structure.timed) sc.time = s.structure.time;
  sc.discrete = discrete_time(s);
  // Parameters first, in order: a state's dimension may be one.
  for (const auto& b : s.params) {
    if (!b.name.arguments.empty()) {
      sc.functions[b.name.name] = &b;
      continue;
    }
    const double v = sc.constant(b.value, b.at, "the parameter '" + b.name.name + "'");
    sc.params.push_back(b.name.name);
    sc.param_values.push_back(v);
  }
  for (const auto& st : s.states) {
    if (st.space.words) continue;
    const std::size_t n = detail::count_of(sc, st.space.exponent, st.at);
    for (const auto& d : st.names) {
      if (sc.declares(d.name)) throw ScriptError(st.at, "'" + d.name + "' is declared twice");
      if (st.space.exponent.empty()) {
        sc.states.push_back(d.name);
      } else {
        auto& parts = sc.vectors[d.name];
        for (std::size_t i = 1; i <= n; i++) parts.push_back(d.name + std::to_string(i)), sc.states.push_back(parts.back());
      }
    }
  }
  for (const auto& in : s.inputs) {
    sc.inputs.push_back(in.name);
    sc.input_values.push_back(in.fallback ? sc.constant(in.fallback, s.at, "the input '" + in.name + "'") : 0.0);
  }
  for (const auto& n : s.noises) {
    const std::string kind = n.source->kind == terms::Kind::call || n.source->kind == terms::Kind::symbol ? n.source->name : std::string();
    if (kind.empty()) throw ScriptError(n.at, "the law of '" + n.name + "' is not a name");
    if (n.law) {
      sc.initial.push_back({n.name, kind, kind == "data" ? std::vector<double>{} : detail::numbers_of(sc, n.source, n.at)});
      continue;
    }
    std::vector<double> numbers = detail::numbers_of(sc, n.source, n.at);
    if (kind == "wiener" && !numbers.empty() && numbers[0] > 1) {
      auto& parts = sc.vectors[n.name];
      for (std::size_t i = 1; i <= static_cast<std::size_t>(numbers[0]); i++) parts.push_back(n.name + std::to_string(i)), sc.sources.push_back({parts.back(), kind, {}});
    } else {
      if (kind == "wiener") numbers.clear();
      sc.sources.push_back({n.name, kind, numbers});
    }
  }
  return sc;
}

// The i-th component of a term over vectors: every vector of the scope is replaced by its i-th part.
inline TermPtr component(const Scope& sc, const TermPtr& t, std::size_t i, const std::string& prefix = "") {
  return terms::rewrite(t, [&](const TermPtr& x) -> TermPtr {
    if (x->kind != terms::Kind::symbol) return nullptr;
    std::string base = x->name;
    const bool differential = !prefix.empty() && base.rfind(prefix, 0) == 0 && sc.vectors.count(base.substr(prefix.size()));
    if (differential) base = base.substr(prefix.size());
    const auto v = sc.vectors.find(base);
    if (v == sc.vectors.end() || i >= v->second.size()) return nullptr;
    return terms::symbol((differential ? prefix : "") + v->second[i], x->primes);
  });
}

}  // namespace kinesis::systems
