// Flows and maps: a declaration of equations made into a system that runs.
//
// A flow is x' = f(t, x; p, u) on ℝⁿ and a map x ↦ f(n, x; p, u), with parameters p and inputs u. dynamics()
// builds one from a declaration:
//
//   - an equation of order k for q (q'' = …) makes the components q, q', …, q⁽ᵏ⁻¹⁾, so that the system is of
//     first order; every declared state needs one equation;
//   - an equation for a state in real^n is n equations, component by component;
//   - parameters with arguments are inlined, so that the field is a term of the states, the parameters, the time
//     and the inputs alone;
//   - the Jacobian ∂f/∂x is derived from the terms, exactly, entry by entry, and only its entries that are not
//     zero are kept;
//   - an event `when a >= b : resets` is a guard g = a − b with the sense of its crossing, and its resets, which
//     are evaluated together on the state before the event.
//
// The values of the parameters and of the inputs can be set after building; the terms do not change.
//
// A Dynamics holds scratch arrays for evaluation: one object is for one thread, and a copy is independent.
#pragma once

#include <kinesis/systems/model.hpp>

#include <calx/linalg/dense.hpp>

#include <string>
#include <utility>
#include <vector>

namespace kinesis::systems {

class Dynamics {
 public:
  std::string name;
  bool discrete = false;
  std::string time = "t";
  std::vector<std::string> states;         // q, q', X1, …
  std::vector<std::string> params;
  std::vector<double> param_values;
  std::vector<std::string> inputs;
  std::vector<double> input_values;
  std::vector<TermPtr> field_terms;        // f, component by component
  std::vector<std::string> observables;
  std::vector<TermPtr> observable_terms;
  std::vector<TermPtr> conserved_terms;
  std::vector<InitialLaw> initial;

  struct Event {
    TermPtr guard;                         // fires when it crosses zero
    bool rising = true;                    // from below (>=, >) or from above (<=, <)
    std::vector<std::pair<std::size_t, TermPtr>> resets;  // (component, its new value)
    std::string text;
  };
  std::vector<Event> events;

  std::size_t dimension() const { return states.size(); }

  // dx = f(t, x), or the next state of a map.
  void field(double t, const State& x, State& dx) const {
    load(t, x);
    for (std::size_t i = 0; i < field_.size(); i++) dx[i] = field_[i].run(values_.data(), stack_.data());
  }
  // J = ∂f/∂x at (t, x), exact.
  void jacobian(double t, const State& x, calx::linalg::Matrix& J) const {
    load(t, x);
    if (J.n != dimension()) J = calx::linalg::Matrix(dimension());
    std::fill(J.a.begin(), J.a.end(), 0.0);
    for (const auto& e : jacobian_) J(e.row, e.column) = e.program.run(values_.data(), stack_.data());
  }
  // out = J(t, x) v, without forming J.
  void jacobian_vector(double t, const State& x, const State& v, State& out) const {
    load(t, x);
    std::fill(out.begin(), out.end(), 0.0);
    for (const auto& e : jacobian_) out[e.row] += e.program.run(values_.data(), stack_.data()) * v[e.column];
  }
  // The terms of the Jacobian that are not zero, as (row, column, term).
  struct Entry {
    std::size_t row, column;
    TermPtr term;
  };
  const std::vector<Entry>& jacobian_terms() const { return jacobian_terms_; }

  double observe(std::size_t i, double t, const State& x) const {
    load(t, x);
    return observable_[i].run(values_.data(), stack_.data());
  }
  double conserved(std::size_t i, double t, const State& x) const {
    load(t, x);
    return conserved_[i].run(values_.data(), stack_.data());
  }
  double guard(std::size_t e, double t, const State& x) const {
    load(t, x);
    return guard_[e].run(values_.data(), stack_.data());
  }
  // The state after event e at (t, x): its resets, evaluated together.
  State reset(std::size_t e, double t, const State& x) const {
    load(t, x);
    State y = x;
    for (std::size_t k = 0; k < events[e].resets.size(); k++) y[events[e].resets[k].first] = reset_[e][k].run(values_.data(), stack_.data());
    return y;
  }

  int state_index(const std::string& n) const { return index_in(states, n); }
  int param_index(const std::string& n) const { return index_in(params, n); }
  int input_index(const std::string& n) const { return index_in(inputs, n); }
  void set_param(const std::string& n, double v) {
    const int i = param_index(n);
    if (i < 0) throw kernel::Refusal("'" + n + "' is not a parameter of " + name);
    param_values[static_cast<std::size_t>(i)] = v;
  }
  double param(const std::string& n) const {
    const int i = param_index(n);
    if (i < 0) throw kernel::Refusal("'" + n + "' is not a parameter of " + name);
    return param_values[static_cast<std::size_t>(i)];
  }

  // A term of the system's names (its states, parameters, time and inputs) compiled for this system's values.
  terms::Program compile(const TermPtr& t) const { return terms::compile(t, frame_.resolver()); }
  double value(const terms::Program& p, double t, const State& x) const {
    load(t, x);
    return p.depth() <= stack_.size() ? p.run(values_.data(), stack_.data()) : p(values_);
  }

 private:
  friend Dynamics dynamics_of(const System&, const Scope&, const std::vector<language::Equation>&);
  static int index_in(const std::vector<std::string>& v, const std::string& n) {
    const auto it = std::find(v.begin(), v.end(), n);
    return it == v.end() ? -1 : static_cast<int>(it - v.begin());
  }
  void load(double t, const State& x) const {
    const std::size_t n = states.size(), p = params.size();
    if (values_.size() != n + p + 1 + inputs.size()) values_.assign(n + p + 1 + inputs.size(), 0.0);
    std::copy(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(n), values_.begin());
    std::copy(param_values.begin(), param_values.end(), values_.begin() + static_cast<std::ptrdiff_t>(n));
    values_[n + p] = t;
    std::copy(input_values.begin(), input_values.end(), values_.begin() + static_cast<std::ptrdiff_t>(n + p + 1));
  }
  void deepen(const terms::Program& p) {
    if (p.depth() > stack_.size()) stack_.assign(p.depth(), 0.0);
  }

  struct Compiled {
    std::size_t row, column;
    terms::Program program;
  };
  terms::Frame frame_;
  std::vector<terms::Program> field_, observable_, conserved_, guard_;
  std::vector<std::vector<terms::Program>> reset_;
  std::vector<Compiled> jacobian_;
  std::vector<Entry> jacobian_terms_;
  mutable std::vector<double> values_, stack_;
};

// The system of the equations given, in the scope given (the equations may be a declaration's own, or derived).
inline Dynamics dynamics_of(const System& s, const Scope& sc, const std::vector<language::Equation>& equations) {
  Dynamics d;
  d.name = s.name, d.discrete = sc.discrete, d.time = sc.time;
  d.params = sc.params, d.param_values = sc.param_values;
  d.inputs = sc.inputs, d.input_values = sc.input_values;
  d.initial = sc.initial;

  // The equation of each scalar state: (order, right side).
  struct Law {
    unsigned order = 0;
    TermPtr value;
    language::Place at;
  };
  std::map<std::string, Law> laws;
  const auto give = [&](const std::string& state, unsigned order, const TermPtr& value, language::Place at) {
    if (!sc.is_state(state)) throw ScriptError(at, "'" + state + "' is not a state of " + s.name);
    if (order == 0) throw ScriptError(at, "an equation gives a derivative or a next value: write " + state + "' = …");
    if (sc.discrete && order != 1) throw ScriptError(at, "in discrete time " + state + "' is the next value; there is no " + state + std::string(order, '\''));
    if (laws.count(state)) throw ScriptError(at, "'" + state + "' has two equations");
    laws[state] = {order, sc.resolve(value), at};
  };
  for (const auto& e : equations) {
    const auto v = sc.vectors.find(e.name);
    if (v != sc.vectors.end()) {
      for (std::size_t i = 0; i < v->second.size(); i++) give(v->second[i], e.primes, component(sc, e.value, i), e.at);
    } else {
      give(e.name, e.primes, e.value, e.at);
    }
  }
  // The components and their slots: states, then parameters, the time, the inputs.
  std::vector<std::pair<std::string, unsigned>> components;
  for (const auto& st : sc.states) {
    const auto l = laws.find(st);
    if (l == laws.end()) throw ScriptError(s.at, "no change is given for the state '" + st + "' of " + s.name);
    for (unsigned k = 0; k < l->second.order; k++) components.push_back({st, k});
  }
  for (const auto& [st, k] : components) d.frame_.add(st, k), d.states.push_back(st + std::string(k, '\''));
  for (const auto& p : d.params) d.frame_.add(p);
  d.frame_.add(d.time);
  for (const auto& in : d.inputs) d.frame_.add(in);

  const auto compiled = [&](const TermPtr& t, language::Place at, const std::string& what) {
    try {
      terms::Program p = terms::compile(t, d.frame_.resolver());
      d.deepen(p);
      return p;
    } catch (const kernel::Refusal& e) {
      throw ScriptError(at, what + ": " + e.what());
    }
  };
  for (const auto& [st, k] : components) {
    const Law& l = laws[st];
    const TermPtr f = k + 1 < l.order ? terms::symbol(st, k + 1) : l.value;
    d.field_terms.push_back(f);
    d.field_.push_back(compiled(f, l.at, "the equation of " + st));
  }
  for (std::size_t i = 0; i < components.size(); i++)
    for (std::size_t j = 0; j < components.size(); j++) {
      TermPtr entry;
      try {
        entry = terms::derivative(d.field_terms[i], components[j].first, components[j].second);
      } catch (const kernel::Refusal& e) {
        throw ScriptError(laws[components[i].first].at, "the equation of " + components[i].first + " cannot be differentiated: " + e.what());
      }
      if (terms::is_number(entry, 0)) continue;
      d.jacobian_terms_.push_back({i, j, entry});
      d.jacobian_.push_back({i, j, compiled(entry, laws[components[i].first].at, "the derivative of the equation of " + components[i].first)});
    }
  for (const auto& e : s.change.events) {
    Dynamics::Event ev;
    ev.guard = sc.resolve(e.left - e.right);
    ev.rising = e.relation == ">=" || e.relation == ">";
    ev.text = terms::text(e.left) + " " + e.relation + " " + terms::text(e.right);
    d.guard_.push_back(compiled(ev.guard, e.at, "the guard of an event"));
    d.reset_.emplace_back();
    for (const auto& a : e.resets) {
      const int i = d.state_index(a.name);
      if (i < 0) throw ScriptError(e.at, "an event resets '" + a.name + "', which is not a state of " + s.name);
      const TermPtr v = sc.resolve(a.value);
      ev.resets.push_back({static_cast<std::size_t>(i), v});
      d.reset_.back().push_back(compiled(v, e.at, "the reset of " + a.name));
    }
    d.events.push_back(std::move(ev));
  }
  for (const auto& ob : s.difference.observations) {
    const TermPtr v = sc.resolve(ob.value);
    d.observables.push_back(ob.name.empty() ? terms::text(ob.value) : ob.name);
    d.observable_terms.push_back(v);
    d.observable_.push_back(compiled(v, s.at, "the observation " + d.observables.back()));
  }
  for (const auto& c : s.constraint.conserved) {
    const TermPtr v = sc.resolve(c);
    d.conserved_terms.push_back(v);
    d.conserved_.push_back(compiled(v, s.at, "the conserved quantity " + terms::text(c)));
  }
  return d;
}

// The flow or the map a declaration of equations is.
inline Dynamics dynamics(const System& s) {
  const Form f = form_of(s);
  if (f != Form::flow && f != Form::map) throw kernel::Refusal(s.name + " is a " + std::string(name(f)) + " system, not a flow or a map");
  return dynamics_of(s, scope(s), s.change.equations);
}

}  // namespace kinesis::systems
