// Systems made of systems: parts of any form, connected at their inputs.
//
// A composite declares parts, each a system of the script, and wires: `connect { a.I = 0.5 * (b.v - a.v) }` gives
// the input I of the part a as a term of the parts' states (written part.state), the composite's parameters and
// the time. A part's input that no wire gives keeps its own value.
//
// Two ways of running one, which say what they are.
//
// Merged. When every part is a flow (or an action, through its stationary paths), the composite is itself a flow:
// merged() builds its declaration, with each part's states renamed part_state, each input replaced by the term
// wired to it, and the parts' parameters inlined as numbers. The result is a system of the language. It is exact:
// no error is added by the composition, and every tool applies to the whole (its exponents, its equilibria, its
// events).
//
// Co-simulated. Parts of different forms (a flow, a map, a stochastic differential) have no common presentation
// to merge into, and are advanced side by side: at each communication step H the wires are evaluated on the
// states of all parts, each part's inputs are held at those values, and each part is advanced over H by its own
// method (a flow by Runge–Kutta steps and its events, a map by one iteration, a stochastic differential by
// Euler–Maruyama steps from its own streams). This is the explicit (Jacobi) scheme with inputs held constant
// (Kübler and Schiehlen 2000, Math. Comput. Model. Dyn. Syst. 6, 93–113; Gomes et al. 2018, ACM Comput. Surv. 51,
// 49). Its error is of the first order in H where a wired signal varies within a step, and it is exact where the
// signals are constant between steps, as when a map drives a flow. A map iterates once in each step H: H is the
// time between its iterations.
//
// The order of the parts does not matter: the wires are all read before any part moves.
//
// Joined by an action. When every part is given by its action, the composite may give the coupling as an action
// of its own, `change { action S = integral(-k / 2 * (a.q - b.q)^2) }` with `histories stationary S`: the whole
// is then the system whose Lagrangian is the sum of the parts' and the interaction's,
//
//     L = Σ_parts L_part + L_int,
//
// and joined() varies it as one (systems/lagrangian.hpp). A coupling given so is reciprocal by construction: each
// part acts on the other through the one term, and what the whole Lagrangian is invariant under is kept by the
// whole (Noether; Arnold 1989, §20): the total energy when L does not hold the time, the total momentum when L_int
// depends on differences of coordinates alone. A wire need not be reciprocal: `b.u = k * (a.q - b.q)` moves b and
// leaves a as it was, and nothing is kept. Couplings between living things are of this kind more often than not
// (what eats is not eaten back), and they are wires; an interaction action is for couplings that are forces.
// *Remark.* Every wired composite of flows is nonetheless the stationarity of an action in doubled variables, the
// response-field action S[x̃, x] = ∫ Σ_a x̃_a (x_a' − f_a(x_a) − g_a(x)) dt of Martin, Siggia and Rose (1973, Phys.
// Rev. A 8, 423–437), stationary in x̃; there a wire is the term −∫ x̃_a g_a(x) dt of an interaction action, and
// noise adds a term quadratic in x̃. That form is recorded here and not built.
#pragma once

#include <kinesis/language/declaration.hpp>
#include <kinesis/systems/dynamics.hpp>
#include <kinesis/systems/lagrangian.hpp>
#include <kinesis/systems/simulate.hpp>
#include <kinesis/systems/stochastic.hpp>

#include <map>
#include <string>
#include <variant>
#include <vector>

namespace kinesis::coupling {

using language::Script;
using language::ScriptError;
using language::System;
using systems::Dynamics;
using systems::Form;
using systems::State;
using systems::Stochastic;
using terms::TermPtr;

namespace detail {
// The declaration a part runs as: its own, or the stationary paths of its action.
inline System runnable(const System& s) { return systems::form_of(s) == Form::lagrangian ? systems::varied(s, s.name) : s; }

inline const System& declaration_of(const Script& script, const System& composite, const language::Part& part) {
  const System* s = script.find(part.system);
  if (!s) throw ScriptError(composite.at, "the part '" + part.name + "' is a '" + part.system + "', which the script does not declare");
  if (s->derived() || !s->parts.empty()) throw ScriptError(composite.at, "the part '" + part.name + "' is itself derived or composed: parts are declared systems");
  return *s;
}
}  // namespace detail

// ---- merged ----------------------------------------------------------------------------------------------

// Whether a composite gives its coupling as an action.
inline bool joined_by_action(const System& composite) { return !composite.parts.empty() && !composite.change.actions.empty(); }

// The system whose Lagrangian is the sum of the parts' and of the composite's interaction: a declaration given by
// an action, its coordinates the parts' named part_coordinate, the parts' parameters as numbers, the composite's
// own kept. The interaction names a part's coordinate part.q and its velocity part.q'.
inline System joined_action(const Script& script, const System& composite, const std::string& as = "") {
  if (!joined_by_action(composite)) throw kernel::Refusal(composite.name + " gives no action that joins its parts");
  if (composite.structure.stationary.empty()) throw ScriptError(composite.at, composite.name + " joins its parts by an action: say which with `histories stationary S`");
  if (!composite.wires.empty()) throw ScriptError(composite.wires.front().at, composite.name + " joins its parts by an action and has wires as well: a coupling is written in the action, or as wires, not both");
  const language::Named* action = nullptr;
  for (const auto& a : composite.change.actions)
    if (a.name == composite.structure.stationary) action = &a;
  if (!action) throw ScriptError(composite.at, "histories are stationary for '" + composite.structure.stationary + "', which is not an action of " + composite.name);
  const systems::Scope outer = systems::scope(composite);
  System whole;
  whole.name = as.empty() ? composite.name + "_joined" : as;
  whole.at = composite.at;
  whole.params = composite.params;
  whole.structure.timed = true, whole.structure.time = outer.time, whole.structure.time_domain = composite.structure.time_domain;
  whole.structure.stationary = action->name;
  whole.difference = composite.difference;
  language::States all;
  all.space.name = "real", all.at = composite.at;
  std::vector<TermPtr> lagrangians;
  for (const auto& part : composite.parts) {
    const System& decl = detail::declaration_of(script, composite, part);
    if (systems::form_of(decl) != Form::lagrangian) throw kernel::Refusal("the part '" + part.name + "' is not given by an action: parts joined by an action are each given by one; join the others by wires");
    systems::Scope sc = systems::scope(decl);
    for (const auto& o : part.overrides) {
      const int i = sc.param_index(o.name.name);
      if (i < 0) throw ScriptError(composite.at, "'" + o.name.name + "' is not a parameter of the part '" + part.name + "'");
      sc.param_values[static_cast<std::size_t>(i)] = outer.constant(o.value, o.at, "the parameter '" + o.name.name + "' of the part '" + part.name + "'");
    }
    const systems::Lagrangian L = systems::vary(decl);
    lagrangians.push_back(terms::rewrite(L.lagrangian, [&](const TermPtr& x) -> TermPtr {
      if (x->kind != terms::Kind::symbol) return nullptr;
      if (sc.is_state(x->name)) return terms::symbol(part.name + "_" + x->name, x->primes);
      if (x->primes) return nullptr;
      const int p = sc.param_index(x->name);
      if (p >= 0) return terms::number(sc.param_values[static_cast<std::size_t>(p)]);
      if (x->name == sc.time) return terms::symbol(outer.time);
      return nullptr;
    }));
    for (const auto& q : L.coordinates) all.names.push_back({part.name + "_" + q, {outer.time}});
  }
  // The interaction in the whole's names: a.q becomes a_q, with its primes.
  lagrangians.push_back(terms::rewrite(systems::detail::integrand(action->value), [&](const TermPtr& x) -> TermPtr {
    if (x->kind != terms::Kind::symbol) return nullptr;
    const auto dot = x->name.find('.');
    if (dot == std::string::npos) return nullptr;
    const std::string name = x->name.substr(0, dot) + "_" + x->name.substr(dot + 1);
    bool known = false;
    for (const auto& d : all.names) known = known || d.name == name;
    if (!known) throw ScriptError(action->at, "the action names '" + x->name + "', which is not a coordinate of a part of " + composite.name);
    return terms::symbol(name, x->primes);
  }));
  // What the composite observes, in the whole's names.
  whole.difference.observations.clear();
  for (const auto& ob : composite.difference.observations)
    whole.difference.observations.push_back({ob.name, terms::rewrite(ob.value, [&](const TermPtr& x) -> TermPtr {
      if (x->kind != terms::Kind::symbol) return nullptr;
      const auto dot = x->name.find('.');
      return dot == std::string::npos ? nullptr : terms::symbol(x->name.substr(0, dot) + "_" + x->name.substr(dot + 1), x->primes);
    })});
  whole.states.push_back(std::move(all));
  whole.change.actions.push_back({action->name, terms::call("integral", {terms::sum(std::move(lagrangians))}), action->at});
  return whole;
}

// The flow a composite of flows is, as a declaration.
inline System merged(const Script& script, const System& composite, const std::string& as = "") {
  if (composite.parts.empty()) throw kernel::Refusal(composite.name + " has no parts");
  if (joined_by_action(composite)) return systems::varied(joined_action(script, composite, composite.name), as.empty() ? composite.name + "_merged" : as);
  const systems::Scope outer = systems::scope(composite);
  System out;
  out.name = as.empty() ? composite.name + "_merged" : as;
  out.at = composite.at;
  out.params = composite.params;
  out.structure = composite.structure;
  language::States all;
  all.space.name = "real";
  all.at = composite.at;
  std::map<std::string, std::map<std::string, TermPtr>> wired;  // part → input → term
  for (const auto& w : composite.wires) wired[w.part][w.input] = w.value;
  // A wire's term in the merged names: a.v becomes a_v.
  const auto global = [&](const TermPtr& t) {
    return terms::rewrite(t, [&](const TermPtr& x) -> TermPtr {
      if (x->kind != terms::Kind::symbol) return nullptr;
      const auto dot = x->name.find('.');
      return dot == std::string::npos ? nullptr : terms::symbol(x->name.substr(0, dot) + "_" + x->name.substr(dot + 1), x->primes);
    });
  };
  for (const auto& part : composite.parts) {
    const System decl = detail::runnable(detail::declaration_of(script, composite, part));
    const Form form = systems::form_of(decl);
    if (form != Form::flow) throw kernel::Refusal("the part '" + part.name + "' is a " + std::string(systems::name(form)) + " system: only flows are merged into one declaration; co-simulate the others");
    systems::Scope sc = systems::scope(decl);
    for (const auto& o : part.overrides) {
      const int i = sc.param_index(o.name.name);
      if (i < 0) throw ScriptError(composite.at, "'" + o.name.name + "' is not a parameter of the part '" + part.name + "'");
      sc.param_values[static_cast<std::size_t>(i)] = outer.constant(o.value, o.at, "the parameter '" + o.name.name + "' of the part '" + part.name + "'");
    }
    // A part's term in the merged names: states prefixed, parameters as numbers, inputs as their wires.
    const auto local = [&](const TermPtr& t) {
      return terms::rewrite(sc.resolve(t), [&](const TermPtr& x) -> TermPtr {
        if (x->kind != terms::Kind::symbol) return nullptr;
        if (sc.is_state(x->name)) return terms::symbol(part.name + "_" + x->name, x->primes);
        if (x->primes) return nullptr;
        const int p = sc.param_index(x->name);
        if (p >= 0) return terms::number(sc.param_values[static_cast<std::size_t>(p)]);
        const auto in = std::find(sc.inputs.begin(), sc.inputs.end(), x->name);
        if (in != sc.inputs.end()) {
          const auto w = wired[part.name].find(x->name);
          return w != wired[part.name].end() ? global(w->second) : terms::number(sc.input_values[static_cast<std::size_t>(in - sc.inputs.begin())]);
        }
        if (x->name == sc.time) return terms::symbol(outer.time);
        return nullptr;
      });
    };
    for (const auto& st : sc.states) all.names.push_back({part.name + "_" + st, {}});
    for (const auto& e : decl.change.equations) {
      const auto v = sc.vectors.find(e.name);
      const std::size_t count = v == sc.vectors.end() ? 1 : v->second.size();
      for (std::size_t k = 0; k < count; k++)
        out.change.equations.push_back({part.name + "_" + (v == sc.vectors.end() ? e.name : v->second[k]), e.primes, {}, local(v == sc.vectors.end() ? e.value : systems::component(sc, e.value, k)), composite.at});
    }
    for (const auto& ev : decl.change.events) {
      language::Event merged_event{local(ev.left), local(ev.right), ev.relation, {}, composite.at};
      for (const auto& a : ev.resets) merged_event.resets.push_back({part.name + "_" + a.name, local(a.value)});
      out.change.events.push_back(std::move(merged_event));
    }
    for (const auto& ob : decl.difference.observations) out.difference.observations.push_back({part.name + "_" + (ob.name.empty() ? "observed" + std::to_string(out.difference.observations.size()) : ob.name), local(ob.value)});
    for (const auto& w : wired[part.name])
      if (std::find(sc.inputs.begin(), sc.inputs.end(), w.first) == sc.inputs.end()) throw ScriptError(composite.at, "'" + w.first + "' is not an input of the part '" + part.name + "'");
  }
  for (const auto& w : composite.wires) {
    bool known = false;
    for (const auto& part : composite.parts) known = known || part.name == w.part;
    if (!known) throw ScriptError(w.at, "'" + w.part + "' is not a part of " + composite.name);
  }
  out.states.push_back(std::move(all));
  for (const auto& ob : composite.difference.observations) out.difference.observations.push_back({ob.name, global(ob.value)});
  return out;
}

// ---- co-simulated ----------------------------------------------------------------------------------------

class Cosimulation {
 public:
  std::string name;
  struct Part {
    std::string name;
    Form form = Form::flow;
    std::variant<Dynamics, Stochastic> system;
    std::size_t offset = 0, dimension = 0;  // its components in the composite's state
  };
  std::vector<Part> parts;
  std::vector<std::string> names;  // of every component: part.state

  std::size_t dimension() const { return names.size(); }
  int index(const std::string& qualified) const {
    const auto it = std::find(names.begin(), names.end(), qualified);
    return it == names.end() ? -1 : static_cast<int>(it - names.begin());
  }

  struct Options {
    double t0 = 0, t1 = 10;
    double communication = 0.01;   // H
    double step = 0.001;           // of the flows and stochastic differentials within H
    std::size_t every = 1;         // communication steps between recorded states
    std::uint64_t seed = 1;
    double event_tolerance = 1e-10;
  };

  systems::Trajectory run(State x, const Options& o) const {
    if (x.size() != dimension()) throw kernel::Refusal(name + " has " + std::to_string(dimension()) + " components; " + std::to_string(x.size()) + " were given");
    std::vector<Part> local = parts;  // each part's inputs are set as the run goes
    std::vector<std::vector<calx::random::Stream>> streams(local.size());
    for (std::size_t p = 0; p < local.size(); p++)
      if (const Stochastic* sde = std::get_if<Stochastic>(&local[p].system))
        for (std::size_t j = 0; j < sde->sources.size(); j++) streams[p].emplace_back(o.seed, 0x636f7570u, static_cast<std::uint32_t>(p), static_cast<std::uint32_t>(j));  // 'coup'
    systems::Trajectory out;
    out.names = names;
    std::vector<double> values(frame_.size(), 0.0), signal(wires_.size());
    std::copy(param_values_.begin(), param_values_.end(), values.begin() + static_cast<std::ptrdiff_t>(dimension()));
    double t = o.t0;
    out.t.push_back(t), out.x.push_back(x);
    const std::size_t count = static_cast<std::size_t>(std::ceil((o.t1 - o.t0) / o.communication - 1e-9));
    systems::Options part_options;
    part_options.event_tolerance = o.event_tolerance;
    for (std::size_t k = 1; k <= count; k++) {
      const double H = std::min(o.communication, o.t1 - t);
      // The wires, all read on the state before any part moves.
      std::copy(x.begin(), x.end(), values.begin());
      values[dimension() + param_values_.size()] = t;
      for (std::size_t w = 0; w < wires_.size(); w++) signal[w] = wires_[w].program(values);
      for (std::size_t w = 0; w < wires_.size(); w++) {
        Part& part = local[wires_[w].part];
        if (Dynamics* d = std::get_if<Dynamics>(&part.system)) d->input_values[wires_[w].input] = signal[w];
      }
      for (std::size_t p = 0; p < local.size(); p++) {
        Part& part = local[p];
        State y(x.begin() + static_cast<std::ptrdiff_t>(part.offset), x.begin() + static_cast<std::ptrdiff_t>(part.offset + part.dimension));
        if (Dynamics* d = std::get_if<Dynamics>(&part.system)) {
          // An event's states are those of its part; they are widened to the composite's. A map's event is
          // recorded by its iteration: in the whole it happens at the end of the communication step.
          const std::size_t before = out.events.size();
          if (d->discrete) {
            double n = static_cast<double>(k - 1);
            systems::iterate(*d, n, y, &out.events);
            for (std::size_t e = before; e < out.events.size(); e++) out.events[e].t = t + H;
          } else {
            const std::size_t sub = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(H / o.step - 1e-9)));
            double s = t;
            for (std::size_t j = 0; j < sub; j++) systems::advance(*d, s, y, H / static_cast<double>(sub), part_options, &out.events);
          }
          for (std::size_t e = before; e < out.events.size(); e++) widen(out.events[e], x, part);
        } else {
          const Stochastic& sde = std::get<Stochastic>(part.system);
          const std::size_t m = sde.sources.size();
          std::vector<double> dW(std::max<std::size_t>(m, 1));
          if (sde.discrete) {
            State next(y.size());
            for (std::size_t j = 0; j < m; j++) dW[j] = systems::draw(sde.sources[j], streams[p][j]);
            sde.drift(static_cast<double>(k - 1), y, next, dW.data());
            y.swap(next);
          } else {
            const std::size_t sub = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(H / o.step - 1e-9)));
            const double h = H / static_cast<double>(sub);
            for (std::size_t j = 0; j < sub; j++) {
              for (std::size_t q = 0; q < m; q++) dW[q] = std::sqrt(h) * streams[p][q].normal();
              sde.step(t + static_cast<double>(j) * h, y, h, dW.data(), false);
            }
          }
        }
        std::copy(y.begin(), y.end(), x.begin() + static_cast<std::ptrdiff_t>(part.offset));
      }
      t += H;
      if (k % std::max<std::size_t>(o.every, 1) == 0 || k == count) out.t.push_back(t), out.x.push_back(x);
      for (double v : x)
        if (!std::isfinite(v)) throw kernel::Refusal(name + " left the finite numbers at t = " + terms::number_text(t));
    }
    out.steps = count;
    return out;
  }

 private:
  friend Cosimulation cosimulation(const Script&, const System&);
  static void widen(systems::Occurrence& e, const State& whole, const Part& part) {
    State before = whole, after = whole;
    std::copy(e.before.begin(), e.before.end(), before.begin() + static_cast<std::ptrdiff_t>(part.offset));
    std::copy(e.after.begin(), e.after.end(), after.begin() + static_cast<std::ptrdiff_t>(part.offset));
    e.before.swap(before), e.after.swap(after);
  }
  struct Wire {
    std::size_t part, input;
    terms::Program program;
  };
  std::vector<Wire> wires_;
  terms::Frame frame_;                 // the components, the composite's parameters, the time
  std::vector<double> param_values_;
};

inline Cosimulation cosimulation(const Script& script, const System& composite) {
  if (composite.parts.empty()) throw kernel::Refusal(composite.name + " has no parts");
  const systems::Scope outer = systems::scope(composite);
  Cosimulation c;
  c.name = composite.name;
  for (const auto& part : composite.parts) {
    const System decl = detail::runnable(detail::declaration_of(script, composite, part));
    const Form form = systems::form_of(decl);
    Cosimulation::Part p;
    p.name = part.name, p.form = form, p.offset = c.names.size();
    std::vector<std::string> states;
    if (form == Form::flow || form == Form::map) {
      Dynamics d = systems::dynamics(decl);
      for (const auto& o : part.overrides) d.set_param(o.name.name, outer.constant(o.value, o.at, "the parameter '" + o.name.name + "' of the part '" + part.name + "'"));
      states = d.states;
      p.system = std::move(d);
    } else if (form == Form::stochastic) {
      Stochastic d = systems::stochastic(decl);
      for (const auto& o : part.overrides) d.set_param(o.name.name, outer.constant(o.value, o.at, "the parameter '" + o.name.name + "' of the part '" + part.name + "'"));
      states = d.states;
      p.system = std::move(d);
    } else {
      throw kernel::Refusal("the part '" + part.name + "' is a " + std::string(systems::name(form)) + " system: flows, maps and stochastic systems are co-simulated");
    }
    p.dimension = states.size();
    for (const auto& st : states) c.names.push_back(part.name + "." + st);
    c.parts.push_back(std::move(p));
  }
  for (const auto& n : c.names) c.frame_.add(n);
  for (const auto& p : outer.params) c.frame_.add(p);
  c.param_values_ = outer.param_values;
  c.frame_.add(outer.time);
  for (const auto& w : composite.wires) {
    Cosimulation::Wire wire{0, 0, {}};
    bool found = false;
    for (std::size_t p = 0; p < c.parts.size() && !found; p++) {
      if (c.parts[p].name != w.part) continue;
      const Dynamics* d = std::get_if<Dynamics>(&c.parts[p].system);
      if (!d) throw ScriptError(w.at, "the part '" + w.part + "' is stochastic: an input of a stochastic part is not wired yet");
      const int i = d->input_index(w.input);
      if (i < 0) throw ScriptError(w.at, "'" + w.input + "' is not an input of the part '" + w.part + "'");
      wire.part = p, wire.input = static_cast<std::size_t>(i), found = true;
    }
    if (!found) throw ScriptError(w.at, "'" + w.part + "' is not a part of " + composite.name);
    try {
      wire.program = terms::compile(outer.resolve(w.value), c.frame_.resolver());
    } catch (const kernel::Refusal& e) {
      throw ScriptError(w.at, "the wire of " + w.part + "." + w.input + ": " + e.what());
    }
    c.wires_.push_back(std::move(wire));
  }
  return c;
}

}  // namespace kinesis::coupling
