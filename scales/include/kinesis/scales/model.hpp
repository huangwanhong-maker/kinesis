// A system across scales: quantities on the levels of a hierarchy, the operators that carry them between levels, and the
// laws that advance them, each at its own pace (decision 0007).
//
// Quantities. A quantity has a value for every unit of its level. It is a state (advanced by a law), a driver (given over
// time from outside), or fixed (a value for each unit that does not change). Its kind says how it is carried up the
// hierarchy: an extensive quantity (people, floor space, money) is summed; an intensive one (a density, a rent, a share) is
// averaged with a declared weight; an ordinal one (a level on a scale) is read through thresholds on a latent value.
// A quantity may be bounded (a share in [0, 1], a stock not below zero); a run keeps it there and says how often it had to.
//
// Derived quantities. Every quantity that crosses a level, or depends on more than one unit, is derived by one of a few
// operators, so that no law names another level by accident:
//
//   value       a term of the level's own quantities (an index: vacant / total)
//   sum, mean, wmean, minimum, maximum, reduce
//               over the units beneath (reduce by a reducer the caller registers: a Gini coefficient, an entropy)
//   parent      a term of the level above, read by each unit beneath it
//   share       x_i / Σ_{j ~ i} x_j among siblings: allocation of a parent's total
//   logit       exp(x_i) / Σ_{j ~ i} exp(x_j) among siblings: a choice among them
//   zscore      (x_i − mean) / sd among siblings (sd with divisor n − 1): the standardisation composite indexes are built from
//   laplacian   (1/V_a) Σ_e c_e (x_b − x_a) on the level's graph (balance laws on graphs; LeVeque 2002, §4.1)
//   neighbours  the conductance-weighted mean of x over a unit's neighbours
//   inflow      what each unit receives when every unit sends its mass O_i to all units in proportion to A_j e^{−β c_ij}:
//               D_j = Σ_i O_i A_j e^{−β c_ij} / Σ_k A_k e^{−β c_ik}, the production-constrained spatial interaction model
//               (Wilson 1971, Environment and Planning 3, 1–32), whose total Σ_j D_j = Σ_i O_i is conserved
//
// Since an aggregate is the sum of its parts and an allocation is a share of its whole, a model in which coarse quantities are
// derived up and fine ones allocated down is consistent across its levels by construction. A ratio index defined as a value
// on summed extensive quantities is the ratio of sums at every level, not a mean of ratios.
//
// Laws. A law advances one state: as a flow, x' = f, integrated over its own step; as a map, x ← g once in each of its periods;
// or as an equilibrium, x = g solved at every instant (a quantity whose own dynamics are so fast that it is always at rest
// given the slow ones: the quasi-steady state of Tikhonov's theorem, Tikhonov 1952, Mat. Sb. 31, 575–586; Fenichel 1979,
// J. Differ. Equ. 31, 53–98). A law is an expression of the kinesis language in its level's quantities, its derived
// quantities, its measures, the model's parameters and the time t.
//
// Pace. Laws with different steps are advanced side by side and exchange their states at the start of each of their steps
// (the explicit, Jacobi scheme of co-simulation, as kinesis's coupling has it): a law with a step of a year reads the fast
// states as they were at the start of the year. When every law shares one step the whole is one system integrated by
// Runge–Kutta's classical step.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "kinesis/kernel/answer.hpp"
#include "kinesis/language/reader.hpp"
#include "kinesis/scales/hierarchy.hpp"
#include "kinesis/terms/program.hpp"
#include "kinesis/terms/term.hpp"

namespace kinesis::scales {

inline constexpr double infinity = std::numeric_limits<double>::infinity();
inline constexpr double missing_value = std::numeric_limits<double>::quiet_NaN();

enum class Kind { extensive, intensive, ordinal, other };
enum class Role { state, driver, fixed };
enum class Op { value, sum, mean, wmean, minimum, maximum, reduce, parent, share, logit, zscore, laplacian, neighbours, inflow };
enum class Form { flow, map, equilibrium };

inline const char* name(Kind k) {
  switch (k) {
    case Kind::extensive: return "extensive";
    case Kind::intensive: return "intensive";
    case Kind::ordinal: return "ordinal";
    default: return "other";
  }
}
inline const char* name(Role r) {
  switch (r) {
    case Role::state: return "state";
    case Role::driver: return "driver";
    default: return "fixed";
  }
}
inline const char* name(Op o) {
  static const char* names[] = {"value", "sum", "mean", "wmean", "minimum", "maximum", "reduce",
                                "parent", "share", "logit", "zscore", "laplacian", "neighbours", "inflow"};
  return names[static_cast<int>(o)];
}
inline const char* name(Form f) {
  switch (f) {
    case Form::flow: return "flow";
    case Form::map: return "map";
    default: return "equilibrium";
  }
}
inline Kind kind_of(const std::string& s) {
  if (s == "extensive") return Kind::extensive;
  if (s == "intensive") return Kind::intensive;
  if (s == "ordinal") return Kind::ordinal;
  if (s.empty() || s == "other") return Kind::other;
  throw Refusal("no kind of quantity is named " + s + " (extensive, intensive, ordinal, other)");
}
inline Role role_of(const std::string& s) {
  if (s == "state") return Role::state;
  if (s == "driver") return Role::driver;
  if (s == "fixed") return Role::fixed;
  throw Refusal("no role of a quantity is named " + s + " (state, driver, fixed)");
}
inline Op op_of(const std::string& s) {
  for (int i = 0; i <= static_cast<int>(Op::inflow); i++)
    if (s == name(static_cast<Op>(i))) return static_cast<Op>(i);
  throw Refusal("no operator is named " + s +
                " (value, sum, mean, wmean, minimum, maximum, reduce, parent, share, logit, zscore, laplacian, neighbours, inflow)");
}
inline Form form_of(const std::string& s) {
  if (s == "flow") return Form::flow;
  if (s == "map") return Form::map;
  if (s == "equilibrium") return Form::equilibrium;
  throw Refusal("no form of law is named " + s + " (flow, map, equilibrium)");
}

struct Quantity {
  std::string name, level;
  Kind kind = Kind::other;
  Role role = Role::state;
  std::string weight;  // an intensive quantity's weight: a quantity or measure of its level
  std::string unit;    // in words; kinesis carries no units
  double lower = -infinity, upper = infinity;
  std::vector<double> thresholds;  // an ordinal quantity's cut points on its latent scale, increasing
};

struct Derived {
  std::string name, level;
  Op op = Op::value;
  std::string of;       // the term reduced, shared or carried (read at the level the operator reads)
  std::string weight;   // wmean and reduce: the weights' term (beneath); inflow: the attractiveness A; laplacian: a measure
  std::string reducer;  // reduce: a registered reducer
  std::string costs;    // inflow: a matrix of costs of the level
  std::string beta;     // inflow: a term of parameters and time
  Kind kind = Kind::other;
  std::string unit;
};

struct Law {
  std::string target;
  Form form = Form::flow;
  std::string expression;
  double step = 0;     // a flow's step or a map's period; 0: the model's step
  double damping = 1;  // an equilibrium's relaxation, in (0, 1]
};

struct Parameter {
  std::string name;
  double value = missing_value;  // its default, when the inputs give none
  double lower = -infinity, upper = infinity;
  bool fit = false;  // an unknown to be fitted to a record
  std::string unit;
};

// A reducer of the values of the units beneath one unit, with their weights.
using Reducer = std::function<double(std::span<const double> values, std::span<const double> weights)>;

struct Model {
  std::vector<Quantity> quantities;
  std::vector<Derived> derived;
  std::vector<Law> laws;
  std::vector<Parameter> parameters;
  double step = 1;  // the base step: every law's step is a whole multiple of it
  std::map<std::string, Reducer> reducers;
};

// ---- compilation -------------------------------------------------------------------------------------------

// Where a compiled expression reads each of its symbols.
struct Ref {
  enum What { column, global, time } what = column;
  std::size_t slot = 0;
};

struct Expression {
  std::size_t level = 0;
  std::string text;
  terms::Program program;
  std::vector<Ref> reads;
  bool global_only() const {
    return std::all_of(reads.begin(), reads.end(), [](const Ref& r) { return r.what != Ref::column; });
  }
};

struct CompiledDerived {
  std::size_t index = 0, level = 0, slot = 0;
  Op op = Op::value;
  Expression of, weight, beta;
  bool has_weight = false;
  Reducer reducer;      // reduce
  std::string measure;  // laplacian
  std::string costs;    // inflow: the name of a matrix of the level's costs
};

struct CompiledLaw {
  std::size_t index = 0, level = 0, slot = 0;
  Form form = Form::flow;
  Expression rhs;
  std::size_t every = 1;  // in base steps
  double damping = 1;
};

class Compiled {
 public:
  Hierarchy hierarchy;
  Model model;
  std::vector<std::vector<std::string>> columns;             // per level: slot → name
  std::vector<std::map<std::string, std::size_t>> slots;     // per level: name → slot
  std::vector<std::string> globals;                          // parameter names, in order
  std::map<std::string, std::size_t> global_slot;
  std::vector<CompiledDerived> derived;                      // in an order in which each depends only on those before
  std::vector<CompiledLaw> laws;
  std::map<std::string, std::size_t> level_of;               // every quantity, derived quantity → level

  const Quantity* quantity(const std::string& n) const {
    for (const auto& q : model.quantities)
      if (q.name == n) return &q;
    return nullptr;
  }
  const Derived* derived_named(const std::string& n) const {
    for (const auto& d : model.derived)
      if (d.name == n) return &d;
    return nullptr;
  }
  const Law* law_of(const std::string& target) const {
    for (const auto& l : model.laws)
      if (l.target == target) return &l;
    return nullptr;
  }
  std::size_t slot(std::size_t level, const std::string& n) const {
    const auto it = slots.at(level).find(n);
    if (it == slots[level].end()) throw Refusal("the level " + hierarchy.level(level).name + " has no quantity " + n);
    return it->second;
  }
  // How a quantity is carried up: its kind, and for an intensive one its weight.
  Kind kind(const std::string& n) const {
    if (const Quantity* q = quantity(n)) return q->kind;
    if (const Derived* d = derived_named(n)) return d->kind;
    throw Refusal("the model has no quantity " + n);
  }
};

namespace detail {

inline Expression compile_at(const Compiled& c, std::size_t level, const std::string& text, const std::string& what) {
  Expression e;
  e.level = level;
  e.text = text;
  if (text.empty()) throw Refusal(what + " has no expression");
  terms::TermPtr t;
  try {
    t = language::read_term(text);
  } catch (const Refusal& r) {
    throw Refusal(what + ": " + r.what());
  }
  std::map<std::string, int> local;
  for (const auto& [name, primes] : terms::symbols(t)) {
    if (primes != 0) throw Refusal(what + " writes " + name + "'; a law's expression reads values, not derivatives");
    Ref r;
    if (name == "t") {
      r.what = Ref::time;
    } else if (const auto it = c.slots[level].find(name); it != c.slots[level].end()) {
      r.what = Ref::column;
      r.slot = it->second;
    } else if (const auto g = c.global_slot.find(name); g != c.global_slot.end()) {
      r.what = Ref::global;
      r.slot = g->second;
    } else {
      const auto lv = c.level_of.find(name);
      if (lv != c.level_of.end())
        throw Refusal(what + " (at " + c.hierarchy.level(level).name + ") names " + name + ", which is at " +
                      c.hierarchy.level(lv->second).name + ": carry it to this level by a derived quantity (parent, sum, mean, …)");
      throw Refusal(what + " names " + name + ", which is no quantity, measure or parameter of the model");
    }
    local[name] = static_cast<int>(e.reads.size());
    e.reads.push_back(r);
  }
  e.program = terms::compile(t, [&](const std::string& n, unsigned primes) -> int {
    if (primes != 0) return -1;
    const auto it = local.find(n);
    return it == local.end() ? -1 : it->second;
  });
  return e;
}

}  // namespace detail

// The model compiled on its hierarchy: every name resolved, every expression compiled for its level, the derived quantities
// ordered so that each depends only on those before it. Refuses, with what to do, a model that cannot run.
inline Compiled compile(const Hierarchy& h, const Model& m) {
  Compiled c;
  c.hierarchy = h;
  c.model = m;
  if (!(m.step > 0)) throw Refusal("the model's step must be positive");
  const std::size_t L = h.levels();
  c.columns.assign(L, {});
  c.slots.assign(L, {});
  std::set<std::string> names;
  const auto claim = [&](const std::string& n, const std::string& what) {
    if (n.empty()) throw Refusal(what + " has no name");
    if (n == "t") throw Refusal("t is the time; " + what + " needs another name");
    if (!names.insert(n).second) throw Refusal("the name " + n + " is given twice");
  };
  const auto add_column = [&](std::size_t l, const std::string& n) {
    c.slots[l][n] = c.columns[l].size();
    c.columns[l].push_back(n);
  };
  for (const auto& q : m.quantities) {
    claim(q.name, "a quantity");
    const std::size_t l = h.index(q.level);
    c.level_of[q.name] = l;
    add_column(l, q.name);
    if (!(q.lower <= q.upper)) throw Refusal("the bounds of " + q.name + " are the wrong way round");
    if (!std::is_sorted(q.thresholds.begin(), q.thresholds.end())) throw Refusal("the thresholds of " + q.name + " must increase");
    if (q.kind == Kind::ordinal && q.thresholds.empty()) throw Refusal("the ordinal quantity " + q.name + " needs its thresholds");
  }
  for (const auto& d : m.derived) {
    claim(d.name, "a derived quantity");
    const std::size_t l = h.index(d.level);
    c.level_of[d.name] = l;
    add_column(l, d.name);
  }
  for (const auto& p : m.parameters) {
    claim(p.name, "a parameter");
    if (!(p.lower <= p.upper)) throw Refusal("the bounds of the parameter " + p.name + " are the wrong way round");
    if (p.fit && !(std::isfinite(p.lower) || std::isfinite(p.upper) || std::isfinite(p.value)))
      throw Refusal("the unknown " + p.name + " needs bounds or a starting value");
    c.global_slot[p.name] = c.globals.size();
    c.globals.push_back(p.name);
  }
  for (std::size_t l = 0; l < L; l++)
    for (const auto& [mname, values] : h.level(l).measures) {
      if (names.count(mname)) throw Refusal("the measure " + mname + " of " + h.level(l).name + " has the name of a quantity or parameter");
      add_column(l, mname);
    }
  for (const auto& q : m.quantities)
    if (q.kind == Kind::intensive && !q.weight.empty() && !c.slots[c.level_of[q.name]].count(q.weight))
      throw Refusal("the weight " + q.weight + " of " + q.name + " is not a quantity or measure of its level");

  // Derived quantities: compile, then order by dependence.
  std::vector<CompiledDerived> made;
  std::vector<std::set<std::string>> needs;
  for (std::size_t i = 0; i < m.derived.size(); i++) {
    const Derived& d = m.derived[i];
    CompiledDerived cd;
    cd.index = i;
    cd.level = c.level_of[d.name];
    cd.slot = c.slots[cd.level][d.name];
    cd.op = d.op;
    const std::string what = "the derived quantity " + d.name + " (" + name(d.op) + ")";
    std::size_t source = cd.level;
    switch (d.op) {
      case Op::sum:
      case Op::mean:
      case Op::wmean:
      case Op::minimum:
      case Op::maximum:
      case Op::reduce:
        if (cd.level == 0) throw Refusal(what + " reduces the units beneath, and " + h.level(0).name + " has none");
        source = cd.level - 1;
        break;
      case Op::parent:
        if (h.top(cd.level)) throw Refusal(what + " reads the level above, and " + h.level(cd.level).name + " is the top");
        source = cd.level + 1;
        break;
      default:
        break;
    }
    cd.of = detail::compile_at(c, source, d.of, what);
    std::set<std::string> uses;
    const auto note_uses = [&](const std::string& text) {
      for (const auto& [n, primes] : terms::symbols(language::read_term(text))) uses.insert(n);
    };
    note_uses(d.of);
    if (d.op == Op::wmean || d.op == Op::reduce) {
      if (d.op == Op::wmean && d.weight.empty()) throw Refusal(what + " needs its weights");
      if (!d.weight.empty()) {
        cd.weight = detail::compile_at(c, source, d.weight, what + "'s weight");
        cd.has_weight = true;
        note_uses(d.weight);
      }
    }
    if (d.op == Op::reduce) {
      const auto it = m.reducers.find(d.reducer);
      if (it == m.reducers.end()) throw Refusal(what + " names the reducer " + d.reducer + ", which the model has not been given");
      cd.reducer = it->second;
    }
    if (d.op == Op::laplacian) {
      cd.measure = d.weight;
      if (!cd.measure.empty()) (void)h.measure(cd.level, cd.measure, 0);
    }
    if (d.op == Op::inflow) {
      if (d.weight.empty() || d.costs.empty() || d.beta.empty()) throw Refusal(what + " needs the attractiveness, the costs and β");
      cd.weight = detail::compile_at(c, cd.level, d.weight, what + "'s attractiveness");
      cd.has_weight = true;
      note_uses(d.weight);
      cd.beta = detail::compile_at(c, cd.level, d.beta, what + "'s β");
      if (!cd.beta.global_only()) throw Refusal(what + "'s β must be a term of parameters and time, one value for every unit");
      note_uses(d.beta);
      (void)c.hierarchy.costs(cd.level, d.costs);
      cd.costs = d.costs;
    }
    made.push_back(std::move(cd));
    needs.push_back(uses);
  }
  // Order: Kahn's algorithm on the dependence of derived quantities on one another.
  std::vector<char> placed(made.size(), 0);
  std::map<std::string, std::size_t> which;
  for (std::size_t i = 0; i < m.derived.size(); i++) which[m.derived[i].name] = i;
  for (std::size_t round = 0; round < made.size(); round++) {
    bool progress = false;
    for (std::size_t i = 0; i < made.size(); i++) {
      if (placed[i]) continue;
      bool ready = true;
      for (const auto& n : needs[i]) {
        const auto it = which.find(n);
        if (it != which.end() && !placed[it->second]) ready = false;
      }
      if (ready) {
        c.derived.push_back(made[i]);
        placed[i] = 1;
        progress = true;
      }
    }
    if (!progress) break;
  }
  if (c.derived.size() != made.size()) {
    std::string cyc;
    for (std::size_t i = 0; i < made.size(); i++)
      if (!placed[i]) cyc += (cyc.empty() ? "" : ", ") + m.derived[i].name;
    throw Refusal("the derived quantities " + cyc + " depend on one another in a circle; break it with a state and an equilibrium law");
  }

  // Laws.
  std::set<std::string> governed;
  for (std::size_t i = 0; i < m.laws.size(); i++) {
    const Law& law = m.laws[i];
    const Quantity* q = c.quantity(law.target);
    if (!q) throw Refusal("a law advances " + law.target + ", which is no quantity of the model");
    if (q->role != Role::state) throw Refusal("a law advances " + law.target + ", which is a " + name(q->role) + ", not a state");
    if (!governed.insert(law.target).second) throw Refusal("the state " + law.target + " has two laws");
    CompiledLaw cl;
    cl.index = i;
    cl.level = c.level_of[law.target];
    cl.slot = c.slots[cl.level][law.target];
    cl.form = law.form;
    cl.rhs = detail::compile_at(c, cl.level, law.expression, "the law of " + law.target);
    if (law.form != Form::equilibrium) {
      const double s = law.step > 0 ? law.step : m.step;
      const double k = s / m.step;
      const double r = std::round(k);
      if (r < 1 || std::abs(k - r) > 1e-9 * std::max(1.0, k))
        throw Refusal("the step of the law of " + law.target + " (" + std::to_string(s) + ") is not a whole multiple of the model's step");
      cl.every = static_cast<std::size_t>(r);
    } else {
      if (!(law.damping > 0 && law.damping <= 1)) throw Refusal("the damping of the equilibrium of " + law.target + " must be in (0, 1]");
      cl.damping = law.damping;
    }
    c.laws.push_back(std::move(cl));
  }
  for (const auto& q : m.quantities)
    if (q.role == Role::state && !governed.count(q.name))
      throw Refusal("the state " + q.name + " has no law: give it one, or make it a driver or fixed");
  return c;
}

}  // namespace kinesis::scales
