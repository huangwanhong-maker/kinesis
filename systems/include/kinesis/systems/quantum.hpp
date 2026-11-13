// Open quantum systems: a declaration of a Hamiltonian and jumps made into a generator that runs.
//
// The state is a density matrix on a product of factors, `state rho : density(2, 12)` being a spin one half and
// an oscillator cut at 12 levels (`hilbert(2, 12)` names the same space by its vectors; the state kept is still
// the density matrix). The Hamiltonian and the jump operators are terms of named operators, of the
// parameters and of the imaginary unit `i`:
//
//     a spin (a factor of dimension 2)       sx, sy, sz, sm, sp   (σ_−, σ_+: sm lowers up to down)
//     an oscillator (any other dimension)    a, adag, num         (num = a†a)
//     every space                            id
//
// A name refers to the first factor of its kind; with several, sx1, a2 name the factor by its place.
//
// Products. Terms of the language are commutative, and operators are not. Operators on different factors commute,
// so a product of them means what it says whatever its order. A product of two operators on one factor does not,
// and is refused: its order is written with `ordered(adag, adag, a, a)`, an application, whose arguments keep
// their order. A power of one operator is that operator multiplied by itself.
//
// The generator is physis's (physis/quantum/lindblad.hpp): dρ/dt = −i[H, ρ] + Σ (LρL† − ½{L†L, ρ}), with ħ = 1
// and each rate inside its jump operator, L = sqrt(gamma) * sm. Several Hamiltonians are summed.
//
// What is observed. The declaration's observations are operators, read as their expectations tr(ρA). With none
// declared, the populations of the basis states are observed.
#pragma once

#include <kinesis/systems/model.hpp>

#include <physis/quantum/lindblad.hpp>

#include <set>
#include <string>
#include <vector>

namespace kinesis::systems {

using physis::quantum::Complex;
using physis::quantum::Operator;

class Quantum {
 public:
  std::string name;
  std::vector<std::size_t> factors;       // the dimension of each factor
  std::size_t dimension = 1;
  physis::quantum::Lindbladian generator;
  std::vector<std::string> observables;
  std::vector<Operator> observed;

  // The basis state with each factor at the level given: for a spin 0 is up and 1 down; for an oscillator the number.
  Operator basis(const std::vector<std::size_t>& levels) const {
    if (levels.size() != factors.size()) throw kernel::Refusal(name + " has " + std::to_string(factors.size()) + " factor(s); a level is given for each");
    std::size_t index = 0;
    for (std::size_t f = 0; f < factors.size(); f++) {
      if (levels[f] >= factors[f]) throw kernel::Refusal("a factor of dimension " + std::to_string(factors[f]) + " has no level " + std::to_string(levels[f]));
      index = index * factors[f] + levels[f];
    }
    return physis::quantum::basis_state(dimension, index);
  }
  std::vector<double> expectations(const Operator& rho) const {
    std::vector<double> out;
    for (const auto& A : observed) out.push_back(physis::quantum::expectation(rho, A).real());
    return out;
  }
};

namespace detail {

// The value of a term of operators: a number, or an operator with the factors it acts on.
struct OperatorValue {
  bool is_operator = false;
  Complex scalar = 0;
  Operator op;
  std::set<std::size_t> support;
};

struct OperatorSpace {
  const Scope& sc;
  std::vector<std::size_t> factors;
  std::size_t dimension = 1;
  language::Place at;

  // The operator `local` on factor f, as an operator on the whole space.
  Operator embedded(const Operator& local, std::size_t f) const {
    Operator out = Operator::identity(1);
    for (std::size_t k = 0; k < factors.size(); k++) out = physis::quantum::kron(out, k == f ? local : Operator::identity(factors[k]));
    return out;
  }
  // The operator a name stands for, or nothing.
  bool named(const std::string& name, OperatorValue& v) const {
    if (name == "id") {
      v.is_operator = true, v.op = Operator::identity(dimension);
      return true;
    }
    // A stem and, perhaps, the place of the factor.
    std::size_t digits = name.size();
    while (digits > 0 && std::isdigit(static_cast<unsigned char>(name[digits - 1]))) digits--;
    const std::string stem = name.substr(0, digits);
    const bool spin = stem == "sx" || stem == "sy" || stem == "sz" || stem == "sm" || stem == "sp";
    const bool mode = stem == "a" || stem == "adag" || stem == "num";
    if (!spin && !mode) return false;
    std::size_t factor = factors.size();
    if (digits < name.size()) {
      const std::size_t place = static_cast<std::size_t>(std::stoul(name.substr(digits)));
      if (place < 1 || place > factors.size() || (factors[place - 1] == 2) != spin) throw ScriptError(at, "'" + name + "' names factor " + std::to_string(place) + ", which is not " + (spin ? "a spin" : "an oscillator"));
      factor = place - 1;
    } else {
      for (std::size_t k = 0; k < factors.size() && factor == factors.size(); k++)
        if ((factors[k] == 2) == spin) factor = k;
      if (factor == factors.size()) throw ScriptError(at, "'" + name + "' is an operator of " + (spin ? "a spin" : "an oscillator") + ", and the state has none");
    }
    Operator local;
    if (stem == "sx") local = physis::quantum::pauli_x();
    else if (stem == "sy") local = physis::quantum::pauli_y();
    else if (stem == "sz") local = physis::quantum::pauli_z();
    else if (stem == "sm") local = physis::quantum::lowering();
    else if (stem == "sp") local = physis::quantum::dagger(physis::quantum::lowering());
    else {
      const Operator a = physis::quantum::annihilation(factors[factor]);
      local = stem == "a" ? a : stem == "adag" ? physis::quantum::dagger(a) : physis::quantum::dagger(a) * a;
    }
    v.is_operator = true, v.op = embedded(local, factor), v.support = {factor};
    return true;
  }
  Operator as_operator(const OperatorValue& v) const { return v.is_operator ? v.op : v.scalar * Operator::identity(dimension); }

  OperatorValue value(const TermPtr& t) const {
    OperatorValue v;
    switch (t->kind) {
      case terms::Kind::number:
        v.scalar = t->value;
        return v;
      case terms::Kind::symbol: {
        const int p = sc.param_index(t->name);
        if (p >= 0 && t->primes == 0) {
          v.scalar = sc.param_values[static_cast<std::size_t>(p)];
          return v;
        }
        if (t->name == "i" && t->primes == 0) {
          v.scalar = Complex(0, 1);
          return v;
        }
        if (t->primes == 0 && named(t->name, v)) return v;
        throw ScriptError(at, "'" + t->name + std::string(t->primes, '\'') + "' is not a parameter or an operator of this state");
      }
      case terms::Kind::sum: {
        bool any = false;
        for (const auto& a : t->args) any = any || value(a).is_operator;
        if (!any) {
          for (const auto& a : t->args) v.scalar += value(a).scalar;
          return v;
        }
        v.is_operator = true, v.op = Operator(dimension);
        for (const auto& a : t->args) {
          const OperatorValue x = value(a);
          v.op = v.op + as_operator(x);
          v.support.insert(x.support.begin(), x.support.end());
        }
        return v;
      }
      case terms::Kind::product: {
        v.scalar = 1;
        for (const auto& a : t->args) {
          const OperatorValue x = value(a);
          if (!x.is_operator) {
            v.scalar *= x.scalar;
            continue;
          }
          if (!v.is_operator) {
            v.is_operator = true, v.op = x.op, v.support = x.support;
            continue;
          }
          for (std::size_t f : x.support)
            if (v.support.count(f)) throw ScriptError(at, "a product of two operators on one factor has an order, which a product of the language does not keep: write it ordered(…, …), as in ordered(adag, a)");
          v.op = v.op * x.op;
          v.support.insert(x.support.begin(), x.support.end());
        }
        if (v.is_operator) v.op = v.scalar * v.op, v.scalar = 0;
        return v;
      }
      case terms::Kind::power: {
        const OperatorValue base = value(t->args[0]), exponent = value(t->args[1]);
        if (exponent.is_operator) throw ScriptError(at, "an exponent is a number");
        if (!base.is_operator) {
          v.scalar = std::pow(base.scalar, exponent.scalar);
          return v;
        }
        const double k = exponent.scalar.real();
        if (exponent.scalar.imag() != 0 || k < 0 || std::floor(k) != k) throw ScriptError(at, "an operator is raised to a whole power that is not negative");
        v.is_operator = true, v.op = physis::quantum::power(base.op, static_cast<unsigned>(k)), v.support = base.support;
        return v;
      }
      case terms::Kind::call: {
        if (t->name == "ordered") {
          v.is_operator = true, v.op = Operator::identity(dimension);
          for (const auto& a : t->args) {
            const OperatorValue x = value(a);
            v.op = v.op * as_operator(x);
            v.support.insert(x.support.begin(), x.support.end());
          }
          return v;
        }
        std::vector<Complex> args;
        for (const auto& a : t->args) {
          const OperatorValue x = value(a);
          if (x.is_operator) throw ScriptError(at, t->name + "(…) of an operator is not known: only numbers are its arguments");
          args.push_back(x.scalar);
        }
        if (t->name == "exp" && args.size() == 1) v.scalar = std::exp(args[0]);
        else if (t->name == "sqrt" && args.size() == 1) v.scalar = std::sqrt(args[0]);
        else if (t->name == "conj" && args.size() == 1) v.scalar = std::conj(args[0]);
        else {
          const int f = terms::function_index(t->name, args.size());
          if (f < 0) throw ScriptError(at, "'" + t->name + "' is not a function of " + std::to_string(args.size()) + " argument(s)");
          double real[2] = {0, 0};
          for (std::size_t k = 0; k < args.size(); k++) {
            if (args[k].imag() != 0) throw ScriptError(at, t->name + "(…) of a number that is not real is not known");
            real[k] = args[k].real();
          }
          v.scalar = terms::functions()[static_cast<std::size_t>(f)].value(real);
        }
        return v;
      }
      case terms::Kind::hole: throw ScriptError(at, "the unknown ?" + t->name + " has no value");
    }
    return v;
  }
};

}  // namespace detail

inline Quantum quantum(const System& s) {
  if (form_of(s) != Form::quantum) throw kernel::Refusal(s.name + " is a " + std::string(name(form_of(s))) + " system, not a quantum one");
  const Scope sc = scope(s);
  Quantum q;
  q.name = s.name;
  for (const auto& st : s.states) {
    if (st.space.name != "density" && st.space.name != "hilbert") continue;
    for (const auto& a : st.space.args) {
      const double d = sc.constant(a, st.at, "the dimension of a factor");
      if (d < 2 || std::floor(d) != d || d > 64) throw ScriptError(st.at, "a factor has a dimension from 2 to 64");
      q.factors.push_back(static_cast<std::size_t>(d));
    }
  }
  if (q.factors.empty()) throw ScriptError(s.at, s.name + " has a Hamiltonian or jumps and no state of the form `state rho : density(2)`");
  for (std::size_t d : q.factors) q.dimension *= d;
  if (q.dimension > 64) throw ScriptError(s.at, "the state has dimension " + std::to_string(q.dimension) + ": the generator is a matrix of order 2n², and is built here for n up to 64");
  detail::OperatorSpace space{sc, q.factors, q.dimension, s.at};
  Operator H(q.dimension);
  for (const auto& h : s.change.hamiltonians) {
    space.at = h.at;
    const Operator part = space.as_operator(space.value(sc.resolve(h.value)));
    const Operator anti = part - physis::quantum::dagger(part);
    for (const auto& v : anti.a)
      if (std::abs(v) > 1e-12) throw ScriptError(h.at, "the Hamiltonian " + h.name + " is not Hermitian");
    H = H + part;
  }
  std::vector<Operator> jumps;
  for (const auto& j : s.change.jumps) {
    space.at = j.at;
    jumps.push_back(space.as_operator(space.value(sc.resolve(j.value))));
  }
  q.generator = physis::quantum::lindbladian(H, jumps);
  space.at = s.at;
  for (const auto& ob : s.difference.observations) {
    q.observables.push_back(ob.name.empty() ? terms::text(ob.value) : ob.name);
    q.observed.push_back(space.as_operator(space.value(sc.resolve(ob.value))));
  }
  if (q.observed.empty())
    for (std::size_t k = 0; k < q.dimension; k++) q.observables.push_back("p" + std::to_string(k)), q.observed.push_back(physis::quantum::basis_state(q.dimension, k));
  return q;
}

struct QuantumTrajectory {
  std::vector<std::string> names;                  // of the observations
  std::vector<double> t;
  std::vector<std::vector<double>> expectations;   // at each time, of each observation
  std::vector<double> purity;                      // tr ρ²
  Operator final_state;
};

// The state propagated from ρ₀ at equal steps, by the exact propagator of one step.
inline QuantumTrajectory evolve(const Quantum& q, const Operator& rho0, double t1, double step) {
  if (!(step > 0) || !(t1 >= 0)) throw kernel::Refusal("a positive step and a time that is not negative");
  const calx::linalg::Matrix P = physis::quantum::propagator(q.generator, step);
  QuantumTrajectory out;
  out.names = q.observables;
  Operator rho = rho0;
  const std::size_t count = static_cast<std::size_t>(std::llround(t1 / step));
  for (std::size_t k = 0; k <= count; k++) {
    out.t.push_back(static_cast<double>(k) * step);
    out.expectations.push_back(q.expectations(rho));
    out.purity.push_back(physis::quantum::purity(rho));
    if (k < count) rho = physis::quantum::stepped(P, rho);
  }
  out.final_state = rho;
  return out;
}

}  // namespace kinesis::systems
