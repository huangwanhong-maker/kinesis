// Terms: the expressions of the language of dynamical systems, as trees that are built in a normal form,
// differentiated, substituted into, evaluated and written out again.
//
// A term is a number, a symbol (a name with the order of its derivative in time: q, q', q''), a sum, a product,
// a power, an application f(a, …) or an unknown ?f(a, …). There is no difference and no quotient: a − b is
// a + (−1)·b and a / b is a·b⁻¹, so that one form stands for each.
//
// The normal form is what the constructors make, not a decision procedure for equality of expressions:
//   - sums and products are flat and sorted by a fixed total order (compare);
//   - a sum collects equal terms into one with a numerical coefficient, and its constant; a number times a sum
//     is distributed inside a sum (x − (y + z) is x − y − z) and kept as a product on its own (2·(x + y));
//   - a product collects equal bases into one power when the exponents are numbers, and its numerical factor;
//   - x⁰ = 1, x¹ = x, numbers are folded, (xᵃ)ᵇ = xᵃᵇ for an integer b.
// Two terms that differ after this may still be equal as functions (sin² + cos² is not folded to 1).
//
// Derivatives are symbolic and exact (derivative): the usual rules, with the derivatives of the named functions
// below. The total derivative in time (total_derivative) is D_t f = ∂f/∂t + Σ ∂f/∂s⁽ᵏ⁾ · s⁽ᵏ⁺¹⁾ over the dynamical
// symbols s, which is how the Euler–Lagrange equations of an action are derived as terms (systems/lagrangian.hpp).
//
// Written out (text), a term is a word of the language again: it parses back to the same term (the check). latex
// writes it for a page.
#pragma once

#include <kinesis/kernel/answer.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace kinesis::terms {

using kernel::Refusal;

enum class Kind { number, symbol, sum, product, power, call, hole };

struct Term;
using TermPtr = std::shared_ptr<const Term>;

struct Term {
  Kind kind = Kind::number;
  double value = 0;           // number
  std::string name;           // symbol, call, hole
  unsigned primes = 0;        // symbol: the order of its derivative in time
  std::vector<TermPtr> args;  // sum, product: the terms; power: base, exponent; call, hole: the arguments
};

// ---- the order of terms ------------------------------------------------------------------------------------

inline int compare(const TermPtr& a, const TermPtr& b) {
  if (a.get() == b.get()) return 0;
  if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
  if (a->kind == Kind::number) return a->value < b->value ? -1 : a->value > b->value ? 1 : 0;
  if (int c = a->name.compare(b->name)) return c < 0 ? -1 : 1;
  if (a->primes != b->primes) return a->primes < b->primes ? -1 : 1;
  for (std::size_t i = 0; i < std::min(a->args.size(), b->args.size()); i++)
    if (int c = compare(a->args[i], b->args[i])) return c;
  return a->args.size() < b->args.size() ? -1 : a->args.size() > b->args.size() ? 1 : 0;
}
inline bool same(const TermPtr& a, const TermPtr& b) { return compare(a, b) == 0; }
struct TermLess {
  bool operator()(const TermPtr& a, const TermPtr& b) const { return compare(a, b) < 0; }
};

// ---- constructors ------------------------------------------------------------------------------------------

inline TermPtr number(double v) {
  auto t = std::make_shared<Term>();
  t->kind = Kind::number, t->value = v;
  return t;
}
inline TermPtr symbol(std::string name, unsigned primes = 0) {
  auto t = std::make_shared<Term>();
  t->kind = Kind::symbol, t->name = std::move(name), t->primes = primes;
  return t;
}
inline TermPtr call(std::string name, std::vector<TermPtr> args) {
  auto t = std::make_shared<Term>();
  t->kind = Kind::call, t->name = std::move(name), t->args = std::move(args);
  return t;
}
inline TermPtr hole(std::string name, std::vector<TermPtr> args = {}) {
  auto t = std::make_shared<Term>();
  t->kind = Kind::hole, t->name = std::move(name), t->args = std::move(args);
  return t;
}
inline bool is_number(const TermPtr& t) { return t->kind == Kind::number; }
inline bool is_number(const TermPtr& t, double v) { return t->kind == Kind::number && t->value == v; }
inline bool is_integer(const TermPtr& t) { return is_number(t) && std::floor(t->value) == t->value && std::fabs(t->value) < 1e15; }

inline TermPtr power(TermPtr base, TermPtr exponent);
inline TermPtr product(std::vector<TermPtr> factors);

// A term as coefficient × core: 3·x·y is (3, x·y); x is (1, x); 5 is (5, 1).
inline std::pair<double, TermPtr> split_coefficient(const TermPtr& t) {
  if (is_number(t)) return {t->value, number(1)};
  if (t->kind == Kind::product && !t->args.empty() && is_number(t->args.front())) {
    std::vector<TermPtr> rest(t->args.begin() + 1, t->args.end());
    if (rest.size() == 1) return {t->args.front()->value, rest.front()};
    auto p = std::make_shared<Term>();
    p->kind = Kind::product, p->args = std::move(rest);
    return {t->args.front()->value, p};
  }
  return {1.0, t};
}

inline TermPtr sum(std::vector<TermPtr> terms) {
  std::map<TermPtr, double, TermLess> collected;
  double constant = 0;
  std::function<void(const TermPtr&)> take = [&](const TermPtr& t) {
    if (t->kind == Kind::sum) {
      for (const auto& a : t->args) take(a);
      return;
    }
    const auto [c, core] = split_coefficient(t);
    if (is_number(core)) {
      constant += c * core->value;
    } else if (core->kind == Kind::sum) {  // a number times a sum, inside a sum, is distributed: x − (y + z) is x − y − z
      for (const auto& a : core->args) take(product({number(c), a}));
    } else {
      collected[core] += c;
    }
  };
  for (const auto& t : terms) take(t);
  std::vector<TermPtr> out;
  for (const auto& [core, c] : collected) {
    if (c == 0) continue;
    out.push_back(c == 1 ? core : product({number(c), core}));
  }
  if (constant != 0 || out.empty()) out.insert(out.begin(), number(constant));
  if (out.size() == 1) return out.front();
  auto s = std::make_shared<Term>();
  s->kind = Kind::sum, s->args = std::move(out);
  return s;
}

inline TermPtr product(std::vector<TermPtr> factors) {
  // base → the sum of its numerical exponents; bases with other exponents are kept as they are.
  std::map<TermPtr, double, TermLess> powers;
  std::vector<TermPtr> others;
  double coefficient = 1;
  std::function<void(const TermPtr&)> take = [&](const TermPtr& t) {
    if (t->kind == Kind::product) {
      for (const auto& a : t->args) take(a);
    } else if (is_number(t)) {
      coefficient *= t->value;
    } else if (t->kind == Kind::power && is_number(t->args[1])) {
      powers[t->args[0]] += t->args[1]->value;
    } else if (t->kind == Kind::power) {
      others.push_back(t);
    } else {
      powers[t] += 1;
    }
  };
  for (const auto& t : factors) take(t);
  if (coefficient == 0) return number(0);
  std::vector<TermPtr> out;
  for (const auto& [base, e] : powers) {
    if (e == 0) continue;
    out.push_back(e == 1 ? base : power(base, number(e)));
  }
  for (auto& t : others) out.push_back(std::move(t));
  std::sort(out.begin(), out.end(), TermLess{});
  if (out.empty()) return number(coefficient);
  if (coefficient != 1) out.insert(out.begin(), number(coefficient));
  if (out.size() == 1) return out.front();
  auto p = std::make_shared<Term>();
  p->kind = Kind::product, p->args = std::move(out);
  return p;
}

inline TermPtr power(TermPtr base, TermPtr exponent) {
  if (is_number(exponent, 0)) return number(1);
  if (is_number(exponent, 1)) return base;
  if (is_number(base) && is_number(exponent)) {
    const double v = std::pow(base->value, exponent->value);
    if (std::isfinite(v)) return number(v);
  }
  if (is_number(base, 1)) return number(1);
  if (base->kind == Kind::power && is_integer(exponent) && is_number(base->args[1]))
    return power(base->args[0], number(base->args[1]->value * exponent->value));
  if (base->kind == Kind::product && is_integer(exponent)) {
    std::vector<TermPtr> parts;
    for (const auto& a : base->args) parts.push_back(power(a, exponent));
    return product(std::move(parts));
  }
  auto p = std::make_shared<Term>();
  p->kind = Kind::power, p->args = {std::move(base), std::move(exponent)};
  return p;
}

inline TermPtr operator+(const TermPtr& a, const TermPtr& b) { return sum({a, b}); }
inline TermPtr operator*(const TermPtr& a, const TermPtr& b) { return product({a, b}); }
inline TermPtr operator-(const TermPtr& a) { return product({number(-1), a}); }
inline TermPtr operator-(const TermPtr& a, const TermPtr& b) { return sum({a, -b}); }
inline TermPtr operator/(const TermPtr& a, const TermPtr& b) { return product({a, power(b, number(-1))}); }

// ---- the named functions -----------------------------------------------------------------------------------

struct Function {
  const char* name;
  int arity;
  double (*value)(const double*);
};

inline const std::vector<Function>& functions() {
  static const std::vector<Function> table = {
      {"sin", 1, [](const double* a) { return std::sin(a[0]); }},
      {"cos", 1, [](const double* a) { return std::cos(a[0]); }},
      {"tan", 1, [](const double* a) { return std::tan(a[0]); }},
      {"exp", 1, [](const double* a) { return std::exp(a[0]); }},
      {"log", 1, [](const double* a) { return std::log(a[0]); }},
      {"sqrt", 1, [](const double* a) { return std::sqrt(a[0]); }},
      {"sinh", 1, [](const double* a) { return std::sinh(a[0]); }},
      {"cosh", 1, [](const double* a) { return std::cosh(a[0]); }},
      {"tanh", 1, [](const double* a) { return std::tanh(a[0]); }},
      {"asin", 1, [](const double* a) { return std::asin(a[0]); }},
      {"acos", 1, [](const double* a) { return std::acos(a[0]); }},
      {"atan", 1, [](const double* a) { return std::atan(a[0]); }},
      {"abs", 1, [](const double* a) { return std::fabs(a[0]); }},
      {"sign", 1, [](const double* a) { return a[0] > 0 ? 1.0 : a[0] < 0 ? -1.0 : 0.0; }},
      {"step", 1, [](const double* a) { return a[0] >= 0 ? 1.0 : 0.0; }},
      {"floor", 1, [](const double* a) { return std::floor(a[0]); }},
      {"atan2", 2, [](const double* a) { return std::atan2(a[0], a[1]); }},
      {"min", 2, [](const double* a) { return std::min(a[0], a[1]); }},
      {"max", 2, [](const double* a) { return std::max(a[0], a[1]); }},
      {"mod", 2, [](const double* a) { return a[0] - a[1] * std::floor(a[0] / a[1]); }},
  };
  return table;
}
inline int function_index(const std::string& name, std::size_t arity) {
  const auto& t = functions();
  for (std::size_t i = 0; i < t.size(); i++)
    if (name == t[i].name && static_cast<std::size_t>(t[i].arity) == arity) return static_cast<int>(i);
  return -1;
}

// ---- reading a term ----------------------------------------------------------------------------------------

// Every symbol of a term, with the orders of its derivative that occur.
inline void symbols(const TermPtr& t, std::set<std::pair<std::string, unsigned>>& out) {
  if (t->kind == Kind::symbol) out.insert({t->name, t->primes});
  for (const auto& a : t->args) symbols(a, out);
}
inline std::set<std::pair<std::string, unsigned>> symbols(const TermPtr& t) {
  std::set<std::pair<std::string, unsigned>> out;
  symbols(t, out);
  return out;
}
inline bool depends_on(const TermPtr& t, const std::string& name, unsigned primes) {
  if (t->kind == Kind::symbol) return t->name == name && t->primes == primes;
  for (const auto& a : t->args)
    if (depends_on(a, name, primes)) return true;
  return false;
}
// The unknowns of a term and the applications of names that are not functions.
inline void holes(const TermPtr& t, std::set<std::string>& out) {
  if (t->kind == Kind::hole) out.insert(t->name);
  for (const auto& a : t->args) holes(a, out);
}
inline void applications(const TermPtr& t, std::set<std::pair<std::string, std::size_t>>& out) {
  if (t->kind == Kind::call && function_index(t->name, t->args.size()) < 0) out.insert({t->name, t->args.size()});
  for (const auto& a : t->args) applications(a, out);
}

// ---- substitution ------------------------------------------------------------------------------------------

// Each node is offered to `replace`; where it answers a term, that term stands for the node, and otherwise the
// node is rebuilt from its parts, through the constructors, so that the result is in normal form.
inline TermPtr rewrite(const TermPtr& t, const std::function<TermPtr(const TermPtr&)>& replace) {
  if (TermPtr r = replace(t)) return r;
  if (t->args.empty()) return t;
  std::vector<TermPtr> args;
  args.reserve(t->args.size());
  for (const auto& a : t->args) args.push_back(rewrite(a, replace));
  switch (t->kind) {
    case Kind::sum: return sum(std::move(args));
    case Kind::product: return product(std::move(args));
    case Kind::power: return power(args[0], args[1]);
    case Kind::call: return call(t->name, std::move(args));
    case Kind::hole: return hole(t->name, std::move(args));
    default: return t;
  }
}
inline TermPtr substitute(const TermPtr& t, const std::string& name, unsigned primes, const TermPtr& by) {
  return rewrite(t, [&](const TermPtr& x) -> TermPtr { return x->kind == Kind::symbol && x->name == name && x->primes == primes ? by : nullptr; });
}

// ---- derivatives -------------------------------------------------------------------------------------------

inline TermPtr derivative(const TermPtr& t, const std::string& name, unsigned primes = 0);

namespace detail {
// The derivative of a named function of one argument at u, as a term in u; null when it is not known.
inline TermPtr outer_derivative(const std::string& f, const TermPtr& u) {
  if (f == "sin") return call("cos", {u});
  if (f == "cos") return -call("sin", {u});
  if (f == "tan") return power(call("cos", {u}), number(-2));
  if (f == "exp") return call("exp", {u});
  if (f == "log") return power(u, number(-1));
  if (f == "sqrt") return product({number(0.5), power(call("sqrt", {u}), number(-1))});
  if (f == "sinh") return call("cosh", {u});
  if (f == "cosh") return call("sinh", {u});
  if (f == "tanh") return number(1) - power(call("tanh", {u}), number(2));
  if (f == "asin") return power(number(1) - power(u, number(2)), number(-0.5));
  if (f == "acos") return -power(number(1) - power(u, number(2)), number(-0.5));
  if (f == "atan") return power(number(1) + power(u, number(2)), number(-1));
  if (f == "abs") return call("sign", {u});
  if (f == "sign" || f == "step" || f == "floor") return number(0);
  return nullptr;
}
}  // namespace detail

inline TermPtr derivative(const TermPtr& t, const std::string& name, unsigned primes) {
  switch (t->kind) {
    case Kind::number: return number(0);
    case Kind::symbol: return number(t->name == name && t->primes == primes ? 1 : 0);
    case Kind::sum: {
      std::vector<TermPtr> parts;
      for (const auto& a : t->args) parts.push_back(derivative(a, name, primes));
      return sum(std::move(parts));
    }
    case Kind::product: {
      std::vector<TermPtr> parts;
      for (std::size_t i = 0; i < t->args.size(); i++) {
        const TermPtr d = derivative(t->args[i], name, primes);
        if (is_number(d, 0)) continue;
        std::vector<TermPtr> f{d};
        for (std::size_t j = 0; j < t->args.size(); j++)
          if (j != i) f.push_back(t->args[j]);
        parts.push_back(product(std::move(f)));
      }
      return sum(std::move(parts));
    }
    case Kind::power: {
      const TermPtr &b = t->args[0], &e = t->args[1];
      const TermPtr db = derivative(b, name, primes);
      if (is_number(e)) return product({e, power(b, number(e->value - 1)), db});
      const TermPtr de = derivative(e, name, primes);
      return product({t, sum({product({de, call("log", {b})}), product({e, db, power(b, number(-1))})})});
    }
    case Kind::call: {
      if (!depends_on(t, name, primes)) return number(0);
      if (t->args.size() == 1) {
        if (TermPtr outer = detail::outer_derivative(t->name, t->args[0])) return product({outer, derivative(t->args[0], name, primes)});
      }
      if (t->name == "atan2" && t->args.size() == 2) {
        const TermPtr &y = t->args[0], &x = t->args[1];
        const TermPtr r2 = power(x, number(2)) + power(y, number(2));
        return (x * derivative(y, name, primes) - y * derivative(x, name, primes)) / r2;
      }
      if ((t->name == "min" || t->name == "max") && t->args.size() == 2) {
        const TermPtr &a = t->args[0], &b = t->args[1];
        const TermPtr first = call("step", {t->name == "min" ? b - a : a - b});  // 1 where the first argument is taken
        return first * derivative(a, name, primes) + (number(1) - first) * derivative(b, name, primes);
      }
      if (t->name == "mod" && t->args.size() == 2) {  // a − b floor(a / b), where it is differentiable
        const TermPtr &a = t->args[0], &b = t->args[1];
        return derivative(a, name, primes) - derivative(b, name, primes) * call("floor", {a / b});
      }
      throw Refusal("the derivative of " + t->name + "(…) is not known");
    }
    case Kind::hole:
      if (!depends_on(t, name, primes)) return number(0);
      throw Refusal("the derivative of the unknown ?" + t->name + " is not known");
  }
  return number(0);
}

// D_t f: the partial derivative in `time`, and for every symbol s⁽ᵏ⁾ whose name is dynamical, ∂f/∂s⁽ᵏ⁾ · s⁽ᵏ⁺¹⁾.
inline TermPtr total_derivative(const TermPtr& t, const std::string& time, const std::set<std::string>& dynamical) {
  std::vector<TermPtr> parts{derivative(t, time, 0)};
  for (const auto& [name, primes] : symbols(t)) {
    if (!dynamical.count(name)) continue;
    parts.push_back(product({derivative(t, name, primes), symbol(name, primes + 1)}));
  }
  return sum(std::move(parts));
}

// ---- values ------------------------------------------------------------------------------------------------

// The value of a term where `lookup` gives each symbol's; a slow reading, for small uses (programs are the fast one).
inline double evaluate(const TermPtr& t, const std::function<double(const std::string&, unsigned)>& lookup) {
  switch (t->kind) {
    case Kind::number: return t->value;
    case Kind::symbol: return lookup(t->name, t->primes);
    case Kind::sum: {
      double s = 0;
      for (const auto& a : t->args) s += evaluate(a, lookup);
      return s;
    }
    case Kind::product: {
      double p = 1;
      for (const auto& a : t->args) p *= evaluate(a, lookup);
      return p;
    }
    case Kind::power: return std::pow(evaluate(t->args[0], lookup), evaluate(t->args[1], lookup));
    case Kind::call: {
      const int f = function_index(t->name, t->args.size());
      if (f < 0) throw Refusal(t->name + " is not a function of " + std::to_string(t->args.size()) + " argument(s)");
      double a[2] = {0, 0};
      for (std::size_t i = 0; i < t->args.size(); i++) a[i] = evaluate(t->args[i], lookup);
      return functions()[static_cast<std::size_t>(f)].value(a);
    }
    case Kind::hole: throw Refusal("the unknown ?" + t->name + " has no value");
  }
  return 0;
}

// ---- writing a term ----------------------------------------------------------------------------------------

inline std::string number_text(double v) {
  if (!std::isfinite(v)) throw Refusal("a number that is not finite cannot be written");
  if (std::floor(v) == v && std::fabs(v) < 1e15) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.0f", v);
    return buf;
  }
  char buf[48];
  for (int p = 1; p <= 17; p++) {
    std::snprintf(buf, sizeof buf, "%.*g", p, v);
    if (std::strtod(buf, nullptr) == v) break;
  }
  std::string s(buf);
  // The language writes a fraction with a digit before the point and an exponent as e±dd; %g already does.
  return s;
}

namespace detail {
// Precedences: 1 a sum, 2 a product, 3 a negation, 4 a power, 5 an atom.
struct Written {
  std::string text;
  int precedence = 5;
};
inline std::string wrap(const Written& w, int needed) { return w.precedence < needed ? "(" + w.text + ")" : w.text; }

inline Written write(const TermPtr& t);

// A product as numerator / denominator, with its sign taken out: {negative, numerator, denominator}.
struct Fraction {
  bool negative = false;
  std::vector<TermPtr> numerator, denominator;
  double coefficient = 1;
};
inline Fraction fraction(const TermPtr& t) {
  Fraction f;
  const auto take = [&](const TermPtr& a) {
    if (is_number(a)) {
      f.coefficient *= a->value;
    } else if (a->kind == Kind::power && is_number(a->args[1]) && a->args[1]->value < 0) {
      f.denominator.push_back(power(a->args[0], number(-a->args[1]->value)));
    } else {
      f.numerator.push_back(a);
    }
  };
  if (t->kind == Kind::product)
    for (const auto& a : t->args) take(a);
  else
    take(t);
  if (f.coefficient < 0) f.negative = true, f.coefficient = -f.coefficient;
  return f;
}
inline Written write_product(const TermPtr& t) {
  const Fraction f = fraction(t);
  std::string top;
  const bool unit = f.coefficient == 1;
  if (!unit || f.numerator.empty()) top = number_text(f.coefficient);
  for (const auto& a : f.numerator) top += (top.empty() ? "" : " * ") + wrap(write(a), 3);
  std::string out = top;
  for (const auto& a : f.denominator) out += " / " + wrap(write(a), 4);
  Written w{out, (f.numerator.size() + (unit ? 0 : 1) > 1 || !f.denominator.empty()) ? 2 : 5};
  if (w.precedence == 5 && !f.numerator.empty()) w.precedence = write(f.numerator.front()).precedence;
  if (f.negative) w = {"-" + wrap(w, 2), 3};
  return w;
}
inline Written write(const TermPtr& t) {
  switch (t->kind) {
    case Kind::number:
      return t->value < 0 ? Written{"-" + number_text(-t->value), 3} : Written{number_text(t->value), 5};
    case Kind::symbol: return {t->name + std::string(t->primes, '\''), 5};
    case Kind::sum: {
      std::string out;
      for (std::size_t i = 0; i < t->args.size(); i++) {
        Written w = write(t->args[i]);
        if (i == 0) {
          out = wrap(w, 2);
        } else if (w.precedence == 3 && w.text[0] == '-') {
          out += " - " + w.text.substr(1);
        } else {
          out += " + " + wrap(w, 2);
        }
      }
      return {out, 1};
    }
    case Kind::product: return write_product(t);
    case Kind::power: {
      if (is_number(t->args[1]) && t->args[1]->value < 0) return write_product(t);
      const Written e = write(t->args[1]);
      return {wrap(write(t->args[0]), 5) + "^" + (e.precedence >= 5 ? e.text : "(" + e.text + ")"), 4};
    }
    case Kind::call:
    case Kind::hole: {
      std::string out = (t->kind == Kind::hole ? "?" : "") + t->name;
      if (t->kind == Kind::call || !t->args.empty()) {
        out += "(";
        for (std::size_t i = 0; i < t->args.size(); i++) out += (i ? ", " : "") + write(t->args[i]).text;
        out += ")";
      }
      return {out, 5};
    }
  }
  return {"0", 5};
}
}  // namespace detail

// The term as the language writes it.
inline std::string text(const TermPtr& t) { return detail::write(t).text; }

namespace detail {
inline std::string latex_name(const std::string& name) {
  static const char* greek[] = {"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu",
                                "nu", "xi", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega",
                                "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi", "Psi", "Omega"};
  // A name is a stem and an index: x1 is x_{1}, k_on is k_{on}, sigma is \sigma.
  std::string stem = name, index;
  const auto under = name.find('_');
  if (under != std::string::npos && under > 0) {
    stem = name.substr(0, under), index = name.substr(under + 1);
  } else {
    std::size_t digits = name.size();
    while (digits > 0 && std::isdigit(static_cast<unsigned char>(name[digits - 1]))) digits--;
    if (digits > 0 && digits < name.size()) stem = name.substr(0, digits), index = name.substr(digits);
  }
  std::string out;
  bool is_greek = false;
  for (const char* g : greek)
    if (stem == g) is_greek = true;
  if (is_greek) out = "\\" + stem;
  else if (stem.size() == 1) out = stem;
  else out = "\\mathrm{" + stem + "}";
  if (!index.empty()) out += "_{" + (index.size() == 1 || std::isdigit(static_cast<unsigned char>(index[0])) ? index : "\\mathrm{" + index + "}") + "}";
  return out;
}
inline Written latex(const TermPtr& t);
// A number as mathematics writes it. A number of few digits is written as it is; one with many that is a
// fraction of small parts to rounding (1/24, 8/3) is written as that fraction; any other to six significant
// digits. This is for reading: text() keeps every digit, so that a term reads back as it was.
inline std::string latex_number(double v) {
  if (v < 0) return "-" + latex_number(-v);
  const std::string plain = number_text(v);
  std::size_t digits = 0;
  for (char c : plain) {
    if (c == 'e' || c == 'E') break;
    if (c >= '0' && c <= '9') digits++;
  }
  if (digits <= 7) return plain;
  for (long d = 2; d <= 1000; d++) {
    const double n = v * static_cast<double>(d), r = std::round(n);
    if (r != 0 && std::fabs(n - r) <= 1e-12 * std::max(1.0, std::fabs(r))) return "\\frac{" + number_text(r) + "}{" + std::to_string(d) + "}";
  }
  char buf[48];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}
inline std::string latex_wrap(const Written& w, int needed) { return w.precedence < needed ? "\\left(" + w.text + "\\right)" : w.text; }
inline Written latex_product(const TermPtr& t) {
  const Fraction f = fraction(t);
  auto join = [&](const std::vector<TermPtr>& parts, double coefficient) {
    if (coefficient == 1 && parts.size() == 1) return latex(parts.front()).text;  // alone above or below a bar
    std::string out;
    if (coefficient != 1 || parts.empty()) out = latex_number(coefficient);
    for (const auto& a : parts) out += (out.empty() ? "" : "\\,") + latex_wrap(latex(a), 3);
    return out;
  };
  Written w;
  if (f.denominator.empty()) {
    w = {join(f.numerator, f.coefficient), f.numerator.size() + (f.coefficient == 1 ? 0 : 1) > 1 ? 2 : 5};
    if (w.precedence == 5 && !f.numerator.empty()) w.precedence = latex(f.numerator.front()).precedence;
  } else {
    w = {"\\frac{" + join(f.numerator, f.coefficient) + "}{" + join(f.denominator, 1) + "}", 5};
  }
  if (f.negative) w = {"-" + latex_wrap(w, 2), 3};
  return w;
}
inline Written latex(const TermPtr& t) {
  switch (t->kind) {
    case Kind::number: return t->value < 0 ? Written{"-" + latex_number(-t->value), 3} : Written{latex_number(t->value), 5};
    case Kind::symbol: {
      const std::string n = latex_name(t->name);
      if (t->primes == 0) return {n, 5};
      if (t->primes == 1) return {"\\dot{" + n + "}", 5};
      if (t->primes == 2) return {"\\ddot{" + n + "}", 5};
      return {n + "^{(" + std::to_string(t->primes) + ")}", 5};
    }
    case Kind::sum: {
      std::string out;
      for (std::size_t i = 0; i < t->args.size(); i++) {
        const Written w = latex(t->args[i]);
        if (i == 0) out = latex_wrap(w, 2);
        else if (w.precedence == 3 && w.text[0] == '-') out += " - " + w.text.substr(1);
        else out += " + " + latex_wrap(w, 2);
      }
      return {out, 1};
    }
    case Kind::product: return latex_product(t);
    case Kind::power: {
      if (is_number(t->args[1]) && t->args[1]->value < 0) return latex_product(t);
      if (is_number(t->args[1], 0.5)) return {"\\sqrt{" + latex(t->args[0]).text + "}", 5};
      return {latex_wrap(latex(t->args[0]), 5) + "^{" + latex(t->args[1]).text + "}", 4};
    }
    case Kind::call:
    case Kind::hole: {
      if (t->kind == Kind::call && t->name == "sqrt" && t->args.size() == 1) return {"\\sqrt{" + latex(t->args[0]).text + "}", 5};
      if (t->kind == Kind::call && t->name == "exp" && t->args.size() == 1) return {"e^{" + latex(t->args[0]).text + "}", 4};
      if (t->kind == Kind::call && t->name == "abs" && t->args.size() == 1) return {"\\left|" + latex(t->args[0]).text + "\\right|", 5};
      static const char* known[] = {"sin", "cos", "tan", "log", "sinh", "cosh", "tanh", "min", "max"};
      std::string head = t->kind == Kind::hole ? "{?}" + latex_name(t->name) : "\\operatorname{" + t->name + "}";
      for (const char* k : known)
        if (t->kind == Kind::call && t->name == k) head = "\\" + t->name;
      if (t->kind == Kind::hole && t->args.empty()) return {head, 5};
      std::string out = head + "\\left(";
      for (std::size_t i = 0; i < t->args.size(); i++) out += (i ? ", " : "") + latex(t->args[i]).text;
      return {out + "\\right)", 5};
    }
  }
  return {"0", 5};
}
}  // namespace detail

// The term in LaTeX, for a page.
inline std::string latex(const TermPtr& t) { return detail::latex(t).text; }

}  // namespace kinesis::terms
