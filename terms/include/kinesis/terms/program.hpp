// Programs: terms compiled for evaluation many times.
//
// A term is compiled once into a straight-line program for a stack machine, with each symbol resolved to a slot
// of an array of values (a state's components, the parameters, the time, the inputs, the noise). Evaluation then
// touches no name and allocates nothing. A product of a few integer powers is evaluated by multiplication, not by
// pow. The value of a program is the value of its term (the check holds the two to rounding).
//
// A symbol that the resolver does not know, an application of a name that is not a function, and an unknown are
// refused at compilation with the name: a program never evaluates what has no meaning.
#pragma once

#include <kinesis/terms/term.hpp>

#include <cstdint>
#include <span>

namespace kinesis::terms {

enum class Op : std::uint8_t { constant, variable, add, multiply, power, integer_power, function };

struct Instruction {
  Op op;
  std::int32_t a = 0;
};

// The slot of a symbol in the array of values, or −1 when it has none.
using Resolver = std::function<int(const std::string& name, unsigned primes)>;

class Program {
 public:
  Program() = default;

  // The value on `values`; `stack` must hold depth() doubles.
  double run(const double* values, double* stack) const {
    double* top = stack;
    for (const Instruction& i : code_) {
      switch (i.op) {
        case Op::constant: *top++ = constants_[static_cast<std::size_t>(i.a)]; break;
        case Op::variable: *top++ = values[i.a]; break;
        case Op::add: top[-2] += top[-1], top--; break;
        case Op::multiply: top[-2] *= top[-1], top--; break;
        case Op::power: top[-2] = std::pow(top[-2], top[-1]), top--; break;
        case Op::integer_power: {
          double b = top[-1], r = 1;
          for (int e = i.a < 0 ? -i.a : i.a; e > 0; e >>= 1, b *= b)
            if (e & 1) r *= b;
          top[-1] = i.a < 0 ? 1 / r : r;
          break;
        }
        case Op::function: {
          const Function& f = functions()[static_cast<std::size_t>(i.a)];
          top -= f.arity;
          *top = f.value(top);
          top++;
          break;
        }
      }
    }
    return top[-1];
  }
  double operator()(std::span<const double> values) const {
    double small[32];
    if (depth_ <= 32) return run(values.data(), small);
    std::vector<double> stack(depth_);
    return run(values.data(), stack.data());
  }
  std::size_t depth() const { return depth_; }
  std::size_t size() const { return code_.size(); }
  bool constant() const { return constant_; }  // it reads no slot

 private:
  friend Program compile(const TermPtr&, const Resolver&);
  std::vector<Instruction> code_;
  std::vector<double> constants_;
  std::size_t depth_ = 1;
  bool constant_ = true;
};

inline Program compile(const TermPtr& term, const Resolver& slot) {
  Program p;
  std::size_t depth = 0, deepest = 0;
  const auto push = [&](Op op, std::int32_t a, int grows) {
    p.code_.push_back({op, a});
    depth = static_cast<std::size_t>(static_cast<long>(depth) + grows);
    deepest = std::max(deepest, depth);
  };
  std::function<void(const TermPtr&)> emit = [&](const TermPtr& t) {
    switch (t->kind) {
      case Kind::number:
        p.constants_.push_back(t->value);
        push(Op::constant, static_cast<std::int32_t>(p.constants_.size() - 1), +1);
        break;
      case Kind::symbol: {
        const int s = slot(t->name, t->primes);
        if (s < 0) throw Refusal("'" + t->name + std::string(t->primes, '\'') + "' is not declared");
        p.constant_ = false;
        push(Op::variable, s, +1);
        break;
      }
      case Kind::sum:
      case Kind::product:
        for (std::size_t i = 0; i < t->args.size(); i++) {
          emit(t->args[i]);
          if (i > 0) push(t->kind == Kind::sum ? Op::add : Op::multiply, 0, -1);
        }
        break;
      case Kind::power:
        emit(t->args[0]);
        if (is_integer(t->args[1]) && std::fabs(t->args[1]->value) <= 64) {
          push(Op::integer_power, static_cast<std::int32_t>(t->args[1]->value), 0);
        } else {
          emit(t->args[1]);
          push(Op::power, 0, -1);
        }
        break;
      case Kind::call: {
        const int f = function_index(t->name, t->args.size());
        if (f < 0) throw Refusal("'" + t->name + "' is not a function of " + std::to_string(t->args.size()) + " argument(s)");
        for (const auto& a : t->args) emit(a);
        push(Op::function, f, 1 - static_cast<int>(t->args.size()));
        break;
      }
      case Kind::hole: throw Refusal("the unknown ?" + t->name + " has no value yet");
    }
  };
  emit(term);
  p.depth_ = std::max<std::size_t>(deepest, 1);
  return p;
}

// Slots by name: the frame of a system's values.
class Frame {
 public:
  int add(const std::string& name, unsigned primes = 0) {
    const auto key = std::make_pair(name, primes);
    const auto it = slots_.find(key);
    if (it != slots_.end()) return it->second;
    const int s = static_cast<int>(slots_.size());
    slots_.emplace(key, s);
    names_.push_back(name + std::string(primes, '\''));
    return s;
  }
  int find(const std::string& name, unsigned primes = 0) const {
    const auto it = slots_.find({name, primes});
    return it == slots_.end() ? -1 : it->second;
  }
  std::size_t size() const { return slots_.size(); }
  const std::vector<std::string>& names() const { return names_; }
  Resolver resolver() const {
    return [this](const std::string& name, unsigned primes) { return find(name, primes); };
  }

 private:
  std::map<std::pair<std::string, unsigned>, int> slots_;
  std::vector<std::string> names_;
};

}  // namespace kinesis::terms
