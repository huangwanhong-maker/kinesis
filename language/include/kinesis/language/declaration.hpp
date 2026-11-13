// A script as it was written: systems, each a declaration in parts or a derivation from another.
//
// These structures hold what the reader found, in the order and under the names of the grammar, with every
// expression as a term. They decide nothing: which kind of system a declaration is, and whether its names
// resolve, is decided by kinesis::systems (model.hpp: form_of and scope).
#pragma once

#include <kinesis/terms/term.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace kinesis::language {

using terms::TermPtr;

// A place in the text of a script, counted from 1.
struct Place {
  std::size_t line = 0, column = 0;
};

// An error in a script, with its place: of syntax (from the parser's chart) or of meaning (from elaboration).
struct ScriptError : kernel::Refusal {
  Place at;
  ScriptError(Place p, const std::string& what)
      : kernel::Refusal((p.line ? "line " + std::to_string(p.line) + ", column " + std::to_string(p.column) + ": " : std::string()) + what), at(p) {}
};

// Where a variable ranges: the reals, the naturals, an interval, a group, or a named or applied geometry.
struct Domain {
  enum class Kind { real, natural, integer, interval, group, named, applied } kind = Kind::real;
  std::string name;            // group, named, applied
  std::vector<TermPtr> args;   // applied: lattice(L, L)
  TermPtr lower, upper;        // interval
};

struct Declarator {
  std::string name;
  std::vector<std::string> arguments;  // a path q(t), a field phi(x), a function beta(t)
};

struct Space {
  bool words = false;                 // words over an alphabet
  std::vector<std::string> alphabet;
  std::string name;                   // real, complex, circle, density, …
  std::vector<TermPtr> args;          // density(2)
  std::string exponent;               // real^3, real^d: the text of the exponent, or empty
};

struct Ranging {
  std::vector<std::string> names;
  Domain domain;
};

struct States {
  std::vector<Declarator> names;
  Space space;
  std::vector<Ranging> ranges;
  std::vector<std::string> start;  // of a word
  Place at;
};

struct Binding {
  Declarator name;
  TermPtr value;
  Place at;
};

struct Noise {
  std::string name;
  std::vector<TermPtr> arguments;  // X(0)
  bool law = false;                // '~': the variable has this law; ':' the source is this process
  TermPtr source;
  Place at;
};

struct Unknown {
  enum class Kind { number, function, sparse } kind = Kind::number;
  std::string name;
  TermPtr library;  // sparse
};

struct Input {
  std::string name;
  TermPtr fallback;  // its value when nothing is connected, or null
};

struct Equation {
  std::string name;
  unsigned primes = 0;
  std::vector<TermPtr> arguments;
  TermPtr value;
  Place at;
};

struct Named {
  std::string name;
  TermPtr value;
  Place at;
};

struct Rule {
  std::string name;
  std::vector<std::string> left, right;
  TermPtr weight;  // or null
  Place at;
};

struct Assignment {
  std::string name;  // with its primes: v, q'
  TermPtr value;
};

// when left ⋈ right : resets. The event happens when left − right crosses zero into the relation.
struct Event {
  TermPtr left, right;
  std::string relation;  // ">=", "<=", ">", "<"
  std::vector<Assignment> resets;
  Place at;
};

struct Change {
  std::vector<Equation> equations;
  std::vector<Named> actions, hamiltonians, jumps;
  std::vector<Rule> rules;
  std::vector<Event> events;
  bool empty() const { return equations.empty() && actions.empty() && hamiltonians.empty() && jumps.empty() && rules.empty() && events.empty(); }
};

struct Structure {
  bool timed = false;
  std::string time;      // the name of the time
  Domain time_domain;
  bool derivations = false;
  std::string stationary;  // histories stationary S
  TermPtr weight;          // histories weight exp(-S)
  std::string weights;     // boolean, probability, real, complex, operators, or empty
  std::string calculus;    // ito, stratonovich, or empty
};

struct Constraint {
  std::vector<Ranging> backgrounds;
  std::vector<TermPtr> invariances, conserved;
  std::vector<Equation> identities;
  bool empty() const { return backgrounds.empty() && invariances.empty() && conserved.empty() && identities.empty(); }
};

struct Observation {
  std::string name;  // or empty
  TermPtr value;
};

struct Difference {
  std::vector<Observation> observations;
  std::vector<TermPtr> up_to;
};

struct Part {
  std::string name, system;
  std::vector<Binding> overrides;
};

struct Wire {
  std::string part, input;
  TermPtr value;
  Place at;
};

struct System {
  std::string name;
  Place at;
  // A derivation: system B = vary A.
  std::string derivation, source;
  // A declaration.
  std::vector<States> states;
  std::vector<Binding> params;
  std::vector<Noise> noises;
  std::vector<Unknown> unknowns;
  std::vector<Input> inputs;
  Change change;
  Structure structure;
  Constraint constraint;
  Difference difference;
  std::vector<Part> parts;
  std::vector<Wire> wires;

  bool derived() const { return !derivation.empty(); }
};

struct Script {
  std::vector<System> systems;
  const System* find(const std::string& name) const {
    for (const auto& s : systems)
      if (s.name == name) return &s;
    return nullptr;
  }
};

}  // namespace kinesis::language
