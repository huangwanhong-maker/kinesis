// The reader: the text of a script, parsed by the grammar and built into declarations.
//
// The grammar (grammar.hpp) is compiled once by morphe's reader of BNF and every script is parsed by Earley's
// parser. A text that is not a script is refused with the place as far as which it is right and what may follow
// there, from the parser's chart (morphe/formal/bnf.hpp). A text that is a script is turned, rule by rule, into the
// declarations of declaration.hpp, with each expression built as a term in normal form: precedence is that of the
// grammar (a power binds tighter than a negation, a negation tighter than a product, a product tighter than a
// sum; a power associates to the right).
//
// write() is the inverse on what it covers: a system written out is a script that reads back to the same
// declarations, so that a system a tool induces is a word of the language.
#pragma once

#include <kinesis/language/declaration.hpp>
#include <kinesis/language/grammar.hpp>

#include <morphe/formal/bnf.hpp>

#include <cstdlib>
#include <string>
#include <string_view>

namespace kinesis::language {

namespace detail {

using Node = morphe::formal::ParseNode;

inline const morphe::formal::Bnf& grammar() {
  static const morphe::formal::Bnf g = morphe::formal::compile_bnf(grammar_text);
  return g;
}

inline Place place(std::string_view text, std::size_t offset) {
  Place p{1, 1};
  for (std::size_t i = 0; i < offset && i < text.size(); i++) {
    if (text[i] == '\n') p.line++, p.column = 1;
    else p.column++;
  }
  return p;
}

struct Builder {
  std::string_view text;

  Place at(const Node& n) const { return place(text, n.begin); }
  static std::string name_of(const Node& n) {
    const Node* c = n.child("name");
    return c ? c->text : std::string();
  }
  static std::vector<std::string> names_in(const Node& n) {  // every direct child called name
    std::vector<std::string> out;
    for (const Node* c : n.all("name")) out.push_back(c->text);
    return out;
  }

  // ---- terms ----
  TermPtr operand(const Node& n) const {
    if (n.name == "number") return terms::number(std::strtod(n.text.c_str(), nullptr));
    if (n.name == "expr") return term(n);
    if (n.name == "paren") return term(*n.child("expr"));
    if (n.name == "neg") return -operand(n.children.at(0));
    if (n.name == "power") return terms::power(operand(n.children.at(0)), operand(n.children.at(1)));
    if (n.name == "evar") {
      const Node* p = n.child("primes");
      return terms::symbol(name_of(n), p ? static_cast<unsigned>(p->text.size()) : 0u);
    }
    if (n.name == "qvar") {  // a part's state, with its primes: a.q, a.q'
      const auto parts = names_in(*n.child("qname"));
      const Node* p = n.child("primes");
      return terms::symbol(parts.at(0) + "." + parts.at(1), p ? static_cast<unsigned>(p->text.size()) : 0u);
    }
    if (n.name == "qname") {
      const auto parts = names_in(n);
      return terms::symbol(parts.at(0) + "." + parts.at(1));
    }
    if (n.name == "eapp") return terms::call(name_of(n), arguments(n.child("args")));
    if (n.name == "hole") return terms::hole(name_of(n), arguments(n.child("args")));
    throw ScriptError(at(n), "the reader does not know the node '" + n.name + "'");
  }
  std::vector<TermPtr> arguments(const Node* args) const {
    std::vector<TermPtr> out;
    if (args)
      for (const Node* e : args->all("expr")) out.push_back(term(*e));
    return out;
  }
  // An expression: operands joined by '*' and '/' into products, and those by '+' and '-' into a sum.
  TermPtr term(const Node& expr) const {
    std::vector<TermPtr> parts;
    TermPtr current;
    char joining = '*';
    bool negative = false;
    const auto close = [&] {
      if (current) parts.push_back(negative ? -current : current);
      current = nullptr;
    };
    for (const Node& c : expr.children) {
      if (c.name == "addop") {
        close();
        negative = c.text == "-";
      } else if (c.name == "mulop") {
        joining = c.text[0];
      } else {
        const TermPtr v = operand(c);
        current = !current ? v : joining == '*' ? current * v : current / v;
      }
    }
    close();
    return parts.size() == 1 ? parts.front() : terms::sum(std::move(parts));
  }

  // ---- the parts of a declaration ----
  Declarator declarator(const Node& decl) const {
    Declarator d;
    d.name = name_of(decl);
    if (const Node* ns = decl.child("names")) d.arguments = names_in(*ns);
    return d;
  }
  Domain domain(const Node& n) const {
    Domain d;
    if (const Node* i = n.child("interval")) {
      const auto es = i->all("expr");
      d.kind = Domain::Kind::interval, d.lower = term(*es.at(0)), d.upper = term(*es.at(1));
    } else if (const Node* g = n.child("group")) {
      d.kind = Domain::Kind::group, d.name = name_of(*g);
    } else if (const Node* a = n.child("eapp")) {
      d.kind = Domain::Kind::applied, d.name = name_of(*a), d.args = arguments(a->child("args"));
    } else {
      d.name = name_of(n);
      d.kind = d.name == "real" ? Domain::Kind::real : d.name == "natural" ? Domain::Kind::natural : d.name == "integer" ? Domain::Kind::integer : Domain::Kind::named;
    }
    return d;
  }
  Ranging ranging(const Node& n) const { return Ranging{names_in(n), domain(*n.child("domain"))}; }
  static std::vector<std::string> word(const Node* w) { return w ? names_in(*w) : std::vector<std::string>(); }
  Binding binding(const Node& n) const { return Binding{declarator(*n.child("decl")), term(*n.child("expr")), at(n)}; }
  Equation equation(const Node& n) const {
    const Node& lhs = *n.child("lhs");
    Equation e;
    e.name = name_of(lhs);
    if (const Node* p = lhs.child("primes")) e.primes = static_cast<unsigned>(p->text.size());
    e.arguments = arguments(lhs.child("args"));
    e.value = term(*n.child("expr"));
    e.at = at(n);
    return e;
  }

  void states(const Node& n, System& s) const {
    States st;
    st.at = at(n);
    for (const Node* d : n.child("decls")->all("decl")) st.names.push_back(declarator(*d));
    const Node& sp = *n.child("space");
    if (const Node* w = sp.child("wordspace")) {
      st.space.words = true;
      st.space.alphabet = names_in(*w->child("set")->child("names"));
    } else {
      st.space.name = name_of(sp);
      st.space.args = arguments(sp.child("args"));
      if (const Node* e = sp.child("exponent")) st.space.exponent = e->text;
    }
    for (const Node* q : n.all("qualifier")) {
      if (const Node* r = q->child("ranging")) st.ranges.push_back(ranging(*r));
      if (const Node* b = q->child("start")) st.start = word(b->child("word"));
    }
    s.states.push_back(std::move(st));
  }
  void change(const Node& n, System& s) const {
    for (const Node& c : n.children) {
      if (c.name == "equation") {
        s.change.equations.push_back(equation(c));
      } else if (c.name == "action" || c.name == "hamiltonian" || c.name == "jump") {
        Named v{name_of(c), term(*c.child("expr")), at(c)};
        (c.name == "action" ? s.change.actions : c.name == "hamiltonian" ? s.change.hamiltonians : s.change.jumps).push_back(std::move(v));
      } else if (c.name == "rule") {
        const auto words = c.all("word");
        Rule r{name_of(c), word(words.at(0)), word(words.at(1)), nullptr, at(c)};
        if (const Node* w = c.child("expr")) r.weight = term(*w);
        s.change.rules.push_back(std::move(r));
      } else if (c.name == "event") {
        const auto es = c.all("expr");
        Event e{term(*es.at(0)), term(*es.at(1)), c.child("relop")->text, {}, at(c)};
        for (const Node* a : c.all("assign")) e.resets.push_back({name_of(*a) + (a->child("primes") ? a->child("primes")->text : std::string()), term(*a->child("expr"))});
        s.change.events.push_back(std::move(e));
      }
    }
  }
  void structure(const Node& n, System& s) const {
    for (const Node& c : n.children) {
      if (c.name == "time") {
        s.structure.timed = true, s.structure.time = name_of(c), s.structure.time_domain = domain(*c.child("domain"));
      } else if (c.name == "derivations") {
        s.structure.derivations = true;
      } else if (c.name == "histories") {
        if (const Node* st = c.child("stationary")) s.structure.stationary = name_of(*st);
        if (const Node* w = c.child("weight")) s.structure.weight = term(*w->child("expr"));
      } else if (c.name == "weights") {
        s.structure.weights = c.child("wkind")->text;
      } else if (c.name == "calculus") {
        s.structure.calculus = c.child("ckind")->text;
      }
    }
  }
  void constraint(const Node& n, System& s) const {
    for (const Node& c : n.children) {
      if (c.name == "background") s.constraint.backgrounds.push_back(ranging(c));
      else if (c.name == "invariance") s.constraint.invariances.push_back(term(*c.child("expr")));
      else if (c.name == "conserved") s.constraint.conserved.push_back(term(*c.child("expr")));
      else if (c.name == "equation") s.constraint.identities.push_back(equation(c));
    }
  }
  void difference(const Node& n, System& s) const {
    for (const Node& c : n.children) {
      if (c.name == "observe") {
        for (const Node* o : c.all("obs")) s.difference.observations.push_back({name_of(*o), term(*o->child("expr"))});
      } else if (c.name == "upto") {
        s.difference.up_to.push_back(term(*c.child("expr")));
      }
    }
  }

  System system(const Node& n) const {
    System s;
    s.name = name_of(n);
    s.at = at(n);
    if (const Node* d = n.child("derivation")) {
      s.derivation = d->child("dop")->text, s.source = name_of(*d);
      return s;
    }
    for (const Node& c : n.child("body")->children) {
      if (c.name == "state") {
        states(c, s);
      } else if (c.name == "param") {
        for (const Node* b : c.all("binding")) s.params.push_back(binding(*b));
      } else if (c.name == "noise") {
        const Node& target = *c.child("target");
        s.noises.push_back({name_of(target), arguments(target.child("args")), c.child("nkind")->text == "~", term(*c.child("expr")), at(c)});
      } else if (c.name == "unknown") {
        const Node& k = *c.child("ukind");
        const Unknown::Kind kind = k.text == "number" ? Unknown::Kind::number : k.text == "function" ? Unknown::Kind::function : Unknown::Kind::sparse;
        for (const auto& name : names_in(*c.child("names"))) s.unknowns.push_back({kind, name, k.child("expr") ? term(*k.child("expr")) : nullptr});
      } else if (c.name == "input") {
        for (const Node* b : c.all("ibinding")) s.inputs.push_back({name_of(*b), b->child("expr") ? term(*b->child("expr")) : nullptr});
      } else if (c.name == "part") {
        std::vector<Binding> overrides;
        for (const Node* b : c.all("binding")) overrides.push_back(binding(*b));
        for (const auto& name : names_in(*c.child("names"))) s.parts.push_back({name, name_of(c), overrides});
      } else if (c.name == "connect") {
        for (const Node* w : c.all("wire")) {
          const auto q = names_in(*w->child("qname"));
          s.wires.push_back({q.at(0), q.at(1), term(*w->child("expr")), at(*w)});
        }
      } else if (c.name == "change") {
        change(c, s);
      } else if (c.name == "structure") {
        structure(c, s);
      } else if (c.name == "constraint") {
        constraint(c, s);
      } else if (c.name == "difference") {
        difference(c, s);
      }
    }
    return s;
  }
};

}  // namespace detail

// The script a text is, or a ScriptError that says how far the text is right and what may follow.
inline Script read(std::string_view text) {
  const auto parse = detail::grammar().parse(text);
  if (!parse.accepted) {
    std::string expected;
    for (const auto& e : parse.expected) {
      if (e.rfind("[ \\t", 0) == 0 || e == "'#'") continue;  // whitespace and comments may always follow
      expected += (expected.empty() ? "" : ", ") + e;
    }
    throw ScriptError(detail::place(text, parse.furthest), "the script is right as far as here" + (expected.empty() ? std::string() : "; what may follow: " + expected));
  }
  const detail::Builder b{text};
  Script script;
  for (const auto* s : parse.tree.all("system")) {
    System sys = b.system(*s);
    if (script.find(sys.name)) throw ScriptError(sys.at, "the system '" + sys.name + "' is declared twice");
    script.systems.push_back(std::move(sys));
  }
  return script;
}

// One expression of the language, as a term.
inline TermPtr read_term(std::string_view text) {
  const auto parse = detail::grammar().parse(text, "expr");
  if (!parse.accepted) throw ScriptError(detail::place(text, parse.furthest), "the expression is right as far as here");
  return detail::Builder{text}.term(parse.tree);
}

// ---- writing ---------------------------------------------------------------------------------------------

namespace detail {
inline std::string join(const std::vector<std::string>& v, const char* sep) {
  std::string out;
  for (std::size_t i = 0; i < v.size(); i++) out += (i ? sep : "") + v[i];
  return out;
}
inline std::string written(const Declarator& d) { return d.arguments.empty() ? d.name : d.name + "(" + join(d.arguments, ", ") + ")"; }
inline std::string written(const std::vector<TermPtr>& args) {
  std::string out;
  for (std::size_t i = 0; i < args.size(); i++) out += (i ? ", " : "") + terms::text(args[i]);
  return out;
}
inline std::string written(const Domain& d) {
  switch (d.kind) {
    case Domain::Kind::real: return "real";
    case Domain::Kind::natural: return "natural";
    case Domain::Kind::integer: return "integer";
    case Domain::Kind::interval: return "[" + terms::text(d.lower) + ", " + terms::text(d.upper) + "]";
    case Domain::Kind::group: return "group " + d.name;
    case Domain::Kind::applied: return d.name + "(" + written(d.args) + ")";
    default: return d.name;
  }
}
inline std::string written(const Equation& e) {
  return e.name + std::string(e.primes, '\'') + (e.arguments.empty() ? "" : "(" + written(e.arguments) + ")") + " = " + terms::text(e.value);
}
}  // namespace detail

// The system as the language writes it.
inline std::string write(const System& s) {
  using namespace detail;
  if (s.derived()) return "system " + s.name + " = " + s.derivation + " " + s.source + "\n";
  std::string o = "system " + s.name + " {\n";
  for (const auto& st : s.states) {
    std::vector<std::string> names;
    for (const auto& d : st.names) names.push_back(written(d));
    o += "  state " + join(names, ", ") + " : ";
    if (st.space.words) o += "words over { " + join(st.space.alphabet, ", ") + " }";
    else o += st.space.name + (st.space.args.empty() ? "" : "(" + written(st.space.args) + ")") + (st.space.exponent.empty() ? "" : "^" + st.space.exponent);
    for (const auto& r : st.ranges) o += ", " + join(r.names, " ") + " in " + written(r.domain);
    if (!st.start.empty()) o += ", start " + join(st.start, " ");
    o += "\n";
  }
  if (!s.params.empty()) {
    o += "  param ";
    for (std::size_t i = 0; i < s.params.size(); i++) o += (i ? ", " : "") + written(s.params[i].name) + " = " + terms::text(s.params[i].value);
    o += "\n";
  }
  for (const auto& n : s.noises) o += "  noise " + n.name + (n.arguments.empty() ? "" : "(" + written(n.arguments) + ")") + (n.law ? " ~ " : " : ") + terms::text(n.source) + "\n";
  for (const auto& u : s.unknowns)
    o += "  unknown " + u.name + " : " + (u.kind == Unknown::Kind::number ? "number" : u.kind == Unknown::Kind::function ? "function" : "sparse " + terms::text(u.library)) + "\n";
  if (!s.inputs.empty()) {
    o += "  input ";
    for (std::size_t i = 0; i < s.inputs.size(); i++) o += (i ? ", " : "") + s.inputs[i].name + (s.inputs[i].fallback ? " = " + terms::text(s.inputs[i].fallback) : "");
    o += "\n";
  }
  for (const auto& p : s.parts) {
    o += "  part " + p.name + " : " + p.system;
    for (std::size_t i = 0; i < p.overrides.size(); i++) o += (i ? ", " : " with ") + written(p.overrides[i].name) + " = " + terms::text(p.overrides[i].value);
    o += "\n";
  }
  if (!s.change.empty()) {
    o += "  change {\n";
    for (const auto& a : s.change.actions) o += "    action " + a.name + " = " + terms::text(a.value) + "\n";
    for (const auto& a : s.change.hamiltonians) o += "    hamiltonian " + a.name + " = " + terms::text(a.value) + "\n";
    for (const auto& a : s.change.jumps) o += "    jump " + a.name + " = " + terms::text(a.value) + "\n";
    for (const auto& e : s.change.equations) o += "    " + written(e) + "\n";
    for (const auto& r : s.change.rules) o += "    " + r.name + " : " + join(r.left, " ") + " -> " + join(r.right, " ") + (r.weight ? " @ " + terms::text(r.weight) : "") + "\n";
    for (const auto& e : s.change.events) {
      o += "    when " + terms::text(e.left) + " " + e.relation + " " + terms::text(e.right) + " : ";
      for (std::size_t i = 0; i < e.resets.size(); i++) o += (i ? ", " : "") + e.resets[i].name + " -> " + terms::text(e.resets[i].value);
      o += "\n";
    }
    o += "  }\n";
  }
  {
    std::vector<std::string> items;
    const auto& st = s.structure;
    if (st.timed) items.push_back("time " + st.time + " in " + written(st.time_domain));
    if (st.derivations) items.push_back("derivations");
    if (!st.stationary.empty()) items.push_back("histories stationary " + st.stationary);
    if (st.weight) items.push_back("histories weight " + terms::text(st.weight));
    if (!st.weights.empty()) items.push_back("weights " + st.weights);
    if (!st.calculus.empty()) items.push_back("calculus " + st.calculus);
    if (!items.empty()) o += "  structure { " + join(items, " ; ") + " }\n";
  }
  if (!s.constraint.empty()) {
    std::vector<std::string> items;
    for (const auto& b : s.constraint.backgrounds) items.push_back(join(b.names, " ") + " on " + written(b.domain));
    for (const auto& t : s.constraint.invariances) items.push_back("invariant under " + terms::text(t));
    for (const auto& t : s.constraint.conserved) items.push_back("conserves " + terms::text(t));
    for (const auto& e : s.constraint.identities) items.push_back(written(e));
    o += "  constraint { " + join(items, " ; ") + " }\n";
  }
  if (!s.difference.observations.empty() || !s.difference.up_to.empty()) {
    std::vector<std::string> items;
    if (!s.difference.observations.empty()) {
      std::string line = "observe ";
      for (std::size_t i = 0; i < s.difference.observations.size(); i++) {
        const auto& ob = s.difference.observations[i];
        line += (i ? ", " : "") + (ob.name.empty() ? "" : ob.name + " = ") + terms::text(ob.value);
      }
      items.push_back(line);
    }
    for (const auto& t : s.difference.up_to) items.push_back("up to " + terms::text(t));
    o += "  difference { " + join(items, " ; ") + " }\n";
  }
  if (!s.wires.empty()) {
    o += "  connect {\n";
    for (const auto& w : s.wires) o += "    " + w.part + "." + w.input + " = " + terms::text(w.value) + "\n";
    o += "  }\n";
  }
  return o + "}\n";
}

inline std::string write(const Script& script) {
  std::string out;
  for (std::size_t i = 0; i < script.systems.size(); i++) out += (i ? "\n" : "") + write(script.systems[i]);
  return out;
}

}  // namespace kinesis::language
