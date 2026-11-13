// Checks of kinesis/language: the grammar, the reader and the writer.
//
// - The grammar compiles in morphe's reader, and the three tests the language was set (a generative grammar, a
//   group field theory with no background, a diffusion model with noise of a stated law) are read, each part
//   landing where it was written.
// - Precedence and association, on expressions whose value decides them: −2² is −4, 2³² … right association of
//   powers, a − b − c, a / b / c, a / b * c.
// - Terms round the trip: a term written and read is the same term, for terms built at random.
// - Systems round the trip: a script written from its declarations reads back to declarations that write the same.
// - Events, inputs, parts and connections, named observations, unknowns and derivations are read.
// - Refusals with their place: a product by juxtaposition, an unclosed block, a system declared twice.

#include <kinesis/language/reader.hpp>

#include <calx/verification/check.hpp>

#include <random>
#include <string>

using namespace kinesis::language;
using kinesis::terms::TermPtr;
namespace kt = kinesis::terms;

static const char* three_tests = R"(
# Test 1: a generative grammar.
system anbn {
  state w : words over { S, a, b }, start S
  change {
    grow : S -> a S b
    stop : S -> a b
  }
  structure { derivations ; weights boolean }
  difference { observe terminal }
}

system anbn_weighted {
  state w : words over { S, a, b }, start S
  change {
    grow : S -> a S b @ 0.7
    stop : S -> a b   @ 0.3
  }
  structure { derivations ; weights probability }
}

# Test 2: a dynamics with no background.
system boulatov {
  state phi(g1, g2, g3) : real, g1 g2 g3 in group G
  param lambda = 0.1
  constraint {
    phi(g1 * h, g2 * h, g3 * h) = phi(g1, g2, g3)
  }
  change {
    action S = 1/2 * integral(phi(g1, g2, g3)^2)
             + lambda/24 * integral(phi(g1, g2, g3) * phi(g1, g4, g5) * phi(g2, g5, g6) * phi(g3, g6, g4))
  }
  structure { histories weight exp(-S) ; weights real }
}

system tetrahedra = expand boulatov

system phi4 {
  state phi(x) : real
  param m = 1, lambda = 0.5
  constraint { x on lattice(L, L) }
  change { action S = sum(1/2 * grad(phi)^2 + 1/2 * m^2 * phi(x)^2 + lambda/24 * phi(x)^4) }
  structure { histories weight exp(-S) ; weights real }
}

# Test 3: a random model.
system forward {
  state X(t) : real^d
  param d = 2, T = 1, beta(t) = 0.1 + 19.9 * t
  noise W : wiener(d)
  noise X(0) ~ data
  change { dX = -beta(t)/2 * X * dt + sqrt(beta(t)) * dW }
  structure { time t in [0, T] ; weights probability ; calculus ito }
}

system generate = reverse forward
)";

static const char* further = R"(
system lorenz {
  state x, y, z : real
  param sigma = 10, rho = 28, beta = 8/3
  change {
    x' = sigma * (y - x)
    y' = x * (rho - z) - y
    z' = x * y - beta * z
  }
  structure { time t in real }
  constraint { invariant under reflect(x, y) }
  difference { observe x, r = sqrt(x^2 + y^2) }
}

system neuron {
  state v, w : real
  param a = 0.02, b = 0.2, c = -65, d = 8
  input I = 10
  change {
    v' = 0.04 * v^2 + 5 * v + 140 - w + I
    w' = a * (b * v - w)
    when v >= 30 : v -> c, w -> w + d
  }
  structure { time t in real }
}

system pair {
  part left, right : neuron with d = 6
  connect {
    left.I = 10 + 0.5 * (right.v - left.v)
    right.I = 9
  }
}

system learned {
  state q(t) : circle
  unknown V : function
  unknown k : number
  change { action S = integral(q'^2 / 2 - ?V(q) - ?k * q^2) }
  structure { time t in real ; histories stationary S }
  constraint { conserves q'^2 / 2 + ?V(q) }
  difference { up to shift(q) }
}

system qubit {
  state rho : density(2)
  param omega = 1, gamma = 0.1
  change {
    hamiltonian H = omega/2 * sz
    jump L = sqrt(gamma) * sm
  }
  structure { time t in real ; weights operators }
}
)";

int main() {
  calx::verification::Suite s;

  // ---- the grammar and the three tests ----
  Script script;
  try {
    script = read(three_tests);
    s.require("the grammar compiles and the three tests are read: 7 systems", script.systems.size() == 7, std::to_string(script.systems.size()) + " systems");
  } catch (const std::exception& e) {
    s.require("the grammar compiles and the three tests are read", false, e.what());
    return s.done();
  }
  {
    const System& g = *script.find("anbn");
    s.require("a grammar: its alphabet, its start and its two rules", g.states.size() == 1 && g.states[0].space.words && g.states[0].space.alphabet.size() == 3 &&
                                                                          g.states[0].start == std::vector<std::string>{"S"} && g.change.rules.size() == 2 &&
                                                                          g.change.rules[0].name == "grow" && g.change.rules[0].left == std::vector<std::string>{"S"} &&
                                                                          g.change.rules[0].right == std::vector<std::string>{"a", "S", "b"});
    s.require("a grammar: derivations with Boolean weights, observed at its terminal words", g.structure.derivations && g.structure.weights == "boolean" && g.difference.observations.size() == 1 && !g.structure.timed);
    const System& w = *script.find("anbn_weighted");
    s.require("a stochastic grammar: the weights 0.7 and 0.3 of its rules", w.change.rules[0].weight && kt::is_number(w.change.rules[0].weight, 0.7) && kt::is_number(w.change.rules[1].weight, 0.3) && w.structure.weights == "probability");
  }
  {
    const System& b = *script.find("boulatov");
    const auto& st = b.states.at(0);
    s.require("a group field theory: a field of three arguments ranging in a group, and no background", st.names[0].name == "phi" && st.names[0].arguments.size() == 3 && st.ranges.size() == 1 &&
                                                                                                        st.ranges[0].names.size() == 3 && st.ranges[0].domain.kind == Domain::Kind::group && st.ranges[0].domain.name == "G" &&
                                                                                                        b.constraint.backgrounds.empty() && !b.structure.timed);
    s.require("a group field theory: its invariance as an identity of three arguments", b.constraint.identities.size() == 1 && b.constraint.identities[0].arguments.size() == 3);
    std::set<std::pair<std::string, std::size_t>> apps;
    kt::applications(b.change.actions.at(0).value, apps);
    s.require("a group field theory: its action applies phi to three arguments and integrates, read across two lines", apps.count({"phi", 3}) && apps.count({"integral", 1}) && apps.size() == 2);
    s.require("a group field theory: histories weighted by exp(-S)", b.structure.weight && kt::text(b.structure.weight) == "exp(-S)" && b.structure.weights == "real");
    const System& t = *script.find("tetrahedra");
    s.require("its expansion is a derivation", t.derived() && t.derivation == "expand" && t.source == "boulatov");
    const System& p = *script.find("phi4");
    s.require("a field on a background declares it as a constraint", p.constraint.backgrounds.size() == 1 && p.constraint.backgrounds[0].names == std::vector<std::string>{"x"} &&
                                                                         p.constraint.backgrounds[0].domain.kind == Domain::Kind::applied && p.constraint.backgrounds[0].domain.name == "lattice" && p.constraint.backgrounds[0].domain.args.size() == 2);
  }
  {
    const System& f = *script.find("forward");
    s.require("a diffusion: a path in real^d, a schedule beta(t), two declared sources of noise", f.states[0].names[0].arguments == std::vector<std::string>{"t"} && f.states[0].space.exponent == "d" && f.params.size() == 3 &&
                                                                                                     f.params[2].name.arguments.size() == 1 && f.noises.size() == 2 && !f.noises[0].law && f.noises[1].law && f.noises[1].arguments.size() == 1);
    s.require("a diffusion: its stochastic differential, its time in [0, T], Ito's calculus", f.change.equations.size() == 1 && f.change.equations[0].name == "dX" && f.structure.timed && f.structure.time_domain.kind == Domain::Kind::interval &&
                                                                                                 f.structure.calculus == "ito" && f.structure.weights == "probability");
    s.require("a diffusion: the generating process is its reversal", script.find("generate")->derivation == "reverse");
  }

  // ---- precedence ----
  const auto value = [&](const std::string& text) {
    return kt::evaluate(read_term(text), [](const std::string& n, unsigned) { return n == "a" ? 7.0 : n == "b" ? 3.0 : 2.0; });
  };
  const auto valued = [&](const std::string& text, double want) { s.check("the value of " + text, value(text) - want, 1e-13); };
  valued("-2^2", -4);
  valued("2^3^2", 512);
  valued("2^-1", 0.5);
  valued("a - b - c", 2);
  valued("a / b / c", 7.0 / 6);
  valued("a / b * c", 14.0 / 3);
  valued("a - b * c^2", -5);
  valued("-(a - b)^2 / c", -8);
  valued("2 * -a", -14);
  valued("1.5e2 + 2.5E-1", 150.25);
  valued("a\n  + b\n  * c", 13);

  // ---- terms round the trip ----
  {
    std::mt19937_64 rng(7);
    const std::vector<TermPtr> leaves = {kt::symbol("x"), kt::symbol("y"), kt::symbol("q", 1), kt::symbol("a.v"), kt::number(2), kt::number(-3), kt::number(0.25), kt::number(1e-9), kt::hole("f", {kt::symbol("x")})};
    std::function<TermPtr(int)> random_term = [&](int depth) -> TermPtr {
      if (depth == 0 || rng() % 5 == 0) return leaves[rng() % leaves.size()];
      switch (rng() % 6) {
        case 0: return random_term(depth - 1) + random_term(depth - 1);
        case 1: return random_term(depth - 1) * random_term(depth - 1);
        case 2: return random_term(depth - 1) - random_term(depth - 1);
        case 3: return random_term(depth - 1) / random_term(depth - 1);
        case 4: return kt::power(random_term(depth - 1), rng() % 2 ? kt::number(static_cast<double>(rng() % 5) - 2) : random_term(depth - 1));
        default: return kt::call(rng() % 2 ? "sin" : "exp", {random_term(depth - 1)});
      }
    };
    int failed = 0;
    std::string first;
    for (int k = 0; k < 400; k++) {
      const TermPtr t = random_term(4);
      const std::string w = kt::text(t);
      try {
        const TermPtr back = read_term(w);
        if (!kt::same(back, t)) {
          if (!failed++) first = w + "  read back as  " + kt::text(back);
        }
      } catch (const std::exception& e) {
        if (!failed++) first = w + ": " + e.what();
      }
    }
    s.require("400 random terms, written and read, are the same terms", failed == 0, failed ? std::to_string(failed) + " differ; the first: " + first : "");
  }

  // ---- systems round the trip, and the further forms ----
  Script more;
  try {
    more = read(further);
    s.require("flows, events, parts, unknowns and a qubit are read: 5 systems", more.systems.size() == 5);
  } catch (const std::exception& e) {
    s.require("flows, events, parts, unknowns and a qubit are read", false, e.what());
    return s.done();
  }
  for (const Script* sc : {&script, &more}) {
    const std::string once = write(*sc);
    try {
      const std::string twice = write(read(once));
      s.require(std::string("a script of ") + std::to_string(sc->systems.size()) + " systems, written and read, writes the same", once == twice, once == twice ? "" : "the two texts differ");
    } catch (const std::exception& e) {
      s.require("a written script is read", false, std::string(e.what()) + "\n" + once);
    }
  }
  {
    const System& l = *more.find("lorenz");
    s.require("a flow: three equations of first order, a real time, an invariance, a named observation", l.change.equations.size() == 3 && l.change.equations[1].name == "y" && l.change.equations[1].primes == 1 && l.structure.timed &&
                                                                                                            l.structure.time_domain.kind == Domain::Kind::real && l.constraint.invariances.size() == 1 &&
                                                                                                            l.difference.observations.size() == 2 && l.difference.observations[1].name == "r");
    s.require("the equation of y, as a term: x (rho - z) - y", kt::same(l.change.equations[1].value, kt::symbol("x") * (kt::symbol("rho") - kt::symbol("z")) - kt::symbol("y")));
    const System& n = *more.find("neuron");
    s.require("an event: a guard v >= 30 and two resets; an input with its value when alone", n.change.events.size() == 1 && n.change.events[0].relation == ">=" && n.change.events[0].resets.size() == 2 &&
                                                                                                 n.change.events[0].resets[1].name == "w" && n.inputs.size() == 1 && kt::is_number(n.inputs[0].fallback, 10));
    const System& p = *more.find("pair");
    s.require("a system of parts: two parts of one system with a parameter set, two inputs given", p.parts.size() == 2 && p.parts[1].name == "right" && p.parts[1].system == "neuron" && p.parts[0].overrides.size() == 1 &&
                                                                                                      p.wires.size() == 2 && p.wires[0].part == "left" && p.wires[0].input == "I" &&
                                                                                                      kt::depends_on(p.wires[0].value, "right.v", 0));
    const System& u = *more.find("learned");
    std::set<std::string> hs;
    kt::holes(u.change.actions.at(0).value, hs);
    s.require("a template: two unknowns, a function and a number, inside an action; a conserved quantity; states up to a shift", u.unknowns.size() == 2 && u.unknowns[0].kind == Unknown::Kind::function &&
                                                                                                                                      u.unknowns[1].kind == Unknown::Kind::number && hs.size() == 2 && u.structure.stationary == "S" &&
                                                                                                                                      u.constraint.conserved.size() == 1 && u.difference.up_to.size() == 1);
    const System& q = *more.find("qubit");
    s.require("an open quantum system: a Hamiltonian and a jump, weights that are operators", q.change.hamiltonians.size() == 1 && q.change.jumps.size() == 1 && q.structure.weights == "operators" && q.states[0].space.name == "density");
  }

  {
    // A part's state in a term may carry its primes: a.q' is the velocity of the part a's q.
    const TermPtr v = read_term("a.q'^2 / 2 - k * (a.q - b.q)^2");
    bool speed = false, place = false;
    for (const auto& [name, primes] : kt::symbols(v)) speed = speed || (name == "a.q" && primes == 1), place = place || (name == "b.q" && primes == 0);
    s.require("a part's state is written a.q, and its velocity a.q'", speed && place, kt::text(v));
    s.require("and is written out as it was read", kt::same(read_term(kt::text(v)), v), kt::text(v));
  }

  // ---- refusals ----
  const auto refused = [&](const std::string& what, const std::string& text, std::size_t line, std::size_t column) {
    try {
      read(text);
      s.require(what, false, "read");
    } catch (const ScriptError& e) {
      s.require(what, e.at.line == line && e.at.column == column, e.what());
    }
  };
  refused("a product by juxtaposition is refused where it goes wrong", "system w {\n  state x, y : real\n  change { x' = sigma (y - x) }\n}\n", 3, 23);
  refused("a block that is not closed is refused at the end", "system w {\n  state x : real\n", 3, 1);
  refused("a system declared twice is refused at the second", "system a {\n}\nsystem a {\n}\n", 3, 1);
  refused("an equation without a right side is refused", "system w {\n  change { x' = }\n}\n", 2, 17);

  return s.done();
}
