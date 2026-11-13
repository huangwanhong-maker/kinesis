// A catalogue of models: laws that many fields have written down, each a declaration in the language of dynamical
// systems, with where it comes from.
//
// What it is for. A model of the catalogue is a starting point for a study (it is read, run and analysed as any
// declaration is), and a hypothesis for a record: with its parameters freed (induction's as_template) it is a
// template, and a record is matched against the catalogue by fitting each template whose shape it could have and
// ranking them (induction's match). candidates() gives those templates for a record of given components.
//
// What an entry is. Its text declares one system, named as the entry's id, with the values of its parameters that
// are in common use. Nothing in an entry names a place or a measurement: a model of the catalogue is a form of
// law, and a domain that owns the objects it describes (a neuron, an epidemic, a reaction) keeps their values.
//
// Shape and the order of components. A record's components have no names a catalogue knows, so a model's states
// are named anew after the record's, in each order that is tried: all orders for three components or fewer (six
// at most), the record's own order above that. A parameter of the model that shares a name with a component of
// the record is renamed first. *Caveat.* Only plain models are offered as candidates: flows and maps given by
// first-order equations, with no events, noise, inputs or actions. A ranking is among the candidates tried and
// says nothing of a law that is not in the catalogue.
//
// Sources are cited by author, year and journal; the values of parameters are those in common use, which are
// not always those of the first paper, and an entry says so where it matters.
#pragma once

#include <kinesis/catalogue/turbulence.hpp>
#include <kinesis/kernel/answer.hpp>
#include <kinesis/language/declaration.hpp>
#include <kinesis/language/reader.hpp>
#include <kinesis/terms/term.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace kinesis::catalogue {

using language::System;

struct Entry {
  std::string id;      // the name of the system its text declares
  std::string title;
  std::string field;   // where the model is at home: a word for a list, not a classification
  std::string about;   // what it is, and what is known of it in closed form
  std::string source;
  std::string text;    // a script of one system
};

inline const std::vector<Entry>& entries() {
  static const std::vector<Entry> all = {
      // ---- flows of three components ----
      {"lorenz", "Lorenz's system", "fluids",
       "Three modes of convection in a layer heated from below. Its divergence is -(sigma + 1 + beta) everywhere, so its exponents sum to it; at rho = 28 it is chaotic.",
       "Lorenz 1963, J. Atmos. Sci. 20",
       R"(system lorenz {
  state x, y, z : real
  param sigma = 10, rho = 28, beta = 8/3
  change {
    x' = sigma * (y - x)
    y' = x * (rho - z) - y
    z' = x * y - beta * z
  }
}
)"},
      {"rossler", "Roessler's system", "chemistry",
       "One nonlinear term, a chaotic band. Its equilibria are at x = (c +- sqrt(c^2 - 4 a b)) / 2, y = -x / a, z = x / a.",
       "Roessler 1976, Phys. Lett. A 57",
       R"(system rossler {
  state x, y, z : real
  param a = 0.2, b = 0.2, c = 5.7
  change {
    x' = -y - z
    y' = x + a * y
    z' = b + z * (x - c)
  }
}
)"},
      {"chen", "Chen's system", "circuits",
       "Lorenz's terms with another linear part. Its divergence is -(a + b - c) everywhere.",
       "Chen and Ueta 1999, Int. J. Bifurcation Chaos 9",
       R"(system chen {
  state x, y, z : real
  param a = 35, b = 3, c = 28
  change {
    x' = a * (y - x)
    y' = (c - a) * x - x * z + c * y
    z' = x * y - b * z
  }
}
)"},
      {"thomas", "Thomas's cyclically symmetric system", "feedback circuits",
       "Three variables, each damped and driven by the sine of the next; unchanged by the cyclic exchange of its variables. Its divergence is -3 b everywhere.",
       "Thomas 1999, Int. J. Bifurcation Chaos 9",
       R"(system thomas {
  state x, y, z : real
  param b = 0.208186
  change {
    x' = sin(y) - b * x
    y' = sin(z) - b * y
    z' = sin(x) - b * z
  }
  constraint { invariant under cycle(x, y, z) }
}
)"},
      {"sir", "The epidemic of Kermack and McKendrick", "epidemics",
       "Susceptible, infectious, removed, as shares of a closed population: s + i + r is kept. When no one was removed at the start, the share never infected solves ln(s / s0) = -(beta / gamma) (1 - s).",
       "Kermack and McKendrick 1927, Proc. R. Soc. A 115",
       R"(system sir {
  state s, i, r : real
  param beta = 0.5, gamma = 0.2
  change {
    s' = -beta * s * i
    i' = beta * s * i - gamma * i
    r' = gamma * i
  }
  constraint { conserves s + i + r }
}
)"},
      {"hindmarsh_rose", "Hindmarsh and Rose's neuron", "neurons",
       "A fast pair that spikes and a slow variable that ends each burst. The values are those in common use for bursting.",
       "Hindmarsh and Rose 1984, Proc. R. Soc. B 221",
       R"(system hindmarsh_rose {
  state x, y, z : real
  param a = 1, b = 3, c = 1, d = 5, r = 0.006, s = 4, x0 = -1.6, I = 3.25
  change {
    x' = y - a * x^3 + b * x^2 - z + I
    y' = c - d * x^2 - y
    z' = r * (s * (x - x0) - z)
  }
}
)"},
      {"sabra", "A shell model of turbulence", "fluids",
       "The Sabra model on fourteen shells: a cascade of energy from the forced shell to those where viscosity removes it. Without viscosity and force it keeps its energy; in a steady state the power put in, the power dissipated and the flux through every shell between are equal on average. It is not a discretisation of the equations of a fluid.",
       "L'vov, Podivilov, Pomyalov, Procaccia and Vandembroucq 1998, Phys. Rev. E 58",
       shell_model()},
      // ---- flows of two components ----
      {"van_der_pol", "Van der Pol's oscillator", "circuits",
       "A damping that is negative at small amplitudes and positive at large ones: one limit cycle. For small mu its amplitude is 2 and its period 2 pi (1 + mu^2 / 16).",
       "van der Pol 1926, Philos. Mag. 2",
       R"(system van_der_pol {
  state x, y : real
  param mu = 1
  change {
    x' = y
    y' = mu * (1 - x^2) * y - x
  }
}
)"},
      {"fitzhugh_nagumo", "FitzHugh and Nagumo's excitable unit", "neurons",
       "Van der Pol's oscillator made excitable: a fast potential and a slow recovery. With the current given it fires periodically.",
       "FitzHugh 1961, Biophys. J. 1; Nagumo, Arimoto and Yoshizawa 1962, Proc. IRE 50",
       R"(system fitzhugh_nagumo {
  state v, w : real
  param a = 0.7, b = 0.8, eps = 0.08, I = 0.5
  change {
    v' = v - v^3 / 3 - w + I
    w' = eps * (v + a - b * w)
  }
}
)"},
      {"lotka_volterra", "Lotka and Volterra's prey and predator", "populations",
       "Closed orbits about (c / d, a / b): d x - c ln x + b y - a ln y is kept along every path.",
       "Lotka 1925, Elements of Physical Biology; Volterra 1926, Nature 118",
       R"(system lotka_volterra {
  state x, y : real
  param a = 1, b = 0.5, c = 0.75, d = 0.25
  change {
    x' = a * x - b * x * y
    y' = d * x * y - c * y
  }
  constraint { conserves d * x - c * log(x) + b * y - a * log(y) }
}
)"},
      {"brusselator", "The Brusselator", "chemistry",
       "A scheme of four reactions with an autocatalytic step. Its equilibrium (a, b / a) loses its stability at b = 1 + a^2, where a limit cycle is born.",
       "Prigogine and Lefever 1968, J. Chem. Phys. 48",
       R"(system brusselator {
  state x, y : real
  param a = 1, b = 3
  change {
    x' = a - (b + 1) * x + x^2 * y
    y' = b * x - x^2 * y
  }
}
)"},
      {"selkov", "Sel'kov's glycolysis", "biochemistry",
       "A substrate and a product that activates its own making. Its equilibrium is (b, b / (a + b^2)).",
       "Sel'kov 1968, Eur. J. Biochem. 4",
       R"(system selkov {
  state x, y : real
  param a = 0.08, b = 0.6
  change {
    x' = -x + a * y + x^2 * y
    y' = b - a * y - x^2 * y
  }
}
)"},
      {"duffing", "Duffing's oscillator, damped and unforced", "mechanics",
       "A particle in the potential alpha x^2 / 2 + beta x^4 / 4 with friction: with alpha < 0 two wells at +-sqrt(-alpha / beta) and a saddle between them.",
       "Duffing 1918, Erzwungene Schwingungen bei veraenderlicher Eigenfrequenz",
       R"(system duffing {
  state x, y : real
  param alpha = -1, beta = 1, delta = 0.3
  change {
    x' = y
    y' = -delta * y - alpha * x - beta * x^3
  }
}
)"},
      {"stuart_landau", "The normal form of Hopf's bifurcation", "oscillations",
       "The form every system takes near the birth of a small oscillation. For mu > 0 its limit cycle is the circle of radius sqrt(mu), turned at the rate omega.",
       "Landau 1944, Dokl. Akad. Nauk SSSR 44; Stuart 1960, J. Fluid Mech. 9",
       R"(system stuart_landau {
  state x, y : real
  param mu = 0.5, omega = 2
  change {
    x' = mu * x - omega * y - (x^2 + y^2) * x
    y' = omega * x + mu * y - (x^2 + y^2) * y
  }
  constraint { invariant under rotate(x, y) }
}
)"},
      // ---- a flow of one component ----
      {"logistic_growth", "Verhulst's growth", "populations",
       "Growth that saturates: x(t) = K / (1 + (K / x0 - 1) exp(-r t)).",
       "Verhulst 1838, Corresp. Math. Phys. 10",
       R"(system logistic_growth {
  state x : real
  param r = 1, K = 10
  change { x' = r * x * (1 - x / K) }
}
)"},
      // ---- maps ----
      {"logistic_map", "The logistic map", "populations",
       "Generations that do not overlap. It doubles its period on the way to chaos; at r = 4 its exponent is ln 2.",
       "May 1976, Nature 261",
       R"(system logistic_map {
  state x : real
  param r = 3.9
  change { x' = r * x * (1 - x) }
  structure { time n in natural }
}
)"},
      {"henon", "Henon's map", "maps of the plane",
       "A fold and a contraction: the determinant of its derivative is -b everywhere, so its exponents sum to ln b.",
       "Henon 1976, Commun. Math. Phys. 50",
       R"(system henon {
  state x, y : real
  param a = 1.4, b = 0.3
  change {
    x' = 1 - a * x^2 + y
    y' = b * x
  }
  structure { time n in natural }
}
)"},
      {"standard_map", "Chirikov's standard map", "mechanics",
       "A rotor kicked once a period. It keeps area: the determinant of its derivative is 1 everywhere, so its exponents sum to zero.",
       "Chirikov 1979, Phys. Rep. 52",
       R"(system standard_map {
  state p, q : real
  param K = 0.97
  change {
    p' = p + K * sin(q)
    q' = q + p + K * sin(q)
  }
  structure { time n in natural }
}
)"},
      // ---- flows with events ----
      {"integrate_and_fire", "The leaky integrate-and-fire neuron", "neurons",
       "A leak toward rest, a constant drive, and a reset at the threshold. It fires with the period tau ln((R I + E - vr) / (R I + E - theta)) when the drive carries it past the threshold.",
       "Lapicque 1907, J. Physiol. Pathol. Gen. 9",
       R"(system integrate_and_fire {
  state v : real
  param tau = 10, E = -65, R = 10, I = 2, theta = -50, vr = -65
  change {
    v' = (E - v + R * I) / tau
    when v >= theta : v -> vr
  }
}
)"},
      {"izhikevich", "Izhikevich's neuron", "neurons",
       "A quadratic potential with a recovery variable and a reset: with other values of a, b, c and d it spikes regularly, bursts or chatters.",
       "Izhikevich 2003, IEEE Trans. Neural Netw. 14",
       R"(system izhikevich {
  state v, w : real
  param a = 0.02, b = 0.2, c = -65, d = 8, I = 10
  change {
    v' = 0.04 * v^2 + 5 * v + 140 - w + I
    w' = a * (b * v - w)
    when v >= 30 : v -> c, w -> w + d
  }
}
)"},
      // ---- systems given by an action ----
      {"pendulum", "The pendulum", "mechanics",
       "Given by its action. Its small swings have the period 2 pi sqrt(l / g), and its energy is kept.",
       "Lagrange 1788, Mecanique analytique",
       R"(system pendulum {
  state q(t) : circle
  param m = 1, l = 1, g = 9.81
  change { action S = integral(m * l^2 * q'^2 / 2 + m * g * l * cos(q)) }
  structure { time t in real ; histories stationary S }
}
)"},
      {"kepler", "Kepler's problem in the plane", "mechanics",
       "A body drawn to a centre as the inverse square, given by its action. Energy and angular momentum are kept; a circular orbit of radius a has the period 2 pi a^(3/2) / sqrt(k).",
       "Newton 1687, Principia",
       R"(system kepler {
  state x(t), y(t) : real
  param k = 1
  change { action S = integral((x'^2 + y'^2) / 2 + k / sqrt(x^2 + y^2)) }
  structure { time t in real ; histories stationary S }
}
)"},
      // ---- systems with noise ----
      {"ornstein_uhlenbeck", "Ornstein and Uhlenbeck's process", "noise",
       "A pull toward zero and a white noise. Its stationary law is normal with the variance sigma^2 / (2 theta).",
       "Uhlenbeck and Ornstein 1930, Phys. Rev. 36",
       R"(system ornstein_uhlenbeck {
  state X : real
  param theta = 1, sigma = 0.5
  noise W : wiener
  change { dX = -theta * X * dt + sigma * dW }
}
)"},
      {"double_well", "A particle in two wells with noise", "noise",
       "Kramers's problem: the particle stays in a well and at long intervals crosses the barrier, at a rate that falls as exp(-2 barrier / s^2).",
       "Kramers 1940, Physica 7",
       R"(system double_well {
  state X : real
  param s = 0.45
  noise W : wiener
  change { dX = (X - X^3) * dt + s * dW }
}
)"},
      {"geometric_brownian", "Geometric Brownian motion", "finance",
       "Growth and noise both in proportion to the state, read in Ito's calculus: the mean grows as exp(mu t), and the logarithm drifts at mu - sigma^2 / 2.",
       "Samuelson 1965, Ind. Manag. Rev. 6",
       R"(system geometric_brownian {
  state X : real
  param mu = 0.05, sigma = 0.2
  noise W : wiener
  change { dX = mu * X * dt + sigma * X * dW }
  structure { calculus ito }
}
)"},
      {"square_root_diffusion", "The square-root diffusion", "finance",
       "A pull toward a level with a noise that vanishes at zero. Its stationary law is a gamma law of shape 2 kappa level / sigma^2 and rate 2 kappa / sigma^2.",
       "Cox, Ingersoll and Ross 1985, Econometrica 53",
       R"(system square_root_diffusion {
  state X : real
  param kappa = 1, level = 1, sigma = 0.5
  noise W : wiener
  change { dX = kappa * (level - X) * dt + sigma * sqrt(X) * dW }
  structure { calculus ito }
}
)"},
      // ---- open quantum systems ----
      {"driven_atom", "A driven atom that decays", "quantum optics",
       "Two levels driven on resonance and decaying by emission. Its stationary excited population is Omega^2 / (gamma^2 + 2 Omega^2).",
       "Breuer and Petruccione 2002, The Theory of Open Quantum Systems",
       R"(system driven_atom {
  state rho : density(2)
  param Omega = 1.5, gamma = 0.5
  change {
    hamiltonian H = Omega / 2 * sx
    jump L = sqrt(gamma) * sm
  }
  difference { observe excited = (id + sz) / 2, coherence = sx }
}
)"},
      {"jaynes_cummings", "Jaynes and Cummings's model with a lossy cavity", "quantum optics",
       "Two levels and one mode of a cavity that exchange a quantum; the cavity loses it. The oscillator is cut at eight levels.",
       "Jaynes and Cummings 1963, Proc. IEEE 51",
       R"(system jaynes_cummings {
  state rho : density(2, 8)
  param g = 0.7, kappa = 0.08
  change {
    hamiltonian H = g * (a * sp + adag * sm)
    jump L = sqrt(kappa) * a
  }
  difference { observe excited = (id + sz) / 2, photons = num }
}
)"},
      // ---- rules ----
      {"anbn", "The words a^n b^n", "grammars",
       "The first language that no finite automaton accepts, as a system whose histories are derivations.",
       "Chomsky 1956, IRE Trans. Inf. Theory 2",
       R"(system anbn {
  state w : words over { S, a, b }, start S
  change {
    grow : S -> a S b
    stop : S -> a b
  }
  structure { derivations ; weights boolean }
}
)"},
      {"balanced", "Balanced brackets", "grammars",
       "Dyck's words on one pair of brackets, by a grammar in which each word has one history: there are as many of 2 n letters as the n-th Catalan number.",
       "Chomsky and Schuetzenberger 1963, in Computer Programming and Formal Systems",
       R"(system balanced {
  state w : words over { S, open, close }, start S
  change {
    nest : S -> open S close
    pair : S -> open close
    nest_then : S -> open S close S
    pair_then : S -> open close S
  }
  structure { derivations ; weights boolean }
}
)"},
  };
  return all;
}

inline const Entry* find(const std::string& id) {
  for (const auto& e : entries())
    if (e.id == id) return &e;
  return nullptr;
}

// The system an entry declares.
inline System system(const Entry& e) {
  const language::Script script = language::read(e.text);
  const System* s = script.find(e.id);
  if (!s) throw kernel::Refusal("the catalogue's entry '" + e.id + "' does not declare a system of its name");
  return *s;
}
inline System system(const std::string& id) {
  const Entry* e = find(id);
  if (!e) throw kernel::Refusal("the catalogue holds no model named '" + id + "'");
  return system(*e);
}

// The names of a model's states, in the order declared.
inline std::vector<std::string> states_of(const System& s) {
  std::vector<std::string> out;
  for (const auto& group : s.states)
    for (const auto& d : group.names) out.push_back(d.name);
  return out;
}

// Whether a model is plain: a flow or a map given by first-order equations alone, which can be named anew.
inline bool plain(const System& s) {
  if (s.derived() || !s.noises.empty() || !s.inputs.empty() || !s.unknowns.empty() || !s.parts.empty() || !s.wires.empty()) return false;
  if (s.change.equations.empty() || !s.change.actions.empty() || !s.change.hamiltonians.empty() || !s.change.jumps.empty() || !s.change.rules.empty() || !s.change.events.empty()) return false;
  for (const auto& group : s.states) {
    if (group.space.words) return false;
    for (const auto& d : group.names)
      if (!d.arguments.empty()) return false;
  }
  for (const auto& e : s.change.equations)
    if (e.primes != 1 || !e.arguments.empty()) return false;
  return true;
}
inline bool discrete(const System& s) {
  return s.structure.timed && (s.structure.time_domain.kind == language::Domain::Kind::natural || s.structure.time_domain.kind == language::Domain::Kind::integer);
}

// A plain model with its states named anew: names[i] stands for its i-th state. A parameter that shares a name
// with a new state is given another first. What the model declared of its symmetries and of what it keeps is
// dropped, since a template with its numbers freed need not keep them.
inline System renamed(const System& model, const std::vector<std::string>& names) {
  const std::vector<std::string> old = states_of(model);
  if (!plain(model)) throw kernel::Refusal(model.name + " is not a plain flow or map: its states are not named anew");
  if (old.size() != names.size()) throw kernel::Refusal(model.name + " has " + std::to_string(old.size()) + " components, and " + std::to_string(names.size()) + " names were given");
  std::map<std::string, std::string> to;
  for (std::size_t i = 0; i < old.size(); i++) to[old[i]] = names[i];
  const auto taken = [&](const std::string& n) {
    if (std::find(names.begin(), names.end(), n) != names.end()) return true;
    for (const auto& p : model.params)
      if (p.name.name == n) return true;
    return false;
  };
  for (const auto& p : model.params) {
    if (std::find(names.begin(), names.end(), p.name.name) == names.end() || to.count(p.name.name)) continue;
    std::string fresh = p.name.name + "_";
    while (taken(fresh)) fresh += "_";
    to[p.name.name] = fresh;
  }
  const auto move = [&](const terms::TermPtr& t) {
    return terms::rewrite(t, [&](const terms::TermPtr& x) -> terms::TermPtr {
      if (x->kind != terms::Kind::symbol) return nullptr;
      const auto it = to.find(x->name);
      return it == to.end() ? nullptr : terms::symbol(it->second, x->primes);
    });
  };
  const auto named = [&](const std::string& n) {
    const auto it = to.find(n);
    return it == to.end() ? n : it->second;
  };
  System out = model;
  for (auto& group : out.states)
    for (auto& d : group.names) d.name = named(d.name);
  for (auto& p : out.params) p.name.name = named(p.name.name), p.value = move(p.value);
  for (auto& e : out.change.equations) e.name = named(e.name), e.value = move(e.value);
  out.constraint = {};
  out.difference = {};
  return out;
}

// The templates a record of these components could be a path of: each plain model of the catalogue with as many
// components and the same kind of time, its states named after the record's in every order (three components or
// fewer) or in the record's own, each still under the model's name. Their parameters are not yet freed.
inline std::vector<System> candidates(const std::vector<std::string>& names, bool maps) {
  std::vector<System> out;
  std::vector<std::size_t> order(names.size());
  for (const auto& e : entries()) {
    const System model = system(e);
    if (!plain(model) || discrete(model) != maps || states_of(model).size() != names.size()) continue;
    std::iota(order.begin(), order.end(), std::size_t{0});
    do {
      std::vector<std::string> as;
      for (std::size_t i : order) as.push_back(names[i]);
      out.push_back(renamed(model, as));
    } while (names.size() <= 3 && std::next_permutation(order.begin(), order.end()));
  }
  return out;
}

}  // namespace kinesis::catalogue
