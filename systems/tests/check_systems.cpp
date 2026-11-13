// Checks of kinesis/systems: declarations made into systems that run, in each form.
//
// Flows and maps
// - Lorenz's field and its Jacobian as terms: 8 entries that are not zero, equal to central differences.
// - An equation of the second order becomes two components; the harmonic oscillator's error falls by 16 when the
//   step is halved (fourth order), and its energy is kept to that order.
// - Parameters with arguments are inlined; pi is a number unless declared.
// - The logistic map at r = 4 follows x_n = sin²(2ⁿθ) (Ulam and von Neumann's conjugacy).
// Events
// - A ball dropped from h₀ strikes at √(2h₀/g) and then at intervals 2eᵏv/g, to the tolerance of the search.
// - A leaky integrator that fires at 1 has the period τ ln(I/(I − 1)).
// Stochastic systems
// - The coefficients of a stochastic differential are read off its right side; one that is not linear in the
//   differentials is refused.
// - A Stratonovich equation dX = σX∘dW becomes Itô's with the drift σ²X/2, as a term.
// - On geometric Brownian motion, against the exact solution on the same Brownian path: the strong error of
//   Euler–Maruyama falls as h^½ and Milstein's as h (Kloeden and Platen 1992, §10.2–10.3).
// - The Ornstein–Uhlenbeck process keeps the variance σ²/(2θ); a path depends on (seed, path) alone.
// - A map with a draw at each step: an autoregression's variance 1/(1 − a²).
// Actions
// - The pendulum's action gives q'' = −(g/l) sin q as a term; the derived system, written out, reads back and
//   runs; its energy is kept. The double pendulum's two equations keep its energy. A central force keeps the
//   angular momentum. A degenerate action is refused.
// Open quantum systems
// - A spin and an oscillator at resonance (Jaynes and Cummings 1963): from the excited spin and no photon the
//   populations are cos²(gt) and sin²(gt). A driven, decaying atom reaches Ω²/(γ² + 2Ω²). An ordered product
//   ordered(adag, adag, a, a) is num² − num. A product of two operators on one factor is refused.
// Rules
// - aⁿbⁿ: the words within a bound, one history each; the bound is said. With weights, P(aⁿbⁿ) = 0.7ⁿ⁻¹·0.3.
// - An ambiguous grammar: x+x+x has 2 histories and x+x+x+x has 5 (Catalan's numbers).
// - aⁿbⁿcⁿ by rules with context: abc and aabbcc within length 6.
// Refusals with their place: a state with no equation, a name that is not declared, a state with two equations.

#include <kinesis/language/reader.hpp>
#include <kinesis/systems/lagrangian.hpp>
#include <kinesis/systems/quantum.hpp>
#include <kinesis/systems/rewriting.hpp>
#include <kinesis/systems/simulate.hpp>
#include <kinesis/systems/stochastic.hpp>

#include <calx/verification/check.hpp>
#include <physis/dynamics/flows.hpp>

#include <cmath>
#include <numbers>
#include <string>

using namespace kinesis;
using systems::State;
namespace kt = kinesis::terms;

static const char* text = R"(
system lorenz {
  state x, y, z : real
  param sigma = 10, rho = 28, beta = 8/3
  change {
    x' = sigma * (y - x)
    y' = x * (rho - z) - y
    z' = x * y - beta * z
  }
  structure { time t in real }
  difference { observe x, r = sqrt(x^2 + y^2) }
}

system oscillator {
  state q : real
  param omega = 2
  change { q'' = -omega^2 * q }
  constraint { conserves q'^2 / 2 + omega^2 * q^2 / 2 }
}

system driven {
  state x : real
  param a = 0.5, f(s) = a * sin(2 * pi * s), g(u, v) = u * f(v)
  change { x' = g(x, t) }
}

system logistic {
  state x : real
  param r = 4
  change { x' = r * x * (1 - x) }
  structure { time n in natural }
}

system ball {
  state h : real
  param g = 9.81, e = 0.8
  change {
    h'' = -g
    when h <= 0 : h' -> -e * h'
  }
}

system integrator {
  state v : real
  param tau = 0.02, I = 1.5
  change {
    v' = (I - v) / tau
    when v >= 1 : v -> 0
  }
}

system gbm {
  state X : real
  param mu = 0.5, sigma = 0.8
  noise W : wiener
  change { dX = mu * X * dt + sigma * X * dW }
  structure { time t in real ; calculus ito }
}

system gbm_stratonovich {
  state X : real
  param sigma = 0.8
  noise W : wiener
  change { dX = sigma * X * dW }
  structure { time t in real ; calculus stratonovich }
}

system ou {
  state X : real^2
  param theta = 1.5, sigma = 0.6
  noise W : wiener(2)
  noise X(0) ~ normal(0, 1)
  change { dX = -theta * X * dt + sigma * dW }
}

system autoregression {
  state x : real
  param a = 0.6
  noise xi : normal(0, 1)
  change { x' = a * x + xi }
  structure { time n in natural }
}

system pendulum {
  state q(t) : circle
  param m = 2, l = 0.5, g = 9.81
  change { action S = integral(m * l^2 * q'^2 / 2 + m * g * l * cos(q)) }
  structure { time t in real ; histories stationary S }
}

system double_pendulum {
  state a(t), b(t) : circle
  param g = 9.81
  change { action S = integral(a'^2 + b'^2 / 2 + a' * b' * cos(a - b) + 2 * g * cos(a) + g * cos(b)) }
  structure { time t in real ; histories stationary S }
}

system central {
  state r(t), phi(t) : real
  change { action S = integral((r'^2 + r^2 * phi'^2) / 2 + 1 / r) }
  structure { time t in real ; histories stationary S }
  constraint { conserves r^2 * phi' }
}

system degenerate {
  state q(t) : real
  change { action S = integral(q' * q - q^2) }
  structure { time t in real ; histories stationary S }
}

system jaynes_cummings {
  state rho : density(2, 8)
  param g = 0.7
  change { hamiltonian H = g * (a * sp + adag * sm) }
  difference { observe excited = (id + sz) / 2, photons = num }
}

system driven_atom {
  state rho : density(2)
  param Omega = 0.5, gamma = 1
  change {
    hamiltonian H = Omega / 2 * sx
    jump L = sqrt(gamma) * sm
  }
  difference { observe excited = (id + sz) / 2 }
}

system kerr {
  state rho : density(10)
  param U = 1.5
  change { hamiltonian H = U / 2 * ordered(adag, adag, a, a) }
}

system kerr_by_numbers {
  state rho : density(10)
  param U = 1.5
  change { hamiltonian H = U / 2 * (num^2 - num) }
}

system anbn {
  state w : words over { S, a, b }, start S
  change {
    grow : S -> a S b @ 0.7
    stop : S -> a b   @ 0.3
  }
  structure { derivations ; weights probability }
}

system sums {
  state w : words over { E, plus, x }, start E
  change {
    add : E -> E plus E
    one : E -> x
  }
  structure { derivations ; weights boolean }
}

system anbncn {
  state w : words over { S, B, C, a, b, c }, start S
  change {
    more : S -> a S B C
    last : S -> a B C
    swap : C B -> B C
    ab : a B -> a b
    bb : b B -> b b
    bc : b C -> b c
    cc : c C -> c c
  }
  structure { derivations ; weights boolean }
}
)";

int main() {
  calx::verification::Suite s;
  language::Script script;
  try {
    script = language::read(text);
  } catch (const std::exception& e) {
    s.require("the script of the checks is read", false, e.what());
    return s.done();
  }
  const auto sys = [&](const char* n) -> const language::System& { return *script.find(n); };

  // ---- flows ----
  {
    const systems::Dynamics d = systems::dynamics(sys("lorenz"));
    const State x{1.5, -2.0, 20.0};
    State dx(3);
    d.field(0, x, dx);
    s.check("Lorenz's field at (1.5, -2, 20)", std::fabs(dx[0] + 35) + std::fabs(dx[1] - 14) + std::fabs(dx[2] + 3 + 160.0 / 3), 1e-13);
    calx::linalg::Matrix J(3);
    d.jacobian(0, x, J);
    const auto f = [&](double t, const State& y, State& dy) { d.field(t, y, dy); };
    const calx::linalg::Matrix N = physis::dynamics::jacobian(f, 0.0, x);
    double worst = 0;
    for (std::size_t i = 0; i < 9; i++) worst = std::max(worst, std::fabs(J.a[i] - N.a[i]));
    s.check("its Jacobian from the terms against central differences", worst, 1e-7);
    s.require("its Jacobian has 8 entries that are not zero, and dz/dz is -beta", d.jacobian_terms().size() == 8 && kt::text(d.jacobian_terms().back().term) == "-beta", kt::text(d.jacobian_terms().back().term));
    State v{1, 2, 3}, jv(3);
    d.jacobian_vector(0, x, v, jv);
    s.check("J v without forming J", std::fabs(jv[0] - (J(0, 0) + 2 * J(0, 1))) + std::fabs(jv[2] - (J(2, 0) + 2 * J(2, 1) + 3 * J(2, 2))), 1e-13);
    s.check("its observation r = sqrt(x^2 + y^2)", d.observe(1, 0, x) - 2.5, 1e-15);
    systems::Dynamics e = d;
    e.set_param("rho", 0.5);
    e.field(0, x, dx);
    s.check("a parameter set after building changes the field, not the terms", dx[1] - (1.5 * (0.5 - 20) + 2), 1e-13);
  }
  {
    const systems::Dynamics d = systems::dynamics(sys("oscillator"));
    s.require("q'' = … makes the components q and q'", d.states == std::vector<std::string>{"q", "q'"});
    const auto error = [&](double h) {
      systems::Options o;
      o.t1 = 10, o.step = h, o.every = 1000000;
      const auto tr = systems::simulate(d, {1, 0}, o);
      return std::fabs(tr.x.back()[0] - std::cos(2 * 10.0));
    };
    const double e1 = error(0.02), e2 = error(0.01);
    s.check("the harmonic oscillator at t = 10, step 0.01", e2, 1e-7);
    s.check("the order of the method from two steps", calx::verification::order(e1, e2) - 4, 0.1);
    systems::Options o;
    o.t1 = 10, o.step = 0.01;
    const auto tr = systems::simulate(d, {1, 0}, o);
    double drift = 0;
    for (std::size_t k = 0; k < tr.x.size(); k++) drift = std::max(drift, std::fabs(d.conserved(0, tr.t[k], tr.x[k]) - 2.0));
    s.check("its energy over the run", drift, 1e-7);
  }
  {
    const systems::Dynamics d = systems::dynamics(sys("driven"));
    State dx(1);
    d.field(0.25, {3.0}, dx);
    s.check("functions of the parameters are inlined: g(x, t) = x a sin(2 pi t) at t = 1/4", dx[0] - 1.5, 1e-14);
    s.require("the field is a term of the state, the parameters and the time", kt::text(d.field_terms[0]).find("sin(") != std::string::npos && kt::text(d.field_terms[0]).find("g(") == std::string::npos, kt::text(d.field_terms[0]));
  }
  {
    const systems::Dynamics d = systems::dynamics(sys("logistic"));
    const double theta = 0.3;
    systems::Options o;
    o.t1 = 12;
    const auto tr = systems::simulate(d, {std::pow(std::sin(theta), 2)}, o);
    double worst = 0;
    for (std::size_t k = 0; k <= 12; k++) worst = std::max(worst, std::fabs(tr.x[k][0] - std::pow(std::sin(std::ldexp(theta, static_cast<int>(k))), 2)));
    s.require("a map in discrete time is iterated: 13 states", d.discrete && tr.x.size() == 13);
    s.check("the logistic map at r = 4 against sin^2(2^n theta), 12 iterations", worst, 1e-9);
  }

  // ---- events ----
  {
    const systems::Dynamics d = systems::dynamics(sys("ball"));
    systems::Options o;
    o.t1 = 3, o.step = 0.001;
    const auto tr = systems::simulate(d, {1.0, 0.0}, o);
    const double g = 9.81, e = 0.8, t1 = std::sqrt(2 / g), v = g * t1;
    double worst = 0, t = t1, speed = v;
    std::size_t counted = 0;
    for (const auto& ev : tr.events) {
      if (counted == 4) break;
      worst = std::max({worst, std::fabs(ev.t - t), std::fabs(ev.before[1] + speed), std::fabs(ev.after[1] - e * speed)});
      speed *= e, t += 2 * speed / g, counted++;
    }
    s.require("a dropped ball strikes at least four times in 3 s", tr.events.size() >= 4, std::to_string(tr.events.size()) + " events");
    s.check("the times of four impacts and the velocities before and after them", worst, 1e-8);
  }
  {
    const systems::Dynamics d = systems::dynamics(sys("integrator"));
    systems::Options o;
    o.t1 = 0.2, o.step = 0.0005;
    const auto tr = systems::simulate(d, {0.0}, o);
    const double T = 0.02 * std::log(1.5 / 0.5);
    double worst = 0;
    for (std::size_t k = 0; k < tr.events.size(); k++) worst = std::max(worst, std::fabs(tr.events[k].t - static_cast<double>(k + 1) * T));
    s.require("an integrator that fires at 1 fires 9 times in 0.2 s", tr.events.size() == 9, std::to_string(tr.events.size()));
    s.check("its firing times against tau ln(I / (I - 1))", worst, 1e-8);
  }

  // ---- stochastic systems ----
  {
    const systems::Stochastic d = systems::stochastic(sys("gbm"));
    s.require("the coefficients of dX = mu X dt + sigma X dW are read off", kt::text(d.drift_terms[0]) == "X * mu" && kt::text(d.diffusion_terms[0][0]) == "X * sigma" && d.diagonal,
              kt::text(d.drift_terms[0]) + "; " + kt::text(d.diffusion_terms[0][0]));
    // Strong errors at T = 1 on shared Brownian paths, from 2^10 fine increments.
    const std::size_t fine = 1024, paths = 400;
    const double mu = 0.5, sigma = 0.8;
    std::vector<double> em, mil;
    for (std::size_t level : {32u, 64u, 128u, 256u}) {
      double e_em = 0, e_mil = 0;
      for (std::size_t p = 0; p < paths; p++) {
        calx::random::Stream r(11, 1, static_cast<std::uint32_t>(p), 0);
        std::vector<double> dw(fine);
        double W = 0;
        for (auto& x : dw) x = std::sqrt(1.0 / static_cast<double>(fine)) * r.normal(), W += x;
        const double exact = std::exp((mu - sigma * sigma / 2) + sigma * W);
        std::vector<std::vector<double>> inc(level, std::vector<double>(1, 0.0));
        for (std::size_t k = 0; k < fine; k++) inc[k * level / fine][0] += dw[k];
        systems::StochasticOptions o;
        o.t1 = 1, o.step = 1.0 / static_cast<double>(level), o.every = level;
        e_em += std::fabs(systems::simulate(d, {1.0}, o, &inc).x.back()[0] - exact);
        o.milstein = true;
        e_mil += std::fabs(systems::simulate(d, {1.0}, o, &inc).x.back()[0] - exact);
      }
      em.push_back(e_em / static_cast<double>(paths)), mil.push_back(e_mil / static_cast<double>(paths));
    }
    const double order_em = std::log(em.front() / em.back()) / std::log(8.0), order_mil = std::log(mil.front() / mil.back()) / std::log(8.0);
    s.check("the strong order of Euler-Maruyama on geometric Brownian motion", order_em - 0.5, 0.12);
    s.check("the strong order of Milstein's scheme", order_mil - 1.0, 0.12);
    s.check("Milstein's strong error at h = 1/256", mil.back(), 5e-3);
  }
  {
    const systems::Stochastic d = systems::stochastic(sys("gbm_stratonovich"));
    const kt::TermPtr want = kt::number(0.5) * kt::power(kt::symbol("sigma"), kt::number(2)) * kt::symbol("X");
    s.require("Stratonovich's dX = sigma X dW is Ito's with the drift sigma^2 X / 2", kt::same(d.drift_terms[0], want) && d.calculus == "stratonovich", kt::text(d.drift_terms[0]));
  }
  {
    const systems::Stochastic d = systems::stochastic(sys("ou"));
    s.require("a state in real^2 driven by wiener(2): two components, two sources, diagonal noise, an initial law", d.states == std::vector<std::string>{"X1", "X2"} && d.sources.size() == 2 && d.diagonal &&
                                                                                                                       d.initial.size() == 1 && d.initial[0].kind == "normal");
    std::vector<double> B;
    d.diffusion(0, {0.3, -0.2}, B);
    s.check("its diffusion matrix is sigma times the identity", std::fabs(B[0] - 0.6) + std::fabs(B[1]) + std::fabs(B[2]) + std::fabs(B[3] - 0.6), 1e-15);
    systems::StochasticOptions o;
    o.t1 = 4000, o.step = 0.01, o.every = 10, o.seed = 5;
    const auto tr = systems::simulate(d, {0, 0}, o);
    double v1 = 0, v2 = 0, c = 0;
    for (const auto& x : tr.x) v1 += x[0] * x[0], v2 += x[1] * x[1], c += x[0] * x[1];
    const double n = static_cast<double>(tr.x.size()), want = 0.36 / 3.0;
    s.check("the Ornstein-Uhlenbeck variance sigma^2 / (2 theta), first component", v1 / n / want - 1, 0.05);
    s.check("the same, second component", v2 / n / want - 1, 0.05);
    s.check("the two components are not correlated", c / n / want, 0.05);
    systems::StochasticOptions a = o, b = o;
    a.t1 = b.t1 = 1, b.path = 1;
    const auto p0 = systems::simulate(d, {0, 0}, a), p0again = systems::simulate(d, {0, 0}, a), p1 = systems::simulate(d, {0, 0}, b);
    s.require("a path depends on (seed, path) alone", p0.x.back() == p0again.x.back() && p0.x.back() != p1.x.back());
  }
  {
    const systems::Stochastic d = systems::stochastic(sys("autoregression"));
    systems::StochasticOptions o;
    o.t1 = 200000, o.seed = 3;
    const auto tr = systems::simulate(d, {0}, o);
    double v = 0;
    for (const auto& x : tr.x) v += x[0] * x[0];
    s.check("a map with a draw at each step: the variance 1 / (1 - a^2) of an autoregression", v / static_cast<double>(tr.x.size()) * (1 - 0.36) - 1, 0.02);
  }
  try {
    systems::stochastic(*language::read("system w {\n  state X : real\n  noise W : wiener\n  change { dX = X * dt * dW }\n}\n").find("w"));
    s.require("a right side that is not linear in the differentials is refused", false);
  } catch (const language::ScriptError& e) {
    s.require("a right side that is not linear in the differentials is refused", e.at.line == 4, e.what());
  }

  // ---- actions ----
  {
    const systems::Lagrangian L = systems::vary(sys("pendulum"));
    const kt::TermPtr want = -(kt::symbol("g") / kt::symbol("l")) * kt::call("sin", {kt::symbol("q")});
    s.require("the pendulum's action gives q'' = -(g / l) sin q", L.acceleration.size() == 1 && kt::same(L.acceleration[0], want), kt::text(L.acceleration[0]));
    s.require("its momentum is m l^2 q' and its Lagrangian does not hold the time", kt::text(L.momenta[0]) == "m * q' * l^2" && L.autonomous, kt::text(L.momenta[0]));
    const language::System derived = systems::varied(sys("pendulum"));
    const std::string written = language::write(derived);
    const language::Script back = language::read(written);
    const systems::Dynamics d = systems::dynamics(*back.find("pendulum_paths"));
    systems::Options o;
    o.t1 = 20, o.step = 0.002;
    const auto tr = systems::simulate(d, {2.0, 0.0}, o);
    double drift = 0;
    for (std::size_t k = 0; k < tr.x.size(); k++) drift = std::max(drift, std::fabs(d.conserved(0, tr.t[k], tr.x[k]) - d.conserved(0, 0, tr.x[0])));
    s.require("the derived system, written out, is a flow of the language with its energy declared", d.states == std::vector<std::string>{"q", "q'"} && d.conserved_terms.size() == 1, written);
    s.check("its energy over 20 s from q = 2", drift, 1e-8);
  }
  {
    const language::System derived = systems::varied(sys("double_pendulum"));
    const systems::Dynamics d = systems::dynamics(derived);
    systems::Options o;
    o.t1 = 10, o.step = 0.0005;
    const auto tr = systems::simulate(d, {2.0, 0.0, 2.5, 0.0}, o);
    double drift = 0;
    for (std::size_t k = 0; k < tr.x.size(); k++) drift = std::max(drift, std::fabs(d.conserved(0, tr.t[k], tr.x[k]) - d.conserved(0, 0, tr.x[0])));
    s.require("the double pendulum: two coupled equations solved as terms, four components", d.dimension() == 4);
    s.check("its energy over 10 s of a large swing", drift, 1e-6);
  }
  {
    const systems::Dynamics d = systems::dynamics(systems::varied(sys("central")));
    systems::Options o;
    o.t1 = 20, o.step = 0.001;
    const auto tr = systems::simulate(d, {1.0, 0.0, 0.0, 1.2}, o);
    double drift_l = 0, drift_e = 0;
    for (std::size_t k = 0; k < tr.x.size(); k++) {
      drift_l = std::max(drift_l, std::fabs(d.conserved(0, tr.t[k], tr.x[k]) - 1.2));
      drift_e = std::max(drift_e, std::fabs(d.conserved(1, tr.t[k], tr.x[k]) - d.conserved(1, 0, tr.x[0])));
    }
    s.check("a central force keeps the angular momentum r^2 phi'", drift_l, 1e-9);
    s.check("and the energy", drift_e, 1e-9);
  }
  try {
    systems::vary(sys("degenerate"));
    s.require("an action that does not determine the accelerations is refused", false);
  } catch (const language::ScriptError& e) {
    s.require("an action that does not determine the accelerations is refused", std::string(e.what()).find("degenerate") != std::string::npos, e.what());
  }

  // ---- open quantum systems ----
  {
    const systems::Quantum q = systems::quantum(sys("jaynes_cummings"));
    const auto tr = systems::evolve(q, q.basis({0, 0}), 6, 0.05);
    double worst = 0;
    for (std::size_t k = 0; k < tr.t.size(); k++) worst = std::max({worst, std::fabs(tr.expectations[k][0] - std::pow(std::cos(0.7 * tr.t[k]), 2)), std::fabs(tr.expectations[k][1] - std::pow(std::sin(0.7 * tr.t[k]), 2))});
    s.require("a spin and an oscillator: a state of dimension 16 with two observations", q.dimension == 16 && q.factors == std::vector<std::size_t>{2, 8} && q.observables == std::vector<std::string>{"excited", "photons"});
    s.check("the vacuum Rabi oscillation: the excited population cos^2(g t), the photons sin^2(g t)", worst, 1e-12);
    s.check("a closed system stays pure", tr.purity.back() - 1, 1e-12);
    const systems::Quantum atom = systems::quantum(sys("driven_atom"));
    const auto st = physis::quantum::stationary(atom.generator);
    s.check("a driven, decaying atom: the stationary excited population Omega^2 / (gamma^2 + 2 Omega^2)", atom.expectations(st.state)[0] - 0.25 / 1.5, 1e-12);
    const auto relax = systems::evolve(atom, atom.basis({1}), 40, 0.1);
    s.check("the state propagated for a long time reaches it, and is mixed", std::fabs(relax.expectations.back()[0] - 0.25 / 1.5) + std::max(0.0, relax.purity.back() - 0.99), 1e-9);
    const systems::Quantum k1 = systems::quantum(sys("kerr")), k2 = systems::quantum(sys("kerr_by_numbers"));
    double gap = 0;
    for (std::size_t k = 0; k < k1.generator.hamiltonian.a.size(); k++) gap = std::max(gap, std::abs(k1.generator.hamiltonian.a[k] - k2.generator.hamiltonian.a[k]));
    s.check("an ordered product: ordered(adag, adag, a, a) is num^2 - num", gap, 1e-12);
    s.require("with no observation declared the populations are observed", k1.observables.size() == 10 && k1.observables[3] == "p3");
  }
  {
    const auto refused_quantum = [&](const std::string& what, const std::string& source, const std::string& names) {
      try {
        systems::quantum(*language::read(source).find("w"));
        s.require(what, false, "built");
      } catch (const language::ScriptError& e) {
        s.require(what, std::string(e.what()).find(names) != std::string::npos, e.what());
      }
    };
    refused_quantum("a product of two operators on one factor is refused, and the remedy named", "system w {\n  state rho : density(2)\n  change { hamiltonian H = sx * sy }\n}\n", "ordered");
    refused_quantum("a Hamiltonian that is not Hermitian is refused", "system w {\n  state rho : density(2)\n  change { hamiltonian H = sm }\n}\n", "not Hermitian");
    refused_quantum("an operator of an oscillator on a state that has none is refused", "system w {\n  state rho : density(2)\n  change { hamiltonian H = num }\n}\n", "oscillator");
  }

  // ---- rules ----
  {
    const systems::Rewriting r = systems::rewriting(sys("anbn"));
    systems::Bounds b;
    b.max_length = 8;
    const auto h = systems::explore(r, b);
    const auto words = systems::language(r, h);
    bool right = words.size() == 4;
    double worst = 0;
    for (std::size_t k = 0; k < words.size() && right; k++) {
      right = words[k].word.size() == 2 * (k + 1) && words[k].histories == 1;
      worst = std::max(worst, std::fabs(words[k].probability - std::pow(0.7, static_cast<double>(k)) * 0.3));
    }
    s.require("a^n b^n: four words within eight letters, one history each; the bound is said", right && !h.complete, h.cut);
    s.check("their probabilities 0.7^(n-1) 0.3", worst, 1e-15);
    s.require("membership by Earley's parser: aabb is a word, aab is not", r.accepts({"a", "a", "b", "b"}) && !r.accepts({"a", "a", "b"}) && r.context_free && r.monotone);
    calx::random::Stream stream(1, 2, 3, 4);
    bool all = true;
    double mean = 0;
    for (int k = 0; k < 2000; k++) {
      const auto history = systems::sample(r, stream);
      all = all && r.terminal(history.back()) && r.accepts(history.back());
      mean += static_cast<double>(history.back().size()) / 2;
    }
    s.require("2000 sampled histories each end in a word of the language", all);
    s.check("the mean n of a sampled a^n b^n against 1 / 0.3", mean / 2000 * 0.3 - 1, 0.06);
  }
  {
    const systems::Rewriting r = systems::rewriting(sys("sums"));
    systems::Bounds b;
    b.max_length = 7;
    const auto words = systems::language(r, systems::explore(r, b));
    const std::vector<double> catalan = {1, 1, 2, 5};
    bool right = words.size() == 4;
    for (std::size_t k = 0; k < words.size() && right; k++) right = words[k].word.size() == 2 * k + 1 && words[k].histories == catalan[k];
    s.require("an ambiguous grammar: x, x+x, x+x+x, x+x+x+x have 1, 1, 2, 5 histories", right, words.size() == 4 ? std::to_string(words[2].histories) + ", " + std::to_string(words[3].histories) : std::to_string(words.size()) + " words");
  }
  {
    const systems::Rewriting r = systems::rewriting(sys("anbncn"));
    systems::Bounds b;
    b.max_length = 6;
    const auto words = systems::language(r, systems::explore(r, b));
    s.require("rules with context: a^n b^n c^n within six letters is abc and aabbcc", !r.context_free && words.size() == 2 && words[0].word == systems::Word{"a", "b", "c"} &&
                                                                                           words[1].word == systems::Word{"a", "a", "b", "b", "c", "c"},
              std::to_string(words.size()) + " words");
  }

  // ---- refusals ----
  const auto refused = [&](const std::string& what, const std::string& source, std::size_t line, const std::string& names) {
    try {
      systems::dynamics(*language::read(source).find("w"));
      s.require(what, false, "built");
    } catch (const language::ScriptError& e) {
      s.require(what, e.at.line == line && std::string(e.what()).find(names) != std::string::npos, e.what());
    }
  };
  {
    // Milstein's scheme asks that no coefficient depend on another driven component, whichever source is declared first.
    const language::Script noises = language::read(R"(
system later {
  state X, Y : real
  noise U : wiener
  noise V : wiener
  change {
    dX = Y * dU
    dY = dV
  }
}
system earlier {
  state X, Y : real
  noise U : wiener
  noise V : wiener
  change {
    dX = dU
    dY = X * dV
  }
}
system apart {
  state X, Y, Z : real
  noise U : wiener
  noise V : wiener
  change {
    dX = X * dU
    dY = Z * Y * dV
    dZ = -Z * dt
  }
}
system pulled {
  state X : real
  noise W : wiener
  change { dX = -X * dt + 0.5 * dW }
}
)");
    const auto later = systems::stochastic(*noises.find("later")), earlier = systems::stochastic(*noises.find("earlier")), apart = systems::stochastic(*noises.find("apart"));
    s.require("noise is diagonal for Milstein's scheme only when no coefficient depends on another driven component, whichever source comes first; one that is not driven may enter",
              !later.diagonal && !earlier.diagonal && apart.diagonal);
    systems::StochasticOptions m;
    m.t1 = 1, m.milstein = true;
    bool refusedMilstein = false;
    try {
      systems::simulate(later, {0, 1}, m);
    } catch (const kernel::Refusal&) {
      refusedMilstein = true;
    }
    s.require("Milstein's scheme is refused where the noise does not commute", refusedMilstein);
    // A transient: the path is run before it is recorded, so a start far from the stationary law is forgotten.
    const auto pulled = systems::stochastic(*noises.find("pulled"));
    systems::StochasticOptions o;
    o.t1 = 1, o.step = 0.01, o.seed = 3;
    const auto direct = systems::simulate(pulled, {40.0}, o);
    o.transient = 30;
    const auto settled = systems::simulate(pulled, {40.0}, o);
    s.require("a stochastic run with a transient begins its record after it: from 40 the record begins near the stationary law, at the time asked",
              direct.x.front()[0] == 40 && std::fabs(settled.x.front()[0]) < 2 && settled.t.front() == 0 && settled.t.size() == direct.t.size(), std::to_string(settled.x.front()[0]));
    bool refusedGiven = false;
    try {
      std::vector<std::vector<double>> inc(100, std::vector<double>(1, 0.0));
      systems::simulate(pulled, {1.0}, o, &inc);
    } catch (const kernel::Refusal&) {
      refusedGiven = true;
    }
    s.require("a run that is given its increments is refused a transient", refusedGiven);
  }
  refused("a state with no equation is refused by name", "system w {\n  state x, y : real\n  change { x' = y }\n}\n", 1, "'y'");
  refused("a name that is not declared is refused at its equation", "system w {\n  state x : real\n  change {\n    x' = -k * x\n  }\n}\n", 4, "'k'");
  refused("a state with two equations is refused at the second", "system w {\n  state x : real\n  change {\n    x' = 1\n    x' = 2\n  }\n}\n", 5, "two equations");

  return s.done();
}
