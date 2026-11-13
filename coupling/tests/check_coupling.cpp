// Checks of kinesis/coupling: systems made of systems.
//
// Merged
// - A response driven by Lorenz's x (Pecora and Carroll 1990, Phys. Rev. Lett. 64, 821–824): the composite of the
//   drive and the response (y, z) is one flow of five components; the response falls onto the drive, and the
//   composite's exponents hold the drive's three and the response's two conditional ones, −1.81 and −1.86 (their
//   Table I, the x drive).
// - The merged declaration, written out, reads back and has the same field.
// - Parameters set on a part are those of the part; an input that is not wired keeps its own value.
// - Events of the parts are events of the whole: two integrators that fire, each with its period.
// Co-simulated
// - Two flows side by side against their merged flow: the error falls by 2 when the communication step is
//   halved (first order), and the order of the parts does not matter.
// - A map that drives a flow is exact: a counter n and x' = n give x(KH) = H K (K − 1) / 2.
// - A stochastic differential that drives a flow: an Ornstein–Uhlenbeck process through a filter of rate a has
//   the variance σ²/(2θ) · a/(a + θ).
// Joined by an action
// - Two equal oscillators joined by the action −k/2 (a.q − b.q)²: the whole is one Lagrangian system, whose modes
//   have the frequencies ω and √(ω² + 2k); from a displaced and b at rest, a = ½(cos ωt + cos √(ω² + 2k) t) and
//   b = ½(cos ωt − cos √(ω² + 2k) t), and the energy of the whole is kept.
// - Two free masses joined so keep their total momentum (the interaction depends on the difference alone:
//   Noether); joined by a wire one way, b.f = k (a.q − b.q), they do not.
// - The joined system is a declaration by an action, which reads back; a part's parameters may be set.
// - Refusals: a part not given by an action, wires beside an action, an action that names no coordinate, and
//   more coordinates than are solved as terms.
// Refusals with their reasons: a part the script does not declare, a wire to what is not an input, a map merged.

#include <kinesis/analysis/invariants.hpp>
#include <kinesis/coupling/composite.hpp>
#include <kinesis/language/reader.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

using namespace kinesis;
using systems::State;

static const char* text = R"(
system lorenz {
  state x, y, z : real
  param sigma = 10, rho = 28, beta = 8/3
  change {
    x' = sigma * (y - x)
    y' = x * (rho - z) - y
    z' = x * y - beta * z
  }
}
system response {
  state y, z : real
  input u
  param rho = 28, beta = 8/3
  change {
    y' = u * (rho - z) - y
    z' = u * y - beta * z
  }
}
system synchronised {
  part drive : lorenz
  part follower : response
  connect { follower.u = drive.x }
}

system oscillator {
  state q, p : real
  input f = 0.25
  param omega = 1
  change {
    q' = p
    p' = -omega^2 * q + f
  }
}
system pair {
  param k = 0.3
  part a : oscillator
  part b : oscillator with omega = 1.7
  connect {
    a.f = k * (b.q - a.q)
    b.f = k * (a.q - b.q)
  }
}
system alone {
  part a : oscillator with omega = 2
}

system integrator {
  state v : real
  input I = 1.5
  param tau = 0.02
  change {
    v' = (I - v) / tau
    when v >= 1 : v -> 0
  }
}
system two_integrators {
  part slow : integrator
  part fast : integrator with tau = 0.01
}

system counter {
  state n : real
  change { n' = n + 1 }
  structure { time k in natural }
}
system accumulator {
  state x : real
  input u
  change { x' = u }
}
system sampled {
  part clock : counter
  part plant : accumulator
  connect { plant.u = clock.n }
}
system wrapping {
  state n : real
  change {
    n' = n + 1
    when n >= 3 : n -> 0
  }
  structure { time k in natural }
}
system sawtooth {
  part plant : accumulator
  part clock : wrapping
  connect { plant.u = clock.n }
}

system ou {
  state X : real
  param theta = 1.5, s = 0.6
  noise W : wiener
  change { dX = -theta * X * dt + s * dW }
}
system filter {
  state y : real
  input u
  param a = 4
  change { y' = a * (u - y) }
}
system filtered {
  part source : ou
  part low : filter
  connect { low.u = source.X }
}

system broken_part { part a : nothing }
system broken_wire {
  part a : oscillator
  connect { a.g = 1 }
}
system mixed {
  part clock : counter
  part plant : accumulator
}

system mass {
  state q(t) : real
  param m = 1, w = 1
  change { action S = integral(m * q'^2 / 2 - m * w^2 * q^2 / 2) }
  structure { time t in real ; histories stationary S }
}
system spring_pair {
  param k = 0.5
  part a : mass
  part b : mass
  change { action S = integral(-k / 2 * (a.q - b.q)^2) }
  structure { time t in real ; histories stationary S }
}
system free_pair {
  param k = 0.5
  part a : mass with w = 0
  part b : mass with w = 0, m = 3
  change { action S = integral(-k / 2 * (a.q - b.q)^2) }
  structure { time t in real ; histories stationary S }
}
system puck {
  state q : real
  input f = 0
  change { q'' = f }
}
system led {
  param k = 0.5
  part a : puck
  part b : puck
  connect { b.f = k * (a.q - b.q) }
}
system four {
  part a : mass
  part b : mass
  part c : mass
  part d : mass
  change { action S = integral(-(a.q - b.q)^2 - (c.q - d.q)^2) }
  structure { time t in real ; histories stationary S }
}
system not_an_action {
  part a : mass
  part b : puck
  change { action S = integral(-(a.q - b.q)^2) }
  structure { time t in real ; histories stationary S }
}
system both_ways {
  part a : mass
  part b : mass
  change { action S = integral(-(a.q - b.q)^2) }
  structure { time t in real ; histories stationary S }
  connect { b.f = a.q }
}
system no_such {
  part a : mass
  part b : mass
  change { action S = integral(-(a.q - b.z)^2) }
  structure { time t in real ; histories stationary S }
}
)";

int main() {
  calx::verification::Suite s;
  const language::Script script = language::read(text);
  const auto sys = [&](const char* n) -> const language::System& { return *script.find(n); };

  // ---- merged ----
  {
    const language::System whole = coupling::merged(script, sys("synchronised"));
    const systems::Dynamics d = systems::dynamics(whole);
    s.require("a drive and a response merge into one flow of five components", d.states == std::vector<std::string>{"drive_x", "drive_y", "drive_z", "follower_y", "follower_z"}, language::write(whole));
    systems::Options o;
    o.t1 = 25, o.step = 0.002, o.every = 100000000;
    const State end = systems::simulate(d, {1, 1, 20, -7, 35}, o).x.back();
    s.check("the response falls onto the drive: y and z after 25", std::fabs(end[3] - end[1]) + std::fabs(end[4] - end[2]), 1e-9);
    physis::dynamics::LyapunovOptions lo;
    lo.step = 0.005, lo.transient = 20, lo.time = 300;
    const auto e = analysis::lyapunov(d, {1, 1, 20, -7, 35}, lo);
    s.check("the composite's exponents hold the response's conditional ones, -1.81 and -1.86, within the error of a finite time", std::fabs(e.exponents[2] + 1.81) + std::fabs(e.exponents[3] + 1.86), 0.12);
    s.check("the two sum to the divergence of the response, -(1 + beta)", e.exponents[2] + e.exponents[3] + 11.0 / 3, 0.02);
    s.check("and their sum is the divergence of the whole, -(sigma + 1 + beta) - (1 + beta)", e.sum + 41.0 / 3 + 11.0 / 3, 1e-4);
    const language::Script back = language::read(language::write(whole));
    const systems::Dynamics again = systems::dynamics(*back.find("synchronised_merged"));
    State a(5), b(5);
    d.field(0, {1.5, -2, 20, 0.3, 7}, a), again.field(0, {1.5, -2, 20, 0.3, 7}, b);
    s.require("the merged declaration, written out, reads back to the same field", a == b);
  }
  {
    const systems::Dynamics d = systems::dynamics(coupling::merged(script, sys("alone")));
    State dx(2);
    d.field(0, {1, 0}, dx);
    s.check("a parameter set on a part and an input left alone: p' = -4 q + 0.25", dx[1] + 3.75, 1e-15);
    const systems::Dynamics p = systems::dynamics(coupling::merged(script, sys("pair")));
    p.field(0, {1, 0, -0.5, 0}, dx = State(4));
    s.check("two parts with their own parameters and the composite's k: a_p' = -1 + 0.3 (-1.5), b_p' = 1.7^2 / 2 + 0.3 (1.5)", std::fabs(dx[1] + 1.45) + std::fabs(dx[3] - (1.445 + 0.45)), 1e-14);
  }
  {
    const systems::Dynamics d = systems::dynamics(coupling::merged(script, sys("two_integrators")));
    systems::Options o;
    o.t1 = 0.2, o.step = 0.0005;
    const auto tr = systems::simulate(d, {0, 0}, o);
    std::size_t slow = 0, fast = 0;
    for (const auto& ev : tr.events) (ev.event == 0 ? slow : fast)++;
    s.require("the events of the parts are events of the whole: 9 of the slow integrator and 18 of the fast", slow == 9 && fast == 18, std::to_string(slow) + ", " + std::to_string(fast));
  }

  // ---- joined by an action ----
  {
    const systems::Dynamics d = systems::dynamics(coupling::merged(script, sys("spring_pair")));
    s.require("two oscillators joined by an action are one flow of the components a_q, a_q', b_q, b_q'", d.states == std::vector<std::string>({"a_q", "a_q'", "b_q", "b_q'"}) && d.conserved_terms.size() == 1);
    systems::Options o;
    o.t1 = 10, o.step = 0.001, o.every = 100;
    const auto tr = systems::simulate(d, {1, 0, 0, 0}, o);
    const double w = 1, wa = std::sqrt(2.0), T = 10;
    const State& end = tr.x.back();
    s.check("from a displaced and b at rest: a = (cos wt + cos sqrt(w^2 + 2k) t) / 2, b = (cos wt - cos sqrt(w^2 + 2k) t) / 2",
            std::fabs(end[0] - 0.5 * (std::cos(w * T) + std::cos(wa * T))) + std::fabs(end[2] - 0.5 * (std::cos(w * T) - std::cos(wa * T))), 1e-9);
    double drift = 0;
    for (const auto& x : tr.x) drift = std::max(drift, std::fabs(d.conserved(0, 0, x) - d.conserved(0, 0, tr.x.front())));
    s.check("the energy of the whole is kept", drift, 1e-11);
    s.check("and is the two parts' with the interaction's: 1/2 + k/2 at the start", d.conserved(0, 0, tr.x.front()) - 0.75, 1e-15);
    const language::System whole = coupling::joined_action(script, sys("spring_pair"));
    const std::string written = language::write(whole);
    s.require("the joined system is a declaration by an action, which reads back", systems::form_of(whole) == systems::Form::lagrangian && language::write(*language::read(written).find(whole.name)) == written, written);
  }
  {
    // Free masses 1 and 3: the interaction depends on a.q - b.q alone, so the total momentum q_a' + 3 q_b' is kept.
    const systems::Dynamics d = systems::dynamics(coupling::merged(script, sys("free_pair")));
    systems::Options o;
    o.t1 = 20, o.step = 0.001, o.every = 100;
    const auto tr = systems::simulate(d, {0, 1, 0.3, 0}, o);
    double drift = 0;
    for (const auto& x : tr.x) drift = std::max(drift, std::fabs(x[1] + 3 * x[3] - 1));
    s.check("two free masses joined by an action keep their total momentum", drift, 1e-11);
    const systems::Dynamics wired = systems::dynamics(coupling::merged(script, sys("led")));
    const auto led = systems::simulate(wired, {0, 1, 0.3, 0}, o);
    double change = 0;
    for (const auto& x : led.x) change = std::max(change, std::fabs(x[1] + x[3] - 1));
    s.require("joined by a wire one way they do not: the first is not moved by the second", change > 0.5 && std::fabs(led.x.back()[1] - 1) < 1e-12, std::to_string(change));
  }
  {
    const auto refused_join = [&](const char* name, const char* system, const char* part) {
      bool refused = false;
      try {
        coupling::merged(script, sys(system));
      } catch (const kernel::Refusal& r) {
        refused = std::string(r.what()).find(part) != std::string::npos;
      }
      s.require(name, refused);
    };
    refused_join("a part that is not given by an action is not joined by one", "not_an_action", "'b'");
    refused_join("wires beside an action are refused", "both_ways", "not both");
    refused_join("an action that names what is no coordinate of a part is refused", "no_such", "'b.z'");
    refused_join("four coordinates are not solved as terms, and that is said", "four", "beyond three");
  }

  // ---- co-simulated ----
  {
    const systems::Dynamics exact = systems::dynamics(coupling::merged(script, sys("pair")));
    systems::Options o;
    o.t1 = 10, o.step = 0.0005, o.every = 100000000;
    const State truth = systems::simulate(exact, {1, 0, -0.5, 0}, o).x.back();
    const coupling::Cosimulation c = coupling::cosimulation(script, sys("pair"));
    const auto error = [&](double H) {
      coupling::Cosimulation::Options co;
      co.t1 = 10, co.communication = H, co.step = 0.0005, co.every = 100000000;
      const State end = c.run({1, 0, -0.5, 0}, co).x.back();
      double worst = 0;
      for (std::size_t i = 0; i < 4; i++) worst = std::max(worst, std::fabs(end[i] - truth[i]));
      return worst;
    };
    const double e1 = error(0.02), e2 = error(0.01);
    s.require("the co-simulation names its components by part", c.names == std::vector<std::string>{"a.q", "a.p", "b.q", "b.p"});
    s.check("two flows side by side against their merged flow, at H = 0.01", e2, 0.02);
    s.check("the error is of the first order in the communication step", calx::verification::order(e1, e2) - 1, 0.1);
  }
  {
    const coupling::Cosimulation c = coupling::cosimulation(script, sys("sampled"));
    coupling::Cosimulation::Options co;
    co.t1 = 2, co.communication = 0.1, co.step = 0.01;
    const auto tr = c.run({0, 0}, co);
    double worst = 0;
    for (std::size_t k = 0; k < tr.x.size(); k++) worst = std::max(worst, std::fabs(tr.x[k][1] - 0.1 * static_cast<double>(k * (k - 1)) / 2));
    s.require("a map and a flow side by side: the counter reaches 20 in 20 steps", tr.x.back()[0] == 20 && c.parts[0].form == systems::Form::map && c.parts[1].form == systems::Form::flow);
    s.check("a map that drives a flow is exact: x(KH) = H K (K - 1) / 2", worst, 1e-12);
  }
  {
    // A map that wraps at 3, second among the parts: its events are events of the whole, with the whole's state
    // and at the whole's time.
    const coupling::Cosimulation c = coupling::cosimulation(script, sys("sawtooth"));
    coupling::Cosimulation::Options co;
    co.t1 = 1, co.communication = 0.1, co.step = 0.01;
    const auto tr = c.run({0, 0}, co);
    bool whole = tr.events.size() == 3;
    for (std::size_t e = 0; e < tr.events.size() && whole; e++) {
      const auto& ev = tr.events[e];
      whole = ev.before.size() == 2 && ev.after.size() == 2 && ev.before[1] == 3 && ev.after[1] == 0 && std::fabs(ev.t - 0.3 * static_cast<double>(e + 1)) < 1e-12 && ev.before[0] == ev.after[0];
    }
    s.require("a map's events in a co-simulation are the whole's: three wraps in ten steps, each with both components and at the time of its step", whole, std::to_string(tr.events.size()) + " events");
  }
  {
    const coupling::Cosimulation c = coupling::cosimulation(script, sys("filtered"));
    coupling::Cosimulation::Options co;
    co.t1 = 4000, co.communication = 0.005, co.step = 0.005, co.every = 4, co.seed = 17;
    const auto tr = c.run({0, 0}, co);
    double vx = 0, vy = 0;
    for (const auto& x : tr.x) vx += x[0] * x[0], vy += x[1] * x[1];
    s.check("a stochastic differential that drives a flow: the filtered variance against s^2 / (2 theta) a / (a + theta)", vy / static_cast<double>(tr.x.size()) / (0.36 / 3 * 4 / 5.5) - 1, 0.06);
    s.check("and the source's own variance", vx / static_cast<double>(tr.x.size()) / 0.12 - 1, 0.05);
  }

  // ---- refusals ----
  const auto refused = [&](const std::string& what, const std::function<void()>& attempt, const std::string& names) {
    try {
      attempt();
      s.require(what, false, "accepted");
    } catch (const std::exception& e) {
      s.require(what, std::string(e.what()).find(names) != std::string::npos, e.what());
    }
  };
  refused("a part the script does not declare is refused by name", [&] { coupling::merged(script, sys("broken_part")); }, "'nothing'");
  refused("a wire to what is not an input is refused", [&] { coupling::merged(script, sys("broken_wire")); }, "'g'");
  refused("the same in a co-simulation", [&] { coupling::cosimulation(script, sys("broken_wire")); }, "'g'");
  refused("a map is not merged with a flow: it is co-simulated", [&] { coupling::merged(script, sys("mixed")); }, "co-simulate");

  return s.done();
}
