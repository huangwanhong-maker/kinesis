// Checks of kinesis/catalogue: each model is read as what it is said to be, and does what is known of it.
//
// Every entry
// - declares one system of its own name, of the form it is listed under; it is written out and read back the same;
//   it has a title, a field, what it is and where it comes from; no two share a name.
// Flows
// - Lorenz, Chen and Thomas have a constant divergence, to which their exponents sum whatever the path:
//   -(sigma + 1 + beta), -(a + b - c), -3 b.
// - Roessler's two equilibria are at x = (c +- sqrt(c^2 - 4 a b)) / 2, y = -x / a, z = x / a.
// - The epidemic keeps s + i + r, and the share never infected solves ln(s / s0) = -(beta / gamma) (1 - s).
// - Van der Pol's cycle at mu = 0.1: amplitude 2 and period 2 pi (1 + mu^2 / 16), both to the order mu^4
//   (the Poincare-Lindstedt series; Strogatz 2015, Nonlinear Dynamics and Chaos, 2nd ed., on the oscillator).
// - FitzHugh and Nagumo's unit with its current has one equilibrium, and it is not stable.
// - Lotka and Volterra's quantity is kept along a path.
// - The Brusselator's equilibrium is (a, b / a) with the real part (b - 1 - a^2) / 2 of its eigenvalues.
// - Sel'kov's equilibrium is (b, b / (a + b^2)); Duffing's are a saddle at 0 and two stable ones at +-1.
// - Hopf's normal form ends on the circle of radius sqrt(mu) and turns once in 2 pi / omega.
// - Verhulst's growth against its closed form.
// Maps
// - The logistic map at r = 4 has the exponent ln 2; Henon's exponents sum to ln b; the standard map's to zero.
// Events
// - The integrate-and-fire neuron fires with the period tau ln((R I + E - vr) / (R I + E - theta)); Izhikevich's fires.
// Actions
// - The pendulum's small swings have the period 2 pi sqrt(l / g); Kepler's circular orbit of radius 1 closes after
//   2 pi, keeping its energy and its angular momentum.
// Noise
// - Ornstein and Uhlenbeck's stationary variance sigma^2 / (2 theta), in its standard error; the others are read as
//   stochastic differentials in Ito's calculus.
// Quantum
// - The driven atom's stationary excited population Omega^2 / (gamma^2 + 2 Omega^2); the cavity's sixteen levels.
// Rules
// - a^n b^n; balanced brackets, as many of 2 n letters as the n-th Catalan number, each with one history.
// Turbulence
// - The shell model's nonlinear term gives its energy and its second quadratic quantity no rate of change, at any
//   state, to rounding: this holds every sign and coefficient of the term to the model's (L'vov et al. 1998). Along
//   a path without viscosity and force the energy is kept to the integrator's error.
// - Forced and viscous, in its steady state: the mean power put in equals the mean power dissipated; the mean flux
//   of energy through a shell is the power put in less what the shells up to it dissipate, so that through the
//   shells of the cascade it is the same; the energy of the shells falls as a power of the wavenumber steeper than
//   Kolmogorov's 2/3; the flatness of a shell's amplitude grows down the cascade (intermittency); and the largest
//   exponent is positive.
// Templates
// - The candidates for a record of three components are the six plain flows of three in each of six orders, for two
//   the seven of two in both, for a map of two the two maps; a model whose parameters share the record's names is
//   named anew and has the same field.
// - A record of Roessler's system given under other names and in another order is matched to Roessler's among all
//   the catalogue's candidates, and its three numbers are recovered.

#include <kinesis/analysis/invariants.hpp>
#include <kinesis/analysis/statistics.hpp>
#include <kinesis/catalogue/models.hpp>
#include <kinesis/induction/match.hpp>
#include <kinesis/language/reader.hpp>
#include <kinesis/systems/dynamics.hpp>
#include <kinesis/systems/lagrangian.hpp>
#include <kinesis/systems/quantum.hpp>
#include <kinesis/systems/rewriting.hpp>
#include <kinesis/systems/simulate.hpp>
#include <kinesis/systems/stochastic.hpp>

#include <calx/verification/check.hpp>

#include <physis/quantum/lindblad.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <string>

using namespace kinesis;
using systems::Form;
using systems::State;

// What a model does over time: a flow or a map as declared; a system given by an action, the flow of its stationary paths.
static systems::Dynamics flow(const std::string& id) {
  const language::System sys = catalogue::system(id);
  return systems::dynamics(systems::form_of(sys) == Form::lagrangian ? systems::varied(sys, sys.name) : sys);
}

// The times at which component i of a recorded path rises through `level`, by linear interpolation.
static std::vector<double> risings(const systems::Trajectory& p, std::size_t i, double level = 0) {
  std::vector<double> out;
  for (std::size_t k = 1; k < p.x.size(); k++) {
    const double a = p.x[k - 1][i] - level, b = p.x[k][i] - level;
    if (a < 0 && b >= 0) out.push_back(p.t[k - 1] + (p.t[k] - p.t[k - 1]) * a / (a - b));
  }
  return out;
}

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;

  // ---- every entry ----
  {
    const std::map<std::string, Form> forms = {
        {"sabra", Form::flow}, {"lorenz", Form::flow}, {"rossler", Form::flow}, {"chen", Form::flow}, {"thomas", Form::flow}, {"sir", Form::flow}, {"hindmarsh_rose", Form::flow},
        {"van_der_pol", Form::flow}, {"fitzhugh_nagumo", Form::flow}, {"lotka_volterra", Form::flow}, {"brusselator", Form::flow}, {"selkov", Form::flow},
        {"duffing", Form::flow}, {"stuart_landau", Form::flow}, {"logistic_growth", Form::flow}, {"logistic_map", Form::map}, {"henon", Form::map},
        {"standard_map", Form::map}, {"integrate_and_fire", Form::flow}, {"izhikevich", Form::flow}, {"pendulum", Form::lagrangian}, {"kepler", Form::lagrangian},
        {"ornstein_uhlenbeck", Form::stochastic}, {"double_well", Form::stochastic}, {"geometric_brownian", Form::stochastic}, {"square_root_diffusion", Form::stochastic},
        {"driven_atom", Form::quantum}, {"jaynes_cummings", Form::quantum}, {"anbn", Form::rules}, {"balanced", Form::rules}};
    std::set<std::string> ids;
    std::string wrong, unlike, bare;
    for (const auto& e : catalogue::entries()) {
      ids.insert(e.id);
      if (e.title.empty() || e.field.empty() || e.about.empty() || e.source.empty()) bare += e.id + " ";
      try {
        const language::System sys = catalogue::system(e);
        const auto it = forms.find(e.id);
        if (it == forms.end() || systems::form_of(sys) != it->second) wrong += e.id + ": " + systems::name(systems::form_of(sys)) + "  ";
        const std::string written = language::write(sys);
        if (language::write(*language::read(written).find(e.id)) != written) unlike += e.id + " ";
      } catch (const kernel::Refusal& r) {
        wrong += e.id + ": " + r.what() + "  ";
      }
    }
    s.require("the catalogue holds thirty models, no two of one name", catalogue::entries().size() == 30 && ids.size() == 30, std::to_string(catalogue::entries().size()));
    s.require("each declares one system of its own name, of the form it is listed under", wrong.empty(), wrong);
    s.require("each is written out and read back the same", unlike.empty(), unlike);
    s.require("each has a title, a field, what it is and where it comes from", bare.empty(), bare);
    bool refused = false;
    try {
      catalogue::system("nobody");
    } catch (const kernel::Refusal&) {
      refused = true;
    }
    s.require("a name the catalogue does not hold is refused", refused && catalogue::find("lorenz") && !catalogue::find("nobody"));
  }

  // ---- flows ----
  {
    physis::dynamics::LyapunovOptions o;
    o.step = 0.005, o.transient = 10, o.time = 60;
    s.check("Lorenz's exponents sum to -(sigma + 1 + beta)", analysis::lyapunov(flow("lorenz"), {1, 1, 20}, o).sum + 41.0 / 3, 5e-5);
    o.step = 0.002;
    s.check("Chen's sum to -(a + b - c) = -10", analysis::lyapunov(flow("chen"), {1, 1, 20}, o).sum + 10, 5e-5);
    o.step = 0.01;
    s.check("Thomas's sum to -3 b", analysis::lyapunov(flow("thomas"), {1, 0.5, -0.5}, o).sum + 3 * 0.208186, 1e-7);
  }
  {
    const auto d = flow("rossler");
    const double a = 0.2, b = 0.2, c = 5.7, root = std::sqrt(c * c - 4 * a * b);
    const auto found = analysis::equilibria(d, analysis::Box{{-10, -40, -10}, {10, 10, 40}}, 400);
    double worst = found.size() == 2 ? 0 : 1;
    for (const auto& q : found) {
      const double x = q.x[0] < 1 ? (c - root) / 2 : (c + root) / 2;
      worst = std::max({worst, std::fabs(q.x[0] - x), std::fabs(q.x[1] + x / a), std::fabs(q.x[2] - x / a)});
    }
    s.check("Roessler's two equilibria at x = (c +- sqrt(c^2 - 4 a b)) / 2, y = -x / a, z = x / a", worst, 1e-9);
  }
  {
    const auto d = flow("sir");
    systems::Options o;
    o.t1 = 300, o.step = 0.01, o.every = 100;
    const auto p = systems::simulate(d, {0.99, 0.01, 0}, o);
    double drift = 0;
    for (const auto& x : p.x) drift = std::max(drift, std::fabs(x[0] + x[1] + x[2] - 1));
    s.check("the epidemic keeps s + i + r", drift, 1e-12);
    const double left = p.x.back()[0];
    s.check("the share never infected solves ln(s / s0) = -(beta / gamma) (1 - s)", std::log(left / 0.99) + 2.5 * (1 - left), 1e-8);
    s.require("and the model declares what it keeps", d.conserved_terms.size() == 1 && std::fabs(d.conserved(0, 0, {0.2, 0.3, 0.5}) - 1) < 1e-15);
  }
  {
    auto d = flow("van_der_pol");
    d.set_param("mu", 0.1);
    systems::Options o;
    o.t1 = 80, o.step = 0.002, o.transient = 300;
    const auto p = systems::simulate(d, {2, 0}, o);
    double top = 0;
    for (const auto& x : p.x) top = std::max(top, x[0]);
    const auto up = risings(p, 0);
    s.check("Van der Pol's cycle at mu = 0.1 has the amplitude 2", top - 2, 2e-3);
    s.check("and the period 2 pi (1 + mu^2 / 16)", (up.back() - up.front()) / static_cast<double>(up.size() - 1) - 2 * pi * (1 + 0.01 / 16), 5e-5);
  }
  {
    const auto found = analysis::equilibria(flow("fitzhugh_nagumo"), analysis::Box{{-3, -3}, {3, 3}}, 200);
    s.require("FitzHugh and Nagumo's unit with its current has one equilibrium, and it is not stable: it fires", found.size() == 1 && !found[0].stable, std::to_string(found.size()));
  }
  {
    const auto d = flow("lotka_volterra");
    systems::Options o;
    o.t1 = 40, o.step = 0.002, o.every = 50;
    const auto p = systems::simulate(d, {2, 1}, o);
    double drift = 0;
    for (const auto& x : p.x) drift = std::max(drift, std::fabs(d.conserved(0, 0, x) - d.conserved(0, 0, p.x.front())));
    s.check("Lotka and Volterra's quantity d x - c ln x + b y - a ln y is kept along a path", drift, 1e-9);
  }
  {
    auto d = flow("brusselator");
    const auto at = analysis::equilibrium_at(d, {1, 3});
    s.check("the Brusselator's equilibrium is (a, b / a), where its eigenvalues have the real part (b - 1 - a^2) / 2", std::fabs(at.eigenvalues[0].real() - 0.5) + at.residual, 1e-12);
    d.set_param("b", 1.5);
    s.require("it is stable below b = 1 + a^2 and not above", analysis::equilibrium_at(d, {1, 1.5}).stable && !at.stable);
    const auto sel = analysis::equilibrium_at(flow("selkov"), {0.6, 0.6 / (0.08 + 0.36)});
    s.check("Sel'kov's equilibrium is (b, b / (a + b^2))", sel.residual, 1e-14);
    const auto wells = analysis::equilibria(flow("duffing"), analysis::Box{{-2, -2}, {2, 2}}, 200);
    s.require("Duffing's oscillator with alpha < 0 has a saddle at 0 and two stable equilibria at +-1",
              wells.size() == 3 && wells[1].kind == "saddle" && std::fabs(wells[1].x[0]) < 1e-9 && wells[0].stable && wells[2].stable && std::fabs(std::fabs(wells[0].x[0]) - 1) < 1e-9, std::to_string(wells.size()));
  }
  {
    const auto d = flow("stuart_landau");
    systems::Options o;
    o.t1 = 20, o.step = 0.002, o.transient = 60;
    const auto p = systems::simulate(d, {0.1, 0}, o);
    const auto up = risings(p, 1);
    s.check("Hopf's normal form ends on the circle of radius sqrt(mu)", std::hypot(p.x.back()[0], p.x.back()[1]) - std::sqrt(0.5), 1e-9);
    s.check("and turns once in 2 pi / omega", (up.back() - up.front()) / static_cast<double>(up.size() - 1) - pi, 1e-6);
    systems::Options g;
    g.t1 = 3, g.step = 0.001;
    s.check("Verhulst's growth against K / (1 + (K / x0 - 1) exp(-r t))", systems::simulate(flow("logistic_growth"), {0.5}, g).x.back()[0] - 10 / (1 + 19 * std::exp(-3.0)), 1e-10);
    s.require("Hindmarsh and Rose's neuron bursts: its potential rises through 1 many times, and stays bounded", [&] {
      systems::Options h;
      h.t1 = 1500, h.step = 0.01, h.every = 5, h.transient = 500;
      const auto q = systems::simulate(flow("hindmarsh_rose"), {-1, 0, 3}, h);
      double worst = 0;
      for (const auto& x : q.x) worst = std::max(worst, std::fabs(x[0]));
      return risings(q, 0, 1).size() > 10 && worst < 5;
    }());
  }

  // ---- turbulence ----
  {
    const catalogue::ShellModel model;
    const auto d = flow("sabra");
    const std::size_t N = model.shells;
    State x0(d.dimension(), 0.0);
    for (const auto& [name, value] : catalogue::start_of(model)) x0[static_cast<std::size_t>(d.state_index(name))] = value;
    const auto index = [&](const std::string& name) {
      for (std::size_t i = 0; i < d.observables.size(); i++)
        if (d.observables[i] == name) return i;
      return d.observables.size();
    };
    s.require("the shell model has two components a shell, and observes each shell's energy, the whole, the second quantity, the power in and out, and two fluxes",
              d.dimension() == 2 * N && d.observables.size() == N + 6 && index("flux6") < d.observables.size());
    {
      // The identity itself, on the field: at any state the nonlinear term gives the energy and the second quantity
      // no rate of change, Σ (x_n x_n' + y_n y_n') = 0 and Σ (a/c)^n (x_n x_n' + y_n y_n') = 0, to rounding.
      auto free = d;
      free.set_param("nu", 0), free.set_param("f", 0);
      State rate(free.dimension());
      double energy = 0, second = 0, size = 0;
      for (int trial = 0; trial < 5; trial++) {
        State x = x0;
        for (std::size_t i = 0; i < x.size(); i++) x[i] *= 1 + 0.3 * std::sin(1.3 * static_cast<double>(i * (trial + 1)));
        free.field(0, x, rate);
        double de = 0, ds = 0, scale = 0;
        for (std::size_t n = 1; n <= N; n++) {
          const double part = x[2 * n - 2] * rate[2 * n - 2] + x[2 * n - 1] * rate[2 * n - 1];
          de += part, ds += std::pow(model.a / model.c, static_cast<double>(n)) * part;
          scale += std::pow(2.0, static_cast<double>(n)) * (std::fabs(x[2 * n - 2] * rate[2 * n - 2]) + std::fabs(x[2 * n - 1] * rate[2 * n - 1]));
        }
        energy = std::max(energy, std::fabs(de)), second = std::max(second, std::fabs(ds)), size = std::max(size, scale);
      }
      s.check("without viscosity and force the nonlinear term gives the energy no rate of change, at five states", energy / size, 1e-14);
      s.check("nor the second quadratic quantity, the sum of (a/c)^n |u_n|^2", second / size, 1e-14);
      // And along a path, to the error of the integrator at the fastest shell.
      systems::Options o;
      o.t1 = 1, o.step = 2e-5, o.every = 5000;
      const auto p = systems::simulate(free, x0, o);
      double drift = 0;
      for (const auto& x : p.x) drift = std::max(drift, std::fabs(free.observe(index("energy"), 0, x) - free.observe(index("energy"), 0, x0)));
      s.check("and the energy is kept along a path", drift, 1e-9);
    }
    analysis::StatisticsOptions o;
    o.step = 5e-4, o.transient = 40, o.time = 200, o.every = 10;
    const auto st = analysis::averages(d, x0, o);
    const auto mean = [&](const std::string& name) { return st.of[index(name)].mean; };
    const double input = mean("input"), dissipated = mean("dissipation");
    s.check("in the steady state the power put in is the power dissipated, within three of their errors", (input - dissipated) / (st.of[index("input")].error + st.of[index("dissipation")].error), 3);
    s.check("which is to three hundredths over this time", (input - dissipated) / input, 0.03);
    for (std::size_t through : {std::size_t{3}, std::size_t{6}}) {
      double below = 0;
      for (std::size_t n = 1; n <= through; n++) below += 2 * model.nu * std::pow(catalogue::wavenumber(model, n), 2) * mean("e" + std::to_string(n));
      s.check("the flux of energy through shell " + std::to_string(through) + " is the power put in less what the shells up to it dissipate, to a hundredth",
              (mean("flux" + std::to_string(through)) + below - input) / input, 0.01);
    }
    s.check("so the flux is the same through shells 3 and 6 of the cascade, to a hundredth", (mean("flux3") - mean("flux6")) / input, 0.01);
    std::vector<double> k, e;
    for (std::size_t n = 3; n <= 11; n++) k.push_back(catalogue::wavenumber(model, n)), e.push_back(mean("e" + std::to_string(n)));
    const auto law = analysis::scaling(k, e);
    s.require("the energy of the shells falls as a power of the wavenumber over eight octaves, steeper than Kolmogorov's 2/3 and less than 0.9",
              -law.exponent > 2.0 / 3 && -law.exponent < 0.9 && law.error < 0.05, std::to_string(-law.exponent) + " +- " + std::to_string(law.error));
    s.require("the flatness of a shell grows down the cascade: about 2, the Gaussian value, at shell 4, and half as much again or more at shell 9",
              st.of[index("e4")].second > 1.6 && st.of[index("e4")].second < 2.4 && st.of[index("e9")].second > 1.5 * st.of[index("e4")].second,
              std::to_string(st.of[index("e4")].second) + ", " + std::to_string(st.of[index("e9")].second));
    physis::dynamics::LyapunovOptions lo;
    lo.count = 1, lo.step = 5e-4, lo.transient = 0, lo.time = 40;
    s.require("its largest exponent is positive: the motion is chaotic", analysis::lyapunov(d, st.final_state, lo).exponents[0] > 1);
    s.require("a path started from zero is refused nothing and stays on the forced shell, which is why a start moves every shell", [&] {
      systems::Options z;
      z.t1 = 0.5, z.step = 5e-4, z.every = 1000;
      const State end = systems::simulate(d, State(d.dimension(), 0.0), z).x.back();
      double others = 0;
      for (std::size_t i = 2; i < end.size(); i++) others = std::max(others, std::fabs(end[i]));
      return others == 0 && end[0] > 0.1;
    }());
  }

  // ---- maps ----
  {
    physis::dynamics::LyapunovOptions o;
    o.transient = 1000, o.time = 100000, o.renormalise_every = 1;
    auto logistic = flow("logistic_map");
    logistic.set_param("r", 4);
    s.check("the logistic map at r = 4 has the exponent ln 2", analysis::lyapunov(logistic, {0.3141}, o).exponents[0] - std::numbers::ln2, 2e-3);
    s.check("Henon's exponents sum to ln b", analysis::lyapunov(flow("henon"), {0.1, 0.1}, o).sum - std::log(0.3), 1e-9);
    o.time = 20000;
    s.check("the standard map keeps area: its exponents sum to zero", analysis::lyapunov(flow("standard_map"), {0.3, 1.2}, o).sum, 1e-9);
  }

  // ---- events ----
  {
    systems::Options o;
    o.t1 = 100, o.step = 0.01, o.every = 100;
    const auto p = systems::simulate(flow("integrate_and_fire"), {-65}, o);
    s.check("the integrate-and-fire neuron fires with the period tau ln((R I + E - vr) / (R I + E - theta)) = 10 ln 4",
            p.events.size() >= 5 ? p.events[4].t - p.events[3].t - 10 * std::log(4.0) : 1.0, 1e-8);
    o.t1 = 300, o.step = 0.05;
    s.require("Izhikevich's neuron fires and is reset", systems::simulate(flow("izhikevich"), {-65, -13}, o).events.size() >= 5);
  }

  // ---- actions ----
  {
    systems::Options o;
    o.t1 = 30, o.step = 0.001;
    const auto p = systems::simulate(flow("pendulum"), {0.01, 0}, o);
    const auto up = risings(p, 0);
    s.check("the pendulum's small swings have the period 2 pi sqrt(l / g)", (up.back() - up.front()) / static_cast<double>(up.size() - 1) - 2 * pi / std::sqrt(9.81), 2e-5);
    const auto d = flow("kepler");
    s.require("Kepler's problem has the components x, x', y, y'", d.states == std::vector<std::string>({"x", "x'", "y", "y'"}));
    systems::Options k;
    k.t1 = 2 * pi, k.step = 2 * pi / 20000;
    const auto orbit = systems::simulate(d, {1, 0, 0, 1}, k);
    const State& end = orbit.x.back();
    s.check("its circular orbit of radius 1 closes after 2 pi", std::hypot(end[0] - 1, end[2]), 1e-9);
    double energy = 0, momentum = 0;
    for (const auto& x : orbit.x) {
      energy = std::max(energy, std::fabs((x[1] * x[1] + x[3] * x[3]) / 2 - 1 / std::hypot(x[0], x[2]) + 0.5));
      momentum = std::max(momentum, std::fabs(x[0] * x[3] - x[2] * x[1] - 1));
    }
    s.check("keeping its energy -1/2 and its angular momentum 1", energy + momentum, 1e-10);
  }

  // ---- noise ----
  {
    const auto d = systems::stochastic(catalogue::system("ornstein_uhlenbeck"));
    systems::StochasticOptions o;
    o.t1 = 4000, o.step = 0.01, o.every = 10, o.transient = 20, o.seed = 7;
    const auto p = systems::simulate(d, {0}, o);
    double sum = 0, squares = 0;
    for (const auto& x : p.x) sum += x[0], squares += x[0] * x[0];
    const double n = static_cast<double>(p.x.size()), variance = squares / n - (sum / n) * (sum / n);
    // Samples a tenth of a time apart are correlated over a time 1 / theta: about n / 20 of them are independent,
    // and the variance of a normal sample of m has the standard error v sqrt(2 / m).
    s.check("Ornstein and Uhlenbeck's stationary variance sigma^2 / (2 theta), in its standard errors", (variance - 0.125) / (0.125 * std::sqrt(2 / (n / 20))), 3);
    std::string unlike;
    for (const char* id : {"double_well", "geometric_brownian", "square_root_diffusion"}) {
      const auto sd = systems::stochastic(catalogue::system(id));
      if (sd.dimension() != 1 || sd.calculus == "stratonovich") unlike += std::string(id) + " ";
    }
    s.require("the other three are differentials of one component, read in Ito's calculus", unlike.empty(), unlike);
  }

  // ---- quantum ----
  {
    const auto atom = systems::quantum(catalogue::system("driven_atom"));
    const auto st = physis::quantum::stationary(atom.generator);
    s.check("the driven atom's stationary excited population Omega^2 / (gamma^2 + 2 Omega^2) = 9/19", atom.expectations(st.state)[0] - 9.0 / 19, 1e-12);
    const auto cavity = systems::quantum(catalogue::system("jaynes_cummings"));
    s.require("Jaynes and Cummings's model has two factors and sixteen levels, and two observables", cavity.dimension == 16 && cavity.factors.size() == 2 && cavity.observables.size() == 2);
  }

  // ---- rules ----
  {
    const auto r = systems::rewriting(catalogue::system("balanced"));
    systems::Bounds b;
    b.max_length = 8;
    const auto words = systems::language(r, systems::explore(r, b));
    std::map<std::size_t, std::size_t> count;
    bool one = true;
    for (const auto& w : words) count[w.word.size()]++, one = one && w.histories == 1;
    s.require("balanced brackets: 1, 2, 5 and 14 words of 2, 4, 6 and 8 letters, the Catalan numbers, each with one history",
              count[2] == 1 && count[4] == 2 && count[6] == 5 && count[8] == 14 && one, std::to_string(words.size()) + " words");
    const auto ab = systems::rewriting(catalogue::system("anbn"));
    systems::Bounds c;
    c.max_length = 10;
    s.require("a^n b^n: five words of at most ten letters", systems::language(ab, systems::explore(ab, c)).size() == 5);
  }

  // ---- templates ----
  {
    const auto three = catalogue::candidates({"u", "v", "w"}, false), two = catalogue::candidates({"p", "q"}, false), maps = catalogue::candidates({"p", "q"}, true);
    std::set<std::string> names3, names2;
    bool permuted = true;
    for (const auto& c : three) {
      names3.insert(c.name);
      auto st = catalogue::states_of(c);
      std::sort(st.begin(), st.end());
      permuted = permuted && st == std::vector<std::string>({"u", "v", "w"});
    }
    for (const auto& c : two) names2.insert(c.name);
    s.require("the candidates for a record of three components are the six plain flows of three, each in six orders", three.size() == 36 && names3.size() == 6 && permuted, std::to_string(three.size()));
    s.require("for two components the seven of two in both orders, and for a map of two the two maps", two.size() == 14 && names2.size() == 7 && maps.size() == 4, std::to_string(two.size()) + ", " + std::to_string(maps.size()));
    s.require("a model with events, noise or an action is not offered: it is not plain", !catalogue::plain(catalogue::system("izhikevich")) && !catalogue::plain(catalogue::system("double_well")) && !catalogue::plain(catalogue::system("pendulum")) && catalogue::plain(catalogue::system("sir")));

    // Roessler's parameters are a, b, c: named after a record whose components are c, a, b, they move aside.
    const auto moved = catalogue::renamed(catalogue::system("rossler"), {"c", "a", "b"});
    const auto original = flow("rossler"), again = systems::dynamics(moved);
    State f(3), g(3);
    original.field(0, {0.3, -1.1, 2.4}, f), again.field(0, {0.3, -1.1, 2.4}, g);
    s.check("a model whose parameters share the record's names is named anew, and has the same field", std::fabs(f[0] - g[0]) + std::fabs(f[1] - g[1]) + std::fabs(f[2] - g[2]), 1e-15);
    s.require("its states are the record's and its parameters have moved aside", again.states == std::vector<std::string>({"c", "a", "b"}) && again.param_index("a_") >= 0 && again.param_index("a") < 0);

    // A record of Roessler's system under other names, in the order (z, x, y).
    systems::Options o;
    o.t1 = 60, o.step = 0.005, o.transient = 50, o.every = 2;
    const auto path = systems::simulate(original, {1, 1, 0}, o);
    systems::Trajectory record;
    record.names = {"u", "v", "w"}, record.t = path.t;
    for (const auto& x : path.x) record.x.push_back({x[2], x[0], x[1]});
    std::vector<language::System> templates;
    for (const auto& c : three) templates.push_back(induction::as_template(c));
    const auto ranking = induction::match(templates, record);
    const auto& best = ranking.front();
    std::string order;
    for (const auto& n : catalogue::states_of(best.fit.system)) order += n;
    s.require("a record of Roessler's system, under other names and in another order, is matched to Roessler's among the catalogue's thirty-six candidates",
              best.fitted && best.name == "rossler" && order == "vwu", best.name + " in the order " + order);
    const auto& found = best.fit.numbers;
    s.check("and its three numbers are recovered", best.fitted ? std::fabs(found.at("a") - 0.2) + std::fabs(found.at("b") - 0.2) + std::fabs(found.at("c") - 5.7) : 1.0, 2e-3);
    std::size_t fitted = 0;
    for (const auto& c : ranking) fitted += c.fitted;
    s.require("the ranking has every candidate, and those that could not be fitted say why", ranking.size() == 36 && fitted >= 6 && std::all_of(ranking.begin(), ranking.end(), [](const auto& c) { return c.fitted || !c.refused.empty(); }),
              std::to_string(fitted) + " fitted");
  }

  return s.done();
}
