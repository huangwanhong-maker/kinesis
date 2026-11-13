// Checks of kinesis/induction: from a record to a written model. Every record here is simulated, so that the
// law to be recovered is known.
//
// - Differences of second order are exact on a parabola; the monomials of n variables up to degree d number
//   C(n + d, d).
// - Discovery with no law proposed: from a record of Lorenz's system, a sparse law of the monomials up to degree
//   2 has Lorenz's seven terms and no other, with its coefficients; written out, it reads back and runs, and
//   follows the record for a time. The same for the logistic map and Hénon's.
// - The weak form: with noise on the record, the law tested against windows keeps the seven terms where
//   differences of the record do not, and its coefficients are closer.
// - Templates. A model's parameters freed and found again (Lorenz's σ, ρ, β); unknowns that are not linear
//   (the capacity of a logistic growth) by the Levenberg–Marquardt iteration; an unknown inside an action
//   (a pendulum's g/l) through its Euler–Lagrange equations.
// - Operators. Dynamic mode decomposition of a linear system gives its rates; on the dictionary {x, y, x²} the
//   system x' = μx, y' = λ(y − x²) has the Koopman eigenvalues μ, λ and 2μ (Brunton, Brunton, Proctor and Kutz
//   2016, PLoS ONE 11, e0150171).
// - A forecaster with no law: one step of Lorenz to a part in a thousand, and a free run that stays near the
//   truth for more than three Lyapunov times.
// - Stochastic: the rate and the noise of an Ornstein–Uhlenbeck process from a path, by a template and in bins.
// - Matching: of three candidate laws for a record of van der Pol's oscillator, van der Pol's ranks first; a
//   candidate whose states the record does not hold is listed with the reason it was refused.

#include <kinesis/induction/fit.hpp>
#include <kinesis/induction/match.hpp>
#include <kinesis/induction/operators.hpp>
#include <kinesis/induction/surrogate.hpp>
#include <kinesis/language/reader.hpp>
#include <kinesis/systems/simulate.hpp>
#include <kinesis/systems/stochastic.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

using namespace kinesis;
using induction::Trajectory;
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
}
system logistic {
  state x : real
  change { x' = 4 * x * (1 - x) }
  structure { time n in natural }
}
system henon {
  state x, y : real
  change {
    x' = 1 - 1.4 * x^2 + y
    y' = 0.3 * x
  }
  structure { time n in natural }
}
system growth {
  state x : real
  param r = 1.5, K = 4
  change { x' = r * x * (1 - x / K) }
}
system growth_template {
  state x : real
  unknown r, K : number
  change { x' = ?r * x * (1 - x / ?K) }
}
system pendulum {
  state q(t) : circle
  change { action S = integral(q'^2 / 2 + 19.62 * cos(q)) }
  structure { time t in real ; histories stationary S }
}
system pendulum_template {
  state q(t) : circle
  unknown k : number
  change { action S = integral(q'^2 / 2 + ?k * cos(q)) }
  structure { time t in real ; histories stationary S }
}
system focus {
  state x, y : real
  change {
    x' = -0.1 * x + 2 * y
    y' = -2 * x - 0.1 * y
  }
}
system slow_manifold {
  state x, y : real
  change {
    x' = -0.05 * x
    y' = -(y - x^2)
  }
}
system ou {
  state X : real
  param theta = 1.5, s = 0.6
  noise W : wiener
  change { dX = -theta * X * dt + s * dW }
}
system ou_template {
  state X : real
  unknown theta, s : number
  noise W : wiener
  change { dX = -?theta * X * dt + ?s * dW }
}
system vdp {
  state x : real
  param mu = 1.2
  change { x'' = mu * (1 - x^2) * x' - x }
}
system spring_template {
  state x : real
  unknown w : number
  change { x'' = -?w * x }
}
system duffing_template {
  state x : real
  unknown d, a, b : number
  change { x'' = -?d * x' - ?a * x - ?b * x^3 }
}
)";

static Trajectory run(const language::System& s, State x0, double t1, double step, std::size_t every = 1, double transient = 0) {
  systems::Options o;
  o.t1 = t1, o.step = step, o.every = every, o.transient = transient;
  return systems::simulate(systems::dynamics(s), std::move(x0), o);
}

// The coefficient of a term in a sparse law, by evaluating the law's derivative: for the checks' polynomials the
// law is read off by its values at chosen points instead, which needs no parsing of its form.
static double value(const kt::TermPtr& t, std::initializer_list<std::pair<const char*, double>> at) {
  return kt::evaluate(t, [&](const std::string& n, unsigned) {
    for (const auto& [name, v] : at)
      if (n == name) return v;
    return 0.0;
  });
}
static std::size_t count_terms(const kt::TermPtr& t) { return t->kind == kt::Kind::sum ? t->args.size() : kt::is_number(t, 0) ? 0 : 1; }

int main() {
  calx::verification::Suite s;
  const language::Script script = language::read(text);
  const auto sys = [&](const char* n) -> const language::System& { return *script.find(n); };

  // ---- differences and libraries ----
  {
    Trajectory p;
    for (int k = 0; k < 9; k++) {
      const double t = 0.3 * k + 0.01 * k * k;  // steps that differ
      p.t.push_back(t), p.x.push_back({2 - 3 * t + 0.5 * t * t});
    }
    double worst = 0;
    const auto d = induction::derivatives(p);
    for (std::size_t k = 0; k < p.t.size(); k++) worst = std::max(worst, std::fabs(d[k][0] - (-3 + p.t[k])));
    s.check("differences of second order are exact on a parabola, at the ends too, for steps that differ", worst, 1e-12);
    s.require("the monomials of 2 variables to degree 3 are 10, of 3 to degree 2 are 10, the constant first", induction::monomials({"x", "y"}, 3).size() == 10 && induction::monomials({"x", "y", "z"}, 2).size() == 10 &&
                                                                                                                 kt::is_number(induction::monomials({"x"}, 2)[0], 1) && induction::monomials({"x", "y"}, 2, false).size() == 5);
    s.require("a declared library poly(x, y, 2) is those monomials", induction::library_of(language::read_term("poly(x, y, 2)")).size() == 6);
  }

  // ---- discovery ----
  const Trajectory lorenz = run(sys("lorenz"), {1, 1, 20}, 20, 0.002, 1, 5);
  {
    induction::FitOptions o;
    o.threshold = 0.1;
    const induction::Fit f = induction::fit(induction::discovery_template(lorenz, 2), lorenz, o);
    const auto &fx = f.terms.at("f1"), &fy = f.terms.at("f2"), &fz = f.terms.at("f3");
    s.require("a sparse law of Lorenz's record has seven terms: 2, 3 and 2", count_terms(fx) == 2 && count_terms(fy) == 3 && count_terms(fz) == 2 && f.parameters == 7,
              kt::text(fx) + " | " + kt::text(fy) + " | " + kt::text(fz));
    const double worst = std::max({std::fabs(value(fx, {{"x", 1}}) + 10), std::fabs(value(fx, {{"y", 1}}) - 10), std::fabs(value(fy, {{"x", 1}}) - 28), std::fabs(value(fy, {{"y", 1}}) + 1),
                                   std::fabs(value(fy, {{"x", 1}, {"z", 1}}) - 27), std::fabs(value(fz, {{"z", 1}}) + 8.0 / 3), std::fabs(value(fz, {{"x", 1}, {"y", 1}}) - 1)});
    s.check("their coefficients: -10, 10; 28, -1, -1; -8/3, 1, to the error of the differences", worst, 0.02);
    s.check("the equations hold on the record to the error of its differences", f.residual, 1e-3);
    const language::Script back = language::read(language::write(f.system));
    const systems::Dynamics found = systems::dynamics(*back.find("discovered"));
    systems::Options so;
    so.t1 = 1, so.step = 0.002;
    const auto again = systems::simulate(found, lorenz.x[0], so);
    double gap = 0;
    for (std::size_t k = 0; k < again.x.size(); k++)
      for (std::size_t i = 0; i < 3; i++) gap = std::max(gap, std::fabs(again.x[k][i] - lorenz.x[k][i]));
    s.check("the law found, written out and read back, follows the record for a unit of time", gap, 0.05);
  }
  {
    systems::Options so;
    so.t1 = 400;
    const Trajectory map = systems::simulate(systems::dynamics(sys("logistic")), {0.2}, so);
    const induction::Fit f = induction::fit(induction::discovery_template(map, 3, "found", true), map);
    s.require("the logistic map from 400 iterates: 4 x - 4 x^2", count_terms(f.terms.at("f1")) == 2 && std::fabs(value(f.terms.at("f1"), {{"x", 0.5}}) - 1) < 1e-9 && std::fabs(value(f.terms.at("f1"), {{"x", 1}})) < 1e-9, kt::text(f.terms.at("f1")));
    const Trajectory h = systems::simulate(systems::dynamics(sys("henon")), {0.1, 0.1}, so);
    const induction::Fit g = induction::fit(induction::discovery_template(h, 2, "found", true), h);
    s.require("Henon's map: 1 - 1.4 x^2 + y and 0.3 x", count_terms(g.terms.at("f1")) == 3 && count_terms(g.terms.at("f2")) == 1 && std::fabs(value(g.terms.at("f1"), {{"x", 1}}) + 0.4) < 1e-9 &&
                                                              std::fabs(value(g.terms.at("f2"), {{"x", 1}}) - 0.3) < 1e-9,
              kt::text(g.terms.at("f1")) + " | " + kt::text(g.terms.at("f2")));
    s.require("a map found is a map of the language: its time is in the naturals", systems::form_of(g.system) == systems::Form::map);
  }

  // ---- the weak form under noise ----
  {
    Trajectory noisy = lorenz;
    calx::random::Stream r(9, 1, 0, 0);
    for (auto& x : noisy.x)
      for (double& v : x) v += 0.05 * r.normal();
    induction::FitOptions plain, weak;
    plain.threshold = weak.threshold = 0.3;
    weak.weak = true, weak.half = 60, weak.stride = 10;
    const auto tmpl = induction::discovery_template(noisy, 2);
    const induction::Fit a = induction::fit(tmpl, noisy, plain), b = induction::fit(tmpl, noisy, weak);
    const auto error = [&](const induction::Fit& f) {
      const auto &fx = f.terms.at("f1"), &fy = f.terms.at("f2"), &fz = f.terms.at("f3");
      return std::max({std::fabs(value(fx, {{"x", 1}}) + 10), std::fabs(value(fx, {{"y", 1}}) - 10), std::fabs(value(fy, {{"x", 1}}) - 28), std::fabs(value(fy, {{"y", 1}}) + 1),
                       std::fabs(value(fz, {{"z", 1}}) + 8.0 / 3), std::fabs(value(fz, {{"x", 1}, {"y", 1}}) - 1)});
    };
    s.require("with noise on the record the weak form keeps the seven terms", b.parameters == 7, std::to_string(b.parameters) + " terms by the weak form, " + std::to_string(a.parameters) + " by differences");
    s.check("its coefficients", error(b), 0.15);
    s.require("and they are closer than those from differences of the record", error(b) < error(a), std::to_string(error(b)) + " < " + std::to_string(error(a)));
  }

  // ---- templates ----
  {
    const induction::Fit f = induction::fit(induction::as_template(sys("lorenz")), lorenz);
    s.check("Lorenz's parameters freed and found: sigma, rho, beta", std::max({std::fabs(f.numbers.at("sigma") - 10), std::fabs(f.numbers.at("rho") - 28), std::fabs(f.numbers.at("beta") - 8.0 / 3)}), 5e-3);
    s.require("the fit is linear, by least squares, and the model returned has parameters no longer", f.linear && f.method == "least squares" && f.system.unknowns.empty() && f.system.params.empty(), f.method);
    const induction::Fit one = induction::fit(induction::as_template(sys("lorenz"), {"rho"}), lorenz);
    s.check("one parameter freed, the others kept", one.numbers.at("rho") - 28, 5e-3);
    s.require("the others are still parameters", one.system.params.size() == 2 && one.numbers.size() == 1);
  }
  {
    const Trajectory g = run(sys("growth"), {0.1}, 8, 0.01);
    const induction::Fit f = induction::fit(sys("growth_template"), g);
    s.require("an unknown that is not linear is found by iteration", !f.linear && f.method == "Levenberg-Marquardt", f.method + ", " + std::to_string(f.iterations) + " iterations");
    s.check("the rate and the capacity of a logistic growth", std::max(std::fabs(f.numbers.at("r") - 1.5), std::fabs(f.numbers.at("K") - 4)), 1e-3);
  }
  {
    const Trajectory swing = run(systems::varied(sys("pendulum")), {2.0, 0.0}, 10, 0.002);
    Trajectory angle;  // the angle alone: its velocity is not recorded
    angle.names = {"q"}, angle.t = swing.t;
    for (const auto& x : swing.x) angle.x.push_back({x[0]});
    const induction::Fit f = induction::fit(sys("pendulum_template"), angle);
    s.check("an unknown inside an action, through its Euler-Lagrange equations: g / l of a pendulum", f.numbers.at("k") - 19.62, 2e-3);
    s.require("it is returned in the action", f.system.change.actions.size() == 1 && f.system.unknowns.empty() && systems::form_of(f.system) == systems::Form::lagrangian, language::write(f.system));
  }

  // ---- operators ----
  {
    const Trajectory lin = run(sys("focus"), {1, 0.5}, 20, 0.001, 50);
    const induction::Linear A = induction::dmd(lin);
    double worst = 0;
    for (const auto& rate : A.rates) worst = std::max(worst, std::fabs(rate.real() + 0.1) + std::fabs(std::fabs(rate.imag()) - 2));
    s.check("dynamic mode decomposition of a linear focus: the rates -0.1 +- 2i", worst, 1e-8);
    const Trajectory slow = run(sys("slow_manifold"), {1.5, -1}, 30, 0.001, 100);
    const induction::Linear K = induction::edmd(slow, {kt::symbol("x"), kt::symbol("y"), kt::power(kt::symbol("x"), kt::number(2))});
    std::vector<double> rates;
    for (const auto& rate : K.rates) rates.push_back(rate.real());
    std::sort(rates.begin(), rates.end());
    s.check("on the dictionary {x, y, x^2}: Koopman's eigenvalues lambda, 2 mu and mu", std::fabs(rates[0] + 1) + std::fabs(rates[1] + 0.1) + std::fabs(rates[2] + 0.05), 1e-6);
    s.check("the dictionary is invariant: a step leaves nothing unexplained", K.residual, 1e-9);
  }

  // ---- a forecaster ----
  {
    const Trajectory coarse = run(sys("lorenz"), {1, 1, 20}, 140, 0.005, 5, 10);  // sampled at 0.025
    Trajectory train;
    train.names = coarse.names;
    const std::size_t cut = 4000;
    train.t.assign(coarse.t.begin(), coarse.t.begin() + cut), train.x.assign(coarse.x.begin(), coarse.x.begin() + cut);
    const induction::Forecaster f = induction::forecaster(train, 2, 1, 1e-6);
    s.require("a forecaster of two delays on three components has 28 features", f.features == 28);
    s.check("one step of Lorenz, against the size of a step", f.training_error, 5e-3);
    // With too small a ridge one step is fitted better and the free run is worse: the delays make a step of
    // several levels, whose spurious modes the ridge damps.
    const induction::Forecaster loose = induction::forecaster(train, 2, 1, 1e-12);
    std::vector<State> past(coarse.x.begin(), coarse.x.begin() + cut), future(coarse.x.begin() + cut, coarse.x.begin() + cut + 400);
    s.require("a ridge too small fits one step better and runs for a shorter time", loose.training_error < f.training_error && induction::horizon(loose.forecast(past, 400), future, 0.3) < induction::horizon(f.forecast(past, 400), future, 0.3));
    std::vector<State> history(coarse.x.begin(), coarse.x.begin() + cut), truth(coarse.x.begin() + cut, coarse.x.begin() + cut + 400);
    const std::size_t held = induction::horizon(f.forecast(history, 400), truth, 0.3);
    s.require("a free run stays near the truth for more than three Lyapunov times (132 samples)", held > 132, std::to_string(held) + " samples, " + std::to_string(static_cast<double>(held) * 0.025 * 0.9056) + " Lyapunov times");
  }

  // ---- stochastic ----
  {
    const systems::Stochastic d = systems::stochastic(sys("ou"));
    systems::StochasticOptions so;
    so.t1 = 3000, so.step = 0.002, so.every = 5, so.seed = 21;
    const auto path = systems::simulate(d, {0}, so);
    Trajectory record;
    record.names = path.names, record.t = path.t, record.x = path.x;
    const induction::Fit f = induction::fit_stochastic(sys("ou_template"), record);
    s.check("the rate of an Ornstein-Uhlenbeck process from a path", f.numbers.at("theta") / 1.5 - 1, 0.06);
    s.check("its noise", f.numbers.at("s") / 0.6 - 1, 0.02);
    const auto column = record.column(0);
    const induction::Moments m = induction::kramers_moyal(column, 0.01, 8, -0.4, 0.4);
    double slope_error = 0, level_error = 0;
    for (std::size_t b = 0; b < 8; b++) {
      if (std::fabs(m.x[b]) > 0.1) slope_error = std::max(slope_error, std::fabs(m.drift[b] / m.x[b] + 1.5));
      level_error = std::max(level_error, std::fabs(m.diffusion[b] / 0.36 - 1));
    }
    s.check("in bins: the drift is -theta x", slope_error, 0.4);
    s.check("in bins: the second coefficient is s^2", level_error, 0.06);
  }

  // ---- matching ----
  {
    const Trajectory osc = run(sys("vdp"), {0.5, 0}, 40, 0.01);
    const auto ranked = induction::match({sys("spring_template"), induction::as_template(sys("vdp")), sys("duffing_template"), induction::as_template(sys("lorenz"))}, osc);
    s.require("of three laws proposed for van der Pol's record, van der Pol's ranks first", ranked.size() == 4 && ranked[0].name == "vdp" && ranked[0].fitted, ranked[0].name);
    s.check("its parameter", ranked[0].fit.numbers.at("mu") - 1.2, 2e-3);
    s.check("its error of one step on the part kept aside", ranked[0].one_step, 1e-3);
    s.require("it follows the part kept aside to its end (1200 samples), and the others do not", ranked[0].horizon == 1200 && ranked[1].horizon < 1200 && ranked[2].horizon < 1200,
              std::to_string(ranked[0].horizon) + ", " + ranked[1].name + " " + std::to_string(ranked[1].horizon) + ", " + ranked[2].name + " " + std::to_string(ranked[2].horizon));
    s.require("a law whose states the record does not hold is listed last, with the reason", ranked[3].name == "lorenz" && !ranked[3].fitted && ranked[3].refused.find("'y'") != std::string::npos, ranked[3].refused);
  }

  return s.done();
}
