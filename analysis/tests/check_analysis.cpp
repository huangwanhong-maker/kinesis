// Checks of kinesis/analysis.
//
// Invariants
// - Lorenz: the exponents sum to −(σ + 1 + β) (the divergence, exact whatever the path); the largest is 0.9056
//   and Kaplan and Yorke's dimension 2.062 (Sprott 2003, Chaos and Time-Series Analysis, Table A.1 values
//   0.9056, 0, −14.5723), within the error of a finite time.
// - Hénon's map: the exponents sum to ln 0.3 exactly; the largest is 0.4192 (Sprott 2003). The logistic map at
//   r = 4: ln 2.
// - Lorenz's equilibria: the origin and C± = (±√(β(ρ − 1)), ±√(β(ρ − 1)), ρ − 1); at ρ = 28 all are saddles, at
//   ρ = 20 the pair are stable foci; the pair lose stability between ρ = 24.7 and 24.8 (Hopf at
//   σ(σ + β + 3)/(σ − β − 1) = 24.7368).
// - A section of the harmonic oscillator at q = 0 is crossed every 2π/ω.
// - Diagrams: the logistic map has 1, 2 and 4 points at r = 2.8, 3.2 and 3.5, the two at 3.2 in closed form;
//   van der Pol's oscillator at μ = 1 has maxima of 2.0086 (its limit cycle's amplitude).
// - Basins of x' = x − x³, y' = −y: two attractors (±1, 0), each with half of a symmetric grid; basin stability
//   ½ within its standard error.
// Criticality
// - An Ornstein–Uhlenbeck process sampled at Δt: variance σ²/(2λ) and lag-1 autocorrelation e^{−λΔt}.
// - A process whose rate falls towards zero: variance and autocorrelation trend upward, significantly; a
//   stationary one does not.
// - DFA: ½ for white noise, 3/2 for a random walk.
// - A power law sampled by inversion is recovered within three standard errors.
// - A branching process seen through a subsample: the slope at lag 1 is biased, the decay of the slopes gives m.
// - Avalanches of the critical Galton–Watson process with Poisson offspring: sizes with an exponent near 3/2.
// Landscapes
// - x' = x − x³: U = x⁴/4 − x²/2 at every grid point, wells at ±1 and a barrier at 0.
// - A gradient flow of two dimensions is recognised, and its potential difference is exact; Lorenz is not one.
// - dX = (X − X³) dt + σ dW: the exact stationary landscape is 2U/σ² up to a constant; the landscape from a long
//   path agrees with it; Kramers's time against the mean of simulated escapes.
// - A linear stochastic system with rotation: −ln p = |x|²/σ² from a long path.
// Operators
// - Ulam's matrix of the logistic map at r = 4: the stationary density is 1/(π√(x(1 − x))); of the doubling map
//   on a Markov partition of 2⁶ cells: the uniform density, reached exactly in six steps.

#include <kinesis/analysis/criticality.hpp>
#include <kinesis/analysis/invariants.hpp>
#include <kinesis/analysis/statistics.hpp>
#include <kinesis/analysis/landscape.hpp>
#include <kinesis/analysis/operators.hpp>
#include <kinesis/language/reader.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>

using namespace kinesis;
using systems::State;
static const double pi = std::numbers::pi;

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
system henon {
  state x, y : real
  param a = 1.4, b = 0.3
  change {
    x' = 1 - a * x^2 + y
    y' = b * x
  }
  structure { time n in natural }
}
system logistic {
  state x : real
  param r = 4
  change { x' = r * x * (1 - x) }
  structure { time n in natural }
}
system doubling {
  state x : real
  change { x' = mod(2 * x, 1) }
  structure { time n in natural }
}
system oscillator {
  state q : real
  param omega = 3
  change { q'' = -omega^2 * q }
}
system vdp {
  state x : real
  param mu = 1
  change { x'' = mu * (1 - x^2) * x' - x }
}
system wells {
  state x, y : real
  change {
    x' = x - x^3
    y' = -y
  }
}
system well {
  state x : real
  change { x' = x - x^3 }
}
system coupled_wells {
  state x, y : real
  change {
    x' = x - x^3 - 0.5 * (x - y)
    y' = y - y^3 - 0.5 * (y - x)
  }
}
system noisy_well {
  state X : real
  param s = 0.45
  noise W : wiener
  change { dX = (X - X^3) * dt + s * dW }
}
system spiral {
  state X, Y : real
  param s = 0.8
  noise W : wiener
  noise V : wiener
  change {
    dX = (-X + 2 * Y) * dt + s * dW
    dY = (-2 * X - Y) * dt + s * dV
  }
}
)";

int main() {
  calx::verification::Suite s;
  const language::Script script = language::read(text);
  const auto flow = [&](const char* n) { return systems::dynamics(*script.find(n)); };

  // ---- invariants ----
  {
    const auto d = flow("lorenz");
    physis::dynamics::LyapunovOptions o;
    o.step = 0.005, o.transient = 20, o.time = 400;
    const auto e = analysis::lyapunov(d, {1, 1, 20}, o);
    s.check("Lorenz: the exponents sum to -(sigma + 1 + beta)", e.sum + 41.0 / 3, 5e-5);
    s.check("Lorenz: the largest exponent against 0.9056", e.exponents[0] - 0.9056, 0.03);
    s.check("Lorenz: the middle exponent against 0", e.exponents[1], 0.02);
    s.check("Lorenz: Kaplan and Yorke's dimension against 2.062", e.dimension - 2.062, 0.01);
    s.require("with no running record no error is given", e.error.empty());
    o.record_every = 40;   // 200 records of the 400 units of time
    const auto r = analysis::lyapunov(d, {1, 1, 20}, o);
    s.require("with a running record, the largest exponent is within four of its errors of 0.9056, and the error is of the size of a finite time's",
              r.error.size() == 3 && r.error[0] > 0.002 && r.error[0] < 0.05 && std::fabs(r.exponents[0] - 0.9056) < 4 * r.error[0],
              "0.9056 - " + std::to_string(r.exponents[0]) + " against an error of " + std::to_string(r.error.empty() ? 0.0 : r.error[0]));
  }
  {
    // A record whose rate over block b of ten is b: the means over blocks 2 to 10 have the standard deviation
    // sqrt(7.5), and the error of their mean is a third of it.
    std::vector<double> times;
    std::vector<std::vector<double>> running;
    double total = 0;
    for (int k = 1; k <= 40; k++) {
      total += 0.5 * static_cast<double>((k - 1) / 4 + 1);
      times.push_back(0.5 * k), running.push_back({total / (0.5 * k), 3.0});
    }
    const auto error = analysis::block_errors(times, running);
    s.check("the error over blocks of a record whose block rates are 1 to 10: sqrt(7.5) / 3", error.at(0) - std::sqrt(7.5) / 3, 1e-13);
    s.check("and of a rate that is constant: none", error.at(1), 1e-13);
    s.require("too short a record gives none", analysis::block_errors({1, 2, 3}, {{1}, {1}, {1}}).empty());
  }
  {
    physis::dynamics::LyapunovOptions o;
    o.transient = 1000, o.time = 200000, o.renormalise_every = 1;
    const auto e = analysis::lyapunov(flow("henon"), {0.1, 0.1}, o);
    s.check("Henon: the exponents sum to ln 0.3", e.sum - std::log(0.3), 1e-9);
    s.check("Henon: the largest exponent against 0.4192", e.exponents[0] - 0.4192, 0.005);
    const auto l = analysis::lyapunov(flow("logistic"), {0.3141}, o);
    s.check("the logistic map at r = 4: ln 2", l.exponents[0] - std::log(2.0), 0.01);
  }
  {
    auto d = flow("lorenz");
    const analysis::Box box{{-30, -30, -5}, {30, 30, 60}};
    const auto found = analysis::equilibria(d, box);
    const double c = std::sqrt(8.0 / 3 * 27);
    bool placed = found.size() == 3;
    if (placed) placed = std::fabs(found[0].x[0] + c) < 1e-9 && std::fabs(found[1].x[0]) < 1e-9 && std::fabs(found[2].x[0] - c) < 1e-9 && std::fabs(found[2].x[2] - 27) < 1e-9;
    s.require("Lorenz at rho = 28: the origin and the pair at (+-sqrt(beta (rho - 1)), ., rho - 1), all saddles", placed && found[0].kind == "saddle" && found[1].kind == "saddle" && found[2].kind == "saddle",
              std::to_string(found.size()) + " equilibria" + (found.empty() ? "" : "; the first " + found[0].kind));
    d.set_param("rho", 20);
    const auto calm = analysis::equilibria(d, box);
    s.require("at rho = 20 the pair are stable foci and the origin a saddle", calm.size() == 3 && calm[0].kind == "stable focus" && calm[2].kind == "stable focus" && calm[1].kind == "saddle");
    const auto leading = [&](double rho) {
      d.set_param("rho", rho);
      State x{std::sqrt(8.0 / 3 * (rho - 1)), std::sqrt(8.0 / 3 * (rho - 1)), rho - 1};
      double worst = -1e300;
      for (const auto& l : analysis::equilibrium_at(d, x).eigenvalues) worst = std::max(worst, l.real());
      return worst;
    };
    s.require("the pair lose stability between rho = 24.7 and 24.8 (Hopf at 24.7368)", leading(24.7) < 0 && leading(24.8) > 0, std::to_string(leading(24.7)) + ", " + std::to_string(leading(24.8)));
    double lo = 24.7, hi = 24.8;
    for (int k = 0; k < 40; k++) (leading(0.5 * (lo + hi)) < 0 ? lo : hi) = 0.5 * (lo + hi);
    s.check("the crossing against sigma (sigma + beta + 3) / (sigma - beta - 1)", lo - 10 * (10 + 8.0 / 3 + 3) / (10 - 8.0 / 3 - 1), 1e-9);
  }
  {
    const auto d = flow("oscillator");
    systems::Options o;
    o.t1 = 20, o.step = 0.005;
    const auto crossings = analysis::section(d, {1, 0}, terms::symbol("q"), true, o);
    double worst = 0;
    for (std::size_t k = 1; k < crossings.size(); k++) worst = std::max(worst, std::fabs(crossings[k].t - crossings[k - 1].t - 2 * pi / 3));
    s.require("the oscillator crosses q = 0 upward nine times in 20", crossings.size() == 9, std::to_string(crossings.size()));
    s.check("the times between crossings against 2 pi / omega", worst, 1e-8);
  }
  {
    analysis::DiagramOptions o;
    o.transient = 2000, o.keep = 64;
    const auto diagram = analysis::bifurcation_diagram(flow("logistic"), "r", {2.8, 3.2, 3.5}, {0.3}, o);
    const auto distinct = [](std::vector<double> v) {
      std::sort(v.begin(), v.end());
      std::vector<double> out;
      for (double x : v)
        if (out.empty() || x - out.back() > 1e-6) out.push_back(x);
      return out;
    };
    const auto p1 = distinct(diagram.points[0]), p2 = distinct(diagram.points[1]), p4 = distinct(diagram.points[2]);
    s.require("the logistic map: 1, 2 and 4 points at r = 2.8, 3.2 and 3.5", p1.size() == 1 && p2.size() == 2 && p4.size() == 4, std::to_string(p1.size()) + ", " + std::to_string(p2.size()) + ", " + std::to_string(p4.size()));
    const double root = std::sqrt((3.2 - 3) * (3.2 + 1));
    s.check("the fixed point 1 - 1/r and the cycle (r + 1 +- sqrt((r - 3)(r + 1))) / (2 r)", std::fabs(p1[0] - (1 - 1 / 2.8)) + std::fabs(p2[0] - (4.2 - root) / 6.4) + std::fabs(p2[1] - (4.2 + root) / 6.4), 1e-9);
    analysis::DiagramOptions f;
    f.transient = 100, f.time = 60, f.step = 0.005, f.keep = 5;
    const auto vdp = analysis::bifurcation_diagram(flow("vdp"), "mu", {1.0}, {0.5, 0}, f);
    double worst = 0;
    for (double m : vdp.points[0]) worst = std::max(worst, std::fabs(m - 2.0086));
    s.require("van der Pol at mu = 1: maxima are kept", !vdp.points[0].empty());
    s.check("its maxima against the amplitude 2.0086 of the limit cycle", worst, 1e-3);
  }
  {
    const auto d = flow("wells");
    std::vector<double> xs, ys;
    for (int i = 0; i < 20; i++) xs.push_back(-1.9 + 0.2 * i), ys.push_back(-1.9 + 0.2 * i);
    analysis::BasinOptions o;
    o.time = 40;
    const auto b = analysis::basins(d, {0, 0}, xs, ys, o);
    s.require("the basins of x' = x - x^3, y' = -y: two attractors, each with half of a symmetric grid", b.attractors.size() == 2 && std::fabs(b.share[0] - 0.5) < 1e-12 && std::fabs(std::fabs(b.attractors[0][0]) - 1) < 1e-6);
    s.require("both ends are at rest and are stable nodes", b.kind.size() == 2 && b.kind[0] == "stable node" && b.kind[1] == "stable node" && b.speed[0] < 1e-9);
    // A grid with a column on x = 0, the boundary between the basins: its paths end at the saddle at the origin,
    // which is named a saddle and comes last, with the share of one column of twenty-one.
    std::vector<double> odd;
    for (int i = 0; i <= 20; i++) odd.push_back(-2.0 + 0.2 * i);
    odd[10] = 0;
    const auto c = analysis::basins(d, {0, 0}, odd, ys, o);
    s.require("a column of the grid on the boundary ends at the saddle, which is named one and has that column's share",
              c.attractors.size() == 3 && c.kind[2] == "saddle" && std::fabs(c.share[2] - 1.0 / 21) < 1e-12 && std::fabs(c.share[0] - 10.0 / 21) < 1e-12 && std::fabs(c.attractors[2][0]) < 1e-9,
              std::to_string(c.attractors.size()) + " ends");
    const auto bs = analysis::basin_stability(d, {1, 0}, analysis::Box{{-2, -2}, {2, 2}}, 400, o);
    s.check("the basin stability of (1, 0) in a symmetric box against 1/2, in standard errors", (bs.share - 0.5) / bs.error, 3);
    s.check("its standard error sqrt(S (1 - S) / N)", bs.error - std::sqrt(bs.share * (1 - bs.share) / 400), 1e-15);
  }

  // ---- averages and scaling ----
  {
    // The mean of a series and its error from blocks: ten blocks whose means are 1 to 10 have the mean 5.5 and the
    // standard error sqrt(55 / 6) / sqrt(10) (the variance of 1..10 with nine degrees of freedom is 55/6).
    std::vector<double> series;
    for (int b = 1; b <= 10; b++)
      for (int k = 0; k < 7; k++) series.push_back(b + (k % 2 ? 0.25 : -0.25) * (k < 6));
    const auto bm = analysis::block_mean(series, 10);
    s.check("the mean of a series by blocks", bm.mean - 5.5, 1e-13);
    s.check("its standard error: that of the ten block means", bm.error - std::sqrt(55.0 / 6 / 10), 1e-13);
    s.require("a series shorter than its blocks has a mean and no error", analysis::block_mean(std::vector<double>{1, 2, 3}, 10).error == 0);
    // A power law recovered exactly, and one with scatter with its error.
    std::vector<double> x, y;
    for (int i = 0; i < 9; i++) x.push_back(std::pow(2.0, i)), y.push_back(3 * std::pow(x.back(), -5.0 / 3));
    const auto law = analysis::scaling(x, y);
    s.check("a power law through exact points: its exponent, its prefactor and no error", std::fabs(law.exponent + 5.0 / 3) + std::fabs(law.prefactor - 3) + law.error, 1e-12);
    for (std::size_t i = 0; i < y.size(); i++) y[i] *= (i % 2 ? 1.1 : 1 / 1.1);
    const auto rough = analysis::scaling(x, y);
    s.require("with scatter the exponent is within three of its errors of the truth, and the error is not zero", rough.error > 0.005 && std::fabs(rough.exponent + 5.0 / 3) < 3 * rough.error,
              std::to_string(rough.exponent) + " +- " + std::to_string(rough.error));
    // What a system does on average: an oscillator over whole periods, and a chaotic map with a known density.
    const language::Script known = language::read(R"(
system turning {
  state x, y : real
  change {
    x' = -y
    y' = x
  }
  difference { observe square = x^2, both = x^2 + y^2 }
}
system quadratic {
  state x : real
  change { x' = 4 * x * (1 - x) }
  structure { time n in natural }
  difference { observe x }
}
system unseen {
  state x : real
  change { x' = -x }
}
)");
    analysis::StatisticsOptions o;
    o.step = 0.001, o.transient = 0, o.time = 20 * std::numbers::pi, o.every = 1;
    const auto circle = analysis::averages(systems::dynamics(*known.find("turning")), {1, 0}, o);
    s.check("the mean of cos^2 over ten turns is 1/2, and its second moment over its mean squared is 3/2", std::fabs(circle.of[0].mean - 0.5) + std::fabs(circle.of[0].second - 1.5), 1e-4);
    s.check("an observation that does not fluctuate has the second moment 1 and no error", std::fabs(circle.of[1].second - 1) + circle.of[1].error, 1e-10);
    analysis::StatisticsOptions m;
    m.transient = 1000, m.time = 400000, m.every = 1;
    const auto chaos = analysis::averages(systems::dynamics(*known.find("quadratic")), {0.3141}, m);
    s.check("the logistic map at r = 4 has the mean 1/2 under its invariant density, in its standard errors", (chaos.of[0].mean - 0.5) / chaos.of[0].error, 3);
    s.check("and the second moment 3/8 over 1/4", chaos.of[0].second - 1.5, 0.01);
    bool refused = false;
    try {
      analysis::averages(systems::dynamics(*known.find("unseen")), {1}, o);
    } catch (const kernel::Refusal&) {
      refused = true;
    }
    s.require("a system that declares nothing to observe is refused its averages", refused);
  }

  // ---- criticality ----
  {
    // An Ornstein-Uhlenbeck process sampled exactly: x_{k+1} = rho x_k + sqrt(v (1 - rho^2)) xi.
    const double lambda = 0.5, sigma = 0.7, dt = 0.4, rho = std::exp(-lambda * dt), v = sigma * sigma / (2 * lambda);
    calx::random::Stream r(4, 1, 0, 0);
    std::vector<double> x(200000);
    double state = 0;
    for (double& value : x) state = rho * state + std::sqrt(v * (1 - rho * rho)) * r.normal(), value = state;
    s.check("an Ornstein-Uhlenbeck process: the variance sigma^2 / (2 lambda)", tyche::descriptive::variance(x) / v - 1, 0.03);
    s.check("its autocorrelation at the sampling lag, exp(-lambda dt)", tyche::time_series::acf(x, 1)[1] - rho, 0.01);
  }
  {
    // The rate falls from 1 to 0.02 over the record: a system brought towards a fold.
    calx::random::Stream r(8, 2, 0, 0);
    const std::size_t n = 700;
    std::vector<double> slowing(n), steady(n);
    double a = 0, b = 0;
    for (std::size_t k = 0; k < n; k++) {
      const double lambda = 1 - 0.98 * static_cast<double>(k) / static_cast<double>(n), rho = std::exp(-lambda * 0.5);
      a = rho * a + 0.1 * r.normal(), slowing[k] = a;
      b = std::exp(-0.5) * b + 0.1 * r.normal(), steady[k] = b;
    }
    analysis::WarningOptions o;
    o.surrogates = 150;
    const auto w = analysis::early_warnings(slowing, o), calm = analysis::early_warnings(steady, o);
    s.require("a process that slows: variance and autocorrelation trend upward", w.tau_variance > 0.6 && w.tau_autocorrelation > 0.6, std::to_string(w.tau_variance) + ", " + std::to_string(w.tau_autocorrelation));
    s.require("significantly, against red-noise surrogates", w.p_variance < 0.05 && w.p_autocorrelation < 0.05, std::to_string(w.p_variance) + ", " + std::to_string(w.p_autocorrelation));
    s.require("a stationary process is not flagged: neither trend is significant at 1 in 100", calm.p_variance > 0.01 && calm.p_autocorrelation > 0.01, std::to_string(calm.p_variance) + ", " + std::to_string(calm.p_autocorrelation));
    s.require("the windows end from the half of the record to its end", w.at.size() == n - n / 2 + 1 && w.at.back() == static_cast<double>(n - 1));
  }
  {
    calx::random::Stream r(15, 3, 0, 0);
    std::vector<double> white(40000), walk(40000);
    double position = 0;
    for (std::size_t k = 0; k < white.size(); k++) white[k] = r.normal(), position += white[k], walk[k] = position;
    s.check("DFA of white noise: 1/2", analysis::dfa(white).exponent - 0.5, 0.04);
    s.check("DFA of a random walk: 3/2", analysis::dfa(walk).exponent - 1.5, 0.06);
  }
  {
    calx::random::Stream r(23, 4, 0, 0);
    std::vector<double> x(20000);
    for (double& v : x) v = 2.0 * std::pow(1 - r.uniform(), -1 / (2.5 - 1));  // p(x) ~ x^-2.5 above 2
    const auto p = analysis::power_law(x, 2.0);
    s.check("a power law of exponent 2.5, in its standard errors", (p.exponent - 2.5) / p.error, 3);
    s.check("its standard error (alpha - 1) / sqrt(n)", p.error - (p.exponent - 1) / std::sqrt(20000.0), 1e-15);
    s.check("the distance between the sample and the fitted law", p.distance, 0.02);
  }
  {
    // A branching process with drive: A' ~ Poisson(m A + h), seen through a subsample of 1 unit in 10.
    const double m = 0.95, h = 5;
    calx::random::Stream r(31, 5, 0, 0);
    std::vector<double> full(120000), seen(120000);
    double a = h / (1 - m);
    for (std::size_t k = 0; k < full.size(); k++) {
      a = static_cast<double>(r.poisson(m * a + h));
      full[k] = a, seen[k] = static_cast<double>(r.binomial(static_cast<std::int64_t>(a), 0.1));
    }
    const auto whole = analysis::branching_ratio(full), part = analysis::branching_ratio(seen);
    s.check("a branching process seen whole: the slope at lag 1 is m", whole.conventional - m, 0.01);
    s.require("seen through a subsample the slope at lag 1 is biased low", part.conventional < m - 0.1, std::to_string(part.conventional));
    s.check("and the decay of the slopes gives m", part.ratio - m, 0.015);
  }
  {
    // The critical Galton-Watson process with Poisson offspring, each avalanche begun by one unit.
    calx::random::Stream r(42, 6, 0, 0);
    std::vector<double> sizes;
    for (int k = 0; k < 40000; k++) {
      double active = 1, size = 0;
      while (active > 0 && size < 1e6) size += active, active = static_cast<double>(r.poisson(active));
      sizes.push_back(size);
    }
    const auto p = analysis::power_law(sizes, 10, true);
    s.check("avalanches of the critical branching process: sizes with the exponent 3/2", p.exponent - 1.5, 0.05);
    const auto a = analysis::avalanches(std::vector<double>{0, 2, 3, 0, 0, 1, 0, 4, 4, 4});
    s.require("avalanches are runs above the threshold: sizes 5, 1, 12 and durations 2, 1, 3", a.sizes == std::vector<double>{5, 1, 12} && a.durations == std::vector<double>{2, 1, 3});
  }

  // ---- landscapes ----
  {
    const auto d = flow("well");
    const auto p = analysis::potential(d, -1.6, 1.6, 320);
    double worst = 0;
    const auto U = [](double x) { return x * x * x * x / 4 - x * x / 2; };
    for (std::size_t i = 0; i < p.x.size(); i++) worst = std::max(worst, std::fabs(p.u[i] - (U(p.x[i]) - U(-1.6))));
    s.check("x' = x - x^3: the potential x^4/4 - x^2/2 on the grid", worst, 1e-12);
    // A parabola through three samples places the extremum of a quartic to the square of the grid's step.
    s.require("its wells at +-1 and its barrier at 0, to the square of the step", p.minima.size() == 2 && p.maxima.size() == 1 && std::fabs(p.minima[0] + 1) < 1e-4 && std::fabs(p.minima[1] - 1) < 1e-4 && std::fabs(p.maxima[0]) < 1e-9,
              std::to_string(p.minima.size()) + " wells, " + std::to_string(p.maxima.size()) + " barrier(s)");
    const auto c = flow("coupled_wells");
    s.check("a coupled pair of wells is a gradient: its Jacobian is symmetric", analysis::asymmetry(c, {-2, -2}, {2, 2}), 1e-14);
    const auto V = [](double x, double y) { return x * x * x * x / 4 - x * x / 2 + y * y * y * y / 4 - y * y / 2 + 0.25 * (x - y) * (x - y); };
    s.check("its potential between two points", analysis::potential_between(c, {0.3, -0.4}, {1.1, 0.9}) - (V(1.1, 0.9) - V(0.3, -0.4)), 1e-10);
    s.require("Lorenz is not a gradient", analysis::asymmetry(flow("lorenz"), {-10, -10, 0}, {10, 10, 30}) > 0.1);
  }
  {
    const systems::Stochastic d = systems::stochastic(*script.find("noisy_well"));
    const double sigma = 0.45;
    const auto exact = analysis::stationary_landscape(d, -1.6, 1.6, 320);
    const auto U = [](double x) { return x * x * x * x / 4 - x * x / 2; };
    double worst = 0;
    for (std::size_t i = 0; i < exact.x.size(); i++) worst = std::max(worst, std::fabs(exact.u[i] - 2 * (U(exact.x[i]) + 0.25) / (sigma * sigma)));
    s.check("dX = (X - X^3) dt + s dW: the stationary landscape is 2 U / s^2 up to a constant", worst, 1e-9);
    analysis::DensityOptions o;
    o.x0 = -1.6, o.x1 = 1.6, o.nx = 32, o.ny = 1, o.time = 40000, o.step = 0.005, o.paths = 8;
    const auto L = analysis::density_landscape(d, {1.0}, o);
    // The depth of the barrier at 0 below... above the wells, from the path: against 2 (1/4) / s^2.
    const double barrier = 0.5 * (L.u[15] + L.u[16]), well = std::min(*std::min_element(L.u.begin(), L.u.begin() + 12), *std::min_element(L.u.begin() + 20, L.u.end()));
    s.check("the landscape from a long path: the barrier above the wells against 1 / (2 s^2)", (barrier - well) / (0.5 / (sigma * sigma)) - 1, 0.1);
    const double kramers = analysis::kramers_time(flow("well"), 1.0, 0.0, sigma);
    s.check("Kramers's time: 2 pi / sqrt(2) exp(1 / (2 s^2))", kramers - 2 * pi / std::sqrt(2.0) * std::exp(0.5 / (sigma * sigma)), 1e-6);
    // The mean time of 300 escapes from the well at 1 to beyond -0.5.
    double total = 0;
    const double h = 0.004;
    for (std::uint32_t path = 0; path < 300; path++) {
      calx::random::Stream r(77, 9, path, 0);
      State x{1.0};
      double t = 0;
      while (x[0] > -0.5 && t < 5000) {
        const double dw = std::sqrt(h) * r.normal();
        d.step(t, x, h, &dw, false);
        t += h;
      }
      total += t;
    }
    s.check("the mean time of 300 simulated escapes against Kramers's, whose error at this noise is of some tenths", total / 300 / kramers - 1, 0.35);
  }
  {
    const systems::Stochastic d = systems::stochastic(*script.find("spiral"));
    analysis::DensityOptions o;
    o.x0 = -1.5, o.x1 = 1.5, o.y0 = -1.5, o.y1 = 1.5, o.nx = 15, o.ny = 15, o.time = 30000, o.step = 0.005, o.paths = 8;
    const auto L = analysis::density_landscape(d, {0, 0}, o);
    // -ln p = |x|^2 / s^2 (the covariance is s^2 / 2 times the identity, whatever the rotation).
    double worst = 0;
    for (std::size_t iy = 0; iy < 15; iy++)
      for (std::size_t ix = 0; ix < 15; ix++) {
        const double r2 = L.xs[ix] * L.xs[ix] + L.ys[iy] * L.ys[iy];
        if (r2 > 1.0) continue;
        worst = std::max(worst, std::fabs(L.u[iy * 15 + ix] - L.u[7 * 15 + 7] - r2 / 0.64));
      }
    s.check("a linear system with rotation: -ln p = |x|^2 / s^2 within the unit disc, from a long path", worst, 0.15);
  }

  // ---- operators ----
  {
    const auto T = analysis::transfer(flow("logistic"), 0, 1, 200, 100);
    double error = 0;
    for (std::size_t i = 10; i < 190; i++) error += std::fabs(T.stationary[i] - 1 / (pi * std::sqrt(T.centres[i] * (1 - T.centres[i])))) / 180;
    s.check("Ulam's matrix of the logistic map at r = 4: the density 1 / (pi sqrt(x (1 - x))), mean error away from the ends", error, 0.06);
    double rows = 0;
    for (std::size_t i = 0; i < 200; i++) {
      double sum = 0;
      for (std::size_t j = 0; j < 200; j++) sum += T.matrix(i, j);
      rows = std::max(rows, std::fabs(sum - 1));
    }
    s.check("its rows sum to 1: no point leaves the interval", rows + T.lost, 1e-12);
    const auto D = analysis::transfer(flow("doubling"), 0, 1, 64, 64);
    double flat = 0;
    for (double v : D.stationary) flat = std::max(flat, std::fabs(v - 1));
    s.check("of the doubling map: the uniform density", flat, 1e-9);
    // 64 = 2^6 cells are a Markov partition of the doubling map: after six steps every cell is spread evenly over
    // all of them, so P^6 is the matrix whose every entry is 1/64, and P has the eigenvalue 1 and no other but 0.
    calx::linalg::Matrix power = D.matrix;
    for (int k = 1; k < 6; k++) power = power * D.matrix;
    double spread = 0;
    for (double v : power.a) spread = std::max(spread, std::fabs(v - 1.0 / 64));
    s.check("on its Markov partition of 64 cells, six steps spread every cell evenly: P^6 is uniform", spread, 1e-15);
    s.check("its leading eigenvalue is 1", std::abs(D.eigenvalues[0]) - 1, 1e-9);
  }

  return s.done();
}
