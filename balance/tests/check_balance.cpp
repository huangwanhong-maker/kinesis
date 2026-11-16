// Checks of kinesis/balance: balance laws on graphs, run and learned from records of the field. Every record here is
// simulated from a law the check declares, or has a closed form, so that what is to be found is known.
//
// - The grid and its operators: cells and faces counted with a hole in the mask; the Laplacian zero on a constant and
//   exactly 4 on x² + y² inside (the five-point stencil is exact on quadratics); diffusion, crowding and drift each
//   conserving the content to rounding; the rate's moving part equal to minus the divergence of the edges' fluxes;
//   a drift up a potential that rises eastward carrying the current east and nowhere north.
// - Runs: one cell with logistic growth against Verhulst's solution, and with Gompertz's against his; a diffusion
//   keeping the content; a drift outrunning the diffusion made to clip, the content rising by exactly what was clipped;
//   Fisher's front on a fine strip running at 2√(rD) less Bramson's delay.
// - Learning: the integral form recovering a law of growth, spread and drift from a record taken often; from a record
//   taken five times, the integral form biased and shooting recovering the law; thresholded least squares on the whole
//   library keeping exactly the terms of Fisher's law; the standard errors covering the truth under noise.
// - Scores: carrying the last field forward has no skill, the truth has all; a rate of its own per cell reproducing
//   exponential growth exactly.
// - Modes: from ten snapshots of eighty cells moved by a linear law of four modes, the eigenvalues, a forecast four
//   steps past the record and a residual of rounding; times not equally spaced refused.

#include <kinesis/balance/balance.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <complex>
#include <set>
#include <string>
#include <vector>

using namespace kinesis::balance;
using calx::verification::relative;

static Grid full_grid(std::size_t rows, std::size_t cols, double h) { return grid(rows, cols, {}, h); }

int main() {
  calx::verification::Suite s;

  // ---- the grid and its operators ----
  {
    std::vector<char> mask(5 * 6, 1);
    mask[2 * 6 + 3] = 0;  // a hole
    const Grid g = grid(5, 6, mask, 1.0);
    // 29 cells; faces: 5 rows × 5 + 4 × 6 = 49 on the full grid, less the hole's four.
    s.require("a 5 × 6 grid with a hole has 29 cells and 45 faces", g.graph.size() == 29 && g.graph.edges.size() == 45,
              std::to_string(g.graph.size()) + " cells, " + std::to_string(g.graph.edges.size()) + " faces");
    const Grid q = full_grid(9, 9, 2.0);
    std::vector<double> one(q.graph.size(), 7.0), quad(q.graph.size());
    for (std::size_t v = 0; v < quad.size(); v++) {
      const double x = 2.0 * q.col[v], y = 2.0 * q.row[v];
      quad[v] = x * x + y * y;
    }
    double worst_c = 0, worst_q = 0;
    const auto lc = laplacian(q.graph, one), lq = laplacian(q.graph, quad);
    for (std::size_t v = 0; v < quad.size(); v++) {
      worst_c = std::max(worst_c, std::abs(lc[v]));
      if (q.row[v] > 0 && q.row[v] < 8 && q.col[v] > 0 && q.col[v] < 8) worst_q = std::max(worst_q, std::abs(lq[v] - 4.0));
    }
    s.check("the Laplacian of a constant is zero", worst_c, 1e-14);
    s.check("the Laplacian of x² + y² is 4 at every inner cell, on a grid of spacing 2", worst_q, 1e-12);

    calx::random::Stream rng(20261116, 1, 0, 0);
    std::vector<double> u(g.graph.size()), phi(g.graph.size());
    for (std::size_t v = 0; v < u.size(); v++) u[v] = 10 * rng.uniform(), phi[v] = rng.uniform();
    const double scale = content(g.graph, u);
    s.check("diffusion keeps the content", content(g.graph, laplacian(g.graph, u)) / scale, 1e-14);
    s.check("crowding keeps the content", content(g.graph, crowding(g.graph, u)) / (scale * 10), 1e-14);
    s.check("a drift keeps the content", content(g.graph, drift(g.graph, u, phi)) / scale, 1e-14);

    Law law{{{Kind::diffusion}, {Kind::crowding}, {Kind::drift, 0}}, {0.7, 0.02, 1.3}};
    const auto r = rate(g.graph, law, u, {phi});
    const auto f = flux(g.graph, law, u, {phi});
    std::vector<double> div(u.size(), 0.0);
    for (std::size_t e = 0; e < f.size(); e++) div[g.graph.edges[e].a] += f[e], div[g.graph.edges[e].b] -= f[e];
    double worst = 0;
    for (std::size_t v = 0; v < u.size(); v++) worst = std::max(worst, std::abs(r[v] + div[v] / g.graph.volume[v]));
    s.check("the rate of the moving terms is minus the divergence of the edges' fluxes", worst, 1e-12);

    const Grid e = full_grid(4, 4, 1.0);
    std::vector<double> flat(e.graph.size(), 3.0), east(e.graph.size());
    for (std::size_t v = 0; v < east.size(); v++) east[v] = static_cast<double>(e.col[v]);
    const Currents c = currents(e, flux(e.graph, Law{{{Kind::drift, 0}}, {2.0}}, flat, {east}));
    double least_east = 1e300, most_north = 0;
    for (std::size_t v = 0; v < flat.size(); v++) least_east = std::min(least_east, c.east[v]), most_north = std::max(most_north, std::abs(c.north[v]));
    s.require("a drift up a potential rising eastward carries the current east, and nowhere north", least_east > 0 && most_north == 0,
              "least east " + std::to_string(least_east) + ", most north " + std::to_string(most_north));
  }

  // ---- runs ----
  {
    Graph one{{1.0}, {}};
    const double r = 0.3, K = 100, u0 = 5;
    std::vector<double> times;
    for (int k = 1; k <= 30; k++) times.push_back(k);
    RunOptions fine;
    fine.most_step = 0.01;
    const Run lg = simulate(one, Law{{{Kind::linear}, {Kind::quadratic}}, {r, -r / K}}, {u0}, 0, times, {}, {}, fine);
    double worst = 0;
    for (std::size_t k = 0; k < times.size(); k++) worst = std::max(worst, relative(lg.u[k][0], K / (1 + (K / u0 - 1) * std::exp(-r * times[k]))));
    s.check("one cell's logistic growth follows Verhulst's solution", worst, 1e-9);
    const Run gz = simulate(one, Law{{{Kind::linear}, {Kind::gompertz}}, {r * std::log(K), -r}}, {u0}, 0, times, {}, {}, fine);
    worst = 0;
    for (std::size_t k = 0; k < times.size(); k++) worst = std::max(worst, relative(gz.u[k][0], K * std::exp(std::log(u0 / K) * std::exp(-r * times[k]))));
    s.check("and Gompertz's growth his", worst, 1e-9);

    const Grid g = full_grid(12, 12, 1.0);
    std::vector<double> bump(g.graph.size(), 0.0);
    for (std::size_t v = 0; v < bump.size(); v++)
      if (g.row[v] >= 5 && g.row[v] <= 6 && g.col[v] >= 5 && g.col[v] <= 6) bump[v] = 100;
    const Run d = simulate(g.graph, Law{{{Kind::diffusion}}, {0.8}}, bump, 0, {5, 10, 20});
    s.check("a diffusion keeps the content over its run", relative(content(g.graph, d.u.back()), content(g.graph, bump)), 1e-12);
    std::vector<double> west(g.graph.size());
    for (std::size_t v = 0; v < west.size(); v++) west[v] = -3.0 * g.col[v];
    const Run dr = simulate(g.graph, Law{{{Kind::drift, 0}, {Kind::diffusion}}, {1.0, 0.01}}, bump, 0, {2}, {west});
    s.require("a drift that outruns its diffusion is clipped, and says so", dr.clipped > 0, std::to_string(dr.clipped));
    s.check("the content then rises by exactly what was clipped", relative(content(g.graph, dr.u.back()), content(g.graph, bump) + dr.clipped), 1e-12);

    // Fisher's front on a strip of 2400 cells of a quarter: D = 1, r = 1, K = 1.
    const double h = 0.25;
    const Grid strip = full_grid(1, 2400, h);
    std::vector<double> start(strip.graph.size(), 0.0);
    for (std::size_t v = 0; v < 40; v++) start[v] = 1.0;
    const Run fr = simulate(strip.graph, Law{{{Kind::diffusion}, {Kind::linear}, {Kind::quadratic}}, {1.0, 1.0, -1.0}}, start, 0, {100, 150});
    auto front = [&](const std::vector<double>& u) {
      for (std::size_t v = 1; v < u.size(); v++)
        if (u[v - 1] >= 0.5 && u[v] < 0.5) return h * (static_cast<double>(v - 1) + (u[v - 1] - 0.5) / (u[v - 1] - u[v]));
      return -1.0;
    };
    const double speed = (front(fr.u[1]) - front(fr.u[0])) / 50;
    // 2√(rD) = 2, less Bramson's 3/(2t) ln-correction between t = 100 and 150, plus the lattice's h²/12: about 1.993.
    s.check("Fisher's front runs at 2√(rD), less Bramson's delay", speed - (2 - 1.5 * std::log(1.5) / 50 + h * h / 12), 0.01);
  }

  // ---- learning ----
  {
    const Grid g = full_grid(16, 16, 1.0);
    std::vector<double> u0(g.graph.size()), centre(g.graph.size());
    for (std::size_t v = 0; v < u0.size(); v++) {
      const double x = g.col[v] - 7.5, y = g.row[v] - 7.5;
      u0[v] = 5 + 60 * std::exp(-(x * x + y * y) / 6);
      centre[v] = -std::sqrt(x * x + y * y) / 4;
    }
    const std::vector<Term> terms{{Kind::linear}, {Kind::quadratic}, {Kind::diffusion}, {Kind::drift, 0}};
    const Law truth{terms, {0.2, -0.002, 0.5, 0.3}};
    // Often: every 0.05 for 2 time units.
    std::vector<double> often;
    for (int k = 1; k <= 40; k++) often.push_back(0.05 * k);
    const Run ro = simulate(g.graph, truth, u0, 0, often, {centre});
    Record rec{{0.0}, {u0}};
    for (std::size_t k = 0; k < often.size(); k++) rec.t.push_back(often[k]), rec.u.push_back(ro.u[k]);
    const Estimate ie = integral(g.graph, terms, rec, {centre});
    double worst = 0;
    for (std::size_t j = 0; j < terms.size(); j++) worst = std::max(worst, relative(ie.law.coefficients[j], truth.coefficients[j]));
    s.check("from a record taken often, the integral form recovers growth, its ceiling, spread and drift", worst, 2e-3);

    // Five times, five apart.
    const Run rc = simulate(g.graph, truth, u0, 0, {5, 10, 15, 20}, {centre});
    Record coarse{{0.0}, {u0}};
    for (std::size_t k = 0; k < 4; k++) coarse.t.push_back(5.0 * (k + 1)), coarse.u.push_back(rc.u[k]);
    const Estimate ce = integral(g.graph, terms, coarse, {centre});
    double bias = 0;
    for (std::size_t j = 0; j < terms.size(); j++) bias = std::max(bias, relative(ce.law.coefficients[j], truth.coefficients[j]));
    s.require("from a record taken five times, five apart, the integral form is biased by the trapezoid", bias > 0.01, std::to_string(bias));
    const Estimate se = shooting(g.graph, ce, coarse, {centre});
    worst = 0;
    for (std::size_t j = 0; j < terms.size(); j++) worst = std::max(worst, relative(se.law.coefficients[j], truth.coefficients[j]));
    s.check("and shooting from it recovers the law", worst, 1e-4);
    s.require("shooting converged", se.converged, std::to_string(se.iterations) + " iterations");

    // Sparse: Fisher's law found among the whole library.
    const Law fisher{{{Kind::linear}, {Kind::quadratic}, {Kind::diffusion}}, {0.2, -0.002, 0.5}};
    const Run rf = simulate(g.graph, fisher, u0, 0, often);
    Record rfr{{0.0}, {u0}};
    for (std::size_t k = 0; k < often.size(); k++) rfr.t.push_back(often[k]), rfr.u.push_back(rf.u[k]);
    const std::vector<Term> library{{Kind::constant}, {Kind::linear}, {Kind::quadratic}, {Kind::gompertz}, {Kind::diffusion}, {Kind::crowding}, {Kind::drift, 0}};
    LearnOptions so;
    so.threshold = 0.05;
    const Estimate sp = integral(g.graph, library, rfr, {centre}, so);
    std::set<std::string> kept;
    for (std::size_t j = 0; j < library.size(); j++)
      if (sp.law.coefficients[j] != 0) kept.insert(name(library[j].kind));
    std::string said;
    for (const auto& k : kept) said += k + " ";
    s.require("thresholded on the whole library, the integral form keeps exactly Fisher's three terms", kept == std::set<std::string>{"linear", "quadratic", "diffusion"}, said);

    // Noise: the truth within three standard errors of shooting's estimates.
    calx::random::Stream rng(7, 2, 0, 0);
    Record noisy = coarse;
    for (std::size_t k = 1; k < noisy.u.size(); k++)
      for (double& v : noisy.u[k]) v *= 1 + 0.02 * rng.normal();
    const Estimate ne = shooting(g.graph, integral(g.graph, terms, noisy, {centre}), noisy, {centre});
    double z = 0;
    for (std::size_t j = 0; j < terms.size(); j++) z = std::max(z, std::abs(ne.law.coefficients[j] - truth.coefficients[j]) / ne.standard_errors[j]);
    s.check("with 2% noise on the record, the truth lies within three of shooting's standard errors", z, 3.0, "standard errors");

    // Well-posedness: a record that sharpens, as a backward diffusion would, gives a negative diffusivity, which is set aside.
    std::vector<double> sharp = u0;
    const auto lap = laplacian(g.graph, u0);
    for (std::size_t v = 0; v < sharp.size(); v++) sharp[v] -= 0.1 * lap[v];
    const Record back{{0.0, 1.0}, {u0, sharp}};
    const Estimate wp = integral(g.graph, {{Kind::diffusion}, {Kind::linear}}, back);
    LearnOptions loose;
    loose.well_posed = false;
    const Estimate ill = integral(g.graph, {{Kind::diffusion}, {Kind::linear}}, back, {}, loose);
    s.require("a record that sharpens gives a negative diffusivity when allowed, and none when the law must be well posed, said so",
              ill.law.coefficients[0] < 0 && wp.law.coefficients[0] == 0 && wp.dropped.size() == 1, std::to_string(ill.law.coefficients[0]));

    // Scores.
    std::vector<std::vector<double>> scaled{coarse.u[3], coarse.u[4]};
    for (auto& u : scaled) for (double& v : u) v *= 1.2;
    const Score sc = score(g.graph, scaled, coarse, 2);
    s.check("a forecast with the right pattern and a fifth too many has all the pattern's skill", 1 - sc.pattern_skill, 1e-12);
    s.require("and less than all the skill", sc.skill < 0.99, std::to_string(sc.skill));
    const Score ps = score(g.graph, {coarse.u[2], coarse.u[2]}, coarse, 2);
    s.check("carrying the last field forward has no skill", ps.skill, 1e-15);
    const Score ts = score(g.graph, {coarse.u[3], coarse.u[4]}, coarse, 2);
    s.check("the truth has all of it", 1 - ts.skill, 1e-15);
    Record ex;
    for (int k = 0; k < 4; k++) {
      ex.t.push_back(10.0 * k);
      std::vector<double> u(5);
      for (std::size_t i = 0; i < 5; i++) u[i] = (1 + i) * std::exp(0.01 * (i + 1) * 10.0 * k);
      ex.u.push_back(u);
    }
    const OwnRates own = own_rates(first(ex, 3), 0.0);
    double ow = 0;
    const auto at = own.at(30.0);
    for (std::size_t i = 0; i < 5; i++) ow = std::max(ow, relative(at[i], ex.u[3][i]));
    s.check("a rate of its own per cell reproduces exponential growth", ow, 1e-12);
  }

  // ---- modes ----
  {
    const std::size_t n = 80;
    calx::random::Stream rng(11, 3, 0, 0);
    std::vector<std::vector<double>> v(4, std::vector<double>(n));
    for (auto& vec : v)
      for (double& x : vec) x = rng.normal();
    const double a1 = 1.05, a2 = 0.8, rho = 0.9, om = 0.3;
    Record r;
    auto truth = [&](double k) {
      std::vector<double> u(n);
      for (std::size_t i = 0; i < n; i++)
        u[i] = 3 * std::pow(a1, k) * v[0][i] + 2 * std::pow(a2, k) * v[1][i] + std::pow(rho, k) * (std::cos(om * k) * v[2][i] - std::sin(om * k) * v[3][i]);
      return u;
    };
    for (int k = 0; k < 10; k++) r.t.push_back(2000 + 5.0 * k), r.u.push_back(truth(k));
    const Modes m = modes(r);
    std::vector<std::complex<double>> want{{a1, 0}, {a2, 0}, std::polar(rho, om), std::polar(rho, -om)};
    double worst = 0;
    for (const auto& w : want) {
      double best = 1e300;
      for (const Mode& md : m.modes) best = std::min(best, std::abs(md.eigenvalue - w));
      worst = std::max(worst, best);
    }
    s.require("ten snapshots of a linear law of four modes have rank four", m.rank == 4, std::to_string(m.rank));
    s.check("their eigenvalues: 1.05, 0.8 and 0.9 e^{±0.3i}", worst, 1e-9);
    const auto f = m.at(2000 + 5.0 * 14), t = truth(14);
    double fe = 0, tn = 0;
    for (std::size_t i = 0; i < n; i++) fe += (f[i] - t[i]) * (f[i] - t[i]), tn += t[i] * t[i];
    s.check("the modes forecast four steps past the record", std::sqrt(fe / tn), 1e-8);
    s.check("and carry each snapshot to the next", m.residual, 1e-10);
    Record uneven = r;
    uneven.t[3] += 1;
    bool refused = false;
    try { (void)modes(uneven); } catch (const Refusal&) { refused = true; }
    s.require("times not equally spaced are refused", refused);
  }

  return s.done();
}
