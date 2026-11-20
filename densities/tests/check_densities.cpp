// Checks of kinesis/densities: a generator's rows sum to zero and its implicit flow keeps mass and positivity; under a
// linear drift the mean moves exactly as the implicit step says it must; the stationary density of Bouchaud and Mézard's
// law (drift J(1 − w), diffusion σ²w²) tends to their inverse gamma law at first order in the step; an intervention and a
// flow composed in both orders differ in their mean by δ(1 − q^N), exactly; two translations of all the mass commute while
// two transfers to those below a line do not; carrying mass keeps mass and mean; the Gini coefficient of two equal masses
// at 0 and 1 is one half. The value equation: with one regime its control is the permanent-income line κ(a + z/r), more
// nearly as the grid is refined; with two regimes and r < ρ the lower regime holds at the limit, with a mass there, the
// higher grows there, the regimes' masses are the switching chain's exactly, and the stationary density solves Aᵀg = 0.

#include <kinesis/densities/generator.hpp>
#include <kinesis/densities/value.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace kinesis::densities;

namespace {

// Bouchaud and Mézard's law (2000, "Wealth condensation in a simple model of economy", Physica A 282, 536–545), mean-field:
// dW = J(1 − W) dt + √2 σ W dB, stationary P(w) = β^α/Γ(α) w^{−α−1} e^{−β/w} with α = 1 + J/σ², β = J/σ².
Generator mean_field(std::size_t cells, double J, double s2, double top) {
  Grid g{0, top, cells};
  std::vector<std::vector<double>> mu(1, std::vector<double>(cells)), D(1, std::vector<double>(cells));
  for (std::size_t i = 0; i < cells; i++) {
    const double w = g.at(i);
    mu[0][i] = J * (1 - w), D[0][i] = s2 * w * w;
  }
  return generator(g, mu, D);
}

double inverse_gamma(double w, double a, double b) { return w > 0 ? std::exp(a * std::log(b) - std::lgamma(a) - (a + 1) * std::log(w) - b / w) : 0; }

double l1_to_law(std::size_t cells) {
  const double J = 1, s2 = 1.0 / 3;
  const auto G = mean_field(cells, J, s2, 40);
  const auto g = stationary(G);
  const double h = G.grid.step();
  double e = 0;
  for (std::size_t i = 0; i < cells; i++) e += std::fabs(g[i] - inverse_gamma(G.grid.at(i), 1 + J / s2, J / s2) * h);
  return e;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  {
    const auto G = mean_field(801, 1, 1.0 / 3, 40);
    s.check("a generator's rows sum to zero", row_defect(G), 1e-9);
    std::vector<double> g(G.states(), 0);
    g[G.index(40, 0)] = 1;  // all the mass at w = 2
    const auto F = flow(G, 0.01);
    double m = mean(G, g);
    bool positive = true;
    double worst = 0;
    for (int k = 0; k < 300; k++) {
      g = F.step(g);
      for (double v : g) positive = positive && v >= 0;
      // (m' − 1)(1 + J dt) − (m − 1) = −dt g'_top D_top / h: the drift's part exactly, and the diffusion's rate that
      // reflection removed at the top node.
      const std::size_t top = G.grid.cells - 1;
      const double Dtop = (1.0 / 3) * G.grid.at(top) * G.grid.at(top);
      const double lhs = (mean(G, g) - 1) * (1 + 0.01) - (m - 1), rhs = -0.01 * g[G.index(top, 0)] * Dtop / G.grid.step();
      worst = std::max(worst, std::fabs(lhs - rhs));
      m = mean(G, g);
    }
    s.check("the implicit flow keeps the mass over 300 steps", total(g) - 1, 1e-12);
    s.require("... and keeps it non-negative", positive);
    s.check("under a linear drift the mean moves exactly as the implicit step says, the rate reflection removed counted", worst, 1e-13);
  }
  {
    const double coarse = l1_to_law(801), fine = l1_to_law(1601);
    s.check("Bouchaud and Mézard's stationary law: the inverse gamma (α = 4, β = 3), in L1 on a step of 0.025", fine, 0.02);
    s.check("... at first order in the step: halving it halves the error", coarse / fine - 2, 0.25);
  }
  {
    // An intervention and a flow in both orders: everyone given δ, then N steps, against N steps, then δ.
    const auto G = mean_field(801, 1, 1.0 / 3, 40);
    auto g = stationary(G);
    const auto F = flow(G, 0.05);
    const double d = 0.4;
    const std::size_t N = 20;
    const auto give = [&](const std::vector<double>& x) { return carry(G, x, [d](double w, std::size_t) { return w + d; }).density; };
    const auto a = F.run(give(g), N), b = give(F.run(g, N));
    const double q = 1 / (1 + 0.05);
    s.check("the order of a transfer and a flow: the means differ by δ(1 − q^N), exactly", (mean(G, b) - mean(G, a)) - d * (1 - std::pow(q, static_cast<double>(N))), 2e-5);
    const auto t1 = [&](const std::vector<double>& x) { return carry(G, x, [](double w, std::size_t) { return w + 0.3; }).density; };
    const auto t2 = [&](const std::vector<double>& x) { return carry(G, x, [](double w, std::size_t) { return w + 0.5; }).density; };
    const auto u = t2(t1(g)), v = t1(t2(g));
    double diff = 0;
    for (std::size_t i = 0; i < u.size(); i++) diff = std::max(diff, std::fabs(u[i] - v[i]));
    s.check("two translations of all the mass commute", diff, 1e-12);
    const auto lo = [&](const std::vector<double>& x) { return transfer(G, x, 0.5, 0.8).density; };
    const auto hi = [&](const std::vector<double>& x) { return transfer(G, x, 0.3, 1.2).density; };
    const double first = below(G, F.run(hi(F.run(lo(g), N)), N), 1.0), second = below(G, F.run(lo(F.run(hi(g), N)), N), 1.0);
    s.require("two transfers to those below a line, with flows between, do not: the mass below 1 differs", std::fabs(first - second) > 1e-4,
              std::to_string(first) + " against " + std::to_string(second));
    const auto moved = transfer(G, g, 0.7, 1.5);
    s.check("carrying mass within the grid keeps the mass", total(moved.density) - total(g), 1e-13);
    s.check("... and adds exactly the transfer to the mean", mean(G, moved.density) - mean(G, g) - 0.7 * below(G, g, 1.5), 1e-6);
  }
  {
    Grid g{0, 1, 3};
    const auto G = generator(g, {{0, 0, 0}}, {{1, 1, 1}});
    s.check("the Gini coefficient of equal masses at 0 and 1 is one half", gini(G, {0.5, 0, 0.5}) - 0.5, 1e-15);
  }
  // ---- the value equation -----------------------------------------------------------------------------------------------
  {
    auto error = [](std::size_t cells) {
      Stock m;
      m.rho = 0.05, m.r = 0.03, m.gamma = 2, m.inflow = {0.1};
      m.grid = {-0.1 / 0.03 + 0.01, 10, cells};
      const auto S = solve(m);
      const double kappa = (m.rho - m.r * (1 - m.gamma)) / m.gamma;
      double worst = 0;
      for (std::size_t i = 0; i < cells; i++) {
        const double a = m.grid.at(i);
        if (a < -2 || a > 8) continue;
        const double c = kappa * (a + 0.1 / 0.03);
        worst = std::max(worst, std::fabs(S.control[0][i] - c) / c);
      }
      return std::make_pair(worst, S.converged);
    };
    const auto c1 = error(500), c2 = error(1000);
    s.require("one regime: the value equation converges", c1.second && c2.second);
    s.check("... and its control is the permanent-income line κ(a + z/r) within 1%", c2.first, 0.01);
    s.require("... more nearly as the grid is refined", c2.first < c1.first, std::to_string(c1.first) + " then " + std::to_string(c2.first));
  }
  {
    Stock m;
    m.rho = 0.05, m.r = 0.03, m.gamma = 2, m.inflow = {0.1, 0.2};
    m.switching = {{0, 0.5}, {1.5, 0}};
    m.grid = {-0.15, 5, 600};
    const auto S = solve(m);
    s.require("two regimes: the value equation converges", S.converged);
    s.require("in the lower regime the stock holds at the limit and falls just above it", S.drift[0][0] == 0 && S.drift[0][1] < 0);
    s.require("in the higher regime it grows at the limit", S.drift[1][0] > 0);
    s.check("its generator's rows sum to zero", row_defect(S.G), 1e-9);
    const auto g = stationary(S.G);
    s.check("the regimes' masses are the switching chain's: 1.5/(0.5 + 1.5) in the lower", type_mass(S.G, g, 0) - 0.75, 1e-9);
    // A mass at the limit stays as the grid is refined; the density beside it, singular as (a − a_min)^(−1/2), puts on
    // the next node a part that shrinks with the step.
    Stock fine = m;
    fine.grid.cells = 2400;
    const auto Sf = solve(fine);
    const auto gf = stationary(Sf.G);
    const double at = g[S.G.index(0, 0)], at_f = gf[Sf.G.index(0, 0)], next = g[S.G.index(1, 0)], next_f = gf[Sf.G.index(1, 0)];
    // The first node holds the mass M and the singular density's part, which scales as √h: quartering the step halves it,
    // so M = 2 g_fine − g_coarse (Richardson's extrapolation in √h).
    const double M = 2 * at_f - at, Mh = 2 * gf[Sf.G.index(0, 1)] - g[S.G.index(0, 1)];
    s.require("the lower regime has a mass at the limit, by extrapolation in √h, while the node beside it empties",
              M > 0.03 && next_f < next / 1.5, std::to_string(M) + "; next " + std::to_string(next) + " → " + std::to_string(next_f));
    s.check("... and the higher regime none", Mh, 0.005);
    const auto r = calx::linalg::multiply(calx::linalg::transpose(S.G.A), g);
    double worst = 0;
    for (double v : r) worst = std::max(worst, std::fabs(v));
    s.check("the stationary density solves Aᵀg = 0: the forward equation is the value equation's generator transposed", worst, 1e-10);
  }
  return s.done();
}
