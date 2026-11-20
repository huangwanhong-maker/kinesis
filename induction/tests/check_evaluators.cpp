// Checks of kinesis/induction/evaluators: the true law of an AR(1) record beats the random walk and the mean and is not beaten
// by the fitted autoregression; a law that is the fitted autoregression is not shown better than it, which is said; a law with the wrong coefficient is beaten by every baseline, which is
// said; over laws x' = a x + u the claim that one order of two inputs is better holds, fails, is arguable or undecided as the
// laws' a say (its difference a³(1 − a) exactly), is undecided when the order is immaterial (a = 1), and undecided when an
// input or a state leaves the record's support, named at its step; a law with a threshold makes order matter where the linear
// integrator does not; bagged sparse fits keep the true terms of a record nearly always and the others seldom.

#include <kinesis/induction/evaluators.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace kinesis::induction;

namespace {
Member linear(double a, double w = 1) {
  return {"a = " + std::to_string(a), w, [a](const std::vector<double>& x, const std::vector<double>& u) { return std::vector<double>{a * x[0] + u[0]}; }};
}
bool says(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }
}  // namespace

int main() {
  calx::verification::Suite s;
  {
    calx::random::Stream r(3, 1, 0, 0);
    std::vector<double> y(700);
    double x = 0;
    for (int t = 0; t < 100; t++) x = 0.6 * x + r.normal();
    for (double& v : y) x = 0.6 * x + r.normal(), v = x;
    const auto truth = held_to_baselines(y, 100, [](std::span<const double> h) { return 0.6 * h.back(); });
    s.require("the true law beats the random walk and the historical mean",
              truth.against[0].verdict == "the law forecasts better" && truth.against[1].verdict == "the law forecasts better");
    s.require("... and is not beaten by the fitted autoregression (which must estimate what the true law knows)", !truth.beaten);
    const auto same = held_to_baselines(y, 100, [](std::span<const double> h) { return tyche::time_series::ar1(h); });
    s.require("a law that is the fitted autoregression is not shown better than it, and the verdict says so",
              same.against[2].verdict == "no difference shown" && same.verdict == "not shown better than the first-order autoregression", same.verdict);
    const auto wrong = held_to_baselines(y, 100, [](std::span<const double> h) { return -0.3 * h.back(); });
    s.require("a law with the wrong coefficient is beaten, by every baseline, and the verdict says so",
              wrong.beaten && wrong.verdict == "beaten by the random walk, historical mean, first-order autoregression", wrong.verdict);
  }
  {
    const std::vector<double> x0{0};
    const std::vector<std::vector<double>> first{{1}, {2}, {0}, {0}, {0}}, second{{2}, {1}, {0}, {0}, {0}};
    const auto final = [](const Run& r) { return r.states.back()[0]; };
    const auto h = rank({linear(0.8), linear(0.9)}, x0, first, second, final);
    s.require("laws with a < 1 all prefer the larger input later: the claim holds", h.standing == "holds" && h.share_first == 1);
    s.check("... its difference under a = 0.8 is a³(1 − a) exactly", h.difference[0] - 0.8 * 0.8 * 0.8 * 0.2, 1e-14);
    s.require("laws with a > 1 all say the opposite: it fails", rank({linear(1.1), linear(1.2)}, x0, first, second, final).standing == "fails");
    const auto ar = rank({linear(0.8), linear(0.9), linear(1.1)}, x0, first, second, final);
    s.require("two of three equal laws for it: arguable, with the range of the difference across them",
              ar.standing == "arguable" && ar.lower < 0 && ar.upper > 0 && std::fabs(ar.share_first - 2.0 / 3) < 1e-12);
    s.require("one law each way: undecided, the laws disagree", rank({linear(0.9), linear(1.1)}, x0, first, second, final).standing == "undecided");
    const auto im = rank({linear(1.0)}, x0, first, second, final);
    s.require("under the linear integrator the order is immaterial, and the claim is undecided for that reason", im.immaterial && im.standing == "undecided" && says(im.reason, "immaterial"));
    // Each unit of input adds one while the state is below 1 and a half once it is not: the effect depends on the state.
    const Member threshold{"a threshold", 1, [](const std::vector<double>& x, const std::vector<double>& u) {
                             return std::vector<double>{x[0] + u[0] * (x[0] < 1 ? 1.0 : 0.5)};
                           }};
    const auto th = rank({linear(1.0), threshold}, x0, first, second, final);
    s.require("a law with a threshold makes the order matter where the linear integrator does not", th.difference[0] == 0 && std::fabs(th.difference[1] + 0.5) < 1e-12);
    // The record held inputs in [0, 2] and states in [−5, 5].
    const Support inputs = support_of({{0}, {2}, {1}}, {"u"}), states = support_of({{-5}, {5}}, {"x"});
    const auto in = rank({linear(0.8), linear(0.9)}, x0, first, second, final, &states, &inputs);
    s.require("within the record's support nothing is flagged", in.flags.empty() && in.standing == "holds");
    const std::vector<std::vector<double>> bold{{1}, {3}, {0}, {0}, {0}};
    const auto out = rank({linear(0.8), linear(0.9)}, x0, bold, second, final, &states, &inputs);
    s.require("an input beyond what the record varied: undecided, named at its step and variable",
              out.standing == "undecided" && says(out.reason, "out of the record's support") && says(out.flags.front(), "step 1: input u = 3"), out.flags.empty() ? "" : out.flags.front());
    const std::vector<std::vector<double>> drive{{4}, {4}, {0}, {0}, {0}};
    const auto st = rank({linear(0.9)}, x0, drive, second, final, &states, nullptr);
    s.require("a state driven beyond what the record held is flagged too", !st.flags.empty() && says(st.flags.front(), "state x = "));
  }
  {
    // y = 1.5 x − 0.8 x³ + noise, on a library of monomials of degree 0 … 5.
    calx::random::Stream r(17, 4, 0, 0);
    const std::size_t N = 200, K = 6;
    Rectangular theta(N, K);
    std::vector<double> y(N);
    for (std::size_t i = 0; i < N; i++) {
      const double x = r.uniform(-2, 2);
      for (std::size_t k = 0; k < K; k++) theta(i, k) = std::pow(x, static_cast<double>(k));
      y[i] = 1.5 * x - 0.8 * x * x * x + 0.05 * r.normal();
    }
    SparseOptions o;
    o.threshold = 0.1;
    const auto B = bagged(theta, y, 100, 23, o);
    s.require("bagged sparse fits keep the true terms (x and x³) nearly always", B.inclusion[1] >= 0.95 && B.inclusion[3] >= 0.95);
    double others = 0;
    for (std::size_t k : {0u, 2u, 4u, 5u}) others = std::max(others, B.inclusion[k]);
    s.check("... and every other term seldom (the largest inclusion among them)", others, 0.3);
    s.check("... their median coefficients near the truth", std::fabs(B.median[1] - 1.5) + std::fabs(B.median[3] + 0.8), 0.05);
  }
  return s.done();
}
