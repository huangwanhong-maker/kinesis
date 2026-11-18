// Checks of kinesis/scales: systems across scales on a hierarchy of units.
//
// - The hierarchy: children, groups of siblings, ancestors; refusals of a parent that does not exist, a coarse unit with nothing
//   beneath it, costs of the wrong size.
// - Carrying up: sums across three levels; a ratio index as a value on summed quantities is the ratio of sums (and the mean of
//   ratios is not), equal to the weighted mean of the ratios by the denominator.
// - Runs held to closed forms: exponential growth in every block, its city total Σ x₀ e^{rt}; a logistic total allocated to
//   blocks by fixed shares, each block's path w_i/Σw times the logistic and the sum equal to the total to rounding; a map with a
//   period of two steps, a^k after k periods; diffusion between two cells, the difference falling as e^{−2Dt} and the content kept.
// - Pace: two laws at two steps against the exact exponential of their linear system, the error halving with the slow step (first
//   order, Jacobi), and of the order of Runge–Kutta's when they share one step.
// - A fast quantity at rest: x' = −y with y = x/2 as an equilibrium gives e^{−t/2} exactly; written as a fast flow ε y' = x/2 − y
//   the error falls in proportion to ε (Tikhonov).
// - Spatial interaction: the inflow against the formula by hand, its total conserved, and β = 0 sharing by attractiveness alone;
//   Harris and Wilson's centres W' = ε(D − κW) reaching D_j = κW_j with Σ W = Σ O/κ.
// - Logit shares summing to one in each group; z-scores among siblings; a registered reducer (a weighted median) and an ordinal level read through steps.
// - What is missing: initial values of two blocks, a parameter, a driver, each refused with its name, level and count; an
//   equilibrium needs no initial value.
// - Ensembles: a growth rate declared on [0.1, 0.3] gives x(1) with median e^{0.2} and 5% and 95% quantiles e^{0.11} and e^{0.29};
//   a quantity's units drawn at one quantile; Morris's screening ranks the wider factor first.
// - Fitting: a logistic's r and K recovered from its exact record; judged on years held aside against persistence; noisy records
//   within four standard errors; two unknowns that enter only as a product reported as not identified.
// - Perfect aggregation: a lumpable linear law and its coarse matrix, and a perturbed one that is not.
// - Refusals with their reasons: a law naming another level, derived quantities in a circle, a state with no law, a step that is
//   not a whole multiple, a reducer not given.

#include <kinesis/scales/aggregation.hpp>
#include <kinesis/scales/ensemble.hpp>
#include <kinesis/scales/fit.hpp>
#include <kinesis/scales/hierarchy.hpp>
#include <kinesis/scales/model.hpp>
#include <kinesis/scales/run.hpp>

#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace kinesis::scales;

static bool refused(const std::function<void()>& f, const std::string& part, std::string* said = nullptr) {
  try {
    f();
  } catch (const Refusal& r) {
    if (said) *said = r.what();
    return std::string(r.what()).find(part) != std::string::npos;
  }
  return false;
}

static Hierarchy three() {
  Level b{"block", {"b0", "b1", "b2", "b3", "b4", "b5"}, {0, 0, 0, 1, 1, 1}, {}, {}, {}};
  Level a{"area", {"a0", "a1"}, {0, 0}, {}, {}, {}};
  Level c{"city", {"c"}, {}, {}, {}, {}};
  return Hierarchy({b, a, c});
}

static Hierarchy flat(std::size_t n, const std::string& level = "unit") {
  Level u{level, {}, {}, {}, {}, {}};
  for (std::size_t i = 0; i < n; i++) u.units.push_back("u" + std::to_string(i));
  return Hierarchy({u});
}

static std::vector<Value> vals(const std::vector<double>& x, Source s = Source::observed) {
  std::vector<Value> v;
  for (double e : x) v.push_back(Value{e, s, 0, missing_value, missing_value});
  return v;
}

static Quantity q(const std::string& n, const std::string& level, Role r, Kind k = Kind::extensive) {
  Quantity x;
  x.name = n, x.level = level, x.role = r, x.kind = k;
  return x;
}
static Derived d(const std::string& n, const std::string& level, Op op, const std::string& of, const std::string& weight = "") {
  Derived x;
  x.name = n, x.level = level, x.op = op, x.of = of, x.weight = weight;
  return x;
}
static Law law(const std::string& target, Form f, const std::string& e, double step = 0) {
  Law x;
  x.target = target, x.form = f, x.expression = e, x.step = step;
  return x;
}
static Parameter par(const std::string& n, double v, bool fit = false, double lo = -infinity, double hi = infinity) {
  Parameter p;
  p.name = n, p.value = v, p.fit = fit, p.lower = lo, p.upper = hi;
  return p;
}

int main() {
  calx::verification::Suite s;

  // ---- the hierarchy ----
  {
    const Hierarchy h = three();
    s.require("the children of area a0 are blocks 0, 1, 2", h.children(1, 0) == std::vector<std::size_t>{0, 1, 2});
    s.require("block b4 is in the second group of siblings", h.group(0, 4) == 1);
    s.require("block b4's ancestor in the city is the city", h.ancestor(0, 4, 2) == 0);
    s.require("the top level is one group of all its units", h.groups(2).size() == 1 && h.groups(2)[0].size() == 1);
    s.require("a parent that does not exist is refused",
              refused([] { Hierarchy({Level{"b", {"x", "y"}, {0, 3}, {}, {}, {}}, Level{"a", {"p"}, {}, {}, {}, {}}}); }, "names a parent"));
    s.require("a coarse unit with nothing beneath it is refused",
              refused([] { Hierarchy({Level{"b", {"x"}, {0}, {}, {}, {}}, Level{"a", {"p", "q"}, {}, {}, {}, {}}}); }, "has no units of"));
    s.require("costs of the wrong size are refused",
              refused([] { Hierarchy({Level{"b", {"x", "y"}, {}, {}, {}, {{"walk", {0, 1, 1}}}}}); }, "n × n"));
  }

  // ---- carrying up: sums and ratio indexes ----
  {
    Model m;
    m.quantities = {q("vacant", "block", Role::fixed), q("total", "block", Role::fixed)};
    m.derived = {d("vacant_a", "area", Op::sum, "vacant"), d("total_a", "area", Op::sum, "total"),
                 d("rate_a", "area", Op::value, "vacant_a / total_a"), d("mean_of_rates", "area", Op::mean, "vacant / total"),
                 d("weighted_rates", "area", Op::wmean, "vacant / total", "total"),
                 d("vacant_c", "city", Op::sum, "vacant_a"), d("total_c", "city", Op::sum, "total_a"),
                 d("rate_c", "city", Op::value, "vacant_c / total_c")};
    const Compiled c = compile(three(), m);
    Inputs in;
    in.values["vacant"] = vals({1, 2, 3, 0, 5, 5});
    in.values["total"] = vals({10, 10, 40, 20, 10, 20});
    Options o;
    o.t0 = 0, o.t1 = 1;
    const Run r = simulate(c, in, o);
    s.check("the area's ratio index is the ratio of its sums, 6/60", r.at("rate_a", 0)[0] - 0.1, 1e-15);
    s.check("... and in the second area 10/50", r.at("rate_a", 0)[1] - 0.2, 1e-15);
    s.check("the mean of the blocks' ratios is another number, 0.125", r.at("mean_of_rates", 0)[0] - 0.125, 1e-15);
    s.check("the ratios weighted by their denominators are the ratio of sums again", r.at("weighted_rates", 0)[1] - 0.2, 1e-15);
    s.check("the city's index from its own sums, 16/110", r.at("rate_c", 0)[0] - 16.0 / 110.0, 1e-15);
  }

  // ---- exponential growth in every block, the city's total ----
  {
    Model m;
    m.quantities = {q("x", "block", Role::state)};
    m.derived = {d("xa", "area", Op::sum, "x"), d("xc", "city", Op::sum, "xa")};
    m.laws = {law("x", Form::flow, "r * x")};
    m.parameters = {par("r", 0.3)};
    m.step = 0.01;
    const Compiled c = compile(three(), m);
    Inputs in;
    in.values["x"] = vals({1, 2, 3, 4, 5, 6});
    Options o;
    o.t0 = 0, o.t1 = 2;
    o.report = {2};
    const Run r = simulate(c, in, o);
    s.check("each block grows as x0 e^{rt}: block 3 at t = 2", r.at("x", 0)[3] - 4 * std::exp(0.6), 1e-10);
    s.check("the city's total is 21 e^{0.6}", r.at("xc", 0)[0] - 21 * std::exp(0.6), 1e-9);
  }

  // ---- a logistic total allocated by fixed shares ----
  {
    Level b{"block", {"b0", "b1", "b2", "b3", "b4", "b5"}, {0, 0, 0, 0, 0, 0}, {}, {}, {}};
    Level c0{"city", {"c"}, {}, {}, {}, {}};
    const Hierarchy h({b, c0});
    Model m;
    m.quantities = {q("P", "city", Role::state), q("w", "block", Role::fixed, Kind::other)};
    m.derived = {d("s", "block", Op::share, "w"), d("Pc", "block", Op::parent, "P"), d("x", "block", Op::value, "s * Pc"),
                 d("sum_x", "city", Op::sum, "x")};
    m.laws = {law("P", Form::flow, "r * P * (1 - P / K)")};
    m.parameters = {par("r", 0.5), par("K", 100)};
    m.step = 0.01;
    const Compiled c = compile(h, m);
    Inputs in;
    in.values["P"] = vals({10});
    in.values["w"] = vals({1, 2, 3, 4, 5, 5});
    Options o;
    o.t0 = 0, o.t1 = 6;
    o.report = {6};
    const Run r = simulate(c, in, o);
    const double logistic = 100 / (1 + 9 * std::exp(-0.5 * 6));
    s.check("the total follows the logistic", r.at("P", 0)[0] - logistic, 1e-8);
    s.check("block 4 holds 5/20 of it", r.at("x", 0)[4] - 0.25 * logistic, 1e-8);
    s.check("what is allocated sums to the total, to rounding", r.at("sum_x", 0)[0] - r.at("P", 0)[0], 1e-12);
  }

  // ---- a map with a period of two steps ----
  {
    Model m;
    m.quantities = {q("x", "unit", Role::state)};
    m.laws = {law("x", Form::map, "a * x", 2)};
    m.parameters = {par("a", 1.5)};
    m.step = 1;
    const Compiled c = compile(flat(1), m);
    Inputs in;
    in.values["x"] = vals({2});
    Options o;
    o.t0 = 0, o.t1 = 10;
    const Run r = simulate(c, in, o);
    s.check("a map of period 2 has not moved at t = 1", r.at("x", 1)[0] - 2, 1e-15);
    s.check("... has moved once at t = 2", r.at("x", 2)[0] - 3, 1e-15);
    s.check("... and five times at t = 10: 2 · 1.5⁵", r.at("x", 10)[0] - 2 * std::pow(1.5, 5), 1e-12);
  }

  // ---- diffusion between two cells ----
  {
    Level u{"cell", {"p", "q"}, {}, {{"area", {1, 1}}}, {{0, 1, 1}}, {}};
    const Hierarchy h({u});
    Model m;
    m.quantities = {q("x", "cell", Role::state)};
    m.derived = {d("lap", "cell", Op::laplacian, "x", "area")};
    m.laws = {law("x", Form::flow, "D * lap")};
    m.parameters = {par("D", 0.5)};
    m.step = 0.01;
    const Compiled c = compile(h, m);
    Inputs in;
    in.values["x"] = vals({3, 1});
    Options o;
    o.t0 = 0, o.t1 = 1;
    o.report = {1};
    const Run r = simulate(c, in, o);
    const auto& x = r.at("x", 0);
    s.check("the difference falls as 2 e^{−2Dt}", (x[0] - x[1]) - 2 * std::exp(-1.0), 1e-10);
    s.check("the content is kept", x[0] + x[1] - 4, 1e-12);
  }

  // ---- two paces ----
  {
    const auto error = [&](double h, std::size_t slow) {
      Model m;
      m.quantities = {q("x", "unit", Role::state), q("y", "unit", Role::state)};
      m.laws = {law("x", Form::flow, "-x + y", h * static_cast<double>(slow)), law("y", Form::flow, "-10 * y + x", h)};
      m.step = h;
      const Compiled c = compile(flat(1), m);
      Inputs in;
      in.values["x"] = vals({1});
      in.values["y"] = vals({0});
      Options o;
      o.t0 = 0, o.t1 = 1;
      o.report = {1};
      const Run r = simulate(c, in, o);
      calx::linalg::Matrix M(2);
      M(0, 0) = -1, M(0, 1) = 1, M(1, 0) = 1, M(1, 1) = -10;
      const calx::linalg::Matrix E = calx::linalg::expm(M);
      return std::hypot(r.at("x", 0)[0] - E(0, 0), r.at("y", 0)[0] - E(1, 0));
    };
    const double e1 = error(0.01, 10), e2 = error(0.005, 10), e0 = error(0.01, 1);
    s.require("two paces: the error halves with the slow step (first order)", e1 / e2 > 1.6 && e1 / e2 < 2.4,
              "ratio " + std::to_string(e1 / e2));
    s.require("one pace: Runge–Kutta's error, below 1e-8", e0 < 1e-8, std::to_string(e0));
  }

  // ---- a fast quantity at rest ----
  {
    Model m;
    m.quantities = {q("x", "unit", Role::state), q("y", "unit", Role::state)};
    m.laws = {law("x", Form::flow, "-y"), law("y", Form::equilibrium, "x / 2")};
    m.step = 0.01;
    const Compiled c = compile(flat(1), m);
    Inputs in;
    in.values["x"] = vals({1});
    Options o;
    o.t0 = 0, o.t1 = 2;
    o.report = {2};
    o.tolerance = 1e-14;
    const Run r = simulate(c, in, o);
    s.check("with y at rest, x(2) = e^{−1}", r.at("x", 0)[0] - std::exp(-1.0), 1e-10);
    s.require("... and the equilibrium needed no initial value", r.counts.unconverged == 0);
    const auto fast = [&](double eps) {
      Model f;
      f.quantities = {q("x", "unit", Role::state), q("y", "unit", Role::state)};
      f.laws = {law("x", Form::flow, "-y"), law("y", Form::flow, "(x / 2 - y) / eps")};
      f.parameters = {par("eps", eps)};
      f.step = 0.001;
      const Compiled cf = compile(flat(1), f);
      Inputs fi;
      fi.values["x"] = vals({1});
      fi.values["y"] = vals({0.5});
      Options fo;
      fo.t0 = 0, fo.t1 = 2;
      fo.report = {2};
      return std::abs(simulate(cf, fi, fo).at("x", 0)[0] - std::exp(-1.0));
    };
    const double a = fast(0.1), b = fast(0.05);
    s.require("as a fast flow, the error falls in proportion to ε (Tikhonov)", a / b > 1.6 && a / b < 2.4, "ratio " + std::to_string(a / b));
  }

  // ---- spatial interaction ----
  {
    Level u{"centre", {"p", "q", "r"}, {}, {}, {}, {{"walk", {0, 1, 2, 1, 0, 1, 2, 1, 0}}}};
    const Hierarchy h({u});
    Model m;
    m.quantities = {q("O", "centre", Role::fixed), q("A", "centre", Role::fixed, Kind::other)};
    Derived inflow = d("D", "centre", Op::inflow, "O", "A");
    inflow.costs = "walk", inflow.beta = "beta";
    m.derived = {inflow};
    m.parameters = {par("beta", 1)};
    const Compiled c = compile(h, m);
    Inputs in;
    in.values["O"] = vals({10, 0, 5});
    in.values["A"] = vals({1, 2, 1});
    Options o;
    o.t0 = 0, o.t1 = 1;
    o.report = {0};
    const Run r = simulate(c, in, o);
    const double O[3] = {10, 0, 5}, A[3] = {1, 2, 1}, C[3][3] = {{0, 1, 2}, {1, 0, 1}, {2, 1, 0}};
    double worst = 0, total = 0;
    for (int j = 0; j < 3; j++) {
      double D = 0;
      for (int i = 0; i < 3; i++) {
        double den = 0;
        for (int k = 0; k < 3; k++) den += A[k] * std::exp(-C[i][k]);
        D += O[i] * A[j] * std::exp(-C[i][j]) / den;
      }
      worst = std::max(worst, std::abs(r.at("D", 0)[j] - D));
      total += r.at("D", 0)[j];
    }
    s.check("the inflow against the formula by hand", worst, 1e-12);
    s.check("its total is the mass sent, 15", total - 15, 1e-12);
    in.parameters["beta"] = Value{0, Source::declared, 0, missing_value, missing_value};
    const Run r0 = simulate(c, in, o);
    s.check("with β = 0 a centre receives by its attractiveness alone: 15 · 2/4", r0.at("D", 0)[1] - 7.5, 1e-12);

    // Harris and Wilson: the centres grow by what they receive less what they cost.
    Model hw;
    hw.quantities = {q("O", "centre", Role::fixed), q("W", "centre", Role::state)};
    Derived dd = d("D", "centre", Op::inflow, "O", "W^alpha");
    dd.costs = "walk", dd.beta = "beta";
    hw.derived = {dd};
    hw.laws = {law("W", Form::flow, "eps * (D - kappa * W)")};
    hw.parameters = {par("alpha", 0.8), par("beta", 0.5), par("kappa", 2), par("eps", 1)};
    hw.step = 0.05;
    const Compiled ch = compile(h, hw);
    Inputs hi;
    hi.values["O"] = vals({10, 5, 8});
    hi.values["W"] = vals({1, 1, 1});
    Options ho;
    ho.t0 = 0, ho.t1 = 40;
    ho.report = {40};
    const Run rh = simulate(ch, hi, ho);
    const auto& W = rh.at("W", 0);
    const auto& D = rh.at("D", 0);
    s.check("Harris–Wilson: the centres together reach Σ O / κ = 11.5", W[0] + W[1] + W[2] - 11.5, 1e-9);
    double gap = 0;
    for (int j = 0; j < 3; j++) gap = std::max(gap, std::abs(D[j] - 2 * W[j]));
    s.check("... each at rest, D_j = κ W_j", gap, 1e-7);
  }

  // ---- logit shares, a registered reducer, an ordinal level ----
  {
    Model m;
    m.quantities = {q("u", "block", Role::fixed, Kind::other), q("z", "block", Role::fixed, Kind::other), q("pop", "block", Role::fixed)};
    m.derived = {d("p", "block", Op::logit, "u"), d("p_sum", "area", Op::sum, "p"), d("zu", "block", Op::zscore, "u"),
                 d("level", "block", Op::value, "step(z - 1) + step(z - 2)")};
    Derived med = d("median_level", "area", Op::reduce, "level", "pop");
    med.reducer = "weighted_median";
    m.derived.push_back(med);
    m.reducers["weighted_median"] = [](std::span<const double> v, std::span<const double> w) {
      std::vector<std::size_t> idx(v.size());
      for (std::size_t i = 0; i < idx.size(); i++) idx[i] = i;
      std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });
      double total = 0, run = 0;
      for (double x : w) total += x;
      for (std::size_t i : idx)
        if ((run += w[i]) >= total / 2) return v[i];
      return v[idx.back()];
    };
    const Compiled c = compile(three(), m);
    Inputs in;
    in.values["u"] = vals({0, 1, 2, -1, 0, 3});
    in.values["z"] = vals({0.5, 1.5, 2.5, 0.2, 0.4, 3});
    in.values["pop"] = vals({1, 1, 5, 4, 4, 1});
    Options o;
    o.t0 = 0, o.t1 = 1;
    o.report = {0};
    const Run r = simulate(c, in, o);
    s.check("logit shares sum to one in each group of siblings", std::abs(r.at("p_sum", 0)[0] - 1) + std::abs(r.at("p_sum", 0)[1] - 1), 1e-14);
    s.check("... each exp(u_i) / Σ exp(u_j)", r.at("p", 0)[1] - std::exp(1) / (1 + std::exp(1) + std::exp(2)), 1e-14);
    s.check("z-scores among siblings: (2 − 1)/1 in the first area", r.at("zu", 0)[2] - 1, 1e-14);
    s.check("... and (−1 − 2/3)/√(13/3) in the second", r.at("zu", 0)[3] - (-1 - 2.0 / 3.0) / std::sqrt(13.0 / 3.0), 1e-14);
    s.require("an ordinal level read through steps on its thresholds", r.at("level", 0)[0] == 0 && r.at("level", 0)[1] == 1 && r.at("level", 0)[2] == 2);
    s.require("the weighted median of the levels: 2 in the first area, 0 in the second",
              r.at("median_level", 0)[0] == 2 && r.at("median_level", 0)[1] == 0);
  }

  // ---- what is missing ----
  {
    Model m;
    m.quantities = {q("x", "block", Role::state), q("y", "block", Role::state), q("rain", "area", Role::driver)};
    m.derived = {d("rain_b", "block", Op::parent, "rain")};
    m.laws = {law("x", Form::flow, "r * x + rain_b"), law("y", Form::equilibrium, "x / 2")};
    m.parameters = {par("r", missing_value)};
    const Compiled c = compile(three(), m);
    Inputs in;
    in.values["x"] = vals({1, 1, missing_value, 1, missing_value, 1});
    const auto ms = missing(c, in, 0, 1);
    bool x_found = false, r_found = false, rain_found = false, y_found = false;
    for (const auto& e : ms) {
      if (e.quantity == "x" && e.units == 2 && e.of == 6 && e.level == "block") x_found = true;
      if (e.quantity == "r" && e.what == "a value of the parameter") r_found = true;
      if (e.quantity == "rain" && e.units == 2) rain_found = true;
      if (e.quantity == "y") y_found = true;
    }
    s.require("two blocks without x's initial value are named, 2 of 6", x_found);
    s.require("the parameter r without a value is named", r_found);
    s.require("the driver rain without a series is named at both areas", rain_found);
    s.require("the equilibrium y needs no initial value", !y_found);
    std::string said;
    Options o;
    o.t0 = 0, o.t1 = 1;
    s.require("the run refuses, saying what to give", refused([&] { simulate(c, in, o); }, "x at block (2 of 6 units)", &said), said);
    s.require("... and how", said.find("declare it with a range") != std::string::npos);
  }

  // ---- ensembles and screening ----
  {
    Model m;
    m.quantities = {q("x", "unit", Role::state)};
    m.laws = {law("x", Form::flow, "r * x")};
    m.parameters = {par("r", missing_value)};
    m.step = 0.01;
    const Compiled c = compile(flat(1), m);
    Inputs in;
    in.values["x"] = vals({1});
    in.parameters["r"] = Value{missing_value, Source::declared, 0, 0.1, 0.3};
    Options o;
    o.t0 = 0, o.t1 = 1;
    o.report = {1};
    EnsembleOptions eo;
    eo.draws = 200;
    const Ensemble e = ensemble(c, in, o, eo);
    const Band& b = e.bands.at("x");
    s.check("a declared range on r: median x(1) = e^{0.2}", b.middle[0][0] - std::exp(0.2), 3e-3);
    s.check("... its 5% quantile e^{0.11}", b.low[0][0] - std::exp(0.11), 3e-3);
    s.check("... its 95% quantile e^{0.29}", b.high[0][0] - std::exp(0.29), 3e-3);
    s.require("one factor, the declared parameter", e.factors.size() == 1 && e.factors[0].name == "r" && e.factors[0].source == Source::declared);

    Inputs two;
    two.values["x"] = {Value{missing_value, Source::declared, 0, 1, 2}, Value{missing_value, Source::declared, 0, 10, 20}};
    const auto fs = factors(two);
    const std::vector<double> u = {0.37};
    const Inputs at = at_quantiles(two, fs, u);
    s.check("a quantity's units are drawn at one quantile", (at.values.at("x")[0].value - 1) - (at.values.at("x")[1].value - 10) / 10, 1e-14);

    Inputs sc = in;
    sc.values["x"] = {Value{missing_value, Source::declared, 0, 0.95, 1.05}};
    const Screening scr = screen(c, sc, o, "x");
    std::size_t ir = 0, ix = 0;
    for (std::size_t k = 0; k < scr.factors.size(); k++) (scr.factors[k].name == "r" ? ir : ix) = k;
    s.require("screening ranks the growth rate above the narrower initial value", scr.mu_star[ir] > scr.mu_star[ix],
              std::to_string(scr.mu_star[ir]) + " vs " + std::to_string(scr.mu_star[ix]));
  }

  // ---- fitting ----
  {
    Model m;
    m.quantities = {q("P", "unit", Role::state)};
    m.laws = {law("P", Form::flow, "r * P * (1 - P / K)")};
    m.parameters = {par("r", 0.2, true, 0.01, 2), par("K", 100, true, 1, 500)};
    m.step = 0.05;
    const Compiled c = compile(flat(1), m);
    Inputs in;
    in.values["P"] = vals({5});
    const auto truth = [](double t) { return 50 / (1 + 9 * std::exp(-0.4 * t)); };
    std::vector<Observation> obs;
    for (int t = 0; t <= 20; t++) obs.push_back({"P", "u0", static_cast<double>(t), truth(t), 1});
    Options o;
    o.t0 = 0, o.t1 = 20;
    const Fitted f = fit(c, in, o, obs);
    s.check("r recovered from the exact record", f.values[0] / 0.4 - 1, 1e-4);
    s.check("K recovered from the exact record", f.values[1] / 50 - 1, 1e-4);
    s.require("... and both identified", f.not_identified.empty());
    FitOptions held;
    held.held_after = 12;
    const Fitted g = fit(c, in, o, obs, held);
    s.require("judged on the years held aside, far better than persistence", g.skill > 0.99, std::to_string(g.skill));
    s.require("eight observations held aside", g.held_count == 8);

    calx::random::Stream rs(20261118, 7, 0, 0);
    std::vector<Observation> noisy = obs;
    for (auto& e : noisy) e.value += 0.5 * rs.normal(), e.sigma = 0.5;
    const Fitted n = fit(c, in, o, noisy);
    s.require("a noisy record: r within four standard errors", std::abs(n.values[0] - 0.4) < 4 * n.standard_errors[0]);
    s.require("... and K", std::abs(n.values[1] - 50) < 4 * n.standard_errors[1]);

    Model prod;
    prod.quantities = {q("x", "unit", Role::state)};
    prod.laws = {law("x", Form::flow, "a * b * x")};
    prod.parameters = {par("a", 0.5, true, 0.01, 5), par("b", 0.5, true, 0.01, 5)};
    prod.step = 0.05;
    const Compiled cp = compile(flat(1), prod);
    Inputs pi;
    pi.values["x"] = vals({1});
    std::vector<Observation> po;
    for (int t = 0; t <= 10; t++) po.push_back({"x", "u0", static_cast<double>(t), std::exp(0.3 * t), 1});
    Options oo;
    oo.t0 = 0, oo.t1 = 10;
    const Fitted p = fit(cp, pi, oo, po);
    s.require("two unknowns that enter only as a product are reported as not identified", !p.not_identified.empty());
  }

  // ---- perfect aggregation ----
  {
    calx::linalg::Rectangular M(4, 4), A(2, 4);
    const double cols[4][4] = {{-3, 1, 1.5, 0.5}, {0.5, -2.5, 1, 1}, {0.3, 0.7, -2, 1}, {0.6, 0.4, 0.5, -1.5}};
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 4; i++) M(i, j) = cols[j][i];
    A(0, 0) = A(0, 1) = A(1, 2) = A(1, 3) = 1;
    const Aggregation g = perfect_aggregation(M, A);
    s.require("a lumpable law aggregates perfectly", g.perfect, std::to_string(g.residual));
    s.check("... into M̃ = [[−2, 1], [2, −1]]", std::abs(g.coarse(0, 0) + 2) + std::abs(g.coarse(0, 1) - 1) +
                                                   std::abs(g.coarse(1, 0) - 2) + std::abs(g.coarse(1, 1) + 1), 1e-12);
    M(0, 1) += 0.3;
    const Aggregation bad = perfect_aggregation(M, A);
    s.require("perturbed, it does not, and says by how much", !bad.perfect && bad.residual > 1e-3, std::to_string(bad.residual));
  }

  // ---- refusals ----
  {
    Model m;
    m.quantities = {q("x", "block", Role::state), q("rent", "area", Role::fixed)};
    m.laws = {law("x", Form::flow, "rent * x")};
    s.require("a law that names another level is refused, with the way to bring it down",
              refused([&] { compile(three(), m); }, "carry it to this level"));
    Model cyc;
    cyc.quantities = {q("x", "block", Role::fixed)};
    cyc.derived = {d("a", "block", Op::value, "b + x"), d("b", "block", Op::value, "a")};
    s.require("derived quantities in a circle are refused", refused([&] { compile(three(), cyc); }, "in a circle"));
    Model lawless;
    lawless.quantities = {q("x", "block", Role::state)};
    s.require("a state with no law is refused", refused([&] { compile(three(), lawless); }, "has no law"));
    Model odd;
    odd.quantities = {q("x", "block", Role::state)};
    odd.laws = {law("x", Form::map, "x", 1.5)};
    s.require("a period that is not a whole number of steps is refused", refused([&] { compile(three(), odd); }, "whole multiple"));
    Model red;
    red.quantities = {q("x", "block", Role::fixed)};
    Derived r = d("g", "area", Op::reduce, "x");
    r.reducer = "gini";
    red.derived = {r};
    s.require("a reducer the model was not given is refused", refused([&] { compile(three(), red); }, "has not been given"));
  }

  return s.done();
}
