// Checks of kinesis/transitions: critical transitions read from records of structured states.
//
// - Charts: an identity chart copies; a composition chart's distance is Aitchison's, with the log total first when asked; a
//   histogram of a normal sample is read at its median and at one standard deviation; a covariance chart is zero at its
//   reference and its norm is the affine-invariant distance from it; a zero part, a matrix that is not positive definite and
//   a state of the wrong size are refused.
// - The linear core: the discrete and continuous Lyapunov equations at their closed forms; a VAR(1) of known propagator
//   recovered with its leading eigenvalue and both eigenvectors; its covariance at the Lyapunov solution; λ_max(C₀⁻¹C) of 2C₀
//   is 2 and unchanged by a linear change of coordinates; shrinkage making a short window's covariance well conditioned.
// - Rates: an OU process's κ by AR(1), by increments and by the drift's slope; its variance σ²/(2κ); the red-noise
//   autocorrelation symmetric in its two rates and continuous where they meet, and an OU process driven by OU noise giving
//   back the pair of rates, which the AR(1) estimate does not.
// - Kendall's τ of exchangeable values with variance 2(2n + 5)/(9n(n − 1)).
// - Surrogates: random phases keep the periodogram, IAAFT keeps the values exactly, the block bootstrap draws only observed values.
// - Detection: a fold approached slowly is a signal, its eigenvalue going to +1; a stationary OU process gives false alarms at
//   about the stated level; rate-induced tipping at constant r and noise that reddens with κ unchanged give no slowing; a
//   double well is read as two wells with switches; a map approaching a period-doubling has its autocorrelation falling and its
//   eigenvalue going to −1; the Hopf normal form's eigenvalue goes to the circle; a drifting rate is preferred by likelihood
//   with a bounded tipping time, a constant one gives none.
// - Structured states: on a network OU process whose rates fall alike the leading mode is uniform and the spatial variance
//   trends less than the largest covariance eigenvalue, while when one unit's rate falls the mode shows there; a critical mode
//   the noise does not excite gives no rise; the log total of a composition carries a slowing the shares do not; on a hierarchy
//   the totals of each group carry a slowing that the shares within them do not.
// - Refusals: an irregular record, a record too short.

#include <kinesis/transitions/transitions.hpp>

#include <calx/random/stream.hpp>
#include <calx/transforms/fft.hpp>
#include <calx/verification/check.hpp>
#include <tyche/compositions/simplex.hpp>

#include <cmath>
#include <functional>
#include <numbers>
#include <string>
#include <vector>

using namespace kinesis::transitions;
using kinesis::kernel::Refusal;

static bool refused(const std::function<void()>& f, const std::string& part) {
  try {
    f();
  } catch (const Refusal& r) {
    return std::string(r.what()).find(part) != std::string::npos;
  }
  return false;
}

static bool near(double a, double b, double relative) { return std::abs(a - b) <= relative * std::abs(b); }

static Options quick(std::uint64_t seed = 1) {
  Options o;
  o.trend.null_draws = 79;
  o.trend.alternative_draws = 40;
  o.trend.surrogate_draws = 39;
  o.trend.robustness = {Null::phase};
  o.trend.grid = false;
  o.trend.max_windows = 60;
  o.trend.seed = seed;
  o.likelihood_options.replicates = 40;
  o.channel_null_draws = 59;
  o.channel_alternative_draws = 30;
  return o;
}

static const Finding* find(const Assessment& a, const std::string& indicator, const std::string& channel = "the record") {
  for (const auto& f : a.findings)
    if (f.indicator == indicator && f.channel == channel) return &f;
  return nullptr;
}

static Record vector_record(const std::vector<double>& t, const std::vector<std::vector<double>>& states) {
  Record r;
  r.times = t;
  r.states = states;
  return r;
}

int main() {
  calx::verification::Suite s;

  // ---- charts ----
  {
    StateSpace id;
    id.units = 2;
    id.parts = 1;
    const Record r = vector_record({0, 1, 2}, {{1, 2}, {3, 4}, {5, 6}});
    const Charted c = chart(id, r);
    s.require("an identity chart copies the states", c.x(2, 1) == 6 && c.unit_of[1] == 1 && c.x.cols == 2);

    StateSpace comp;
    comp.parts = 3;
    comp.chart = Chart::composition;
    const std::vector<double> x = {0.2, 0.3, 0.5}, y = {0.6, 0.1, 0.3};
    const Charted cc = chart(comp, vector_record({0, 1}, {x, y}));
    const double euclid = std::hypot(cc.x(0, 0) - cc.x(1, 0), cc.x(0, 1) - cc.x(1, 1));
    s.require("a composition chart's Euclidean distance is Aitchison's", near(euclid, tyche::compositions::distance(x, y), 1e-12));
    comp.with_total = true;
    const Charted ct = chart(comp, vector_record({0, 1}, {{2, 3, 5}, {6, 1, 3}}));
    s.require("with the total, the first coordinate is its logarithm", near(ct.x(0, 0), std::log(10.0), 1e-12) && ct.x.cols == 3);
    s.require("a zero part of a composition is refused", refused([&] { chart(comp, vector_record({0, 1}, {{0, 3, 5}, {6, 1, 3}})); }, "not positive"));

    StateSpace hist;
    hist.parts = 400;
    hist.chart = Chart::histogram;
    for (int k = 0; k <= 400; k++) hist.bin_edges.push_back(-10 + 20.0 * k / 400);
    hist.probabilities = {0.5, 0.841344746};
    std::vector<double> counts(400);
    for (int k = 0; k < 400; k++) {
      const double a = hist.bin_edges[k], b = hist.bin_edges[k + 1];
      counts[k] = 0.5 * (std::erf((b - 2) / (1.5 * std::sqrt(2.0))) - std::erf((a - 2) / (1.5 * std::sqrt(2.0))));
    }
    const Charted ch = chart(hist, vector_record({0, 1}, {counts, counts}));
    s.require("a histogram of N(2, 1.5²) is read at its median 2", std::abs(ch.x(0, 0) - 2) < 1e-3);
    s.require("and at one standard deviation above it, 3.5", std::abs(ch.x(0, 1) - 3.5) < 1e-3);

    StateSpace cov;
    cov.parts = 4;
    cov.chart = Chart::covariance;
    cov.reference = 2;
    // R = [[2, 0.5], [0.5, 1]]; C = R^{1/2} diag(e^{0.3}, e^{-0.4}) R^{1/2} is at affine-invariant distance √(0.3² + 0.4²) = 0.5.
    const Matrix R = [] { Matrix m(2); m(0, 0) = 2; m(0, 1) = m(1, 0) = 0.5; m(1, 1) = 1; return m; }();
    const Matrix Rh = detail::symmetric_function(R, [](double v) { return std::sqrt(v); });
    Matrix D(2);
    D(0, 0) = std::exp(0.3), D(1, 1) = std::exp(-0.4);
    const Matrix C = Rh * D * Rh;
    const Charted cv = chart(cov, vector_record({0, 1, 2}, {R.a, R.a, C.a}));
    double n0 = 0, n2 = 0;
    for (std::size_t k = 0; k < 3; k++) n0 += cv.x(0, k) * cv.x(0, k), n2 += cv.x(2, k) * cv.x(2, k);
    s.require("a covariance chart is zero at its reference", std::sqrt(n0) < 1e-10);
    s.require("and its norm is the affine-invariant distance from the reference, 0.5", near(std::sqrt(n2), 0.5, 1e-9));
    Matrix bad(2);
    bad(0, 0) = 1, bad(0, 1) = bad(1, 0) = 2, bad(1, 1) = 1;
    s.require("a covariance that is not positive definite is refused", refused([&] { chart(cov, vector_record({0, 1, 2}, {R.a, R.a, bad.a})); }, "positive definite"));
    s.require("a state of the wrong size is refused", refused([&] { chart(id, vector_record({0, 1}, {{1}, {2}})); }, "declares"));
  }

  // ---- the linear core ----
  {
    Matrix J(2), B(2);
    J(0, 0) = -1, J(1, 1) = -2, B(0, 0) = 1, B(1, 1) = 1;
    const Matrix Cc = continuous_lyapunov(J, B);
    s.require("J C + C Jᵀ + I = 0 at diag(1/2, 1/4)", near(Cc(0, 0), 0.5, 1e-12) && near(Cc(1, 1), 0.25, 1e-12) && std::abs(Cc(0, 1)) < 1e-14);
    Matrix A(2);
    A(0, 0) = 0.5, A(1, 1) = 0.5;
    const Matrix Cd = discrete_lyapunov(A, B);
    s.require("C = A C Aᵀ + I at 4/3 for A = I/2", near(Cd(0, 0), 4.0 / 3, 1e-12));

    // A non-normal propagator with eigenvalues 0.9 and 0.5.
    Matrix P(2);
    P(0, 0) = 0.9, P(0, 1) = 0.3, P(1, 0) = 0.0, P(1, 1) = 0.5;
    Matrix Bn(2);
    Bn(0, 0) = 1, Bn(1, 1) = 1;
    calx::random::Stream rng(20261119, 135, 1, 0);
    const auto states = models::var1(40000, [&](double) { return P; }, Bn, rng);
    Rectangular X(states.size(), 2);
    for (std::size_t t = 0; t < states.size(); t++) X(t, 0) = states[t][0], X(t, 1) = states[t][1];
    const Propagator est = propagator(X, 0.0);
    s.require("a VAR(1)'s leading eigenvalue 0.9 is recovered", std::abs(std::abs(est.spectrum[0]) - 0.9) < 0.01, std::to_string(std::abs(est.spectrum[0])));
    s.require("and its second, 0.5", std::abs(std::abs(est.spectrum[1]) - 0.5) < 0.02);
    const Leading L = leading(est, 1.0);
    // Right eigenvector of 0.9: (1, 0); left: (0.4, 0.3)/0.5 ∝ (4, 3) normalised... wᵀP = 0.9 wᵀ gives w ∝ (1, 0.75).
    const double rcos = std::abs(L.right[0]) / std::hypot(std::abs(L.right[0]), std::abs(L.right[1]));
    const double wl = std::hypot(std::abs(L.left[0]), std::abs(L.left[1]));
    const double lcos = (std::abs(L.left[0]) * 1 + std::abs(L.left[1]) * 0.75) / (wl * std::hypot(1.0, 0.75));
    s.require("its right eigenvector (where it shows) is (1, 0)", rcos > 0.999, std::to_string(rcos));
    s.require("and its left eigenvector (what excites it) is along (1, 0.75)", lcos > 0.999, std::to_string(lcos));
    Matrix Q(2);
    Q(0, 0) = 1, Q(1, 1) = 1;
    const Matrix exact = discrete_lyapunov(P, Q);
    const Covariance emp = covariance(X, false);
    s.require("its covariance is the Lyapunov solution", near(emp.c(0, 0), exact(0, 0), 0.05) && near(emp.c(1, 1), exact(1, 1), 0.05) && near(emp.c(0, 1), exact(0, 1), 0.08));

    Matrix C0(2), M(2);
    C0(0, 0) = 2, C0(0, 1) = C0(1, 0) = 0.3, C0(1, 1) = 1;
    M(0, 0) = 1, M(0, 1) = 2, M(1, 0) = -0.5, M(1, 1) = 3;
    Matrix C2 = C0;
    for (double& v : C2.a) v *= 2;
    s.require("λ_max(C₀⁻¹ 2C₀) = 2", near(relative_growth(C0, C2), 2, 1e-12));
    Matrix Ct = C0;
    Ct(0, 0) = 5;
    auto congruent = [&](const Matrix& C) {
      Matrix Mt(2);
      for (std::size_t i = 0; i < 2; i++)
        for (std::size_t j = 0; j < 2; j++) Mt(i, j) = M(j, i);
      return M * C * Mt;
    };
    s.require("and λ_max(C₀⁻¹C) is unchanged by a linear change of coordinates",
              near(relative_growth(congruent(C0), congruent(Ct)), relative_growth(C0, Ct), 1e-9));

    calx::random::Stream rs(20261119, 135, 2, 0);
    Rectangular W(30, 20);
    for (auto& v : W.a) v = rs.normal();
    const auto e_shrunk = calx::linalg::symmetric_eigen(covariance(W, true).c), e_raw = calx::linalg::symmetric_eigen(covariance(W, false).c);
    const double cond_shrunk = e_shrunk.values.back() / e_shrunk.values.front(), cond_raw = e_raw.values.back() / std::max(1e-300, e_raw.values.front());
    s.require("shrinkage makes a window of 30 states of 20 coordinates well conditioned", cond_shrunk < 10 && cond_raw > 50,
              std::to_string(cond_shrunk) + " against " + std::to_string(cond_raw));
  }

  // ---- rates ----
  {
    calx::random::Stream rng(20261119, 135, 3, 0);
    const auto t = models::times(4000, 0.1);
    const auto x = models::ou(t, [](double) { return 0.5; }, 1.0, rng);
    const Rate a = ar1_rate(x, 0.1), inc = increment_rate(x, 0.1), dr = drift_rate(x, 0.1);
    s.require("an OU process's κ = 0.5 by AR(1), within its interval", a.lower <= 0.5 && 0.5 <= a.upper && near(a.kappa, 0.5, 0.2), std::to_string(a.kappa));
    s.require("by increments with AR(1) noise", near(inc.kappa, 0.5, 0.2), std::to_string(inc.kappa));
    s.require("by the drift's slope about the mode", near(dr.kappa, 0.5, 0.3), std::to_string(dr.kappa));
    s.require("its variance σ²/(2κ) = 1", near(tyche::descriptive::variance(x, 1), 1.0, 0.15));
    s.require("the red-noise autocorrelation is symmetric in its two rates", std::abs(red_autocorrelation(0.5, 3, 0.7) - red_autocorrelation(3, 0.5, 0.7)) < 1e-14);
    s.require("is one at lag zero and continuous where the rates meet",
              red_autocorrelation(0.5, 3, 0) == 1 && std::abs(red_autocorrelation(1, 1 + 1e-7, 0.7) - 1.7 * std::exp(-0.7)) < 1e-6);
    {
      const auto yy = models::reddening(40000, 0.05, 0.5, 5.0, 5.0, 1.0, rng);
      const auto emp = tyche::time_series::acf(yy, 60);
      double worst = 0;
      for (int k = 1; k <= 60; k++) worst = std::max(worst, std::abs(emp[k] - red_autocorrelation(0.5, 5.0, k * 0.05)));
      s.require("a simulated OU process driven by OU noise has the red-noise autocorrelation, within its sampling error", worst < 0.08, std::to_string(worst));
    }
    const auto y = models::reddening(40000, 0.05, 0.5, 5.0, 5.0, 1.0, rng);
    const Rate r = red_rate(y, 0.05), ar = ar1_rate(y, 0.05);
    s.require("an OU process driven by OU noise (κ = 0.5, θ = 5) gives back κ under red noise", near(r.kappa, 0.5, 0.3), std::to_string(r.kappa));
    s.require("and the noise's rate θ", near(r.noise_rate, 5.0, 0.4), std::to_string(r.noise_rate));
    s.require("which the AR(1) estimate confounds", std::abs(ar.kappa - 0.5) > 0.3, std::to_string(ar.kappa));
    s.require("the red estimate says the record identifies only the pair", r.assumptions.size() >= 2 && r.assumptions[1].find("pair") != std::string::npos);
  }

  // ---- Kendall's τ under exchangeability ----
  {
    calx::random::Stream rng(20261119, 135, 4, 0);
    const std::size_t n = 30;
    std::vector<double> at(n);
    for (std::size_t i = 0; i < n; i++) at[i] = static_cast<double>(i);
    double s2 = 0;
    const int draws = 4000;
    for (int k = 0; k < draws; k++) {
      std::vector<double> v(n);
      for (auto& z : v) z = rng.normal();
      const double tau = tyche::descriptive::kendall_tau_b(at, v);
      s2 += tau * tau;
    }
    const double expected = 2.0 * (2 * n + 5) / (9.0 * n * (n - 1));
    s.require("Kendall's τ of exchangeable values has variance 2(2n + 5)/(9n(n − 1))", near(s2 / draws, expected, 0.07),
              std::to_string(s2 / draws) + " against " + std::to_string(expected));
  }

  // ---- surrogates ----
  {
    calx::random::Stream rng(20261119, 135, 5, 0);
    const auto t = models::times(256, 1.0);
    const auto x = models::ou(t, [](double) { return 0.2; }, 1.0, rng);
    Rectangular X(256, 1);
    for (std::size_t i = 0; i < 256; i++) X(i, 0) = x[i];
    const Rectangular ph = surrogate(X, Null::phase, rng), ia = surrogate(X, Null::iaaft, rng), bl = surrogate(X, Null::block, rng);
    auto amplitudes = [](const Rectangular& r) {
      std::vector<calx::transforms::Complex> z(r.rows);
      for (std::size_t i = 0; i < r.rows; i++) z[i] = r(i, 0);
      calx::transforms::fft(z);
      std::vector<double> a(z.size());
      for (std::size_t i = 0; i < z.size(); i++) a[i] = std::abs(z[i]);
      return a;
    };
    const auto a0 = amplitudes(X), a1 = amplitudes(ph);
    double worst = 0;
    for (std::size_t i = 0; i < a0.size(); i++) worst = std::max(worst, std::abs(a0[i] - a1[i]));
    s.require("random phases keep the periodogram", worst < 1e-8);
    auto sorted = [](const Rectangular& r) {
      std::vector<double> v(r.a);
      std::sort(v.begin(), v.end());
      return v;
    };
    s.require("IAAFT keeps the values exactly", sorted(ia) == sorted(X));
    bool subset = true;
    const auto values = sorted(X);
    for (double v : bl.a) subset &= std::binary_search(values.begin(), values.end(), v);
    s.require("the block bootstrap draws only observed values", subset);
  }

  // ---- detection on one coordinate ----
  {
    calx::random::Stream rng(20261119, 135, 6, 0);
    const std::size_t T = 1500;
    const double dt = 0.2;
    const auto x = models::fold(T, dt, 1.0, 0.02, 0.05, rng);
    StateSpace one;
    Options o = quick(7);
    o.trend.grid = true;
    o.trend.windows = {0.4, 0.5, 0.6};
    o.trend.bandwidths = {0.1};
    o.trend.grid_draws = 39;
    const Assessment a = assess(one, scalar_record(models::times(T, dt), x), o);
    const Finding* v = find(a, "variance");
    const Finding* k = find(a, "restoring rate (AR(1))");
    s.require("a fold approached slowly: variance gives a signal", v && v->verdict == Verdict::signal, v ? v->reasons.front() : "missing");
    s.require("and the restoring rate falls, a signal", k && k->verdict == Verdict::signal, k ? k->reasons.front() : "missing");
    s.require("the assessment's overall verdict is a signal", a.overall == Verdict::signal);
    s.require("its dominant eigenvalue goes to +1", a.kind.approach == Approach::stationary, a.kind.sentence);
    s.require("every finding says what it is blind to, rate-induced tipping among it",
              v && std::find(v->blind.begin(), v->blind.end(), blind_to_all()[0]) != v->blind.end() && !a.cannot_see.empty());
    s.require("the variance's test detects the declared alternative most of the time", v && v->trend.detection > 0.6, v ? std::to_string(v->trend.detection) : "");
  }
  {
    // False alarms on a stationary OU process at about the stated level.
    const auto ind = univariate_indicators(0, {RateEstimator::ar1});
    std::vector<Indicator> two = {ind[0], ind[2]};  // variance and the AR(1) rate
    TrendOptions to;
    to.null_draws = 99, to.alternative_draws = 0, to.robustness = {}, to.grid = false, to.max_windows = 50;
    int alarms = 0;
    const int records = 40;
    for (int k = 0; k < records; k++) {
      calx::random::Stream rng(20261119, 135, 100 + k, 0);
      const auto t = models::times(300, 0.2);
      const auto x = models::ou(t, [](double) { return 0.5; }, 1.0, rng);
      Rectangular X(300, 1);
      for (std::size_t i = 0; i < 300; i++) X(i, 0) = x[i];
      to.seed = 1000 + k;
      const auto b = trends(X, t, 0.2, two, to);
      alarms += b.trends[0].p < 0.05;
    }
    s.require("a stationary OU process gives false alarms at about the stated level (at most 6 in 40)", alarms <= 6, std::to_string(alarms));
  }
  {
    // Rate-induced tipping at constant r < 1: the stable state's rate is constant, so nothing slows.
    const auto ind = univariate_indicators(0, {RateEstimator::ar1});
    std::vector<Indicator> two = {ind[0], ind[2]};
    TrendOptions to;
    to.null_draws = 99, to.alternative_draws = 0, to.robustness = {}, to.grid = false, to.max_windows = 50;
    int alarms = 0;
    for (int k = 0; k < 10; k++) {
      calx::random::Stream rng(20261119, 135, 200 + k, 0);
      const auto x = models::rate_tipping(400, 0.1, 0.6, 0.05, rng);
      Rectangular X(400, 1);
      for (std::size_t i = 0; i < 400; i++) X(i, 0) = x[i];
      to.seed = 2000 + k;
      const auto b = trends(X, models::times(400, 0.1), 0.1, two, to);
      alarms += b.trends[1].p < 0.05;
    }
    s.require("rate-induced tipping at constant r = 0.6 shows no falling rate (at most 2 alarms in 10)", alarms <= 2, std::to_string(alarms));
    calx::random::Stream rng(20261119, 135, 299, 0);
    const auto tipped = models::rate_tipping(400, 0.1, 1.3, 0.05, rng);
    const double y_end = tipped.back() + 1.3 * 39.9;
    s.require("while at r = 1.3 the system tips (the moving coordinate leaves −1)", y_end > 0 || std::abs(tipped.back()) > 40);
  }
  {
    // Noise that reddens with κ unchanged: variance and autocorrelation rise, the red-noise rate does not fall.
    calx::random::Stream rng(20261119, 135, 7, 0);
    const auto y = models::reddening(24000, 0.05, 0.5, 8.0, 1.0, 1.0, rng);
    std::vector<double> early(y.begin(), y.begin() + 8000), late(y.end() - 8000, y.end());
    const Rate re = red_rate(early, 0.05), rl = red_rate(late, 0.05), ae = ar1_rate(early, 0.05), al = ar1_rate(late, 0.05);
    s.require("reddening noise lowers the AR(1) rate by more than a third", al.kappa < 0.67 * ae.kappa, std::to_string(ae.kappa) + " → " + std::to_string(al.kappa));
    s.require("while the rate under red noise stays within 40% of 0.5", near(re.kappa, 0.5, 0.4) && near(rl.kappa, 0.5, 0.4),
              std::to_string(re.kappa) + " → " + std::to_string(rl.kappa));
  }
  {
    calx::random::Stream rng(20261119, 135, 8, 0);
    const auto w2 = models::double_well(4000, 0.1, 0.45, rng);
    const Wells two = wells(w2, 0.1);
    s.require("a double well at fixed parameters is read as two wells", two.wells == 2, std::to_string(two.wells));
    s.require("with switches between them", two.switches >= 2, std::to_string(two.switches));
    s.require("and a finite Kramers time", std::isfinite(two.escape_time));
    const auto o1 = models::ou(models::times(4000, 0.1), [](double) { return 0.5; }, 1.0, rng);
    s.require("an OU process is read as one well", wells(o1, 0.1).wells == 1);
  }
  {
    calx::random::Stream rng(20261119, 135, 9, 0);
    const std::size_t T = 600;
    const auto x = models::flip(T, 0.6, 0.03, 1.0, rng);
    const Assessment a = assess(StateSpace{}, scalar_record(models::times(T, 1.0), x), quick(9));
    const Finding* ac = find(a, "lag-1 autocorrelation");
    const Finding* dm = find(a, "dominant eigenvalue modulus");
    s.require("before a period-doubling the lag-1 autocorrelation falls", ac && ac->trend.tau < -0.3, ac ? std::to_string(ac->trend.tau) : "");
    s.require("while the dominant eigenvalue's modulus rises, a signal", dm && dm->verdict == Verdict::signal, dm ? dm->reasons.front() : "");
    s.require("and goes to −1", a.kind.approach == Approach::flip, a.kind.sentence);
  }
  {
    calx::random::Stream rng(20261119, 135, 10, 0);
    const std::size_t T = 800;
    {
      const auto zs = models::hopf(3000, 0.25, -0.05, -0.05, 2.0, 0.02, rng);
      Rectangular Z(3000, 2);
      for (std::size_t i = 0; i < 3000; i++) Z(i, 0) = zs[i][0], Z(i, 1) = zs[i][1];
      const auto p = propagator(Z, 1e-4);
      s.require("a stationary Hopf normal form at μ = −0.05 has its propagator's modulus e^{μΔt} = 0.988", std::abs(std::abs(p.spectrum.front()) - std::exp(-0.05 * 0.25)) < 0.006,
                std::to_string(std::abs(p.spectrum.front())));
      s.require("and its argument ωΔt = 0.5", std::abs(std::abs(std::arg(p.spectrum.front())) - 0.5) < 0.02, std::to_string(std::arg(p.spectrum.front())));
    }
    const auto z = models::hopf(T, 0.25, -1.0, -0.05, 2.0, 0.02, rng);
    StateSpace plane;
    plane.parts = 2;
    const Assessment a = assess(plane, vector_record(models::times(T, 0.25), z), quick(10));
    const Finding* lead = find(a, "leading propagator modulus");
    s.require("the Hopf normal form's leading propagator modulus rises, a signal", lead && lead->verdict == Verdict::signal, lead ? lead->reasons.front() : "");
    s.require("and its eigenvalue goes to the circle off the real axis", a.kind.approach == Approach::oscillatory, a.kind.sentence);
    s.require("at a frequency near ω = 2", std::abs(a.kind.frequency - 2) < 0.4, std::to_string(a.kind.frequency));
    std::vector<double> xs;
    for (const auto& v : z) xs.push_back(v[0]);
    const Assessment b = assess(StateSpace{}, scalar_record(models::times(T, 0.25), xs), quick(11));
    s.require("one coordinate of it, read by an AR(2), points to the circle too", b.kind.approach == Approach::oscillatory, b.kind.sentence);
  }
  {
    calx::random::Stream rng(20261119, 135, 11, 0);
    const auto t = models::times(500, 0.1);
    const double duration = t.back();
    const auto x = models::ou(t, [duration](double s) { return 1.0 - 0.9 * s / duration; }, 0.5, rng);
    LikelihoodOptions lo;
    lo.replicates = 60;
    const auto test = likelihood_test(t, x, lo);
    const double truth = duration / 0.9;
    s.require("a linearly falling rate is preferred to a constant one by likelihood", test.p < 0.05, std::to_string(test.p));
    s.require("with an area under its ROC curve above 0.8", test.auc > 0.8, std::to_string(test.auc));
    s.require("and a tipping time, when one is given, with the truth inside its interval",
              !test.tipping_reported || (test.tipping_lower <= truth && truth <= test.tipping_upper),
              test.tipping_sentence);
    const auto flat = models::ou(t, [](double) { return 0.5; }, 0.5, rng);
    const auto none = likelihood_test(t, flat, lo);
    s.require("a constant rate gives no tipping time", !none.tipping_reported, none.tipping_sentence);
  }

  // ---- structured states ----
  {
    const std::size_t n = 8;
    std::vector<std::pair<std::size_t, std::size_t>> ring;
    for (std::size_t i = 0; i < n; i++) ring.push_back({i, (i + 1) % n});
    StateSpace net;
    net.units = n;
    for (const auto& [a, b] : ring) net.edges.push_back({a, b, 1.0});
    const std::size_t T = 700;
    const double dt = 0.2;
    calx::random::Stream rng(20261119, 135, 12, 0);
    const auto alike = models::network_ou(T, dt, n, ring, 1.0, [&](double t) { return std::vector<double>(n, 1.0 - 0.95 * t / ((T - 1) * dt)); }, 0.3, rng);
    const Assessment a = assess(net, vector_record(models::times(T, dt), alike), quick(12));
    const Finding* big = find(a, "largest covariance eigenvalue");
    const Finding* sv = find(a, "spatial variance");
    s.require("rates falling alike on a network: the largest covariance eigenvalue gives a signal", big && big->verdict == Verdict::signal, big ? big->reasons.front() : "");
    double cosine = 0;
    if (a.where)
      for (double v : a.where->loading) cosine += std::sqrt(v) / std::sqrt(static_cast<double>(n));
    s.require("the leading mode is uniform over the units (its loadings' cosine with the uniform mode above 0.95)", a.where && cosine > 0.95, std::to_string(cosine));
    s.require("and the spatial variance changes less than the largest eigenvalue (the mode shows in the mean, not the spread)",
              big && sv && sv->trend.change < 0.6 * big->trend.change, sv && big ? std::to_string(sv->trend.change) + " against " + std::to_string(big->trend.change) : "");
    calx::random::Stream rng2(20261119, 135, 13, 0);
    const std::size_t T2 = 1000;
    const auto one = models::network_ou(T2, dt, n, ring, 0.05, [&](double t) {
      std::vector<double> g(n, 1.0);
      g[2] = 1.0 - 0.99 * t / ((T2 - 1) * dt);
      return g;
    }, 0.3, rng2);
    const Assessment b = assess(net, vector_record(models::times(T2, dt), one), quick(13));
    std::size_t top = 0;
    if (b.where)
      for (std::size_t u = 1; u < n; u++)
        if (b.where->loading[u] > b.where->loading[top]) top = u;
    s.require("when one unit's rate falls, the leading mode shows there", b.where && top == 2, b.where ? b.where->sentence : "");
    const Finding* rg = find(b, "relative variance growth");
    s.require("and the variance's growth relative to the first window gives a signal", rg && rg->verdict == Verdict::signal, rg ? rg->reasons.front() : "");
    const Finding* sv2 = find(b, "spatial variance");
    s.require("while the spatial variance, diluted over eight units, rises less", sv2 && rg && sv2->trend.change > 0 && sv2->trend.change < rg->trend.change,
              sv2 && rg ? std::to_string(sv2->trend.change) + " against " + std::to_string(rg->trend.change) : "");
  }
  {
    // A critical mode the noise does not excite: x₁ = μ(t) x₁ is never pushed, x₂ is driven only through itself.
    calx::random::Stream rng(20261119, 135, 14, 0);
    Matrix Bs(2);
    Bs(1, 1) = 1;
    const auto states = models::var1(600, [](double s) {
      Matrix A(2);
      A(0, 0) = 0.5 + 0.49 * s, A(1, 0) = 0.2, A(1, 1) = 0.5;
      return A;
    }, Bs, rng);
    StateSpace plane;
    plane.parts = 2;
    const Assessment a = assess(plane, vector_record(models::times(600, 1.0), states), quick(14));
    const Finding* big = find(a, "largest covariance eigenvalue");
    s.require("a critical mode the noise does not excite gives no signal in the covariance", big && big->verdict != Verdict::signal, big ? big->reasons.front() : "");
    s.require("and every covariance indicator says it is blind to it",
              big && std::any_of(big->blind.begin(), big->blind.end(), [](const std::string& b) { return b.find("does not excite") != std::string::npos; }));
  }
  {
    // A population of three parts whose total slows while its shares do not.
    calx::random::Stream rng(20261119, 135, 15, 0);
    const std::size_t T = 1200;
    const double dt = 0.2;
    const auto t = models::times(T, dt);
    const double duration = t.back();
    const auto total = models::ou(t, [duration](double s) { return 1.0 - 0.97 * s / duration; }, 0.2, rng);
    const auto b1 = models::ou(t, [](double) { return 1.0; }, 0.2, rng), b2 = models::ou(t, [](double) { return 1.0; }, 0.2, rng);
    std::vector<std::vector<double>> counts;
    for (std::size_t i = 0; i < T; i++) {
      const auto shares = tyche::compositions::ilr_inverse(std::vector<double>{b1[i], b2[i]}, 1.0);
      const double size = 1000 * std::exp(total[i]);
      counts.push_back({size * shares[0], size * shares[1], size * shares[2]});
    }
    StateSpace pop;
    pop.parts = 3;
    pop.chart = Chart::composition;
    pop.with_total = true;
    const Assessment with_total = assess(pop, vector_record(t, counts), quick(15));
    pop.with_total = false;
    const Assessment shares_only = assess(pop, vector_record(t, counts), quick(16));
    auto any_signal = [](const Assessment& a) {
      for (const char* name : {"largest covariance eigenvalue", "relative variance growth", "leading propagator modulus"}) {
        const Finding* f = find(a, name);
        if (f && f->verdict == Verdict::signal) return true;
      }
      return false;
    };
    s.require("the log total of a composition carries a slowing", any_signal(with_total), with_total.summary.back());
    s.require("which the shares alone do not", !any_signal(shares_only), shares_only.summary.back());
  }
  {
    // A hierarchy: eight blocks in two districts; each district's blocks move together with a rate that falls.
    using kinesis::scales::Hierarchy;
    using kinesis::scales::Level;
    Level blocks{"block", {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8"}, {0, 0, 0, 0, 1, 1, 1, 1}, {}, {}, {}};
    Level districts{"district", {"north", "south"}, {}, {}, {}, {}};
    StateSpace city;
    city.units = 8;
    city.kind = Kind::extensive;
    city.hierarchy = Hierarchy({blocks, districts});
    calx::random::Stream rng(20261119, 135, 16, 0);
    const std::size_t T = 600;
    const double dt = 0.2;
    const auto t = models::times(T, dt);
    const double duration = t.back();
    const auto north = models::ou(t, [duration](double s) { return 1.0 - 0.95 * s / duration; }, 0.2, rng);
    const auto south = models::ou(t, [duration](double s) { return 1.0 - 0.95 * s / duration; }, 0.2, rng);
    std::vector<std::vector<double>> states(T, std::vector<double>(8));
    std::vector<std::vector<double>> idio(8);
    for (std::size_t u = 0; u < 8; u++) idio[u] = models::ou(t, [](double) { return 1.0; }, 0.1, rng);
    for (std::size_t i = 0; i < T; i++)
      for (std::size_t u = 0; u < 8; u++) states[i][u] = 100 * std::exp((u < 4 ? north[i] : south[i]) + idio[u][i]);
    const Assessment a = assess(city, vector_record(t, states), quick(17));
    const Finding* totals = find(a, "relative unit variance", "district: totals");
    const Finding* within = find(a, "relative unit variance", "district: shares within");
    const Finding* tl = find(a, "leading propagator modulus", "district: totals");
    s.require("on a hierarchy, the districts' totals carry the slowing", (tl && tl->verdict == Verdict::signal) || (totals && totals->verdict == Verdict::signal),
              tl ? tl->reasons.front() : "missing");
    s.require("and the shares within them do not", within && within->verdict != Verdict::signal, within ? within->reasons.front() : "missing");
  }

  // ---- refusals ----
  {
    Record irregular = scalar_record({0, 1, 3, 4}, {1, 2, 3, 4});
    s.require("an irregular record is refused for slowing down", refused([&] { assess(StateSpace{}, irregular); }, "regularly sampled"));
    std::vector<double> t(20), x(20);
    for (int i = 0; i < 20; i++) t[i] = i, x[i] = std::sin(i);
    s.require("a record too short is refused", refused([&] { assess(StateSpace{}, scalar_record(t, x)); }, "at least 40"));
  }

  return s.done();
}
