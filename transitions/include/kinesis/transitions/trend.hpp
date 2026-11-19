// A change in an indicator, and how often the test of it is wrong.
//
// An indicator is a function of a window of the residual record (the chart's coordinates less their slow trend, removed by
// a Gaussian kernel smoother). Its series over rolling windows is reported with Kendall's τ against time, as the literature
// reports it (Dakos et al. 2012, PLoS ONE 7, e41010). The test, though, is of the change between the first window and the
// last: the logarithm of their ratio for an indicator that is positive by construction (a variance, an eigenvalue of a
// covariance), their difference otherwise. τ measures monotonicity and not magnitude, and with windows of half the record
// the rolling series of a stationary record has only one or two independent pieces, so that it is monotone, and τ near
// ±1, often enough for a test of τ to have little power however large the change (the theory note, §6.1; Boettiger and
// Hastings 2012, J. R. Soc. Interface 9, 2527–2539, found tests of τ weak by their ROC curves).
//
// The change's significance is the share of records of a null with a change at least as large in the indicator's direction
// (the record itself counted among them, so that a significance is never zero). The primary null is the stationary linear
// process fitted to the residuals; under it, and under the declared alternative of nulls.hpp (the leading mode slowing),
// the same pipeline (detrending, windows, change) is run on simulated records, which gives the test's operating
// characteristics at the record's own length: its false-alarm rate at the chosen level, its detection rate against the
// alternative, and the area under its ROC curve (Boettiger and Hastings's principle). Other nulls (random phases, IAAFT,
// the block bootstrap) give significances of their own, a robustness check; and a grid over window and bandwidth gives the
// share of choices under which the change is significant, and whether any choice reverses it, since a window chosen after
// looking inflates a trend. The joint test across the directed indicators takes the largest change, each standardised by
// its spread under the fitted null, against the same largest under the null, so that asking many indicators does not
// multiply false alarms.
//
// Kendall's τ of n distinct values against time has mean 0 and variance 2(2n + 5)/(9n(n − 1)) when they are exchangeable
// (Kendall 1938, Biometrika 30, 81–93); values from overlapping windows are not, which is one more reason significance
// comes from simulated records.
#pragma once

#include <kinesis/analysis/criticality.hpp>
#include <kinesis/transitions/nulls.hpp>

#include <tyche/descriptive/rank.hpp>

#include <functional>

namespace kinesis::transitions {

enum class Direction { rising, falling, either };

inline const char* name(Direction d) {
  switch (d) {
    case Direction::rising: return "rising";
    case Direction::falling: return "falling";
    default: return "changing";
  }
}

// An indicator: a value of the residual record's rows [from, to) at step Δt, with what it estimates, which way it moves
// before a transition it can see, whether it is positive by construction (its change then read on a logarithmic scale),
// and what it sees and is blind to.
struct Indicator {
  std::string name;
  std::string estimates;
  Direction direction = Direction::rising;
  std::function<double(const Rectangular& residual, std::size_t from, std::size_t to, double dt)> value;
  std::vector<std::string> sees, blind, assumes;
  bool positive = false;
};

struct TrendOptions {
  double window = 0.5;            // as a share of the record
  double bandwidth = 0.1;         // the smoother's standard deviation as a share of the record; 0 takes away the mean only
  std::size_t max_windows = 120;  // windows are taken at a stride that keeps at most this many
  std::size_t null_draws = 200, alternative_draws = 100;
  double remaining = 0.1;         // the alternative: the leading rate falls to this share of its value
  std::vector<Null> robustness = {Null::phase, Null::iaaft, Null::block};
  std::size_t surrogate_draws = 100;
  bool grid = true;
  std::vector<double> windows = {0.3, 0.4, 0.5};
  std::vector<double> bandwidths = {0.05, 0.1, 0.2};
  std::size_t grid_draws = 60;
  double alpha = 0.05;
  std::uint64_t seed = 1;
};

struct GridCell {
  double window = 0, bandwidth = 0, tau = 0, change = 0;
  double p = 1;                     // significance in the indicator's direction
  double p_reverse = 1;             // significance the other way
};

struct Trend {
  std::string indicator;
  Direction direction = Direction::rising;
  bool positive = false;
  std::vector<double> at, series;   // the time at which each window ends, and the indicator in it
  double tau = 0;                   // Kendall's τ of the series against time
  double change = 0;                // last window against first: log ratio (positive indicators) or difference
  double p = 1;                     // the change's significance under the fitted null
  double false_alarm = 0;           // the null's rejection rate at the level, as simulated
  double detection = 0;             // the alternative's rejection rate at the level
  double auc = 0.5;                 // the area under the ROC curve of the statistic, alternative against null
  double alternative_share = 1;     // the share of the alternative's records with a change no larger than the observed
  std::vector<std::pair<std::string, double>> robustness;  // significance under each further null
  std::vector<GridCell> grid;
  double grid_significant = 0;      // the share of grid cells significant in the indicator's direction
  bool grid_flips = false;          // some cell is significant the other way
};

struct TrendBattery {
  std::vector<Trend> trends;
  double joint_p = 1;               // the significance of the largest standardised change among the directed indicators
  std::string null_described;
  std::size_t window = 0, stride = 1;
};

namespace detail {

inline Rectangular detrend(const Rectangular& x, double bandwidth) {
  const std::size_t T = x.rows, d = x.cols;
  Rectangular out(T, d);
  for (std::size_t c = 0; c < d; c++) {
    const auto col = column(x, c);
    if (bandwidth > 0) {
      const auto s = analysis::smooth(col, bandwidth * static_cast<double>(T));
      for (std::size_t t = 0; t < T; t++) out(t, c) = col[t] - s[t];
    } else {
      const double m = tyche::descriptive::mean(col);
      for (std::size_t t = 0; t < T; t++) out(t, c) = col[t] - m;
    }
  }
  return out;
}

inline double score(Direction d, double change) {
  switch (d) {
    case Direction::rising: return change;
    case Direction::falling: return -change;
    default: return std::abs(change);
  }
}

struct Rolled {
  std::vector<double> at;
  std::vector<std::vector<double>> series;  // one per indicator
  std::vector<double> tau, change;
};

inline Rolled roll(const Rectangular& residual, const std::vector<double>& times, const std::vector<Indicator>& indicators,
                   std::size_t window, std::size_t stride, double dt) {
  Rolled r;
  r.series.assign(indicators.size(), {});
  for (std::size_t end = window; end <= residual.rows; end += stride) {
    r.at.push_back(times[end - 1]);
    for (std::size_t k = 0; k < indicators.size(); k++) {
      double v = std::numeric_limits<double>::quiet_NaN();
      try {
        v = indicators[k].value(residual, end - window, end, dt);
      } catch (const std::exception&) {  // a window the indicator cannot read (a constant one, a singular one) gives no value
      }
      r.series[k].push_back(v);
    }
  }
  // The last window is the one ending at the record's end, whatever the stride.
  if (r.at.empty() || r.at.back() != times[residual.rows - 1]) {
    r.at.push_back(times[residual.rows - 1]);
    for (std::size_t k = 0; k < indicators.size(); k++) {
      double v = std::numeric_limits<double>::quiet_NaN();
      try {
        v = indicators[k].value(residual, residual.rows - window, residual.rows, dt);
      } catch (const std::exception&) {
      }
      r.series[k].push_back(v);
    }
  }
  r.tau.resize(indicators.size());
  r.change.resize(indicators.size());
  for (std::size_t k = 0; k < indicators.size(); k++) {
    std::vector<double> a, s;
    for (std::size_t i = 0; i < r.at.size(); i++)
      if (std::isfinite(r.series[k][i])) a.push_back(r.at[i]), s.push_back(r.series[k][i]);
    r.tau[k] = s.size() >= 3 ? tyche::descriptive::kendall_tau_b(a, s) : 0.0;
    if (s.size() >= 2) {
      const double first = s.front(), last = s.back();
      r.change[k] = indicators[k].positive ? (first > 0 && last > 0 ? std::log(last / first) : 0.0) : last - first;
    } else {
      r.change[k] = 0;
    }
  }
  return r;
}

inline double upper_share(const std::vector<double>& draws, double value) {
  std::size_t above = 0;
  for (double v : draws) above += v >= value;
  return static_cast<double>(above + 1) / static_cast<double>(draws.size() + 1);
}

inline double lower_share(const std::vector<double>& draws, double value) {
  std::size_t below = 0;
  for (double v : draws) below += v <= value;
  return static_cast<double>(below + 1) / static_cast<double>(draws.size() + 1);
}

inline double quantile_of(std::vector<double> v, double q) {
  if (v.empty()) return std::numeric_limits<double>::infinity();
  std::sort(v.begin(), v.end());
  return tyche::descriptive::quantile_sorted(v, q);
}

inline double auc(const std::vector<double>& alternative, const std::vector<double>& null) {
  if (alternative.empty() || null.empty()) return 0.5;
  double s = 0;
  for (double a : alternative)
    for (double b : null) s += a > b ? 1.0 : (a == b ? 0.5 : 0.0);
  return s / (static_cast<double>(alternative.size()) * static_cast<double>(null.size()));
}

inline double spread(const std::vector<double>& v) {
  if (v.size() < 2) return 1;
  const double m = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
  double s = 0;
  for (double x : v) s += (x - m) * (x - m);
  const double sd = std::sqrt(s / static_cast<double>(v.size() - 1));
  return sd > 0 ? sd : 1;
}

}  // namespace detail

// The changes of a set of indicators on a record in chart coordinates, each with its significance, operating
// characteristics, robustness and grid.
inline TrendBattery trends(const Rectangular& coordinates, const std::vector<double>& times, double dt,
                           const std::vector<Indicator>& indicators, const TrendOptions& o = {}) {
  const std::size_t T = coordinates.rows;
  if (T < 40) throw Refusal("a trend in an indicator needs a record of at least 40 states; this one has " + std::to_string(T));
  if (indicators.empty()) throw Refusal("no indicators were asked for");
  auto window_of = [&](double share) {
    return std::clamp<std::size_t>(static_cast<std::size_t>(std::llround(share * static_cast<double>(T))), 10, T - 3);
  };
  auto stride_of = [&](std::size_t w) { return std::max<std::size_t>(1, (T - w + 1 + o.max_windows - 1) / o.max_windows); };
  TrendBattery out;
  out.window = window_of(o.window);
  out.stride = stride_of(out.window);
  const Rectangular residual = detail::detrend(coordinates, o.bandwidth);
  const auto observed = detail::roll(residual, times, indicators, out.window, out.stride, dt);

  const LinearNull fitted = fit_null(residual);
  out.null_described = fitted.described;
  const std::size_t K = indicators.size();
  std::vector<std::vector<double>> null_scores(K), alternative_scores(K);
  calx::random::Stream rng(o.seed, 0x74726e64u, 1, 0);  // 'trnd'
  for (std::size_t s = 0; s < o.null_draws; s++) {
    const Rectangular sim = detail::detrend(simulate_null(fitted, T, rng), o.bandwidth);
    const auto r = detail::roll(sim, times, indicators, out.window, out.stride, dt);
    for (std::size_t k = 0; k < K; k++) null_scores[k].push_back(detail::score(indicators[k].direction, r.change[k]));
  }
  calx::random::Stream rng_alt(o.seed, 0x74726e64u, 2, 0);
  for (std::size_t s = 0; s < o.alternative_draws; s++) {
    const Rectangular sim = detail::detrend(simulate_null(fitted, T, rng_alt, o.remaining), o.bandwidth);
    const auto r = detail::roll(sim, times, indicators, out.window, out.stride, dt);
    for (std::size_t k = 0; k < K; k++) alternative_scores[k].push_back(detail::score(indicators[k].direction, r.change[k]));
  }
  // The joint test: the largest standardised change among the directed indicators.
  std::vector<double> sd(K, 1), centre(K, 0);
  for (std::size_t k = 0; k < K; k++) {
    sd[k] = detail::spread(null_scores[k]);
    centre[k] = null_scores[k].empty() ? 0 : std::accumulate(null_scores[k].begin(), null_scores[k].end(), 0.0) / static_cast<double>(null_scores[k].size());
  }
  auto joint_of = [&](auto&& score_of) {
    double m = -1e300;
    for (std::size_t k = 0; k < K; k++)
      if (indicators[k].direction != Direction::either) m = std::max(m, (score_of(k) - centre[k]) / sd[k]);
    return m;
  };
  std::vector<double> joint_null;
  for (std::size_t s = 0; s < o.null_draws; s++) joint_null.push_back(joint_of([&](std::size_t k) { return null_scores[k][s]; }));
  out.joint_p = detail::upper_share(joint_null, joint_of([&](std::size_t k) { return detail::score(indicators[k].direction, observed.change[k]); }));

  // Robustness: significances under further nulls, shared across indicators.
  std::vector<std::vector<std::vector<double>>> robust(o.robustness.size(), std::vector<std::vector<double>>(K));
  for (std::size_t q = 0; q < o.robustness.size(); q++) {
    calx::random::Stream rs(o.seed, 0x726f6275u, static_cast<std::uint32_t>(q), 0);  // 'robu'
    for (std::size_t s = 0; s < o.surrogate_draws; s++) {
      const Rectangular sur = detail::detrend(surrogate(residual, o.robustness[q], rs, &fitted), o.bandwidth);
      const auto r = detail::roll(sur, times, indicators, out.window, out.stride, dt);
      for (std::size_t k = 0; k < K; k++) robust[q][k].push_back(detail::score(indicators[k].direction, r.change[k]));
    }
  }

  // The grid of windows and bandwidths, each cell against its own reduced fitted null.
  std::vector<std::vector<GridCell>> cells(K);
  if (o.grid) {
    for (double bw : o.bandwidths)
      for (double ws : o.windows) {
        const std::size_t w = window_of(ws), stride = stride_of(w);
        const Rectangular res = detail::detrend(coordinates, bw);
        const auto obs = detail::roll(res, times, indicators, w, stride, dt);
        const LinearNull fit = fit_null(res);
        std::vector<std::vector<double>> draws(K);
        calx::random::Stream rg(o.seed, 0x67726964u, static_cast<std::uint32_t>(std::llround(ws * 1000)),
                                static_cast<std::uint32_t>(std::llround(bw * 1000)));  // 'grid'
        for (std::size_t s = 0; s < o.grid_draws; s++) {
          const Rectangular sim = detail::detrend(simulate_null(fit, T, rg), bw);
          const auto r = detail::roll(sim, times, indicators, w, stride, dt);
          for (std::size_t k = 0; k < K; k++) draws[k].push_back(detail::score(indicators[k].direction, r.change[k]));
        }
        for (std::size_t k = 0; k < K; k++) {
          const double sc = detail::score(indicators[k].direction, obs.change[k]);
          cells[k].push_back({ws, bw, obs.tau[k], obs.change[k], detail::upper_share(draws[k], sc), detail::lower_share(draws[k], sc)});
        }
      }
  }

  for (std::size_t k = 0; k < K; k++) {
    Trend t;
    t.indicator = indicators[k].name;
    t.direction = indicators[k].direction;
    t.positive = indicators[k].positive;
    t.at = observed.at;
    t.series = observed.series[k];
    t.tau = observed.tau[k];
    t.change = observed.change[k];
    const double s_obs = detail::score(t.direction, t.change);
    t.p = detail::upper_share(null_scores[k], s_obs);
    const double critical = detail::quantile_of(null_scores[k], 1 - o.alpha);
    std::size_t fa = 0, hit = 0;
    for (double v : null_scores[k]) fa += v > critical;
    for (double v : alternative_scores[k]) hit += v > critical;
    t.false_alarm = static_cast<double>(fa) / static_cast<double>(std::max<std::size_t>(1, null_scores[k].size()));
    t.detection = static_cast<double>(hit) / static_cast<double>(std::max<std::size_t>(1, alternative_scores[k].size()));
    t.auc = detail::auc(alternative_scores[k], null_scores[k]);
    t.alternative_share = alternative_scores[k].empty() ? 1.0 : detail::lower_share(alternative_scores[k], s_obs);
    for (std::size_t q = 0; q < o.robustness.size(); q++) t.robustness.push_back({name(o.robustness[q]), detail::upper_share(robust[q][k], s_obs)});
    t.grid = cells[k];
    if (!t.grid.empty()) {
      std::size_t significant = 0;
      for (const auto& c : t.grid) {
        if (c.p < o.alpha) significant++;
        if (t.direction != Direction::either && c.p_reverse < o.alpha) t.grid_flips = true;
      }
      t.grid_significant = static_cast<double>(significant) / static_cast<double>(t.grid.size());
    }
    out.trends.push_back(t);
  }
  return out;
}

}  // namespace kinesis::transitions
