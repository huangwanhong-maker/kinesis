// An assessment of a record: every indicator its space admits, each with what it sees and is blind to, its null, its
// operating characteristics at the record's length, and a verdict of signal, no signal or undecidable with the reasons.
//
// The verdict follows the survey's synthesis (Part I §9.2, Part III §6.3, Part V §10.2):
//   signal       the trend is significant under the fitted null at the stated level, under at least half of the further
//                nulls, and in at least half of the grid of windows and bandwidths, with no choice reversing it;
//   no signal    the change is not significant, the test would have detected the declared alternative (the leading mode
//                slowing to a tenth of its rate over the record) at least as often as the required power, and a change as
//                small as the one observed is rare under that alternative: an absence is reported only where a presence could
//                have been seen and was not;
//   undecidable  otherwise, with each reason said: significant under one null and not another, sensitive to the window,
//                too weak to have seen the alternative, a record whose lag-1 autocorrelation is already near one, a variance
//                that rises while no estimate of the restoring rate falls (a change of noise, not slowing down), or, where the
//                noise is red, a slowing the restoring rate under red noise does not confirm (the survey, Part I §9.2).
//
// Where the state has many coordinates the assessment also says where the slowing shows (the leading mode's loadings over
// the units) and what excites it, which kind of loss of stability the dominant eigenvalue points to, and, on a hierarchy,
// at which level and in which channel a warning lives: the totals of each group (which keep the mode uniform within a
// group) or the contrasts within it (which take that mode away). Every assessment lists what none of its detectors can
// see.
#pragma once

#include <kinesis/transitions/indicators.hpp>
#include <kinesis/transitions/likelihood.hpp>
#include <kinesis/transitions/states.hpp>

#include <set>

namespace kinesis::transitions {

enum class Verdict { signal, no_signal, undecidable };

inline const char* name(Verdict v) {
  switch (v) {
    case Verdict::signal: return "signal";
    case Verdict::no_signal: return "no signal";
    default: return "undecidable";
  }
}

struct Finding {
  std::string channel;                 // "the record", or a level and channel of a hierarchy
  std::string indicator, estimates;
  Direction direction = Direction::rising;
  std::vector<std::string> sees, blind, assumes;
  std::string null;                    // the null in words
  Trend trend;
  Verdict verdict = Verdict::undecidable;
  std::vector<std::string> reasons;
};

struct Where {
  std::vector<std::string> units;
  std::vector<double> loading;         // the leading mode's weight in each unit (where it shows)
  std::vector<double> excitation;      // the left eigenvector's weight (what excites it)
  std::string sentence;
};

struct Channel {
  std::string name;
  double joint_p = 1;
  std::size_t coordinates = 0;
};

struct Assessment {
  std::string space;
  std::size_t states = 0, coordinates = 0;
  double step = 0;
  std::size_t window = 0;
  std::vector<Finding> findings;
  std::vector<Channel> channels;
  KindReading kind;
  std::optional<Wells> wells;
  std::optional<LikelihoodTest> likelihood;
  std::optional<Where> where;
  Verdict overall = Verdict::undecidable;
  std::vector<std::string> summary;
  std::vector<std::string> cannot_see;
  std::vector<std::string> refused;    // detectors the record or the space did not admit, with why
  std::optional<Rate> noise;           // one coordinate: the red-noise fit of the last window, and whether the noise is red
};

struct Options {
  TrendOptions trend;
  std::vector<RateEstimator> rates = {RateEstimator::ar1, RateEstimator::red};
  bool likelihood = true;
  LikelihoodOptions likelihood_options;
  bool wells = true;
  bool spatial = true;
  bool hierarchy = true;
  std::size_t channel_null_draws = 100, channel_alternative_draws = 60;
  std::size_t max_dimension = 60;      // beyond it the covariance indicators are not read at the base level
  double power_required = 0.8;
  double grid_required = 0.5;
};

namespace detail {

inline void judge(Finding& f, double alpha, double power_required, double grid_required, double remaining = 0.1) {
  const Trend& t = f.trend;
  std::size_t robust_significant = 0;
  for (const auto& [n, p] : t.robustness) robust_significant += p < alpha;
  const bool robust = t.robustness.empty() || 2 * robust_significant >= t.robustness.size();
  const bool gridded = t.grid.empty() || (t.grid_significant >= grid_required && !t.grid_flips);
  auto pct = [](double v) { return std::to_string(static_cast<int>(std::lround(100 * v))) + "%"; };
  const std::string alternative = "the declared alternative (the leading rate falling to " + pct(remaining) + " of its fitted value)";
  if (t.p < alpha) {
    if (robust && gridded) {
      f.verdict = Verdict::signal;
      f.reasons.push_back("the change from the first window to the last (" + std::string(t.positive ? "log ratio " : "difference ") +
                          std::to_string(t.change).substr(0, 6) + ", τ = " + std::to_string(t.tau).substr(0, 5) + ") is significant (p = " + std::to_string(t.p).substr(0, 5) +
                          ") under the fitted null" +
                          (t.robustness.empty() ? std::string() : ", under " + std::to_string(robust_significant) + " of " + std::to_string(t.robustness.size()) + " further nulls") +
                          (t.grid.empty() ? std::string() : " and in " + pct(t.grid_significant) + " of the windows and bandwidths"));
    } else {
      f.verdict = Verdict::undecidable;
      if (!robust) {
        std::string which;
        for (const auto& [n, p] : t.robustness)
          if (p >= alpha) which += (which.empty() ? "" : ", ") + n + " (p = " + std::to_string(p).substr(0, 5) + ")";
        f.reasons.push_back("significant under the fitted null but not under " + which);
      }
      if (!gridded)
        f.reasons.push_back(t.grid_flips ? "some choices of window and bandwidth give a significant change the other way"
                                         : "significant in only " + pct(t.grid_significant) + " of the windows and bandwidths");
    }
  } else if (f.direction != Direction::either && t.detection >= power_required && t.alternative_share <= 1 - power_required) {
    f.verdict = Verdict::no_signal;
    f.reasons.push_back("no significant change (p = " + std::to_string(t.p).substr(0, 5) + "), where the test would have detected " + alternative + " " +
                        pct(t.detection) + " of the time, and a change this small comes in only " + pct(t.alternative_share) + " of its records");
  } else if (f.direction != Direction::either && t.detection >= power_required) {
    f.verdict = Verdict::undecidable;
    f.reasons.push_back("no significant change (p = " + std::to_string(t.p).substr(0, 5) + "), but a change this large is common under " + alternative + " (" +
                        pct(t.alternative_share) + " of its records show less), so an absence is not shown");
  } else {
    f.verdict = Verdict::undecidable;
    f.reasons.push_back("no significant change (p = " + std::to_string(t.p).substr(0, 5) + "), but the test detects " + alternative + " only " +
                        pct(t.detection) + " of the time at this record's length, so an absence says little");
  }
}

inline std::vector<Finding> findings_of(const TrendBattery& b, const std::vector<Indicator>& ind, const std::string& channel, const Options& o) {
  std::vector<Finding> out;
  for (std::size_t k = 0; k < ind.size(); k++) {
    Finding f;
    f.channel = channel;
    f.indicator = ind[k].name;
    f.estimates = ind[k].estimates;
    f.direction = ind[k].direction;
    f.sees = ind[k].sees;
    f.blind = ind[k].blind;
    f.assumes = ind[k].assumes;
    f.null = b.null_described;
    f.trend = b.trends[k];
    judge(f, o.trend.alpha, o.power_required, o.grid_required, o.trend.remaining);
    out.push_back(f);
  }
  return out;
}

inline bool is_rate(const std::string& indicator) { return indicator.rfind("restoring rate", 0) == 0; }

// A record of derived coordinates as a charted record, for the channels of a hierarchy.
struct Derived {
  Rectangular x;
  StateSpace space;   // identity chart over the derived units, with their neighbourhood
  std::string name;
};

// The two channels of one level of a hierarchy, from the native states of an identity or composition chart.
inline std::vector<Derived> channels_of(const StateSpace& s, const Record& r, std::size_t level) {
  const auto& h = *s.hierarchy;
  const std::size_t n = s.units, P = s.parts, T = r.times.size();
  const std::size_t groups = h.size(level);
  std::vector<std::size_t> group_of(n);
  for (std::size_t u = 0; u < n; u++) group_of[u] = h.ancestor(0, u, level);
  std::vector<std::size_t> members(groups, 0);
  for (std::size_t u = 0; u < n; u++) members[group_of[u]]++;
  const bool extensive = s.kind == Kind::extensive || s.chart == Chart::composition;
  std::vector<Derived> out;

  // Totals (or means) of each group.
  Derived agg;
  agg.name = h.level(level).name + ": " + (extensive ? "totals" : "means");
  {
    StateSpace a;
    a.units = groups;
    a.parts = P;
    a.chart = s.chart == Chart::composition ? Chart::composition : Chart::identity;
    a.with_total = s.with_total;
    a.kind = s.kind;
    a.unit_names = h.level(level).units;
    a.part_names = s.part_names;
    a.edges = h.level(level).edges;
    Record ar;
    ar.times = r.times;
    for (std::size_t t = 0; t < T; t++) {
      std::vector<double> state(groups * P, 0);
      for (std::size_t u = 0; u < n; u++)
        for (std::size_t k = 0; k < P; k++) state[group_of[u] * P + k] += r.states[t][u * P + k];
      if (!extensive)
        for (std::size_t g = 0; g < groups; g++)
          for (std::size_t k = 0; k < P; k++) state[g * P + k] /= static_cast<double>(members[g]);
      ar.states.push_back(state);
    }
    const Charted c = chart(a, ar);
    agg.x = c.x;
    a.chart = Chart::identity;
    a.parts = a.coordinates_per_unit();
    a.with_total = false;
    StateSpace ident;
    ident.units = groups;
    ident.parts = c.x.cols / groups;
    ident.unit_names = h.level(level).units;
    ident.edges = h.level(level).edges;
    agg.space = ident;
  }
  out.push_back(agg);

  // Contrasts within each group: log shares centred within the group (extensive), or deviations from its mean (intensive).
  Derived within;
  within.name = h.level(level).name + ": " + (extensive ? "shares within" : "deviations within");
  {
    within.x = Rectangular(T, n * P);
    for (std::size_t t = 0; t < T; t++)
      for (std::size_t k = 0; k < P; k++) {
        std::vector<double> sum(groups, 0);
        for (std::size_t u = 0; u < n; u++) {
          const double v = r.states[t][u * P + k];
          if (extensive && !(v > 0))
            throw Refusal("a value of " + s.unit_name(u) + " at time " + std::to_string(r.times[t]) +
                          " is not positive; shares within a group are read in log-ratios, so a zero must be replaced by a declared value first");
          sum[group_of[u]] += extensive ? std::log(v) : v;
        }
        for (std::size_t u = 0; u < n; u++) {
          const double v = r.states[t][u * P + k];
          within.x(t, u * P + k) = (extensive ? std::log(v) : v) - sum[group_of[u]] / static_cast<double>(members[group_of[u]]);
        }
      }
    StateSpace ident;
    ident.units = n;
    ident.parts = P;
    ident.unit_names = s.unit_names;
    ident.edges = s.neighbourhood();
    within.space = ident;
  }
  out.push_back(within);
  return out;
}

inline std::string describe(const StateSpace& s) {
  std::string d = std::to_string(s.units) + (s.units == 1 ? " unit" : " units") + " of " + std::to_string(s.parts) +
                  (s.parts == 1 ? " part" : " parts") + ", read as " + name(s.chart);
  if (s.with_total) d += " with each unit's log total";
  if (!s.neighbourhood().empty()) d += ", with a neighbourhood of " + std::to_string(s.neighbourhood().size()) + " edges";
  if (s.hierarchy) d += ", on a hierarchy of " + std::to_string(s.hierarchy->levels()) + " levels";
  return d;
}

}  // namespace detail

inline Assessment assess(const StateSpace& space, const Record& record, const Options& o = {}) {
  check(space);
  check(record);
  const double dt = require_regular(record, "an assessment of slowing down");
  const Charted c = chart(space, record);
  const std::size_t T = c.x.rows, d = c.x.cols;
  Assessment a;
  a.space = detail::describe(space);
  a.states = T;
  a.coordinates = d;
  a.step = dt;
  if (T < 40) throw Refusal("an assessment needs a record of at least 40 states; this one has " + std::to_string(T));

  std::set<std::string> blind;
  auto run = [&](const Rectangular& x, const std::vector<Indicator>& ind, const std::string& channel, const TrendOptions& to) {
    const TrendBattery b = trends(x, c.times, dt, ind, to);
    a.window = std::max(a.window, b.window);
    a.channels.push_back({channel, b.joint_p, x.cols});
    auto f = detail::findings_of(b, ind, channel, o);
    for (const auto& one : f)
      for (const auto& s : one.blind) blind.insert(s);
    a.findings.insert(a.findings.end(), f.begin(), f.end());
    return b;
  };

  const Rectangular residual = detail::detrend(c.x, o.trend.bandwidth);
  const std::size_t window = std::clamp<std::size_t>(static_cast<std::size_t>(std::llround(o.trend.window * static_cast<double>(T))), 10, T - 3);
  // The record's lag-1 autocorrelation in its first window, averaged over coordinates.
  double first_ac = 0;
  std::size_t varying = 0;
  for (std::size_t j = 0; j < d; j++) {
    const auto w = detail::column_rows(residual, j, 0, window);
    if (tyche::descriptive::variance(w, 1) > 0) first_ac += tyche::time_series::acf(w, 1)[1], varying++;
  }
  first_ac = varying ? first_ac / static_cast<double>(varying) : 0;

  if (d == 1) {
    const auto ind = univariate_indicators(0, o.rates);
    run(c.x, ind, "the record", o.trend);
    // Kind: the dominant eigenvalue of an AR(1) or AR(2) in each window.
    std::vector<double> at;
    std::vector<Complex> path;
    const std::size_t stride = std::max<std::size_t>(1, (T - window + 1 + o.trend.max_windows - 1) / o.trend.max_windows);
    for (std::size_t end = window; end <= T; end += stride) {
      at.push_back(c.times[end - 1]);
      path.push_back(ar_dominant(detail::column_rows(residual, 0, end - window, end)).mu);
    }
    a.kind = classify(at, path, dt);
    if (o.likelihood) {
      try {
        auto lo = o.likelihood_options;
        lo.alpha = o.trend.alpha;
        a.likelihood = likelihood_test(c.times, column(c.x, 0), lo);
      } catch (const Refusal& e) {
        a.refused.push_back(std::string("the likelihood test: ") + e.what());
      }
    }
    if (o.wells) {
      try {
        a.wells = wells(detail::column_rows(residual, 0, T - window, T), dt);
      } catch (const Refusal& e) {
        a.refused.push_back(std::string("counting wells: ") + e.what());
      }
    }
  } else {
    std::vector<Indicator> ind;
    if (d <= o.max_dimension && window >= d + 3) {
      ind = multivariate_indicators();
    } else {
      a.refused.push_back("the covariance and propagator indicators at the base level: " + std::to_string(d) +
                          " coordinates against windows of " + std::to_string(window) + " states" +
                          (d > o.max_dimension ? " (more than " + std::to_string(o.max_dimension) + ")" : "") +
                          "; read the hierarchy's levels or fewer coordinates");
    }
    if (o.spatial && space.units >= 4) {
      const auto sp = spatial_indicators(space);
      ind.insert(ind.end(), sp.begin(), sp.end());
    }
    if (!ind.empty()) run(c.x, ind, "the record", o.trend);
    if (d <= o.max_dimension && window >= d + 3) {
      // Kind and place: the propagator's leading mode in each window and at the end.
      std::vector<double> at;
      std::vector<Complex> path;
      const std::size_t stride = std::max<std::size_t>(1, (T - window + 1 + o.trend.max_windows - 1) / o.trend.max_windows);
      for (std::size_t end = window; end <= T; end += stride) {
        at.push_back(c.times[end - 1]);
        path.push_back(propagator(rows(residual, end - window, end), 1e-4).spectrum.front());
      }
      a.kind = classify(at, path, dt);
      const Propagator P = propagator(rows(residual, T - window, T), 1e-4);
      const Leading L = leading(P, dt);
      Where w;
      w.units.resize(space.units);
      w.loading.assign(space.units, 0);
      w.excitation.assign(space.units, 0);
      for (std::size_t u = 0; u < space.units; u++) w.units[u] = space.unit_name(u);
      for (std::size_t j = 0; j < d; j++) {
        w.loading[c.unit_of[j]] += std::norm(L.right[j]);
        w.excitation[c.unit_of[j]] += std::norm(L.left[j]);
      }
      double sl = 0, se = 0;
      for (std::size_t u = 0; u < space.units; u++) sl += w.loading[u], se += w.excitation[u];
      for (std::size_t u = 0; u < space.units; u++) {
        w.loading[u] = sl > 0 ? w.loading[u] / sl : 0;
        w.excitation[u] = se > 0 ? w.excitation[u] / se : 0;
      }
      std::size_t top = 0;
      for (std::size_t u = 1; u < space.units; u++)
        if (w.loading[u] > w.loading[top]) top = u;
      w.sentence = "in the last window the leading mode (|μ| = " + std::to_string(std::abs(L.mode.mu)).substr(0, 5) + ") shows most in " +
                   space.unit_name(top) + " (" + std::to_string(static_cast<int>(std::lround(100 * w.loading[top]))) + "% of its weight)";
      a.where = w;
      if (o.wells) {
        // The record projected on its leading covariance direction.
        const auto dir = covariance_spectrum(covariance(rows(residual, T - window, T), true).c).direction;
        std::vector<double> proj(window, 0);
        for (std::size_t t = 0; t < window; t++)
          for (std::size_t j = 0; j < d; j++) proj[t] += dir[j] * residual(T - window + t, j);
        try {
          a.wells = wells(proj, dt);
        } catch (const Refusal& e) {
          a.refused.push_back(std::string("counting wells: ") + e.what());
        }
      }
    }
    if (o.hierarchy && space.hierarchy) {
      if (space.chart != Chart::identity && space.chart != Chart::composition) {
        a.refused.push_back("the hierarchy's channels: they are read for identity and composition charts");
      } else {
        TrendOptions to = o.trend;
        to.null_draws = o.channel_null_draws;
        to.alternative_draws = o.channel_alternative_draws;
        to.grid = false;
        to.robustness = {Null::phase};
        for (std::size_t level = 1; level < space.hierarchy->levels(); level++) {
          for (const auto& ch : detail::channels_of(space, record, level)) {
            std::vector<Indicator> ci;
            if (ch.x.cols <= o.max_dimension && window >= ch.x.cols + 3) ci = multivariate_indicators();
            if (ch.x.cols == 1) ci = univariate_indicators(0, {RateEstimator::ar1});
            if (o.spatial && ch.space.units >= 4) {
              const auto sp = spatial_indicators(ch.space);
              ci.insert(ci.end(), sp.begin(), sp.end());
            }
            if (ci.empty()) {
              a.refused.push_back("the channel " + ch.name + ": " + std::to_string(ch.x.cols) + " coordinates against windows of " + std::to_string(window));
              continue;
            }
            run(ch.x, ci, ch.name, to);
          }
        }
      }
    }
  }

  // Conditions across findings.
  const bool near_one = first_ac >= 0.95;
  bool rate_signal = false, rate_no_signal = false;
  for (const auto& f : a.findings)
    if (f.channel == "the record" && detail::is_rate(f.indicator)) {
      rate_signal |= f.verdict == Verdict::signal;
      rate_no_signal |= f.verdict == Verdict::no_signal;
    }
  // Is the noise red? In the last window, the red-noise fit leaving under three quarters of white noise's residual, its noise rate
  // resolved (θΔt < 3): declared bounds, provisional.
  bool red = false, red_rate_signal = false, red_in_battery = false;
  std::string red_said;
  if (d == 1) {
    for (const auto& f : a.findings)
      if (f.channel == "the record" && f.indicator == name(RateEstimator::red)) red_in_battery = true, red_rate_signal = f.verdict == Verdict::signal;
    if (red_in_battery) {
      try {
        const Rate r = red_rate(detail::column_rows(residual, 0, T - window, T), dt);
        a.noise = r;
        red = r.red_noise;
        red_said = "in the last window the noise is red (the red-noise autocorrelation leaves " + std::to_string(r.fit_white > 0 ? r.fit_red / r.fit_white : 1).substr(0, 4) +
                   " of white noise's residual, under the declared bound 0.75, with the noise's rate " + std::to_string(r.noise_rate).substr(0, 6) + " resolved at this sampling)";
      } catch (const std::exception&) {
      }
    }
  }
  for (auto& f : a.findings) {
    if (red && !red_rate_signal && f.verdict == Verdict::signal && f.channel == "the record" &&
        (f.indicator == "variance" || f.indicator == "lag-1 autocorrelation" || f.indicator == name(RateEstimator::ar1) ||
         f.indicator == name(RateEstimator::increment) || f.indicator == name(RateEstimator::drift) || f.indicator == "dominant eigenvalue modulus")) {
      f.verdict = Verdict::undecidable;
      f.reasons.push_back(red_said + ", and the restoring rate under red noise does not confirm the slowing: a reddening of the noise is not excluded");
    }
    if (near_one && f.verdict != Verdict::undecidable) {
      f.verdict = Verdict::undecidable;
      f.reasons.push_back("the record's lag-1 autocorrelation is already " + std::to_string(first_ac).substr(0, 5) +
                          " in its first window, so a window holds few independent states");
    }
    if (d == 1 && (f.indicator == "variance" || f.indicator == "lag-1 autocorrelation") && f.verdict == Verdict::signal && !rate_signal && rate_no_signal) {
      f.verdict = Verdict::undecidable;
      f.reasons.push_back("it rises while no estimate of the restoring rate falls: the noise may have grown or reddened, which is not slowing down");
    }
  }

  // The overall verdict.
  std::size_t signals = 0, directed = 0, absent = 0;
  for (const auto& f : a.findings) {
    if (f.direction == Direction::either) continue;
    directed++;
    signals += f.verdict == Verdict::signal;
    absent += f.verdict == Verdict::no_signal;
  }
  double joint = 1;
  for (const auto& ch : a.channels) joint = std::min(joint, ch.joint_p * static_cast<double>(a.channels.size()));
  if (signals > 0 && joint < o.trend.alpha && !near_one) a.overall = Verdict::signal;
  else if (directed > 0 && absent == directed) a.overall = Verdict::no_signal;
  else a.overall = Verdict::undecidable;

  auto pct = [](double v) { return std::to_string(static_cast<int>(std::lround(100 * v))) + "%"; };
  a.summary.push_back("A record of " + std::to_string(T) + " states of " + a.space + ", " + std::to_string(d) + " coordinates in its chart.");
  for (const auto& f : a.findings)
    if (f.verdict == Verdict::signal) a.summary.push_back(f.channel + ": " + f.indicator + " " + name(f.direction) + ", a signal (detection rate of the test " + pct(f.trend.detection) + ").");
  if (signals == 0) a.summary.push_back("No indicator gives a signal that survives its nulls and its grid.");
  if (!a.kind.path.empty()) a.summary.push_back(std::string("The dominant eigenvalue: ") + a.kind.sentence + "." + (a.kind.shape.empty() ? "" : " " + a.kind.shape + "."));
  if (a.where) a.summary.push_back("Where: " + a.where->sentence + ".");
  if (a.likelihood) a.summary.push_back("Likelihood: " + a.likelihood->tipping_sentence + ".");
  if (a.wells && a.wells->wells >= 2) a.summary.push_back("States: " + a.wells->sentence + ".");
  a.summary.push_back(std::string("Overall: ") + name(a.overall) + ", the joint test across the directed indicators giving p = " + std::to_string(joint).substr(0, 5) + ".");
  a.cannot_see.assign(blind.begin(), blind.end());
  return a;
}

}  // namespace kinesis::transitions
