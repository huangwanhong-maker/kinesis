// Trajectories of flows and maps, with their events.
//
// A flow is stepped by the classical Runge–Kutta method of physis (physis/dynamics/flows.hpp) at a fixed step h,
// its global error O(h⁴) where the field is smooth. A map is iterated.
//
// Events of a flow. Before each step the guards g_e are read. A guard that is negative before a step and not
// after it (for a rising event; the reverse for a falling one) has crossed within the step. The first crossing is
// found by bisection on the length of the step, each trial a Runge–Kutta step from the state before, until the
// bracket is shorter than `event_tolerance`; the state is taken there, the resets are applied, the occurrence is
// recorded with the states before and after, and the rest of the step follows from the new state. The time of an
// event is therefore known to the tolerance and to the integrator's error, not to the step. A guard that is zero
// or beyond before a step does not fire: an event is a crossing, not a level.
//
// Events of a map are levels: after an iteration, an event whose relation holds applies its resets.
#pragma once

#include <kinesis/systems/dynamics.hpp>

#include <physis/dynamics/flows.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

namespace kinesis::systems {

struct Occurrence {
  double t = 0;
  std::size_t event = 0;
  State before, after;
};

struct Trajectory {
  std::vector<std::string> names;
  std::vector<double> t;
  std::vector<State> x;
  std::vector<Occurrence> events;
  std::size_t steps = 0;

  std::vector<double> column(std::size_t i) const {
    std::vector<double> out(x.size());
    for (std::size_t k = 0; k < x.size(); k++) out[k] = x[k][i];
    return out;
  }
};

struct Options {
  double t0 = 0, t1 = 100;       // a flow's times; a map runs from iteration t0 to t1
  double step = 0.01;            // a flow's step
  std::size_t every = 1;         // steps between recorded states
  double transient = 0;          // time (iterations) run before recording begins
  double event_tolerance = 1e-10;
  std::size_t max_events = 1000000;
};

namespace detail {
inline bool crossed(const Dynamics::Event& e, double before, double after) { return e.rising ? before < 0 && after >= 0 : before > 0 && after <= 0; }
inline bool holds(const Dynamics::Event& e, double g) { return e.rising ? g >= 0 : g <= 0; }
}  // namespace detail

// One step of a flow of length h from (t, x), through whatever events occur in it. Returns the events' count.
inline std::size_t advance(const Dynamics& d, double& t, State& x, double h, const Options& o, std::vector<Occurrence>* record = nullptr) {
  const auto f = [&](double s, const State& y, State& dy) { d.field(s, y, dy); };
  const std::size_t ne = d.events.size();
  if (ne == 0) {
    physis::dynamics::rk4_step(f, t, x, h);
    t += h;
    return 0;
  }
  std::size_t fired = 0;
  double remaining = h;
  std::vector<double> g0(ne);
  State y;
  while (remaining > 0) {
    for (std::size_t e = 0; e < ne; e++) g0[e] = d.guard(e, t, x);
    y = x;
    physis::dynamics::rk4_step(f, t, y, remaining);
    const auto any = [&](double s, const State& at) {
      for (std::size_t e = 0; e < ne; e++)
        if (detail::crossed(d.events[e], g0[e], d.guard(e, s, at))) return static_cast<long>(e);
      return -1L;
    };
    if (any(t + remaining, y) < 0) {
      x.swap(y), t += remaining;
      break;
    }
    double lo = 0, hi = remaining;
    State at_hi = y;
    while (hi - lo > o.event_tolerance) {
      const double mid = 0.5 * (lo + hi);
      State m = x;
      physis::dynamics::rk4_step(f, t, m, mid);
      if (any(t + mid, m) >= 0) hi = mid, at_hi.swap(m);
      else lo = mid;
    }
    const std::size_t e = static_cast<std::size_t>(any(t + hi, at_hi));
    Occurrence occ{t + hi, e, at_hi, d.reset(e, t + hi, at_hi)};
    x = occ.after, t += hi, remaining -= hi;
    if (record) record->push_back(std::move(occ));
    if (++fired >= o.max_events) throw kernel::Refusal("more than " + std::to_string(o.max_events) + " events in one step of " + d.name);
  }
  return fired;
}

// One iteration of a map, with its events as levels.
inline void iterate(const Dynamics& d, double& n, State& x, std::vector<Occurrence>* record = nullptr) {
  State next(x.size());
  d.field(n, x, next);
  n += 1;
  for (std::size_t e = 0; e < d.events.size(); e++) {
    if (!detail::holds(d.events[e], d.guard(e, n, next))) continue;
    State after = d.reset(e, n, next);
    if (record) record->push_back({n, e, next, after});
    next.swap(after);
  }
  x.swap(next);
}

// The trajectory from x0: a flow from t0 to t1 at the step given, or a map from iteration t0 to t1.
inline Trajectory simulate(const Dynamics& d, State x, const Options& o = {}) {
  if (x.size() != d.dimension()) throw kernel::Refusal(d.name + " has " + std::to_string(d.dimension()) + " components; " + std::to_string(x.size()) + " were given");
  Trajectory out;
  out.names = d.states;
  double t = o.t0;
  const std::size_t every = std::max<std::size_t>(o.every, 1);
  if (d.discrete) {
    for (double k = 0; k < o.transient; k += 1) iterate(d, t, x);
    t = o.t0;
    out.t.push_back(t), out.x.push_back(x);
    const std::size_t count = static_cast<std::size_t>(std::max(0.0, o.t1 - o.t0));
    for (std::size_t k = 1; k <= count; k++) {
      iterate(d, t, x, &out.events);
      if (k % every == 0 || k == count) out.t.push_back(t), out.x.push_back(x);
    }
    out.steps = count;
    return out;
  }
  if (!(o.step > 0)) throw kernel::Refusal("the step of a flow is positive");
  for (double s = 0; s < o.transient - 0.5 * o.step; s += o.step) advance(d, t, x, o.step, o);
  t = o.t0;
  out.t.push_back(t), out.x.push_back(x);
  const std::size_t count = static_cast<std::size_t>(std::ceil((o.t1 - o.t0) / o.step - 1e-9));
  for (std::size_t k = 1; k <= count; k++) {
    const double h = std::min(o.step, o.t1 - t);
    advance(d, t, x, h, o, &out.events);
    if (k % every == 0 || k == count) out.t.push_back(t), out.x.push_back(x);
    for (double v : x)
      if (!std::isfinite(v)) throw kernel::Refusal(d.name + " left the finite numbers at t = " + terms::number_text(t));
  }
  out.steps = count;
  return out;
}

// The initial state of a system: its named values, zero elsewhere.
inline State initial(const Dynamics& d, const std::vector<std::pair<std::string, double>>& values) {
  State x(d.dimension(), 0.0);
  for (const auto& [n, v] : values) {
    const int i = d.state_index(n);
    if (i < 0) throw kernel::Refusal("'" + n + "' is not a component of " + d.name);
    x[static_cast<std::size_t>(i)] = v;
  }
  return x;
}

}  // namespace kinesis::systems
