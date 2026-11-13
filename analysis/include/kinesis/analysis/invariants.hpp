// What a flow or a map does over a long time: its exponents, its equilibria, its sections, its attractors as a
// parameter moves, and its basins.
//
// Each needs the system as a presentation: the field and, for the exponents and the equilibria, its Jacobian,
// which the terms give exactly (systems/dynamics.hpp).
//
// Lyapunov's exponents. The tangent vectors are carried by the exact Jacobian and orthonormalised as in physis
// (physis/dynamics/flows.hpp: Benettin, Galgani, Giorgilli and Strelcyn 1980). Their sum is the mean divergence
// of the field along the path, which for a field of constant divergence is known exactly: the check of the
// method. Kaplan and Yorke's dimension follows.
//
// Equilibria. Newton's method from seeds spread over a box, with the exact Jacobian; the roots found are told
// apart to a tolerance and classified by the eigenvalues of the Jacobian (for a flow, by the signs of their real
// parts; for a map, by their moduli against 1). What is found is every root the seeds converged to: an
// equilibrium no seed reached is not found, and the answer does not say there is none.
//
// Sections. A Poincaré section is the set of points where a path crosses a surface g(x) = 0 in one sense; a
// crossing is located as an event is (systems/simulate.hpp).
//
// Bifurcation diagrams. For each value of a parameter the system is run past a transient and a summary of the
// attractor is kept: the iterates of a map, or the maxima of an observable of a flow (the points of its section
// g = dφ/dt = 0, falling). The state is carried from one value to the next, so that a branch is followed.
//
// Basins. Each point of a grid of initial states is run; the states it ends near are clustered into attractors,
// and the grid is labelled. The share of a box that leads to an attractor is its basin stability (Menck,
// Heitzig, Marwan and Kurths 2013, Nature Phys. 9, 89–92); from N independent draws its standard error is
// √(S(1 − S)/N), whatever the dimension.
#pragma once

#include <kinesis/systems/simulate.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>
#include <calx/random/stream.hpp>
#include <physis/dynamics/flows.hpp>

#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace kinesis::analysis {

using systems::Dynamics;
using systems::State;

// ---- Lyapunov's exponents --------------------------------------------------------------------------------

struct Exponents {
  std::vector<double> exponents;   // λ_1 ≥ … ≥ λ_k, per unit of time or per iteration
  double sum = 0;
  double dimension = 0;            // Kaplan and Yorke's
  double time = 0;
  std::vector<double> times;
  std::vector<std::vector<double>> running;
  std::vector<double> error;       // the standard error of each, from the means over blocks; empty where not recorded
  State final_state;
};

// The error of a finite time. The time averaged over is cut into blocks, the mean rate over each is read from the
// running record (T_b λ(T_b) − T_{b−1} λ(T_{b−1}) over T_b − T_{b−1}), and the standard error of their mean is
// returned. The first block is left out: in it the frame is still turning toward its directions. The blocks are
// taken as independent, which they are when each is long against the time the path remembers; where they are not,
// the error is understated. Empty if fewer than two records fall in a block.
inline std::vector<double> block_errors(const std::vector<double>& times, const std::vector<std::vector<double>>& running, std::size_t blocks = 10) {
  const std::size_t m = times.size();
  if (blocks < 3 || m < 2 * blocks || running.size() != m) return {};
  const std::size_t k = running.front().size();
  std::vector<double> out(k, 0.0);
  for (std::size_t j = 0; j < k; j++) {
    std::vector<double> means;
    for (std::size_t b = 2; b <= blocks; b++) {
      const std::size_t lo = (b - 1) * m / blocks - 1, hi = b * m / blocks - 1;
      means.push_back((times[hi] * running[hi][j] - times[lo] * running[lo][j]) / (times[hi] - times[lo]));
    }
    double mean = 0, var = 0;
    for (double v : means) mean += v / static_cast<double>(means.size());
    for (double v : means) var += (v - mean) * (v - mean) / static_cast<double>(means.size() - 1);
    out[j] = std::sqrt(var / static_cast<double>(means.size()));
  }
  return out;
}

inline Exponents lyapunov(const Dynamics& d, const State& x0, physis::dynamics::LyapunovOptions o = {}) {
  if (!d.events.empty()) throw kernel::Refusal("Lyapunov's exponents are computed here for systems without events; " + d.name + " has " + std::to_string(d.events.size()));
  physis::dynamics::LyapunovSpectrum spectrum;
  if (d.discrete) {
    const std::size_t n = d.dimension(), k = o.count == 0 ? n : std::min(o.count, n);
    State x = x0, next(n), image(n);
    for (std::size_t s = 0; s < static_cast<std::size_t>(o.transient); s++) d.field(static_cast<double>(s), x, next), x.swap(next);
    std::vector<State> v = physis::dynamics::detail::initial_frame(n, k);
    std::vector<double> sums(k, 0.0);
    const std::size_t steps = static_cast<std::size_t>(o.time), every = std::max<std::size_t>(o.renormalise_every, 1);
    std::size_t renormalisations = 0;
    for (std::size_t s = 1; s <= steps; s++) {
      for (std::size_t j = 0; j < k; j++) d.jacobian_vector(static_cast<double>(s), x, v[j], image), v[j].swap(image);
      d.field(static_cast<double>(s), x, next), x.swap(next);
      if (s % every == 0 || s == steps) {
        const auto r = physis::dynamics::detail::orthonormalise(v);
        for (std::size_t j = 0; j < k; j++) sums[j] += std::log(r[j]);
        if (o.record_every && ++renormalisations % o.record_every == 0) {
          spectrum.times.push_back(static_cast<double>(s));
          spectrum.running.emplace_back(k);
          for (std::size_t j = 0; j < k; j++) spectrum.running.back()[j] = sums[j] / static_cast<double>(s);
        }
      }
    }
    spectrum.time = static_cast<double>(steps);
    spectrum.exponents.resize(k);
    for (std::size_t j = 0; j < k; j++) spectrum.exponents[j] = sums[j] / spectrum.time;
    spectrum.final_state = x;
  } else {
    const auto f = [&](double t, const State& x, State& dx) { d.field(t, x, dx); };
    const auto jv = [&](double t, const State& x, const State& v) {
      State out(v.size());
      d.jacobian_vector(t, x, v, out);
      return out;
    };
    spectrum = physis::dynamics::lyapunov_spectrum(f, jv, x0, o);
  }
  Exponents e;
  e.exponents = spectrum.exponents, e.time = spectrum.time, e.times = spectrum.times, e.running = spectrum.running, e.final_state = spectrum.final_state;
  for (double l : e.exponents) e.sum += l;
  e.error = block_errors(e.times, e.running);
  e.dimension = physis::dynamics::kaplan_yorke(e.exponents);
  return e;
}

// ---- equilibria ------------------------------------------------------------------------------------------

struct Equilibrium {
  State x;
  std::vector<std::complex<double>> eigenvalues;
  bool stable = false, hyperbolic = true;
  std::string kind;   // stable node, stable focus, saddle, unstable node, unstable focus, centre; for a map: attracting, repelling, saddle
  double residual = 0;
};

struct Box {
  std::vector<double> lower, upper;
};

namespace detail {
inline std::string classify(const std::vector<std::complex<double>>& eig, bool discrete, bool& stable, bool& hyperbolic) {
  const double tolerance = 1e-9;
  std::size_t in = 0, out = 0, on = 0;
  bool complex_pair = false;
  for (const auto& l : eig) {
    const double measure = discrete ? std::abs(l) - 1 : l.real();
    if (measure < -tolerance) in++;
    else if (measure > tolerance) out++;
    else on++;
    if (std::fabs(l.imag()) > tolerance) complex_pair = true;
  }
  hyperbolic = on == 0;
  stable = out == 0 && on == 0;
  if (!hyperbolic) return out == 0 && in == 0 ? "centre" : "not hyperbolic";
  if (discrete) return out == 0 ? "attracting" : in == 0 ? "repelling" : "saddle";
  if (out == 0) return complex_pair ? "stable focus" : "stable node";
  if (in == 0) return complex_pair ? "unstable focus" : "unstable node";
  return "saddle";
}
}  // namespace detail

// The equilibrium Newton's method reaches from x, or nothing.
inline bool newton(const Dynamics& d, State& x, double tolerance = 1e-12, std::size_t most = 60) {
  const std::size_t n = d.dimension();
  State f(n);
  calx::linalg::Matrix J(n);
  for (std::size_t it = 0; it < most; it++) {
    d.field(0, x, f);
    if (d.discrete)
      for (std::size_t i = 0; i < n; i++) f[i] -= x[i];
    double norm = 0;
    for (double v : f) norm = std::max(norm, std::fabs(v));
    if (!std::isfinite(norm)) return false;
    if (norm < tolerance) return true;
    d.jacobian(0, x, J);
    if (d.discrete)
      for (std::size_t i = 0; i < n; i++) J(i, i) -= 1;
    calx::linalg::Matrix rhs(n);
    for (std::size_t i = 0; i < n; i++) rhs(i, 0) = -f[i];
    try {
      const calx::linalg::Matrix step = calx::linalg::solve(J, rhs);
      for (std::size_t i = 0; i < n; i++) x[i] += step(i, 0);
    } catch (const std::exception&) {
      return false;
    }
  }
  return false;
}

inline Equilibrium equilibrium_at(const Dynamics& d, const State& x) {
  Equilibrium e;
  e.x = x;
  calx::linalg::Matrix J(d.dimension());
  d.jacobian(0, x, J);
  e.eigenvalues = calx::linalg::eigenvalues(J);
  e.kind = detail::classify(e.eigenvalues, d.discrete, e.stable, e.hyperbolic);
  State f(d.dimension());
  d.field(0, x, f);
  for (std::size_t i = 0; i < f.size(); i++) e.residual = std::max(e.residual, std::fabs(d.discrete ? f[i] - x[i] : f[i]));
  return e;
}

// The equilibria (fixed points, for a map) that Newton's method reaches from `seeds` points of the box.
inline std::vector<Equilibrium> equilibria(const Dynamics& d, const Box& box, std::size_t seeds = 400, std::uint64_t seed = 1, double distinct = 1e-6) {
  const std::size_t n = d.dimension();
  if (box.lower.size() != n || box.upper.size() != n) throw kernel::Refusal("the box has a bound for each of the " + std::to_string(n) + " components of " + d.name);
  calx::random::Stream r(seed, 0x65717569u, 0, 0);  // 'equi'
  std::vector<Equilibrium> found;
  for (std::size_t k = 0; k < seeds; k++) {
    State x(n);
    for (std::size_t i = 0; i < n; i++) x[i] = r.uniform(box.lower[i], box.upper[i]);
    if (!newton(d, x)) continue;
    bool fresh = true;
    for (const auto& e : found) {
      double gap = 0;
      for (std::size_t i = 0; i < n; i++) gap = std::max(gap, std::fabs(e.x[i] - x[i]));
      if (gap < distinct * (1 + std::fabs(x[0]))) fresh = false;
    }
    if (fresh) found.push_back(equilibrium_at(d, x));
  }
  std::sort(found.begin(), found.end(), [](const Equilibrium& a, const Equilibrium& b) { return a.x < b.x; });
  return found;
}

// ---- sections --------------------------------------------------------------------------------------------

struct Crossing {
  double t = 0;
  State x;
};

// The points where a path crosses g(x) = 0 with g rising (or falling), g a term of the system's names.
inline std::vector<Crossing> section(const Dynamics& d, State x, const terms::TermPtr& g, bool rising, const systems::Options& o) {
  if (d.discrete) throw kernel::Refusal("a section is of a flow");
  const terms::Program guard = d.compile(g);
  std::vector<Crossing> out;
  double t = o.t0;
  for (double s = 0; s < o.transient - 0.5 * o.step; s += o.step) systems::advance(d, t, x, o.step, o);
  t = o.t0;
  const auto f = [&](double s, const State& y, State& dy) { d.field(s, y, dy); };
  const std::size_t count = static_cast<std::size_t>(std::ceil((o.t1 - o.t0) / o.step - 1e-9));
  double before = d.value(guard, t, x);
  for (std::size_t k = 0; k < count; k++) {
    State y = x;
    double ty = t;
    systems::advance(d, ty, y, o.step, o);
    const double after = d.value(guard, ty, y);
    if (rising ? before < 0 && after >= 0 : before > 0 && after <= 0) {
      double lo = 0, hi = o.step;
      State at = y;
      while (hi - lo > o.event_tolerance) {
        const double mid = 0.5 * (lo + hi);
        State m = x;
        physis::dynamics::rk4_step(f, t, m, mid);
        const double gm = d.value(guard, t + mid, m);
        if (rising ? gm >= 0 : gm <= 0) hi = mid, at.swap(m);
        else lo = mid;
      }
      out.push_back({t + hi, at});
    }
    x.swap(y), t = ty, before = after;
  }
  return out;
}

// ---- bifurcation diagrams --------------------------------------------------------------------------------

struct Diagram {
  std::string parameter, observable;
  std::vector<double> values;                  // of the parameter
  std::vector<std::vector<double>> points;     // for each value, what was kept of the attractor
};

struct DiagramOptions {
  double transient = 500;      // iterations of a map, or time of a flow
  std::size_t keep = 200;      // iterates kept, or the most maxima kept
  double time = 200;           // a flow's time after the transient
  double step = 0.01;
  std::size_t component = 0;   // which state is kept (a map), or whose maxima (a flow)
  bool carry = true;           // the state at one value starts the next
};

inline Diagram bifurcation_diagram(Dynamics d, const std::string& parameter, const std::vector<double>& values, State x0, const DiagramOptions& o = {}) {
  Diagram out;
  out.parameter = parameter, out.observable = d.states.at(o.component), out.values = values;
  State x = x0;
  for (double v : values) {
    d.set_param(parameter, v);
    if (!o.carry) x = x0;
    std::vector<double> kept;
    if (d.discrete) {
      double n = 0;
      for (double k = 0; k < o.transient; k += 1) systems::iterate(d, n, x);
      for (std::size_t k = 0; k < o.keep; k++) {
        systems::iterate(d, n, x);
        kept.push_back(x[o.component]);
      }
    } else {
      systems::Options so;
      so.step = o.step, so.transient = o.transient, so.t1 = o.time;
      // The maxima of the component: where its derivative, the field's component, falls through zero.
      const auto crossings = section(d, x, d.field_terms[o.component], false, so);
      for (const auto& c : crossings) {
        if (kept.size() >= o.keep) break;
        kept.push_back(c.x[o.component]);
      }
      if (!crossings.empty()) x = crossings.back().x;
      else {
        so.every = 1000000;
        x = systems::simulate(d, x, so).x.back();
        kept.push_back(x[o.component]);  // no maximum: the path has settled
      }
    }
    for (double w : x)
      if (!std::isfinite(w)) x = x0;  // the path escaped at this value: begin again at the next
    out.points.push_back(std::move(kept));
  }
  return out;
}

// ---- basins ----------------------------------------------------------------------------------------------

struct Basins {
  std::size_t nx = 0, ny = 0;
  std::vector<double> xs, ys;                 // the grid of the two components varied
  std::vector<int> label;                     // label[iy * nx + ix]: the attractor reached, or −1 when the path left the box of escape
  std::vector<State> attractors;              // a state on each attractor found
  std::vector<double> share;                  // of the grid, for each attractor; the attractors are in its order
  std::vector<double> speed;                  // how fast each end state still moves: zero at an equilibrium
  std::vector<std::string> kind;              // of an end at rest, the kind of equilibrium it is; else "moving"
};

struct BasinOptions {
  std::size_t component_x = 0, component_y = 1;
  double time = 200, step = 0.01;             // a flow's; a map runs `time` iterations
  double distinct = 1e-2;                     // two end states within it are one attractor
  double escape = 1e6;                        // a path beyond it has left
  double rest = 1e-6;                         // an end state that moves more slowly is at rest
};

// The attractor each state of a grid leads to. The other components start at `base`. End states are compared at
// the final time only, so an attractor that is not a point (a cycle, a strange set) is one label when its points
// are within `distinct` of one already seen at that time, and may be split otherwise: the labels are of where
// paths end, and say so.
inline Basins basins(const Dynamics& d, const State& base, const std::vector<double>& xs, const std::vector<double>& ys, const BasinOptions& o = {}) {
  Basins b;
  b.nx = xs.size(), b.ny = ys.size(), b.xs = xs, b.ys = ys;
  b.label.assign(b.nx * b.ny, -1);
  systems::Options so;
  so.t1 = o.time, so.step = o.step, so.every = 100000000;
  for (std::size_t iy = 0; iy < b.ny; iy++)
    for (std::size_t ix = 0; ix < b.nx; ix++) {
      State x = base;
      x[o.component_x] = xs[ix], x[o.component_y] = ys[iy];
      State end;
      try {
        end = systems::simulate(d, x, so).x.back();
      } catch (const kernel::Refusal&) {
        continue;
      }
      bool escaped = false;
      for (double v : end) escaped = escaped || !(std::fabs(v) < o.escape);
      if (escaped) continue;
      int found = -1;
      for (std::size_t a = 0; a < b.attractors.size() && found < 0; a++) {
        double gap = 0;
        for (std::size_t i = 0; i < end.size(); i++) gap = std::max(gap, std::fabs(end[i] - b.attractors[a][i]));
        if (gap < o.distinct) found = static_cast<int>(a);
      }
      if (found < 0) found = static_cast<int>(b.attractors.size()), b.attractors.push_back(end);
      b.label[iy * b.nx + ix] = found;
    }
  // A path that began near the boundary of a basin is slow to arrive, and ends apart from the others: each end
  // state found is run on for four times as long, and those that then end together are one attractor.
  if (b.attractors.size() > 1) {
    systems::Options longer = so;
    longer.t1 = 4 * o.time;
    std::vector<State> settled;
    std::vector<int> to(b.attractors.size(), -1);
    for (std::size_t a = 0; a < b.attractors.size(); a++) {
      State end = b.attractors[a];
      try {
        end = systems::simulate(d, end, longer).x.back();
      } catch (const kernel::Refusal&) {
      }
      for (std::size_t k = 0; k < settled.size() && to[a] < 0; k++) {
        double gap = 0;
        for (std::size_t i = 0; i < end.size(); i++) gap = std::max(gap, std::fabs(end[i] - settled[k][i]));
        if (gap < o.distinct) to[a] = static_cast<int>(k);
      }
      if (to[a] < 0) to[a] = static_cast<int>(settled.size()), settled.push_back(end);
    }
    for (int& l : b.label)
      if (l >= 0) l = to[static_cast<std::size_t>(l)];
    b.attractors = settled;
  }
  b.share.assign(b.attractors.size(), 0.0);
  for (int l : b.label)
    if (l >= 0) b.share[static_cast<std::size_t>(l)] += 1.0 / static_cast<double>(b.label.size());
  // In the order of their shares, the greatest first, so that the first label is the attractor most paths reach.
  std::vector<std::size_t> order(b.attractors.size());
  for (std::size_t a = 0; a < order.size(); a++) order[a] = a;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t p, std::size_t q) { return b.share[p] > b.share[q]; });
  std::vector<int> place(order.size());
  std::vector<State> attractors;
  std::vector<double> share;
  for (std::size_t k = 0; k < order.size(); k++) place[order[k]] = static_cast<int>(k), attractors.push_back(b.attractors[order[k]]), share.push_back(b.share[order[k]]);
  for (int& l : b.label)
    if (l >= 0) l = place[static_cast<std::size_t>(l)];
  b.attractors = std::move(attractors), b.share = std::move(share);
  // How fast each end state still moves: the size of the field there (a flow), or of one step (a map). At rest it
  // is a point attractor; moving, it is a point of a cycle or a strange set, or a path that has not arrived (one
  // that began on the boundary between two basins creeps toward a saddle and never leaves it).
  State image(d.dimension());
  for (const State& a : b.attractors) {
    d.field(o.time, a, image);
    double speed = 0;
    for (std::size_t i = 0; i < image.size(); i++) speed = std::max(speed, std::fabs(d.discrete ? image[i] - a[i] : image[i]));
    b.speed.push_back(speed);
    // An end at rest is an equilibrium, and is named by its linearisation: an attractor, or a saddle, which only
    // the paths that began on its stable set reach (on a grid, a line of cells at most).
    b.kind.push_back(speed < o.rest ? equilibrium_at(d, a).kind : "moving");
  }
  return b;
}

// Basin stability: the share of N states drawn uniformly from a box that end within `distinct` of `attractor`,
// with its standard error.
struct BasinStability {
  double share = 0, error = 0;
  std::size_t draws = 0;
};
inline BasinStability basin_stability(const Dynamics& d, const State& attractor, const Box& box, std::size_t draws = 500, const BasinOptions& o = {}, std::uint64_t seed = 1) {
  calx::random::Stream r(seed, 0x62617369u, 0, 0);  // 'basi'
  systems::Options so;
  so.t1 = o.time, so.step = o.step, so.every = 100000000;
  std::size_t hits = 0;
  for (std::size_t k = 0; k < draws; k++) {
    State x(d.dimension());
    for (std::size_t i = 0; i < x.size(); i++) x[i] = r.uniform(box.lower[i], box.upper[i]);
    try {
      const State end = systems::simulate(d, x, so).x.back();
      double gap = 0;
      for (std::size_t i = 0; i < end.size(); i++) gap = std::max(gap, std::fabs(end[i] - attractor[i]));
      if (gap < o.distinct) hits++;
    } catch (const kernel::Refusal&) {
    }
  }
  BasinStability s;
  s.draws = draws, s.share = static_cast<double>(hits) / static_cast<double>(draws);
  s.error = std::sqrt(s.share * (1 - s.share) / static_cast<double>(draws));
  return s;
}

}  // namespace kinesis::analysis
