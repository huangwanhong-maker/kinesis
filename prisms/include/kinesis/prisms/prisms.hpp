// Prisms: what a record of positions identifies of a mover's path under a bound on its speed (decision 0006).
//
// A record holds fixes of one mover in the plane: at the instant t_k it was within r_k of (x_k, y_k). Its histories are the
// continuous paths through every fix's disc at the fix's instant that never go faster than a declared speed v. The set of places
// such a path can occupy at a time t between two fixes is the space-time prism of Hägerstrand (1970, What about people in regional
// science?, Papers of the Regional Science Association 24, 7–21), formalised by Lenntorp (1976) and by Miller (2005, A measurement
// theory for time geography, Geographical Analysis 37, 17–45): the intersection of the future cone of the earlier fix and the past
// cone of the later one, two discs at each instant,
//
//     D(x_k, r_k + v (t − t_k)) ∩ D(x_{k+1}, r_{k+1} + v (t_{k+1} − t)),       t_k ≤ t ≤ t_{k+1},
//
// the Alexandrov interval of the cone structure of speed v. What the record identifies is a set (kinesis's invariant 5), and the
// questions asked of it are answered so.
//
// CONSISTENCY. The record admits a history iff points q_k ∈ D(x_k, r_k) can be chosen with |q_{k+1} − q_k| ≤ v (t_{k+1} − t_k).
// It is necessary that every pair of fixes satisfies |x_k − x_j| ≤ v |t_k − t_j| + r_j + r_k; with exact fixes (every r_k = 0)
// the consecutive pairs suffice. With radii the pairwise conditions do not suffice (three discs may meet pairwise and not
// together), and a witness is sought: cyclic projections onto the discs and the bounds on each step (Bregman 1965; Gubin, Polyak
// and Raik 1967), whose limit lies in the intersection when it is not empty; a witness found is checked exactly. So the answer is
// consistent with a witness, inconsistent with the pair that breaks it, or undecided with what was tried; never a guess.
//
// THE SLICE at an instant is the two discs of the fixes that bracket it (one cone before the first fix or after the last). Their
// intersection holds every place a history can be at that instant, exactly when the fixes are exact, and as an outer bound when
// they have radii (the other fixes can only narrow it).
//
// PRESENCE in a region (a disc during an interval) is the set of values in, out over the admissible histories:
//   - "in" is impossible (certified out, an alibi) when no slice within the interval meets the region. Between consecutive fixes
//     the function δ(t) = min_p max_i (|p − c_i| − ρ_i(t)), over the slice's discs and the region, is convex in t (the ρ_i are
//     affine there, and partial minimisation of a jointly convex function is convex), so its least value on each piece is found
//     by Brent's method (calx) and the region is out of reach iff it is positive everywhere. δ at an instant is exact: the
//     minimiser of the largest excess over finitely many discs is a point where one, two or three of them are active (a
//     centre; a point on the segment between two centres; an additively weighted Apollonius point of three), and every such
//     candidate is evaluated;
//   - "in" is possible when a history through the region at the best instant is found and checked (the region taken as a fix);
//   - "out" is admissible when the witness history of the record stays out of the region during the interval, and impossible
//     when the disc of some fix within the interval lies in the region (every history is then inside at that instant).
//   What is not decided is said: an answer that lists one value and is not complete leaves the other open.
//
// THE TIMES at which a mover may be in a region: on each piece, the set {t : δ(t) ≤ 0} is an interval (δ is convex there), found
// about the least value by bisection. With exact fixes they are exact; with radii they are an outer bound.
//
// Units are the caller's: metres and seconds, or any length and time with v in their ratio. Times are absolute numbers (seconds
// since an epoch); the plane is a local frame (a map projection, or a city's metres). The geometry is Euclidean: the shortest
// path between two places is a straight line; a street network or a river that must be crossed by a bridge is not seen.
#pragma once

#include <kinesis/kernel/answer.hpp>

#include <calx/optimisation/minimise.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace kinesis::prisms {

using kernel::Identified;
using kernel::Refusal;

inline constexpr double inf = std::numeric_limits<double>::infinity();

// The mover was within `radius` of (x, y) at the instant t.
struct Fix {
  std::string id;
  double x = 0, y = 0, t = 0, radius = 0;
};

struct Disc {
  double x = 0, y = 0, r = 0;
};

// A disc of the plane during an interval of time.
struct Region {
  double x = 0, y = 0, radius = 0, from = 0, until = 0;
};

// ---- The deepest point of a finite set of discs --------------------------------------------------------------------------------

// min over p of max_i (|p − c_i| − r_i), with the point: the discs share a point iff the depth is ≤ 0, and −depth is then how far
// inside all of them the deepest point lies.
struct Deepest {
  double depth = inf;
  double x = 0, y = 0;
};

namespace detail {

inline double excess(const std::vector<Disc>& d, double x, double y) {
  double m = -inf;
  for (const auto& c : d) m = std::max(m, std::hypot(x - c.x, y - c.y) - c.r);
  return m;
}

inline std::string metres(double v) {
  char b[40];
  std::snprintf(b, sizeof b, "%.4g", v);
  return b;
}

}  // namespace detail

inline Deepest deepest(const std::vector<Disc>& d) {
  if (d.empty()) throw Refusal("kinesis::prisms: the deepest point of no discs");
  for (const auto& c : d)
    if (!std::isfinite(c.x) || !std::isfinite(c.y) || !(c.r >= 0)) throw Refusal("kinesis::prisms: a disc has a finite centre and a radius that is not negative");
  Deepest best;
  const auto consider = [&](double x, double y) {
    if (!std::isfinite(x) || !std::isfinite(y)) return;
    const double v = detail::excess(d, x, y);
    if (v < best.depth) best = {v, x, y};
  };
  const std::size_t n = d.size();
  // One disc active: its centre.
  for (const auto& c : d) consider(c.x, c.y);
  // Two active: on the segment between the centres, where |p − c_i| − r_i = |p − c_j| − r_j.
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      const double dx = d[j].x - d[i].x, dy = d[j].y - d[i].y, dist = std::hypot(dx, dy);
      if (dist == 0) continue;
      const double s = (dist + d[i].r - d[j].r) / 2;
      if (s >= 0 && s <= dist) consider(d[i].x + s * dx / dist, d[i].y + s * dy / dist);
    }
  // Three active: |p − c_m| = r_m + δ for m = i, j, k. Subtracting the squared equation of i from those of j and k gives
  // 2 (c_m − c_i)·p + 2 (r_m − r_i) δ = |c_m|² − |c_i|² − r_m² + r_i², linear in p for each δ; the first equation then fixes δ.
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++)
      for (std::size_t k = j + 1; k < n; k++) {
        const Disc &a = d[i], &b = d[j], &c = d[k];
        const double a11 = 2 * (b.x - a.x), a12 = 2 * (b.y - a.y), a21 = 2 * (c.x - a.x), a22 = 2 * (c.y - a.y);
        const double b1 = b.x * b.x + b.y * b.y - a.x * a.x - a.y * a.y - b.r * b.r + a.r * a.r;
        const double b2 = c.x * c.x + c.y * c.y - a.x * a.x - a.y * a.y - c.r * c.r + a.r * a.r;
        const double e1 = -2 * (b.r - a.r), e2 = -2 * (c.r - a.r);
        const double det = a11 * a22 - a12 * a21;
        const double scale = std::max({std::fabs(a11 * a22), std::fabs(a12 * a21), 1e-300});
        if (std::fabs(det) <= 1e-12 * scale) continue;  // collinear centres: two active discs suffice
        const double Ax = (b1 * a22 - a12 * b2) / det, Bx = (e1 * a22 - a12 * e2) / det;
        const double Ay = (a11 * b2 - b1 * a21) / det, By = (a11 * e2 - e1 * a21) / det;
        const double ux = Ax - a.x, uy = Ay - a.y;
        const double qa = Bx * Bx + By * By - 1, qb = 2 * (ux * Bx + uy * By) - 2 * a.r, qc = ux * ux + uy * uy - a.r * a.r;
        std::vector<double> roots;
        if (std::fabs(qa) < 1e-14) {
          if (qb != 0) roots.push_back(-qc / qb);
        } else {
          const double disc = qb * qb - 4 * qa * qc;
          if (disc >= 0) {
            const double sq = std::sqrt(disc), q = -0.5 * (qb + (qb >= 0 ? sq : -sq));
            roots.push_back(q / qa);
            if (q != 0) roots.push_back(qc / q);
          }
        }
        for (double delta : roots) consider(Ax + Bx * delta, Ay + By * delta);
      }
  return best;
}

// ---- Fixes in order -------------------------------------------------------------------------------------------------------------

namespace detail {

inline std::vector<Fix> ordered(std::vector<Fix> f) {
  for (const auto& x : f)
    if (!std::isfinite(x.x) || !std::isfinite(x.y) || !std::isfinite(x.t) || !(x.radius >= 0)) throw Refusal("kinesis::prisms: the fix " + x.id + " has a finite place and time and a radius that is not negative");
  std::stable_sort(f.begin(), f.end(), [](const Fix& a, const Fix& b) { return a.t < b.t; });
  return f;
}

inline void require_speed(double v) {
  if (!(v > 0) || !std::isfinite(v)) throw Refusal("kinesis::prisms: a bound on speed is positive and finite");
}

inline std::string name_of(const Fix& f, std::size_t k) { return f.id.empty() ? "fix " + std::to_string(k + 1) : f.id; }

}  // namespace detail

// ---- Consistency ----------------------------------------------------------------------------------------------------------------

// A history the record admits, as its places at the fixes' instants: joined by straight steps it never goes faster than v.
struct History {
  std::vector<Fix> at;   // the fixes in order of time, each moved to the place the history passes through (radius 0)
};

struct Consistency {
  std::string status = "undecided";       // consistent, inconsistent, undecided
  double margin = inf;                    // the least slack of the pairwise condition, in the unit of length (negative: broken)
  std::vector<std::pair<std::string, std::string>> broken;   // pairs of fixes no history can join
  std::optional<History> witness;
  std::string reason;
  bool consistent() const { return status == "consistent"; }
  bool inconsistent() const { return status == "inconsistent"; }
};

namespace detail {

// Cyclic projections onto the fixes' discs and onto the bounds on each step; a witness if the result satisfies them all within `tol`.
inline std::optional<History> project(const std::vector<Fix>& f, double v, double tol, int iterations = 600) {
  const std::size_t n = f.size();
  std::vector<double> qx(n), qy(n);
  for (std::size_t k = 0; k < n; k++) qx[k] = f[k].x, qy[k] = f[k].y;
  const auto satisfied = [&]() {
    for (std::size_t k = 0; k < n; k++)
      if (std::hypot(qx[k] - f[k].x, qy[k] - f[k].y) > f[k].radius + tol) return false;
    for (std::size_t k = 0; k + 1 < n; k++)
      if (std::hypot(qx[k + 1] - qx[k], qy[k + 1] - qy[k]) > v * (f[k + 1].t - f[k].t) + tol) return false;
    return true;
  };
  for (int it = 0; it <= iterations; it++) {
    if (satisfied()) {
      History h;
      for (std::size_t k = 0; k < n; k++) h.at.push_back({f[k].id, qx[k], qy[k], f[k].t, 0});
      return h;
    }
    for (std::size_t k = 0; k < n; k++) {
      const double dx = qx[k] - f[k].x, dy = qy[k] - f[k].y, d = std::hypot(dx, dy);
      if (d > f[k].radius) qx[k] = f[k].x + dx * f[k].radius / d, qy[k] = f[k].y + dy * f[k].radius / d;
    }
    for (std::size_t k = 0; k + 1 < n; k++) {
      const double s = v * (f[k + 1].t - f[k].t), dx = qx[k + 1] - qx[k], dy = qy[k + 1] - qy[k], d = std::hypot(dx, dy);
      if (d <= s) continue;
      const double cut = (d - s) / 2 / d;  // each end moves half the excess toward the other
      qx[k] += dx * cut, qy[k] += dy * cut, qx[k + 1] -= dx * cut, qy[k + 1] -= dy * cut;
    }
  }
  return std::nullopt;
}

}  // namespace detail

// Whether the record admits a history under the speed v, with a witness or the pair that breaks it. `tol` is a length.
inline Consistency consistent(const std::vector<Fix>& fixes, double v, double tol = 1e-6) {
  detail::require_speed(v);
  Consistency out;
  const auto f = detail::ordered(fixes);
  if (f.empty()) {
    out.status = "consistent", out.witness = History{}, out.reason = "a record without fixes admits every history";
    return out;
  }
  for (std::size_t j = 0; j < f.size(); j++)
    for (std::size_t k = j + 1; k < f.size(); k++) {
      const double slack = v * (f[k].t - f[j].t) + f[j].radius + f[k].radius - std::hypot(f[k].x - f[j].x, f[k].y - f[j].y);
      out.margin = std::min(out.margin, slack);
      if (slack < -tol) out.broken.push_back({detail::name_of(f[j], j), detail::name_of(f[k], k)});
    }
  if (!out.broken.empty()) {
    out.status = "inconsistent";
    const auto& [a, b] = out.broken.front();
    out.reason = "no path at the speed bound joins " + a + " and " + b + " (short by " + detail::metres(-out.margin) + " at the worst pair)";
    return out;
  }
  out.witness = detail::project(f, v, tol);
  if (out.witness) {
    out.status = "consistent";
    out.reason = "a history through every fix was found and checked; the least slack between two fixes is " + detail::metres(out.margin);
  } else {
    out.reason = "every pair of fixes can be joined, but no history through all of them was found by cyclic projections: it is not decided";
  }
  return out;
}

// ---- The slice at an instant -----------------------------------------------------------------------------------------------------

// The discs whose intersection holds every place a history can occupy at t: the bracketing fixes' cones (exact for exact fixes, an
// outer bound otherwise). Before the first fix and after the last, one cone.
inline std::vector<Disc> slice(const std::vector<Fix>& fixes, double v, double t) {
  detail::require_speed(v);
  const auto f = detail::ordered(fixes);
  if (f.empty()) throw Refusal("kinesis::prisms: a record without fixes places the mover nowhere in particular");
  if (t <= f.front().t) return {{f.front().x, f.front().y, f.front().radius + v * (f.front().t - t)}};
  if (t >= f.back().t) return {{f.back().x, f.back().y, f.back().radius + v * (t - f.back().t)}};
  std::size_t k = 0;
  while (k + 1 < f.size() && f[k + 1].t < t) k++;
  return {{f[k].x, f[k].y, f[k].radius + v * (t - f[k].t)}, {f[k + 1].x, f[k + 1].y, f[k + 1].radius + v * (f[k + 1].t - t)}};
}

// ---- When a mover may be in a region ---------------------------------------------------------------------------------------------

namespace detail {

// The pieces of [from, until] on which the slice's radii are affine: cut at the fixes' instants.
inline std::vector<std::pair<double, double>> pieces(const std::vector<Fix>& f, double from, double until) {
  std::vector<double> cuts{from};
  for (const auto& x : f)
    if (x.t > from && x.t < until) cuts.push_back(x.t);
  cuts.push_back(until);
  std::vector<std::pair<double, double>> out;
  for (std::size_t k = 0; k + 1 < cuts.size(); k++)
    if (cuts[k + 1] >= cuts[k]) out.push_back({cuts[k], cuts[k + 1]});
  return out;
}

inline Deepest depth_at(const std::vector<Fix>& f, double v, const Region& r, double t) {
  auto d = slice(f, v, t);
  d.push_back({r.x, r.y, r.radius});
  return deepest(d);
}

}  // namespace detail

struct Window {
  double from = 0, until = 0;
};

struct Reach {
  std::vector<Window> windows;   // the instants within the region's interval at which the slice meets the region
  double least = inf;            // the least depth over the interval (≤ 0: reachable; > 0: short by that much)
  double when = 0, x = 0, y = 0; // the instant and place of the least depth
};

// When the mover may be in the region during its interval: the instants at which its slice meets the region, and the deepest
// meeting. `tol` is a length (meeting within it counts) and the instants are found to within `tolerance_time`.
inline Reach reach(const std::vector<Fix>& fixes, double v, const Region& region, double tol = 1e-6, double tolerance_time = 1e-3) {
  detail::require_speed(v);
  if (!(region.until >= region.from) || !(region.radius >= 0)) throw Refusal("kinesis::prisms: a region has a radius that is not negative and an interval that ends no earlier than it begins");
  const auto f = detail::ordered(fixes);
  if (f.empty()) throw Refusal("kinesis::prisms: a record without fixes places the mover nowhere in particular");
  Reach out;
  const auto depth = [&](double t) { return detail::depth_at(f, v, region, t).depth; };
  for (const auto& [a, b] : detail::pieces(f, region.from, region.until)) {
    double tm = a, dm = depth(a);
    if (b > a) {
      const auto m = calx::optimisation::brent_minimise(depth, a, b, tolerance_time / 4);
      for (double t : {a, b, m.x}) {
        const double d = depth(t);
        if (d < dm) dm = d, tm = t;
      }
    }
    const auto at = detail::depth_at(f, v, region, tm);
    if (at.depth < out.least) out.least = at.depth, out.when = tm, out.x = at.x, out.y = at.y;
    // Out of reach on this piece only when short by more than Brent's tolerance can hide (the depth moves at most v per unit time).
    if (dm > tol + v * tolerance_time) continue;
    // The sublevel set of a convex function on [a, b] is an interval about tm: bisect for its ends.
    const auto edge = [&](double inside, double outside) {
      if (depth(outside) <= tol + v * tolerance_time) return outside;
      for (int it = 0; it < 200 && std::fabs(outside - inside) > tolerance_time; it++) {
        const double mid = (inside + outside) / 2;
        (depth(mid) <= tol + v * tolerance_time ? inside : outside) = mid;
      }
      return inside;
    };
    const Window w{edge(tm, a), edge(tm, b)};
    if (!out.windows.empty() && w.from <= out.windows.back().until + tolerance_time) out.windows.back().until = std::max(out.windows.back().until, w.until);
    else out.windows.push_back(w);
  }
  return out;
}

// ---- Presence -------------------------------------------------------------------------------------------------------------------

struct Presence {
  Identified<std::string> identified;   // values among "in" and "out"
  Reach reach;
  std::optional<History> in, out;       // a history through the region, and one that stays out of it, when found
  std::string necessary_by;             // the fix whose disc lies in the region, when every history is inside
};

namespace detail {

// Whether the piecewise-linear history stays farther than the region's radius from its centre throughout [from, until].
inline bool stays_out(const History& h, const Region& r, double tol) {
  const auto& a = h.at;
  if (a.empty()) return false;
  const auto outside = [&](double x, double y) { return std::hypot(x - r.x, y - r.y) > r.radius + tol; };
  // Before the first fix and after the last the history rests where it is.
  if (r.from < a.front().t && !outside(a.front().x, a.front().y)) return false;
  if (r.until > a.back().t && !outside(a.back().x, a.back().y)) return false;
  for (std::size_t k = 0; k + 1 < a.size(); k++) {
    const double t0 = std::max(a[k].t, r.from), t1 = std::min(a[k + 1].t, r.until);
    if (t1 < t0) continue;
    const double span = a[k + 1].t - a[k].t;
    const auto at = [&](double t) {
      const double s = span > 0 ? (t - a[k].t) / span : 0;
      return std::pair<double, double>{a[k].x + s * (a[k + 1].x - a[k].x), a[k].y + s * (a[k + 1].y - a[k].y)};
    };
    const auto [x0, y0] = at(t0);
    const auto [x1, y1] = at(t1);
    // The least distance from the centre to the segment from (x0, y0) to (x1, y1).
    const double dx = x1 - x0, dy = y1 - y0, len2 = dx * dx + dy * dy;
    const double s = len2 > 0 ? std::clamp(((r.x - x0) * dx + (r.y - y0) * dy) / len2, 0.0, 1.0) : 0;
    if (!outside(x0 + s * dx, y0 + s * dy)) return false;
  }
  return true;
}

}  // namespace detail

// Whether the mover was in the region during its interval, over every history the record admits under the speed v.
inline Presence presence(const std::vector<Fix>& fixes, double v, const Region& region, double tol = 1e-6) {
  Presence out;
  auto& id = out.identified;
  const std::string bound = "speed at most " + detail::metres(v);
  id.assumptions = {bound, "the plane is Euclidean: paths are straight lines, with no obstacle", "each fix within its stated radius at its stated instant"};
  const auto f = detail::ordered(fixes);
  const auto c = consistent(f, v, tol);
  if (c.inconsistent()) {
    id.reason = "the record admits no history: " + c.reason;
    return out;
  }
  out.reach = reach(f, v, region, tol);
  // "in": impossible when no slice meets the region; possible when a history through it at the deepest instant is found.
  bool in_open = false, out_open = false;
  if (out.reach.windows.empty()) {
    // certified out of reach
  } else {
    auto g = f;
    g.push_back({"the region", region.x, region.y, out.reach.when, region.radius});
    const auto with = consistent(g, v, tol);
    if (with.consistent()) out.in = with.witness, id.add("in");
    else in_open = true;
  }
  // "out": impossible when some fix within the interval lies in the region; admissible when the record's witness stays out.
  for (std::size_t k = 0; k < f.size(); k++)
    if (f[k].t >= region.from && f[k].t <= region.until && std::hypot(f[k].x - region.x, f[k].y - region.y) + f[k].radius <= region.radius + tol) out.necessary_by = detail::name_of(f[k], k);
  if (out.necessary_by.empty()) {
    if (c.witness && detail::stays_out(*c.witness, region, tol)) out.out = c.witness, id.add("out");
    else out_open = true;
  }
  id.complete = !in_open && !out_open;
  const auto when = [&]() {
    std::string s;
    for (const auto& w : out.reach.windows) s += (s.empty() ? "" : ", ") + detail::metres(w.from) + " to " + detail::metres(w.until);
    return s;
  };
  if (out.reach.windows.empty()) {
    id.reason = "no history reaches the region during its interval: at best it falls short by " + detail::metres(out.reach.least);
    if (out_open) id.reason += " (and the record's own witness was not found)";
  } else if (id.certified() && id.values.front() == "in") {
    id.reason = "every history is in the region: the disc of " + out.necessary_by + " lies within it";
  } else if (id.complete) {
    id.reason = "a history through the region and one that stays out of it were both found; it may be there at " + when();
  } else {
    id.reason = std::string(in_open ? "the region is within the outer bound at " + when() + ", but no history through it was found" : "a history through the region was found (it may be there at " + when() + ")") +
                (out_open ? "; whether a history stays out of it was not decided" : "");
  }
  return out;
}

// The place of a history at t, between its fixes in straight steps (resting before the first and after the last).
inline std::pair<double, double> position(const History& h, double t) {
  const auto& a = h.at;
  if (a.empty()) throw Refusal("kinesis::prisms: an empty history has no place");
  if (t <= a.front().t) return {a.front().x, a.front().y};
  if (t >= a.back().t) return {a.back().x, a.back().y};
  std::size_t k = 0;
  while (k + 1 < a.size() && a[k + 1].t < t) k++;
  const double span = a[k + 1].t - a[k].t, s = span > 0 ? (t - a[k].t) / span : 0;
  return {a[k].x + s * (a[k + 1].x - a[k].x), a[k].y + s * (a[k + 1].y - a[k].y)};
}

}  // namespace kinesis::prisms
