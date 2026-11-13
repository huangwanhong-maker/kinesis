// Braids read from records of planar trajectories, and what such a record identifies.
//
// The setting is that of Huang (2026), "Rethinking trajectory-based credibility under symbol inflation:
// possibilities for building epistemic infrastructure through topological methods with braid group theory as an
// example" (a draft; cited below as the draft, by its section and equation numbers). A record holds the positions
// of n labelled points of the plane at finitely many times. A history is a continuous, collision-free motion of
// the points through the recorded positions; its class is a morphism of the fundamental groupoid of the
// configuration space (the draft, §3.1). A record does not determine it: between two samples a point may turn
// about another any number of times (its Proposition 6.6). What a record identifies is relative to assumptions.
//
// The reading. The piecewise-linear interpolant of a record is one history. Its braid is read by projection on
// a line (Thiffeault 2010, Chaos 20, 017516, §II; the draft, §3.2): the points are ordered by their projection,
// and each time two neighbours in that order exchange places a letter σ_i^{±1} is written, i the position of the
// left one, the sign positive when the left one passes below (the counterclockwise half-exchange of
// morphe/algebra/braid.hpp). Configurations whose projections are distinct form contractible cells (the top
// cells of Fox and Neuwirth 1962), so that a path between two such configurations has a braid word that is well
// defined once the line is chosen: the convention for endpoints the draft asks for (its §3.1; Guillemot and
// Lairez 2026, arXiv:2601.23073, use the same). A record whose first or last configuration has two equal
// projections, or in which three points cross at one time, is not generic for the line and is refused.
//
// The certificate. Assume each point i moves no faster than V_i and is observed within η. Let r be the observed
// position of j relative to i at the two ends of an interval of length h. If
//
//     |r_k| + |r_{k+1}| − (V_i + V_j) h − 8η > 0                                                        (*)
//
// for every interval and every pair, then every history the assumptions admit is homotopic, in the
// configuration space, to the interpolant, and the braid read is the braid of every one of them.
// *Proof.* The true relative position stays in the ellipse with foci at the true r_k, r_{k+1} and string
// (V_i + V_j) h, which is convex; (*) with |r̃ − r| ≤ 2η puts the origin outside it, so the straight-line homotopy
// from the true path to the true interpolant never meets the origin, for all pairs at once. The observed and the
// true interpolants differ by at most 2η in each relative position, and (*) keeps the observed chord farther than
// 2η from the origin, so the straight-line homotopy between them is collision-free too. □ With η = 0 and one
// bound V, (*) is implied by the draft's Corollary 6.4 (δ̂ > VΔ) and is weaker; for two points it is sharp: where
// it fails strictly another class is admissible. The endpoints of the two interpolants lie in one cell when the
// projections at the first and last samples are farther apart than 2η, which is asked as well.
// This certificate was derived for this library; the draft's coarser margin is reported beside it.
//
// The identified windings of a pair. For two points with exact samples and a least separation ρ, the classes of
// the relative path on an interval are its total angles Θ = θ₀ + 2πm. The shortest path of total angle Θ that
// keeps |r| ≥ ρ is the chord when it has that angle and clears the disc, and otherwise two tangents and an arc,
//
//     √(a² − ρ²) + √(b² − ρ²) + ρ (|Θ| − arccos(ρ/a) − arccos(ρ/b)),   a = |r_k|, b = |r_{k+1}|,
//
// and a class is admissible exactly when this length is at most (V_i + V_j) h. The winding of a closed relative
// path is the sum over the intervals, so its identified set (the draft's Eq. 20) is a sum of finite sets of
// integers, computed exactly. With ρ = 0 the tangent path has length a + b whatever the angle, so that one
// admissible turn admits every number of turns: the set is then infinite, and is reported as not complete.
#pragma once

#include <kinesis/kernel/answer.hpp>

#include <morphe/algebra/braid.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <set>
#include <string>
#include <vector>

namespace kinesis::braids {

using kernel::Identified;
using kernel::Refusal;
using morphe::algebra::Braid;

namespace detail {
// A number as a sentence writes it: five significant digits, no trailing zeros.
inline std::string written(double v) {
  char buffer[40];
  std::snprintf(buffer, sizeof buffer, "%.5g", v);
  return buffer;
}
}  // namespace detail

struct Point {
  double x = 0, y = 0;
};

// The positions of n labelled points at increasing times: at[k][i] is point i at times[k].
struct Samples {
  std::vector<double> times;
  std::vector<std::vector<Point>> at;
  std::vector<std::string> labels;  // or empty

  std::size_t points() const { return at.empty() ? 0 : at.front().size(); }
  void validate() const {
    if (times.size() != at.size() || times.size() < 2) throw Refusal("a record of trajectories holds at least two times, with the positions at each");
    for (std::size_t k = 0; k < at.size(); k++) {
      if (at[k].size() != points()) throw Refusal("a record of trajectories holds the same points at every time: a point that arrives or leaves begins another segment");
      if (k > 0 && !(times[k] > times[k - 1])) throw Refusal("the times of a record increase");
    }
  }
  // Every `stride`-th sample, and the last.
  Samples thinned(std::size_t stride) const {
    Samples s;
    s.labels = labels;
    for (std::size_t k = 0; k < times.size(); k += stride) s.times.push_back(times[k]), s.at.push_back(at[k]);
    if (s.times.back() != times.back()) s.times.push_back(times.back()), s.at.push_back(at.back());
    return s;
  }
};

struct Crossing {
  double time = 0;
  int letter = 0;             // ±i
  std::size_t left = 0, right = 0;  // the labels of the two points, the left one first
};

struct Reading {
  Braid braid;
  std::vector<Crossing> crossings;
  std::vector<std::size_t> first, last;  // the points in the order of their projection, at the first and last times
};

// The braid of the piecewise-linear interpolant, read on the line at `angle` to the x-axis.
inline Reading read(const Samples& s, double angle = 0) {
  s.validate();
  const std::size_t n = s.points();
  const double c = std::cos(angle), sn = std::sin(angle);
  const auto along = [&](const Point& p) { return c * p.x + sn * p.y; };
  const auto across = [&](const Point& p) { return -sn * p.x + c * p.y; };
  std::vector<std::size_t> order(n);
  for (std::size_t i = 0; i < n; i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return along(s.at[0][a]) < along(s.at[0][b]); });
  for (std::size_t i = 0; i + 1 < n; i++)
    if (along(s.at[0][order[i]]) == along(s.at[0][order[i + 1]])) throw Refusal("two points have the same projection at the first time: the record is not generic for this line");
  Reading r;
  r.braid = Braid(std::max<std::size_t>(n, 1));
  r.first = order;
  struct Pending {
    double fraction;
    std::size_t a, b;
  };
  for (std::size_t k = 0; k + 1 < s.times.size(); k++) {
    const auto& p = s.at[k];
    const auto& q = s.at[k + 1];
    std::vector<Pending> pending;
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t b = a + 1; b < n; b++) {
        const double d0 = along(p[b]) - along(p[a]), d1 = along(q[b]) - along(q[a]);
        // The order of two points at a sample is the sign of the difference of their projections, a difference of
        // exactly zero counted as positive: one rule at both ends of every interval, so that an exchange that
        // falls on a sample is read once (at the end of the interval before it or the start of the one after),
        // and a touch that does not pass is read as two letters that cancel.
        if ((d0 < 0) == (d1 < 0)) continue;  // the order is kept
        pending.push_back({d0 / (d0 - d1), a, b});
      }
    std::sort(pending.begin(), pending.end(), [](const Pending& x, const Pending& y) { return x.fraction < y.fraction; });
    for (const Pending& e : pending) {
      const auto ia = std::find(order.begin(), order.end(), e.a) - order.begin(), ib = std::find(order.begin(), order.end(), e.b) - order.begin();
      const auto lo = std::min(ia, ib), hi = std::max(ia, ib);
      if (hi - lo != 1) throw Refusal("three points cross at one time near t = " + detail::written(s.times[k]) + ": the record is not generic for this line");
      const std::size_t left = order[static_cast<std::size_t>(lo)], right = order[static_cast<std::size_t>(hi)];
      const auto where = [&](std::size_t i) {
        return Point{p[i].x + e.fraction * (q[i].x - p[i].x), p[i].y + e.fraction * (q[i].y - p[i].y)};
      };
      const double dy = across(where(left)) - across(where(right));
      if (dy == 0) throw Refusal("two points meet near t = " + detail::written(s.times[k]) + ": a collision");
      const int letter = (dy < 0 ? 1 : -1) * static_cast<int>(lo + 1);  // the left point below: counterclockwise
      r.braid.word.push_back(letter);
      r.crossings.push_back({s.times[k] + e.fraction * (s.times[k + 1] - s.times[k]), letter, left, right});
      std::swap(order[static_cast<std::size_t>(lo)], order[static_cast<std::size_t>(hi)]);
    }
  }
  r.last = order;
  const auto& end = s.at.back();
  for (std::size_t i = 0; i + 1 < n; i++)
    if (along(end[order[i]]) == along(end[order[i + 1]])) throw Refusal("two points have the same projection at the last time: the record is not generic for this line");
  return r;
}

// What is assumed of the histories a record admits.
struct Assumptions {
  std::vector<double> speed;   // V_i, one for each point, or one for all
  double error = 0;            // η: each observed position is within it of the true one
  double separation = 0;       // ρ: the points are never closer than it
  double angle = 0;            // the line of projection

  double speed_of(std::size_t i) const {
    if (speed.empty()) throw Refusal("a certificate needs a bound on the speed of the points");
    return speed.size() == 1 ? speed[0] : speed.at(i);
  }
  std::vector<std::string> stated() const {
    std::vector<std::string> out;
    if (!speed.empty()) out.push_back(speed.size() == 1 ? "each point moves no faster than " + detail::written(speed[0]) : "each point moves no faster than its stated bound");
    out.push_back(error == 0 ? "the positions are exact" : "each position is within " + detail::written(error) + " of the true one");
    if (separation > 0) out.push_back("two points are never closer than " + detail::written(separation));
    out.push_back("the labels are right");
    return out;
  }
};

struct Certificate {
  bool certified = false;
  double margin = 0;                 // the least of (*) over the intervals and pairs
  std::size_t interval = 0, a = 0, b = 0;  // where it is least
  double coarse_margin = 0;          // the draft's: δ̂ − VΔ − 2η, with V the largest bound and Δ the longest interval
  bool endpoints = true;             // the first and last configurations are in one cell whatever the error
  std::size_t open = 0;              // the (interval, pair)s where (*) fails
  std::string reason;
};

namespace detail {
inline double norm(double x, double y) { return std::sqrt(x * x + y * y); }
// The distance from the origin to the segment from (ax, ay) to (bx, by).
inline double segment_distance(double ax, double ay, double bx, double by) {
  const double dx = bx - ax, dy = by - ay, l2 = dx * dx + dy * dy;
  if (l2 == 0) return norm(ax, ay);
  const double u = std::clamp(-(ax * dx + ay * dy) / l2, 0.0, 1.0);
  return norm(ax + u * dx, ay + u * dy);
}
}  // namespace detail

// Whether every admissible history has the braid of the interpolant.
inline Certificate certify(const Samples& s, const Assumptions& A) {
  s.validate();
  const std::size_t n = s.points();
  Certificate c;
  c.margin = INFINITY;
  double clearance = INFINITY, longest = 0, fastest = 0;
  for (std::size_t i = 0; i < n; i++) fastest = std::max(fastest, A.speed_of(i));
  for (std::size_t k = 0; k + 1 < s.times.size(); k++) {
    const double h = s.times[k + 1] - s.times[k];
    longest = std::max(longest, h);
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t b = a + 1; b < n; b++) {
        const double x0 = s.at[k][b].x - s.at[k][a].x, y0 = s.at[k][b].y - s.at[k][a].y;
        const double x1 = s.at[k + 1][b].x - s.at[k + 1][a].x, y1 = s.at[k + 1][b].y - s.at[k + 1][a].y;
        const double m = detail::norm(x0, y0) + detail::norm(x1, y1) - (A.speed_of(a) + A.speed_of(b)) * h - 8 * A.error;
        if (!(m > 0)) c.open++;
        if (m < c.margin) c.margin = m, c.interval = k, c.a = a, c.b = b;
        clearance = std::min(clearance, detail::segment_distance(x0, y0, x1, y1));
      }
  }
  c.coarse_margin = clearance - fastest * longest - 2 * A.error;
  const double cs = std::cos(A.angle), sn = std::sin(A.angle);
  for (const auto* end : {&s.at.front(), &s.at.back()})
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t b = a + 1; b < n; b++)
        if (std::fabs(cs * ((*end)[b].x - (*end)[a].x) + sn * ((*end)[b].y - (*end)[a].y)) <= 2 * A.error) c.endpoints = false;
  c.certified = c.open == 0 && c.endpoints;
  if (c.certified) c.reason = "every admissible history is homotopic to the interpolant: the least margin is " + detail::written(c.margin);
  else if (c.open) c.reason = std::to_string(c.open) + " interval-and-pair(s) admit another class; the first near t = " + detail::written(s.times[c.interval]) + " between points " + std::to_string(c.a) + " and " + std::to_string(c.b);
  else c.reason = "the order of the points at an end of the record is within the error of observation";
  return c;
}

// The least length of a relative path from (ax, ay) to (bx, by) of total angle `theta` that keeps |r| ≥ rho.
inline double least_length(double ax, double ay, double bx, double by, double theta, double rho) {
  const double a = detail::norm(ax, ay), b = detail::norm(bx, by);
  const double chord_angle = std::atan2(ax * by - ay * bx, ax * bx + ay * by);
  if (std::fabs(theta - chord_angle) < 1e-12 && detail::segment_distance(ax, ay, bx, by) >= rho) return detail::norm(bx - ax, by - ay);
  if (rho <= 0) return a + b;  // an infimum: by the origin
  const double arc = std::fabs(theta) - std::acos(std::min(1.0, rho / a)) - std::acos(std::min(1.0, rho / b));
  return std::sqrt(std::max(0.0, a * a - rho * rho)) + std::sqrt(std::max(0.0, b * b - rho * rho)) + rho * std::max(0.0, arc);
}

// The windings of point b about point a that the record admits, for a closed relative path with exact samples:
// the draft's identified-value set of the pair's winding number (its Eqs. 15 and 20).
inline Identified<long> pair_winding(const Samples& s, std::size_t a, std::size_t b, const Assumptions& A, long most = 1000) {
  s.validate();
  if (A.error > 0) throw Refusal("the windings a record admits are computed for exact positions; with an error of observation only the certificate is given");
  const auto rel = [&](std::size_t k) { return Point{s.at[k][b].x - s.at[k][a].x, s.at[k][b].y - s.at[k][a].y}; };
  const Point first = rel(0), last = rel(s.times.size() - 1);
  if (detail::norm(first.x - last.x, first.y - last.y) > 1e-9 * (1 + detail::norm(first.x, first.y)))
    throw Refusal("a winding number is an integer for a relative path that closes; this one does not");
  Identified<long> out;
  out.assumptions = A.stated();
  const double two_pi = 2 * std::numbers::pi;
  double chords = 0;            // the total angle of the interpolant
  std::set<long> sums{0};       // the sums of the extra turns over the intervals so far
  for (std::size_t k = 0; k + 1 < s.times.size(); k++) {
    const Point p = rel(k), q = rel(k + 1);
    const double budget = (A.speed_of(a) + A.speed_of(b)) * (s.times[k + 1] - s.times[k]);
    const double theta0 = std::atan2(p.x * q.y - p.y * q.x, p.x * q.x + p.y * q.y);
    chords += theta0;
    // The extra turns m admitted on this interval. The least length grows with |m|, so the search stops at the
    // first |m| neither sign of which is admitted. With no least separation only a few are listed.
    std::vector<long> extra;
    const long limit = A.separation > 0 ? most : 2;
    const auto admitted = [&](long m) { return least_length(p.x, p.y, q.x, q.y, theta0 + two_pi * static_cast<double>(m), A.separation) <= budget; };
    if (admitted(0)) extra.push_back(0);
    for (long m = 1; m <= limit; m++) {
      const bool plus = admitted(m), minus = admitted(-m);
      if (plus) extra.push_back(m);
      if (minus) extra.push_back(-m);
      if (!plus && !minus) break;
      if (m == limit) out.complete = false;
    }
    if (extra.empty()) {  // not even the chord: the assumptions contradict the record
      out.values.clear();
      out.reason = "the samples near t = " + detail::written(s.times[k]) + " are farther apart than the speeds allow";
      return out;
    }
    if (A.separation <= 0 && extra.size() > 1) out.complete = false;  // one admissible turn admits every number of turns
    std::set<long> next;
    for (long x : sums)
      for (long m : extra) next.insert(x + m);
    sums.swap(next);
    if (sums.size() > 100000) throw Refusal("the windings admitted are more than 100000: the record is too sparse for its speed bound");
  }
  const long base = std::lround(chords / two_pi);
  for (long x : sums) out.add(base + x);
  if (out.certified()) out.reason = "every admissible history winds " + std::to_string(out.values.front()) + " time(s)";
  else if (!out.complete) out.reason = A.separation <= 0 ? "with no least separation, an admissible turn admits every number of turns: the windings listed are only some of those admitted" : "more windings are admitted than were listed";
  else out.reason = std::to_string(out.values.size()) + " windings are admitted, from " + std::to_string(out.values.front()) + " to " + std::to_string(out.values.back());
  return out;
}

// The record seen in a mirror (y ↦ −y), a homeomorphism of the plane that reverses orientation: it carries path
// classes to path classes (the draft's Proposition 7.2) and inverts every letter of a braid read on the same line.
inline Samples mirrored(Samples s) {
  for (auto& row : s.at)
    for (auto& p : row) p.y = -p.y;
  return s;
}
inline Braid mirrored(const Braid& b) {
  Braid m(b.strands);
  for (int l : b.word) m.word.push_back(-l);
  return m;
}

// ---- records in which points arrive and leave ------------------------------------------------------------

// An observation of some of the points at one time.
struct Observation {
  double time = 0;
  std::vector<std::string> labels;
  std::vector<Point> at;
};

// A maximal run of observations with one set of points, as a record with a fixed number of strands.
struct Segment {
  Samples samples;
  std::vector<std::string> arrived, left;  // against the segment before
};

// The observations cut where the set of points changes. Each segment has its own configuration space; a braid
// is read within a segment and never across the arrival or the departure of a point (the draft, §3.4).
inline std::vector<Segment> segments(std::vector<Observation> observations) {
  std::sort(observations.begin(), observations.end(), [](const Observation& x, const Observation& y) { return x.time < y.time; });
  std::vector<Segment> out;
  std::vector<std::string> current;
  for (const Observation& o : observations) {
    if (o.labels.size() != o.at.size()) throw Refusal("an observation holds a position for each of its labels");
    std::vector<std::string> sorted = o.labels;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) throw Refusal("an observation names a point twice");
    if (out.empty() || sorted != current) {
      Segment seg;
      std::set_difference(sorted.begin(), sorted.end(), current.begin(), current.end(), std::back_inserter(seg.arrived));
      std::set_difference(current.begin(), current.end(), sorted.begin(), sorted.end(), std::back_inserter(seg.left));
      seg.samples.labels = sorted;
      out.push_back(std::move(seg));
      current = sorted;
    }
    std::vector<Point> row(current.size());
    for (std::size_t i = 0; i < o.labels.size(); i++) row[static_cast<std::size_t>(std::lower_bound(current.begin(), current.end(), o.labels[i]) - current.begin())] = o.at[i];
    out.back().samples.times.push_back(o.time), out.back().samples.at.push_back(std::move(row));
  }
  return out;
}

}  // namespace kinesis::braids
