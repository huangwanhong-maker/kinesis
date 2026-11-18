// Checks of kinesis/prisms: what a record of positions identifies of a mover's path under a bound on its speed.
//
// - The deepest point of discs, held to closed forms (two discs apart; three equal discs about the centroid of an equilateral
//   triangle; a disc within a disc) and to a fine grid on 200 random sets of three to five discs: never above the grid's least
//   value, and never below it by more than the grid's half-diagonal (the excess is 1-Lipschitz).
// - Consistency: a walk recorded exactly is consistent at its own speed and not below it, the pair that breaks it named; three
//   fixes whose discs meet pairwise and not together are never said to be consistent, and are when the ends leave room the middle
//   disc reaches. (Where the ends allow exactly one place, cyclic projections only approach it, and the answer stays undecided.)
// - The slice between two exact fixes is the two cones' discs.
// - Presence: a mover at home at both ends of ten minutes cannot reach a place 700 m away at 2 m/s (it falls short by 45 m, at the
//   middle instant, in closed form); a place 500 m away it may reach between 245 s and 355 s, and both "in" and "out" are found
//   with their witnesses; a fix whose disc lies in the region makes "in" necessary.
// - The window of a cyclist seen at an inn and later at home, who may have passed a door between: in closed form, from the earliest
//   arrival (the distance from the inn less the radii over the speed) to the latest departure.
// - A record that breaks the bound answers "inconsistent" with no value.

#include <kinesis/prisms/prisms.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace kinesis::prisms;

static bool has(const kinesis::kernel::Identified<std::string>& s, const std::string& v) { return std::find(s.values.begin(), s.values.end(), v) != s.values.end(); }

int main() {
  calx::verification::Suite s;

  // ---- the deepest point ----
  {
    const auto two = deepest({{0, 0, 2}, {10, 0, 3}});
    s.check("two discs 10 apart, radii 2 and 3: depth (10 - 5) / 2", two.depth - 2.5, 1e-12);
    s.check("... at 4.5 from the first centre, on the segment", two.x - 4.5, 1e-12);
    const double L = 10, r = 4;
    std::vector<Disc> tri;
    for (int k = 0; k < 3; k++) tri.push_back({L * std::cos(2 * M_PI * k / 3), L * std::sin(2 * M_PI * k / 3), r});
    const auto t3 = deepest(tri);
    s.check("three equal discs about an equilateral triangle's centroid: depth L - r", t3.depth - (L - r), 1e-9);
    s.require("... at the centroid", std::hypot(t3.x, t3.y) < 1e-7);
    const auto in = deepest({{0, 0, 5}, {1, 0, 1}});
    s.check("a disc within a disc: depth -1 at the small disc's centre", in.depth + 1, 1e-12);
    calx::random::Stream rs(20261118, 1, 0, 0);
    std::size_t above = 0, below = 0;
    double worst = 0;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 3 + rs.below(3);
      std::vector<Disc> d;
      for (std::size_t k = 0; k < n; k++) d.push_back({rs.uniform(-10, 10), rs.uniform(-10, 10), rs.uniform(0, 8)});
      const auto best = deepest(d);
      const double h = 0.05;
      double grid = INFINITY;
      for (double x = -20; x <= 20; x += h)
        for (double y = -20; y <= 20; y += h) {
          double m = -INFINITY;
          for (const auto& c : d) m = std::max(m, std::hypot(x - c.x, y - c.y) - c.r);
          grid = std::min(grid, m);
        }
      if (best.depth > grid + 1e-9) above++;
      if (best.depth < grid - h * M_SQRT1_2 - 1e-9) below++;
      worst = std::max(worst, grid - best.depth);
    }
    s.require("on 200 random sets of discs the deepest point is never above a fine grid's least value", above == 0, std::to_string(above) + " above");
    s.require("... and never below it by more than the grid's half-diagonal", below == 0, std::to_string(below) + " below; the widest gap " + std::to_string(worst));
  }

  // ---- consistency ----
  {
    std::vector<Fix> walk;
    for (int k = 0; k <= 10; k++) walk.push_back({"w" + std::to_string(k), 1.4 * 60 * k, 0, 60.0 * k, 0});
    const auto ok = consistent(walk, 1.4), slow = consistent(walk, 1.3);
    s.require("a walk at 1.4 m/s recorded exactly is consistent at 1.4 m/s, with a witness", ok.consistent() && ok.witness && ok.witness->at.size() == 11, ok.reason);
    s.require("... and inconsistent at 1.3 m/s, the pair named", slow.inconsistent() && !slow.broken.empty(), slow.reason);
    s.check("... the least slack at 1.3 m/s: the whole walk short by 0.1 m/s for ten minutes", slow.margin + 0.1 * 600, 1e-9);
    const std::vector<Fix> pairwise = {{"a", 0, 0, 0, 0}, {"b", 5, 4, 5, 3.5}, {"c", 10, 0, 10, 0}};
    const auto p = consistent(pairwise, 1);
    s.require("three fixes whose discs meet pairwise and not together are never said to be consistent", !p.consistent(), p.status + ": " + p.reason);
    auto widened = pairwise;
    widened[2].x = 9;  // the ends now leave room: a lens, not one place
    const auto w = consistent(widened, 1);
    s.require("... and are when the ends leave room that the middle disc reaches", w.consistent(), w.reason);
    if (w.witness) {
      const auto& a = w.witness->at;
      bool good = true;
      for (std::size_t k = 0; k < a.size(); k++) good = good && std::hypot(a[k].x - widened[k].x, a[k].y - widened[k].y) <= widened[k].radius + 1e-6;
      for (std::size_t k = 0; k + 1 < a.size(); k++) good = good && std::hypot(a[k + 1].x - a[k].x, a[k + 1].y - a[k].y) <= 1 * (a[k + 1].t - a[k].t) + 1e-6;
      s.require("... the witness passes every disc at its instant and never goes faster than the bound", good);
    }
  }

  // ---- the slice ----
  {
    const auto sl = slice({{"A", 0, 0, 0, 0}, {"B", 100, 0, 100, 0}}, 2, 50);
    s.require("the slice between two exact fixes is two discs", sl.size() == 2);
    s.check("... the earlier fix's cone has grown to v (t - t_A)", sl[0].r - 100, 1e-12);
    s.check("... the later fix's has shrunk to v (t_B - t)", sl[1].r - 100, 1e-12);
    const auto before = slice({{"A", 0, 0, 0, 0}, {"B", 100, 0, 100, 0}}, 2, -10);
    s.require("before the first fix, its past cone alone", before.size() == 1 && std::fabs(before[0].r - 20) < 1e-12);
  }

  // ---- presence ----
  {
    const std::vector<Fix> home = {{"home at the start", 0, 0, 0, 0}, {"home at the end", 0, 0, 600, 0}};
    const auto far = presence(home, 2, {700, 0, 10, 0, 600});
    s.require("at home at both ends of ten minutes, a place 700 m away is out of reach at 2 m/s: certified out", far.identified.certified() && far.identified.value() == "out", far.identified.reason);
    s.check("... it falls short by (700 - 600 - 10) / 2 = 45 m", far.reach.least - 45, 1e-6);
    s.check("... at the middle instant", far.reach.when - 300, 1e-2);
    const auto near = presence(home, 2, {500, 0, 10, 0, 600});
    s.require("a place 500 m away may be reached: both in and out are found, with witnesses", near.identified.complete && has(near.identified, "in") && has(near.identified, "out") && near.in && near.out, near.identified.reason);
    s.require("... in one window", near.reach.windows.size() == 1);
    if (!near.reach.windows.empty()) {
      s.check("... from 490 / 2 = 245 s", near.reach.windows[0].from - 245, 1e-2);
      s.check("... to 600 - 245 = 355 s", near.reach.windows[0].until - 355, 1e-2);
    }
    if (near.in) {
      const auto [x, y] = position(*near.in, near.reach.when);
      s.require("... the witness through the region is in it at the instant found", std::hypot(x - 500, y) <= 10 + 1e-6);
    }
    const std::vector<Fix> seen = {{"home", 0, 0, 0, 0}, {"seen at the door", 500, 0, 300, 5}, {"home again", 0, 0, 600, 0}};
    const auto must = presence(seen, 2, {500, 0, 10, 250, 350});
    s.require("a fix whose disc lies in the region during its interval makes 'in' necessary", must.identified.certified() && must.identified.value() == "in" && must.necessary_by == "seen at the door", must.identified.reason);
    const auto broken = presence({{"a", 0, 0, 0, 0}, {"b", 1000, 0, 10, 0}}, 2, {0, 0, 1, 0, 10});
    s.require("a record that breaks the bound answers inconsistent, with no value", broken.identified.inconsistent(), broken.identified.reason);
  }

  // ---- a window in closed form ----
  {
    // Seen at an inn (within 3 m) and, eight minutes later, at home (within 15 m), on a bicycle at 6 m/s: the door between them.
    const double px = -560, py = 210, hx = 1650, hy = 940, v = 6;
    const std::vector<Fix> f = {{"inn", px, py, 0, 3}, {"home", hx, hy, 480, 15}};
    const auto r = reach(f, v, {0, 0, 5, 0, 480});
    const double earliest = (std::hypot(px, py) - 3 - 5) / v, latest = 480 - (std::hypot(hx, hy) - 15 - 5) / v;
    s.require("a cyclist seen at an inn and later at home may have passed a door between, in one window", r.windows.size() == 1);
    if (!r.windows.empty()) {
      s.check("... from the earliest arrival, (|inn - door| - 3 - 5) / 6", r.windows[0].from - earliest, 1e-2);
      s.check("... to the latest departure, 480 - (|home - door| - 15 - 5) / 6", r.windows[0].until - latest, 1e-2);
    }
    const auto p = presence(f, v, {0, 0, 5, 0, 480});
    s.require("... so 'in' is admissible, with a history that passes the door", has(p.identified, "in") && p.in.has_value(), p.identified.reason);
    const auto slower = presence(f, 5.0, {0, 0, 5, 0, 480});
    s.require("at 5 m/s the record still holds (the inn to home directly) and the door is out of reach: an alibi, certified", slower.identified.certified() && slower.identified.value() == "out", slower.identified.reason);
  }

  return s.done();
}
