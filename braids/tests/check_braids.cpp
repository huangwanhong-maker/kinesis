// Checks of kinesis/braids: braids read from records, what a record identifies, and entropy.
//
// The draft's synthetic archive (Huang 2026, §8, Table 2), state by state. Point a moves once counterclockwise
// round point b at rest, z_a = exp(2πit), z_b = 0, 0 ≤ t ≤ 1; a record holds samples of it.
//   v0  the two ends, V = 4π            ambiguous
//   v1  quarter periods, V = 4π         ambiguous; the windings admitted include 1 and the draft's alias −3
//   v2  sixteenths, V = 4π              certified, winding +1; the draft's margin cos(π/16) − π/4 is reproduced
//   v3  the same in a mirror            transported: winding −1, the braid's letters inverted
//   v4  sixteenths, V = 16π             withdrawn: the windings admitted include 1 and the draft's alias 0
//   v5  sixty-fourths, V = 16π          certified again, winding +1; the draft's margin cos(π/64) − π/4
// and the dispositions between them are those of the table. The sets contract under refinement (its
// Proposition 7.1): I(v2) ⊆ I(v1) ⊆ I(v0).
//
// - The sharper certificate (*) certifies records that the draft's margin does not, and never the reverse; where
//   it fails for two points an alias is exhibited: a path through the same samples, within the speed bound, of
//   another winding.
// - Aliasing (its Proposition 6.6): a point at rest and one that turns N times have the same N + 1 samples.
// - The least length of a path about a disc against a direct construction.
// - Reading: a turn of two points is σ_1² and winds once; three points on a figure of eight (the choreography of
//   Moore 1993 and of Chenciner and Montgomery 2000) braid as (σ_1 σ_2⁻¹)³ up to conjugacy and inversion: a pure
//   braid of exponent sum 0 whose entropy is 3 log((3 + √5)/2). The reading does not depend on the sampling once
//   it is certified, and on another line the braid is conjugate (equal entropy, exponent sum and Alexander
//   polynomial).
// - Entropy: the Burau bound is attained for three strands at t = −1, and lies below the estimate for four;
//   a braid of finite order and a braid of twists have none.
// - Observations with a point that arrives are cut into segments with their own numbers of strands (its §3.4).

#include <kinesis/braids/entropy.hpp>
#include <kinesis/braids/trajectories.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>

using namespace kinesis::braids;
using kinesis::kernel::Disposition;
namespace ma = morphe::algebra;

static const double pi = std::numbers::pi;

// Point 0 at rest at the origin and point 1 turning `turns` times, at `count` + 1 equal times in [0, 1].
static Samples circle(std::size_t count, double turns = 1) {
  Samples s;
  s.labels = {"b", "a"};
  for (std::size_t k = 0; k <= count; k++) {
    const double t = static_cast<double>(k) / static_cast<double>(count), angle = 2 * pi * turns * t;
    s.times.push_back(t);
    // The last sample is the first, exactly: the path closes.
    s.at.push_back({Point{0, 0}, k == count ? Point{1, 0} : Point{std::cos(angle), std::sin(angle)}});
  }
  return s;
}

static bool holds(const kinesis::kernel::Identified<long>& set, long v) { return std::find(set.values.begin(), set.values.end(), v) != set.values.end(); }

int main() {
  calx::verification::Suite s;
  Assumptions slow, fast;
  slow.speed = {4 * pi}, slow.separation = 0.1;
  fast.speed = {16 * pi}, fast.separation = 0.1;

  // ---- the draft's archive ----
  const Samples v0 = circle(1), v1 = circle(4), v2 = circle(16), v5 = circle(64);
  const auto w0 = pair_winding(v0, 0, 1, slow), w1 = pair_winding(v1, 0, 1, slow), w2 = pair_winding(v2, 0, 1, slow);
  const auto w4 = pair_winding(v2, 0, 1, fast), w5 = pair_winding(v5, 0, 1, fast);
  s.require("v0, the two ends: ambiguous, and 1 is among the windings admitted", w0.ambiguous() && w0.complete && holds(w0, 1) && holds(w0, 0) && holds(w0, -2), w0.reason);
  s.require("v1, quarter periods: ambiguous; 1 and the draft's alias -3 are admitted", w1.ambiguous() && holds(w1, 1) && holds(w1, -3), w1.reason);
  s.require("v2, sixteenths: the winding +1 is certified", w2.certified() && w2.value() == 1, w2.reason);
  s.require("v4, the same samples under the speed 16 pi: 1 and the draft's alias 0 are admitted", w4.ambiguous() && holds(w4, 1) && holds(w4, 0), w4.reason);
  s.require("v5, sixty-fourths under 16 pi: the winding +1 is certified again", w5.certified() && w5.value() == 1, w5.reason);
  s.require("the sets contract under refinement: I(v2) in I(v1) in I(v0)", w2.refines(w1) && w1.refines(w0));
  s.require("the dispositions of the table: ambiguous, ambiguous, issued, withdrawn, issued again",
            kinesis::kernel::disposition(kinesis::kernel::Identified<long>{}, w0) == Disposition::ambiguous && kinesis::kernel::disposition(w0, w1) == Disposition::ambiguous &&
                kinesis::kernel::disposition(w1, w2) == Disposition::reissued && kinesis::kernel::disposition(w2, w4) == Disposition::withdrawn &&
                kinesis::kernel::disposition(w4, w5) == Disposition::reissued);
  {
    const Certificate c2 = certify(v2, slow), c1 = certify(v1, slow), c4 = certify(v2, fast), c5 = certify(v5, fast);
    s.require("the certificate agrees: v1 open, v2 certified, v4 open, v5 certified", !c1.certified && c2.certified && !c4.certified && c5.certified, c2.reason + " | " + c4.reason);
    s.check("the draft's margin of v2, cos(pi/16) - pi/4 (its Eq. 25)", c2.coarse_margin - (std::cos(pi / 16) - pi / 4), 1e-14);
    s.check("the draft's margin of v5, cos(pi/64) - pi/4 (its Eq. 27)", c5.coarse_margin - (std::cos(pi / 64) - pi / 4), 1e-14);
    s.check("the sharper margin of v2: 2 - 2 (4 pi) / 16", c2.margin - (2 - pi / 2), 1e-14);
    const Reading r2 = read(v2);
    s.require("the braid read from v2 is s1^2: two points that wind once", r2.braid.word == std::vector<int>{1, 1} && ma::crossings(r2.braid)[0][1] == 2, r2.braid.str());
    const Reading r3 = read(mirrored(v2));
    const auto w3 = pair_winding(mirrored(v2), 0, 1, slow);
    s.require("v3, in a mirror: the winding -1 is certified and the braid is the mirror's", w3.certified() && w3.value() == -1 && r3.braid.word == mirrored(r2.braid).word, r3.braid.str());
  }

  // ---- the sharper certificate, and aliases where it fails ----
  {
    // With each point no faster than pi and N samples of the turn, (*) is 2 - 2 pi / N and the draft's margin is
    // cos(pi / N) - pi / N: at N = 4 the first is positive and the second negative; from N = 5 both are positive.
    std::size_t sharper_only = 0, coarse_only = 0, both = 0;
    for (std::size_t count = 3; count <= 40; count++) {
      Assumptions A;
      A.speed = {pi};  // each point; the relative speed is at most 2 pi, which the turning point attains alone
      const Certificate c = certify(circle(count), A);
      const bool coarse = c.coarse_margin > 0;
      sharper_only += c.certified && !coarse, coarse_only += coarse && !c.certified, both += coarse && c.certified;
    }
    s.require("(*) certifies every record the draft's margin does, and some it does not", coarse_only == 0 && sharper_only > 0 && both > 0,
              std::to_string(sharper_only) + " by (*) alone, " + std::to_string(both) + " by both");
  }
  {
    // Where (*) fails strictly for two points, another class is admitted: an alias is built and measured.
    Assumptions A;
    A.speed = {4 * pi}, A.separation = 0.05;
    const Samples sparse = circle(4);
    const Certificate c = certify(sparse, A);
    const auto set = pair_winding(sparse, 0, 1, A);
    // The alias of the draft: in to the radius 0.1, clockwise by 3 pi / 2 there, and out, in each quarter.
    const double r = 0.1, length = 2 * (1 - r) + 3 * pi * r / 2, budget = 2 * 4 * pi / 4;
    s.require("where (*) fails the set has more than one winding", !c.certified && set.values.size() > 1);
    s.require("the draft's alias path is within the speed bound: its length against the budget of a quarter", length < budget, std::to_string(length) + " < " + std::to_string(budget));
    s.require("the least length of the class of that alias, about a disc of radius 0.1, is below the path's", least_length(1, 0, 0, 1, -3 * pi / 2, r) < length,
              std::to_string(least_length(1, 0, 0, 1, -3 * pi / 2, r)) + " < " + std::to_string(length));
    s.check("the least length about a disc: two tangents and an arc", least_length(1, 0, 0, 1, -3 * pi / 2, r) - (2 * std::sqrt(1 - r * r) + r * (3 * pi / 2 - 2 * std::acos(r))), 1e-14);
    s.check("the chord is the least length of its own class", least_length(1, 0, 0, 1, pi / 2, r) - std::sqrt(2.0), 1e-15);
    s.check("with no disc every other class costs |a| + |b|", least_length(1, 0, 0, 2, -3 * pi / 2, 0) - 3, 1e-15);
  }
  {
    Assumptions A;
    A.speed = {4 * pi};  // no least separation
    const auto set = pair_winding(circle(4), 0, 1, A);
    s.require("with no least separation an admitted turn admits every number: the set is not complete, and says so", !set.complete && !set.certified() && set.ambiguous(), set.reason);
    Assumptions slowest;
    slowest.speed = {0.1}, slowest.separation = 0.1;
    const auto none = pair_winding(circle(4), 0, 1, slowest);
    s.require("a speed bound the samples exceed leaves no history: inconsistent", none.inconsistent(), none.reason);
    Assumptions noisy = A;
    noisy.error = 0.01;
    try {
      pair_winding(circle(16), 0, 1, noisy);
      s.require("windings are refused when positions have an error", false);
    } catch (const kinesis::kernel::Refusal&) {
      // (*) with an error eta: 2 - 2 (4 pi) / 16 - 8 eta, positive at eta = 0.01 and negative at 0.06.
      Assumptions noisier = noisy;
      noisier.error = 0.06;
      s.require("windings are refused when positions have an error; the certificate is given, and an error too large withdraws it", certify(circle(16), noisy).certified && !certify(circle(16), noisier).certified);
      s.check("the margin with an error of 0.01", certify(circle(16), noisy).margin - (2 - pi / 2 - 0.08), 1e-14);
    }
  }
  {
    // Aliasing: N turns sampled N + 1 times is a point at rest.
    const std::size_t N = 7;
    Samples turning = circle(N, static_cast<double>(N)), still = circle(N, 0);
    double worst = 0;
    for (std::size_t k = 0; k <= N; k++) worst = std::max(worst, std::hypot(turning.at[k][1].x - still.at[k][1].x, turning.at[k][1].y - still.at[k][1].y));
    s.check("aliasing: 7 turns at 8 samples are the samples of a point at rest", worst, 1e-14);
    s.require("and the braid read from them is the identity", read(turning).braid.word.empty());
  }

  // ---- three points on a figure of eight ----
  {
    const auto eight = [](std::size_t count) {
      Samples e;
      for (std::size_t k = 0; k <= count; k++) {
        const double t = 2 * pi * static_cast<double>(k) / static_cast<double>(count);
        e.times.push_back(t);
        std::vector<Point> row;
        for (int j = 0; j < 3; j++) {
          const double u = t + 2 * pi * j / 3 + 0.1;  // a phase, so that no two projections are equal at the ends
          row.push_back({std::sin(u), std::sin(2 * u) / 2});
        }
        e.at.push_back(row);
      }
      return e;
    };
    const Reading fine = read(eight(600)), coarse = read(eight(60)), turned = read(eight(600), 0.4);
    const double golden = std::log((3 + std::sqrt(5.0)) / 2);
    const Entropy h = entropy(fine.braid);
    s.require("three points on a figure of eight: a pure braid of six crossings and exponent sum 0", ma::pure(fine.braid) && fine.braid.length() == 6 && ma::exponent_sum(fine.braid) == 0, fine.braid.str());
    s.check("its entropy per period is 3 log((3 + sqrt 5) / 2)", h.estimate - 3 * golden, 1e-9);
    s.require("it is (s1 s2^-1)^3 up to conjugacy and inversion: one of the two is conjugate by a letter or equal",
              [&] {
                const ma::Braid pigtail = ma::power(ma::Braid(3, {1, -2}), 3), mirror = ma::power(ma::Braid(3, {-1, 2}), 3);
                for (const ma::Braid& target : {pigtail, mirror})
                  for (int g : {0, 1, -1, 2, -2}) {
                    const ma::Braid c = g == 0 ? ma::Braid(3) : ma::Braid(3, {g});
                    for (const ma::Braid& d : {c, c * c, c * ma::Braid(3, {g == 1 ? 2 : 1})})
                      if (ma::equal(d * fine.braid * ma::inverse(d), target).value()) return true;
                  }
                return false;
              }(),
              fine.braid.str());
    s.require("a coarser sampling reads the same braid", ma::equal(coarse.braid, fine.braid).value(), coarse.braid.str());
    s.check("on another line the braid has the same entropy", entropy(turned.braid).estimate - h.estimate, 1e-9);
    s.require("and the same exponent sum and Alexander polynomial of its closure", ma::exponent_sum(turned.braid) == 0 && ma::alexander(turned.braid) == ma::alexander(fine.braid), turned.braid.str());
    Assumptions A;
    A.speed = {1.5};  // |z'| <= sqrt(cos^2 + cos^2(2u)) <= sqrt 2
    s.require("with 600 samples the reading is certified, and with 12; with 6 it is not", certify(eight(600), A).certified && certify(eight(12), A).certified && !certify(eight(6), A).certified,
              certify(eight(12), A).reason + " | " + certify(eight(6), A).reason);
    s.check("the braiding exponent of one period over its time", braiding_exponent(fine.braid, 2 * pi) * 2 * pi - std::log(ma::act(fine.braid, ma::generating_loops(3)).norm().to_double()), 1e-12);
  }

  // ---- entropy and its bound ----
  {
    const double golden = std::log((3 + std::sqrt(5.0)) / 2);
    const Entropy three = entropy(ma::Braid(3, {1, -2}));
    s.check("s1 s2^-1: the estimate", three.estimate - golden, 1e-9);
    s.check("s1 s2^-1: the Burau bound is attained", three.lower - golden, 1e-6);
    const Entropy four = entropy(ma::Braid(4, {1, 2, -3}));
    s.require("s1 s2 s3^-1: the bound lies below the estimate, which settled", four.settled && four.lower <= four.estimate + 1e-9 && four.lower > 0.3 && four.estimate > 0.83 && four.estimate < 0.832,
              std::to_string(four.lower) + " <= " + std::to_string(four.estimate));
    const Entropy finite = entropy(ma::Braid(4, {1, 2, 3})), twists = entropy(ma::Braid(3, {1, 1}));
    s.require("a braid of finite order and a braid of twists have no entropy", finite.estimate == 0 && twists.estimate == 0 && finite.lower < 1e-5,
              std::to_string(finite.estimate) + ", " + std::to_string(twists.estimate) + "; the bound " + std::to_string(finite.lower));
    s.require("two strands have none", entropy(ma::Braid(2, {1, 1, 1})).estimate == 0);
  }

  // ---- a point that arrives ----
  {
    std::vector<Observation> obs;
    for (int k = 0; k <= 16; k++) {
      const double t = k / 16.0;
      obs.push_back({t, {"a", "b"}, {Point{std::cos(2 * pi * t), std::sin(2 * pi * t)}, Point{0, 0}}});
    }
    obs.push_back({1.25, {"c", "a", "b"}, {Point{3, 0}, Point{1, 0}, Point{0, 0}}});
    obs.push_back({1.5, {"a", "b", "c"}, {Point{1, 0}, Point{0, 0}, Point{3, 0}}});
    const auto parts = segments(obs);
    s.require("a point that arrives begins another segment: two points at 17 times, then three at 2", parts.size() == 2 && parts[0].samples.points() == 2 && parts[0].samples.times.size() == 17 &&
                                                                                                          parts[1].samples.points() == 3 && parts[1].samples.times.size() == 2 && parts[1].arrived == std::vector<std::string>{"c"} &&
                                                                                                          parts[1].left.empty());
    s.require("each segment is read in its own group: s1^2 on two strands, the identity on three", read(parts[0].samples).braid.word == std::vector<int>{1, 1} && read(parts[1].samples).braid.strands == 3 && read(parts[1].samples).braid.word.empty());
    Samples uneven = parts[0].samples;
    uneven.at[3].push_back(Point{5, 5});
    try {
      read(uneven);
      s.require("a record whose number of points changes is refused", false);
    } catch (const kinesis::kernel::Refusal&) {
      s.require("a record whose number of points changes is refused", true);
    }
  }

  // ---- an exchange that falls exactly on a sample ----
  {
    // Point 0 at rest at the origin; point 1 moves along y = 1, its projection given at three times.
    const auto line = [](double x0, double x1, double x2) {
      Samples p;
      p.labels = {"a", "b"};
      p.times = {0, 1, 2};
      p.at = {{Point{0, 0}, Point{x0, 1}}, {Point{0, 0}, Point{x1, 1}}, {Point{0, 0}, Point{x2, 1}}};
      return read(p);
    };
    const Reading passed = line(-1, 0, 1), back = line(1, 0, -1), touched = line(-1, 0, -1), grazed = line(1, 0, 1);
    s.require("a point that passes over another exactly at a sample exchanges once, from either side, and the two letters are inverse to each other",
              passed.braid.word.size() == 1 && back.braid.word.size() == 1 && passed.braid.word[0] == -back.braid.word[0] && passed.crossings[0].time == 1 && back.crossings[0].time == 1,
              std::to_string(passed.braid.word.size()) + ", " + std::to_string(back.braid.word.size()));
    s.require("a point that touches the other's projection at a sample and turns back does not exchange: two letters that cancel, or none",
              touched.braid.word.size() == 2 && touched.braid.word[0] == -touched.braid.word[1] && grazed.braid.word.empty(),
              std::to_string(touched.braid.word.size()) + ", " + std::to_string(grazed.braid.word.size()));
  }

  return s.done();
}
