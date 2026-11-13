# Braids of Trajectories

Status: **draft** (2026-11-13). Module `kinesis::braids`.

## 1. Concepts

A record holds the positions of n labelled points of the plane at finitely many times. A history is a continuous, collision-free motion of the points through the recorded positions, and its class is a morphism of the fundamental groupoid of the configuration space. A record does not determine it: between two samples a point may turn about another any number of times. What a record identifies is relative to assumptions (a bound on each point's speed, an error of observation, a least separation), and is answered in the kernel's vocabulary: certified, ambiguous or inconsistent.

The piecewise-linear interpolant of a record is one history, and its braid is read by projection on a line: each time two neighbours in the order of projection exchange places a letter σ_i^{±1} is written, positive when the left one passes below. A certificate then says whether every history the assumptions admit has that braid. For a pair of points, the windings the record admits are computed exactly, as a set. A braid, once read, bounds the topological entropy of any flow that moves points along it, which is why braids of trajectories measure mixing. Theory: [docs/theory.md](../docs/theory.md).

## 2. Content

- `trajectories.hpp`:
  - `Point`; `Samples{times, at, labels}` with `points()`, `validate()`, `thinned(stride)`; `Crossing`, `Reading{braid, crossings, first, last}`; `read(s, angle)`, the braid of the interpolant on the line at `angle`.
  - `Assumptions{speed, error, separation, angle}` with `speed_of(i)` and `stated()`; `Certificate{certified, margin, interval, a, b, coarse_margin, endpoints, open, reason}`; `certify(s, A)`, which holds when |r_k| + |r_{k+1}| − (V_i + V_j) h − 8η > 0 for every interval and every pair, and the two ends of the record are in one cell whatever the error.
  - `least_length(ax, ay, bx, by, theta, rho)`, the shortest relative path of a total angle that keeps |r| ≥ ρ; `pair_winding(s, a, b, A, most)`, an `Identified<long>`: the windings of point b about point a that the record admits.
  - `mirrored(Samples)` and `mirrored(Braid)`: the record in a mirror, which inverts every letter of a braid read on the same line.
  - `Observation`, `Segment{samples, arrived, left}`, `segments(observations)`: observations cut where the set of points changes.
- `entropy.hpp`: `Entropy{estimate, lower, applications, settled}`; `entropy(b, tolerance, most)`, the growth of the generating loops in Dynnikov's coordinates under the braid; `burau_bound(b, grid)`, log sup ρ(ψ_t(β)) over a grid of t on the unit circle; `braiding_exponent(b, duration)`, the exponent of a braid read over a finite time.

*Remark.* The certificate was derived for this library; the draft's coarser margin is reported beside it as `coarse_margin`. With η = 0 and one bound V the condition is implied by the draft's Corollary 6.4 and is weaker; for two points it is sharp.

*Caveat.* A record whose first or last configuration has two equal projections, or in which three points cross at one time, is not generic for the line and is refused. A braid is read within a segment and never across the arrival or the departure of a point.

*Caveat.* `pair_winding` is for exact positions and a relative path that closes. With no least separation one admissible turn admits every number of turns: the set is then infinite, and is reported as not complete.

*Caveat.* The entropy estimate is a limit, and says whether it was reached (`settled`); a rate that settles below the tolerance is reported as zero growth. The Burau bound is taken over a grid of t, so it is itself a bound from below on the supremum. The braiding exponent uses the norm of Dynnikov's coordinates in the place of the count of intersections, which has the same exponential rate and not the same value at short times.

## 3. Checks

44: the draft's synthetic archive (Huang 2026, §8, Table 2), state by state, for a point turning once round another at rest: the two ends ambiguous; quarter periods ambiguous, the windings admitted including 1 and the draft's alias −3; sixteenths certified with winding +1 and the draft's margin cos(π/16) − π/4 reproduced; the same in a mirror transported, with winding −1 and the letters inverted; a faster bound withdrawing the certificate, with the alias 0 admitted; sixty-fourths certified again; the dispositions between them those of the table, and the sets contracting under refinement (its Proposition 7.1). The sharper certificate certifies records that the draft's margin does not, and never the reverse; where it fails for two points an alias is exhibited. Aliasing (its Proposition 6.6): a point at rest and one that turns N times have the same N + 1 samples. The least length of a path about a disc against a direct construction. A turn of two points read as σ_1² and winding once; three points on a figure of eight braiding as (σ_1 σ_2⁻¹)³ up to conjugacy and inversion, a pure braid of exponent sum 0 with entropy 3 log((3 + √5)/2); the reading independent of the sampling once certified, and conjugate on another line (equal entropy, exponent sum and Alexander polynomial). The Burau bound attained for three strands at t = −1 and below the estimate for four; no entropy for a braid of finite order or a braid of twists. Observations with a point that arrives cut into segments with their own numbers of strands.

## 4. Dependencies

`kinesis::kernel` (`Identified`, `Refusal`); morphe's `algebra` (`braid.hpp`: `Braid`, the Burau matrix, Dynnikov's coordinates and the action on them) and `arithmetic`; calx's `linalg` (eigenvalues). The checks use calx's `verification`.

## 5. References

Huang (2026), "Rethinking trajectory-based credibility under symbol inflation: possibilities for building epistemic infrastructure through topological methods with braid group theory as an example" (a draft): §3.1, §3.2, §3.4, Corollary 6.4, Proposition 6.6, Proposition 7.2, Eqs. 15 and 20. Reading: Thiffeault (2010), Chaos 20, 017516; Fox and Neuwirth (1962); Guillemot and Lairez (2026), arXiv:2601.23073. Entropy: Fathi, Laudenbach and Poénaru (1979); Boyland (1994); Boyland, Aref and Stremler (2000); Thiffeault (2005); Moussafir (2006); Fried (1986); Kolev (1989); Band and Boyland (2007); Budišić and Thiffeault (2015). In the checks: Moore (1993); Chenciner and Montgomery (2000).
