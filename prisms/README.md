# Prisms

Status: **draft** (2026-11-18). Module `kinesis::prisms` ([decision 0006](../development/decisions/0006-prisms.md)).

## 1. Concepts

A record of a mover holds fixes: at an instant it was within a radius of a place. Its histories are the paths through every fix's disc at the fix's instant that never go faster than a declared speed. Between two fixes the places a history can occupy at an instant are the intersection of two discs, the future cone of the earlier fix and the past cone of the later: Hägerstrand's space-time prism (1970), formalised by Lenntorp (1976) and Miller (2005), the Alexandrov interval of a cone structure of speed v. What the record identifies is a set (the kernel's invariant 5), and each question is answered so.

- **Consistency.** Every pair of fixes must be joinable within the speed (a necessary condition, checked over every pair, the pair that breaks named and how far short). With exact fixes the consecutive pairs suffice; with radii a witness history is sought by cyclic projections onto the discs and the bounds on each step, and checked exactly; otherwise the answer is undecided.
- **The slice** at an instant is the two bracketing cones' discs: exact for exact fixes, an outer bound otherwise.
- **Reach.** On each piece between fixes, the depth of the slice and a region (how far inside all of them the deepest point lies) is convex in time, so the least depth is found by Brent's method and the instants at which the region can be reached form one interval a piece.
- **Presence** in a region during an interval is the set {in, out} over every admissible history: in impossible when the region is out of reach (an alibi, certified), possible with a witness through it; out impossible when a fix's disc lies within the region, admissible with a witness that stays out. What is not decided is said.

The geometry is Euclidean: straight paths, no obstacle. Units and frame are the caller's.

## 2. Use

```cpp
#include <kinesis/prisms/prisms.hpp>
namespace kp = kinesis::prisms;
std::vector<kp::Fix> seen = {{"inn", -560, 210, 780, 3}, {"home", 1650, 940, 1260, 15}};
auto c = kp::consistent(seen, 6.0);                      // consistent, with a witness history
auto r = kp::reach(seen, 6.0, {0, 0, 5, 780, 1260});      // when the door can be reached: one window
auto p = kp::presence(seen, 5.0, {0, 0, 5, 780, 1260});  // at 5 m/s: certified "out"
```

## 3. Checks

`tests/check_prisms.cpp`, 32 checks: the deepest point of discs against closed forms and a fine grid on 200 random sets; a walk consistent at its own speed and not below it; three discs that meet pairwise and not together never said consistent; a place out of reach short by 45 m at the middle instant; a window from 245 s to 355 s; presence made necessary by a fix; a cyclist's window past a door in closed form, and the door out of reach at a lower speed.

## 4. References

- Hägerstrand, T. (1970), What about people in regional science?, *Papers of the Regional Science Association* 24, 7–21.
- Lenntorp, B. (1976), *Paths in Space-Time Environments*, Lund Studies in Geography B 44.
- Miller, H. J. (2005), A measurement theory for time geography, *Geographical Analysis* 37, 17–45.
- Bregman, L. M. (1965), The method of successive projection for finding a common point of convex sets, *Soviet Mathematics Doklady* 6, 688–692.
