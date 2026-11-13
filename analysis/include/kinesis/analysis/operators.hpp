// The lift to operators, computed: a system's transfer operator on a partition of its state space.
//
// A system acts on states; lifted, it acts on functions of the state (Koopman 1931) and, dually, on densities
// (the transfer operator of Perron and Frobenius; Lasota and Mackey 1994, Chaos, Fractals and Noise, 2nd ed.,
// ch. 3). The lift is linear whatever the system is (the idea of 2026-11-12, §4). Ulam's method (Ulam 1960, A
// Collection of Mathematical Problems, §VI.4; Li 1976, J. Approx. Theory 17, 177–186) makes it a matrix: the
// state space is cut into cells, and P_ij is the share of cell i that the map sends into cell j, here estimated
// from points placed evenly in each cell. P is a stochastic matrix, the chain of a walk on the cells.
//
//   - Its stationary vector approximates the invariant density of the map (for a map with one absolutely
//     continuous invariant measure, in the limit of fine cells: Li 1976 for expanding maps of the interval).
//   - Its other eigenvalues approximate the rates at which densities relax, with the caution that Ulam's matrix
//     adds a diffusion of the size of a cell, so that only eigenvalues well inside the unit disc of a map with a
//     spectral gap are meaningful.
//
// Under the lift, Difference is read from the operator: a change is invertible on what the partition can tell
// apart exactly when P is a permutation matrix, and how far it is from one is what the partition forgets.
#pragma once

#include <kinesis/systems/simulate.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace kinesis::analysis {

using systems::Dynamics;
using systems::State;

struct Transfer {
  std::vector<double> centres;              // of the cells
  calx::linalg::Matrix matrix;              // P_ij: from cell i to cell j
  std::vector<double> stationary;           // the invariant density on the cells (it integrates to 1)
  std::vector<std::complex<double>> eigenvalues;  // by decreasing modulus
  double lost = 0;                          // the share of points that left the interval
};

// Ulam's matrix of one component of a system on [lower, upper]: a map's iteration, or a flow's map over `time`.
// The other components are held at `base` (a system of one dimension has none).
inline Transfer transfer(const Dynamics& d, double lower, double upper, std::size_t cells = 200, std::size_t points = 50, double time = 1, std::size_t component = 0, State base = {}) {
  if (base.empty()) base.assign(d.dimension(), 0.0);
  Transfer T;
  T.matrix = calx::linalg::Matrix(cells);
  const double width = (upper - lower) / static_cast<double>(cells);
  systems::Options so;
  so.t1 = time, so.every = 100000000;
  std::size_t lost = 0;
  State x(d.dimension()), next(d.dimension());
  for (std::size_t i = 0; i < cells; i++) {
    T.centres.push_back(lower + (static_cast<double>(i) + 0.5) * width);
    for (std::size_t k = 0; k < points; k++) {
      x = base;
      x[component] = lower + (static_cast<double>(i) + (static_cast<double>(k) + 0.5) / static_cast<double>(points)) * width;
      double image;
      if (d.discrete) {
        d.field(0, x, next);
        image = next[component];
      } else {
        image = systems::simulate(d, x, so).x.back()[component];
      }
      const double f = (image - lower) / (upper - lower);
      if (!(f >= 0) || f > 1) {
        lost++;
        continue;
      }
      T.matrix(i, std::min(cells - 1, static_cast<std::size_t>(f * static_cast<double>(cells)))) += 1.0 / static_cast<double>(points);
    }
  }
  T.lost = static_cast<double>(lost) / static_cast<double>(cells * points);
  // The stationary vector by iteration of v ↦ v P from the uniform density.
  std::vector<double> v(cells, 1.0 / static_cast<double>(cells)), w(cells);
  for (std::size_t it = 0; it < 20000; it++) {
    std::fill(w.begin(), w.end(), 0.0);
    for (std::size_t i = 0; i < cells; i++) {
      if (v[i] == 0) continue;
      for (std::size_t j = 0; j < cells; j++) w[j] += v[i] * T.matrix(i, j);
    }
    double sum = 0, change = 0;
    for (double z : w) sum += z;
    if (sum <= 0) break;
    for (std::size_t j = 0; j < cells; j++) w[j] /= sum, change = std::max(change, std::fabs(w[j] - v[j]));
    v.swap(w);
    if (change < 1e-14) break;
  }
  T.stationary.resize(cells);
  for (std::size_t j = 0; j < cells; j++) T.stationary[j] = v[j] / width;
  T.eigenvalues = calx::linalg::eigenvalues(T.matrix);
  std::sort(T.eigenvalues.begin(), T.eigenvalues.end(), [](const std::complex<double>& a, const std::complex<double>& b) { return std::abs(a) > std::abs(b); });
  return T;
}

}  // namespace kinesis::analysis
