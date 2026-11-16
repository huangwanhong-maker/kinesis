// Balance laws on graphs: a quantity per cell of a graph of cells, and the operators of a balance law on it (decision 0005).
//
// A field here is a density u_i on the cells of a graph: each cell i has an area (a volume) V_i, and two cells that
// share a face are joined by an edge e with a conductance c_e = (width of the face) / (distance between the centres),
// which is 1 on a square grid. A balance law says how u changes: by what crosses the faces and by what is made or lost
// within a cell. Written in finite volumes (LeVeque 2002, Finite Volume Methods for Hyperbolic Problems, §4.1;
// Eymard, Gallouët and Herbin 2000, Finite volume methods, §§9–10), the flux of a density whose flux density is
//
//   q = −D ∇u − D₂ u ∇u + χ u ∇Φ
//
// across the face from cell a to cell b is, in units of the quantity per unit time,
//
//   F_e = c_e [ −D (u_b − u_a) − D₂ (u_b² − u_a²)/2 + χ (Φ_b − Φ_a)(u_a + u_b)/2 ],
//
// and the cell's own change is du_a/dt = −(1/V_a) Σ_{e from a} F_e + (what is made in a). The three parts of the flux are
// diffusion (Fick's law: a density spreads from where it is high), crowding (a diffusivity that grows with the density,
// D₂ u: the porous-medium equation of a gas through a porous solid, Vázquez 2007, The Porous Medium Equation, §2.1, of
// which ∇·(u∇u) = Δ(u²/2)), and drift up a potential Φ (a velocity χ ∇Φ carries the density: advection, or chemotaxis
// in Keller and Segel's sense, 1970, J. Theor. Biol. 26, 399–415). Every part is conservative by construction: what
// leaves a by e enters b, so Σ_i V_i (du_i/dt) is what is made, whatever the law. The drift's density on a face is the
// mean of the two cells' (a centred scheme, second order); it is not monotone when the drift outruns the diffusion
// across a cell (a cell Péclet number |χ ΔΦ| / (2D) above one), and a run says how much it had to clip.
//
// A square grid's cells are joined to their four neighbours, which gives the five-point Laplacian (exact on quadratics)
// on a grid of spacing h with c_e = 1 and V = h². Cells left out of the mask (the sea, a lake) are not cells: no flux
// crosses into them, which is the condition of no flux on the coast.
#pragma once

#include <kinesis/kernel/answer.hpp>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace kinesis::balance {

using kernel::Refusal;

struct Edge {
  std::size_t a = 0, b = 0;
  double conductance = 1;  // face width over the distance between the centres
};

struct Graph {
  std::vector<double> volume;  // each cell's area
  std::vector<Edge> edges;
  std::size_t size() const { return volume.size(); }
};

// A square grid's kept cells (mask non-zero, rows from the north), each joined to its four neighbours.
struct Grid {
  Graph graph;
  std::size_t rows = 0, cols = 0;
  double spacing = 1;
  std::vector<long> vertex;                       // rows * cols entries: the cell's vertex, or −1 when left out
  std::vector<std::size_t> row, col;              // each vertex's row and column
  std::vector<int> direction;                     // each edge's direction from a to b: 0 east, 1 south
};

inline Grid grid(std::size_t rows, std::size_t cols, const std::vector<char>& mask, double spacing) {
  if (rows == 0 || cols == 0) throw Refusal("a grid needs at least one row and one column");
  if (!mask.empty() && mask.size() != rows * cols) throw Refusal("the mask has " + std::to_string(mask.size()) + " entries, not rows × cols = " + std::to_string(rows * cols));
  if (!(spacing > 0)) throw Refusal("a grid's spacing must be positive");
  Grid g;
  g.rows = rows, g.cols = cols, g.spacing = spacing;
  g.vertex.assign(rows * cols, -1);
  for (std::size_t i = 0; i < rows; i++)
    for (std::size_t j = 0; j < cols; j++)
      if (mask.empty() || mask[i * cols + j]) {
        g.vertex[i * cols + j] = static_cast<long>(g.row.size());
        g.row.push_back(i), g.col.push_back(j);
        g.graph.volume.push_back(spacing * spacing);
      }
  if (g.row.empty()) throw Refusal("the mask keeps no cell");
  for (std::size_t i = 0; i < rows; i++)
    for (std::size_t j = 0; j < cols; j++) {
      const long v = g.vertex[i * cols + j];
      if (v < 0) continue;
      if (j + 1 < cols && g.vertex[i * cols + j + 1] >= 0) {
        g.graph.edges.push_back({static_cast<std::size_t>(v), static_cast<std::size_t>(g.vertex[i * cols + j + 1]), 1.0});
        g.direction.push_back(0);
      }
      if (i + 1 < rows && g.vertex[(i + 1) * cols + j] >= 0) {
        g.graph.edges.push_back({static_cast<std::size_t>(v), static_cast<std::size_t>(g.vertex[(i + 1) * cols + j]), 1.0});
        g.direction.push_back(1);
      }
    }
  return g;
}

inline void require_size(const Graph& g, const std::vector<double>& u, const char* what) {
  if (u.size() != g.size()) throw Refusal(std::string(what) + " has " + std::to_string(u.size()) + " values for " + std::to_string(g.size()) + " cells");
}

// Δu: (1/V_a) Σ_e c_e (u_b − u_a), the divergence of the gradient in finite volumes.
inline std::vector<double> laplacian(const Graph& g, const std::vector<double>& u) {
  require_size(g, u, "the field");
  std::vector<double> out(g.size(), 0.0);
  for (const Edge& e : g.edges) {
    const double f = e.conductance * (u[e.b] - u[e.a]);
    out[e.a] += f, out[e.b] -= f;
  }
  for (std::size_t i = 0; i < out.size(); i++) out[i] /= g.volume[i];
  return out;
}

// ∇·(u ∇u) = Δ(u²/2): crowding, a diffusivity proportional to the density.
inline std::vector<double> crowding(const Graph& g, const std::vector<double>& u) {
  require_size(g, u, "the field");
  std::vector<double> half(u.size());
  for (std::size_t i = 0; i < u.size(); i++) half[i] = 0.5 * u[i] * u[i];
  return laplacian(g, half);
}

// −∇·(u ∇Φ): what a drift up the potential Φ brings to each cell, the face's density the mean of its cells'.
inline std::vector<double> drift(const Graph& g, const std::vector<double>& u, const std::vector<double>& potential) {
  require_size(g, u, "the field");
  require_size(g, potential, "the potential");
  std::vector<double> out(g.size(), 0.0);
  for (const Edge& e : g.edges) {
    const double f = e.conductance * (potential[e.b] - potential[e.a]) * 0.5 * (u[e.a] + u[e.b]);  // from a to b
    out[e.a] -= f, out[e.b] += f;
  }
  for (std::size_t i = 0; i < out.size(); i++) out[i] /= g.volume[i];
  return out;
}

// The integral of a density over the graph, Σ V_i u_i: the people, if u is people per unit area.
inline double content(const Graph& g, const std::vector<double>& u) {
  require_size(g, u, "the field");
  double s = 0, c = 0;  // Kahan's summation
  for (std::size_t i = 0; i < u.size(); i++) {
    const double y = g.volume[i] * u[i] - c, t = s + y;
    c = (t - s) - y, s = t;
  }
  return s;
}

// The current through each cell of a grid: half the flux of each of its faces, as a vector (east, north) in units of the
// quantity per unit time, from the fluxes F_e of its edges. It is what an arrow draws of where the density is carried.
struct Currents {
  std::vector<double> east, north;
};

inline Currents currents(const Grid& g, const std::vector<double>& flux) {
  if (flux.size() != g.graph.edges.size()) throw Refusal("the fluxes are not one per edge");
  Currents c{std::vector<double>(g.graph.size(), 0.0), std::vector<double>(g.graph.size(), 0.0)};
  for (std::size_t k = 0; k < flux.size(); k++) {
    const Edge& e = g.graph.edges[k];
    if (g.direction[k] == 0) c.east[e.a] += 0.5 * flux[k], c.east[e.b] += 0.5 * flux[k];
    else c.north[e.a] -= 0.5 * flux[k], c.north[e.b] -= 0.5 * flux[k];  // rows run south
  }
  return c;
}

}  // namespace kinesis::balance
