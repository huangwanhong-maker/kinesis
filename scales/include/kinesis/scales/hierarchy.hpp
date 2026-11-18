// A hierarchy of units on which a system across scales is declared (decision 0007).
//
// The units of a system across scales are arranged in levels, the finest first: cells within districts within a city, or
// patches within a landscape within a region. Each unit of a level below the top has exactly one parent in the level above,
// so the levels are a chain of partitions: the units of a level partition the units of the level beneath it. A level may
// also carry what its units are measured by (an area, a length, a capacity), which units are neighbours (a graph with a
// conductance on each edge, as balance laws on graphs have it), and pairwise costs between its units (a time or a distance
// between every two units, from which interactions are computed).
//
// Units are named; a quantity of the system is a value for every unit of its level. Siblings are the units of a level
// that share a parent; the units of the top level are all siblings of one another. Allocation down the hierarchy (a share
// of a parent's total) is among siblings, so that what is allocated sums to what the parent holds.
#pragma once

#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "kinesis/kernel/answer.hpp"

namespace kinesis::scales {

using kernel::Refusal;

// One level of the hierarchy: its units, each unit's parent in the level above (none at the top), what the units are measured
// by, the graph of neighbours within the level, and named matrices of pairwise costs (row-major, n × n).
struct Level {
  std::string name;
  std::vector<std::string> units;
  std::vector<std::size_t> parent;                       // index into the next level's units; empty at the top
  std::map<std::string, std::vector<double>> measures;   // per unit: area, length, capacity, …
  struct Edge {
    std::size_t a = 0, b = 0;
    double conductance = 1;
  };
  std::vector<Edge> edges;                               // undirected
  std::map<std::string, std::vector<double>> costs;      // name → n × n, row-major
};

class Hierarchy {
 public:
  Hierarchy() = default;
  // The levels, the finest first. Every unit of a level below the top must name a parent in the level above.
  explicit Hierarchy(std::vector<Level> levels) : levels_(std::move(levels)) {
    if (levels_.empty()) throw Refusal("a hierarchy needs at least one level");
    std::map<std::string, std::size_t> seen;
    for (std::size_t l = 0; l < levels_.size(); l++) {
      const Level& v = levels_[l];
      if (v.name.empty()) throw Refusal("level " + std::to_string(l) + " has no name");
      if (!seen.emplace(v.name, l).second) throw Refusal("the level " + v.name + " is named twice");
      if (v.units.empty()) throw Refusal("the level " + v.name + " has no units");
      const std::size_t n = v.units.size();
      if (l + 1 < levels_.size()) {
        if (v.parent.size() != n)
          throw Refusal("every unit of the level " + v.name + " needs a parent in " + levels_[l + 1].name + ": " +
                        std::to_string(v.parent.size()) + " parents for " + std::to_string(n) + " units");
        for (std::size_t u = 0; u < n; u++)
          if (v.parent[u] >= levels_[l + 1].units.size())
            throw Refusal("the unit " + v.units[u] + " of " + v.name + " names a parent that " + levels_[l + 1].name + " does not have");
      } else if (!v.parent.empty()) {
        throw Refusal("the top level " + v.name + " has no level above it for its units' parents");
      }
      for (const auto& [name, m] : v.measures)
        if (m.size() != n) throw Refusal("the measure " + name + " of " + v.name + " needs one value a unit");
      for (const auto& e : v.edges)
        if (e.a >= n || e.b >= n || e.a == e.b || !(e.conductance >= 0))
          throw Refusal("an edge of " + v.name + " joins units that are not two of its own, or has a negative conductance");
      for (const auto& [name, c] : v.costs)
        if (c.size() != n * n) throw Refusal("the costs " + name + " of " + v.name + " need n × n values");
    }
    children_.resize(levels_.size());
    for (std::size_t l = 1; l < levels_.size(); l++) {
      children_[l].assign(levels_[l].units.size(), {});
      for (std::size_t u = 0; u < levels_[l - 1].units.size(); u++) children_[l][levels_[l - 1].parent[u]].push_back(u);
      for (std::size_t u = 0; u < levels_[l].units.size(); u++)
        if (children_[l][u].empty())
          throw Refusal("the unit " + levels_[l].units[u] + " of " + levels_[l].name + " has no units of " + levels_[l - 1].name +
                        " beneath it: the levels must partition one another");
    }
    neighbours_.resize(levels_.size());
    for (std::size_t l = 0; l < levels_.size(); l++) {
      neighbours_[l].assign(levels_[l].units.size(), {});
      for (const auto& e : levels_[l].edges) {
        neighbours_[l][e.a].push_back({e.b, e.conductance});
        neighbours_[l][e.b].push_back({e.a, e.conductance});
      }
    }
    groups_.resize(levels_.size());
    for (std::size_t l = 0; l < levels_.size(); l++) {
      if (l + 1 < levels_.size()) {
        groups_[l] = children_[l + 1];
      } else {
        groups_[l].assign(1, {});
        for (std::size_t u = 0; u < levels_[l].units.size(); u++) groups_[l][0].push_back(u);
      }
    }
  }

  std::size_t levels() const { return levels_.size(); }
  const Level& level(std::size_t l) const { return levels_.at(l); }
  std::size_t size(std::size_t l) const { return levels_.at(l).units.size(); }
  bool top(std::size_t l) const { return l + 1 == levels_.size(); }

  std::size_t index(const std::string& name) const {
    for (std::size_t l = 0; l < levels_.size(); l++)
      if (levels_[l].name == name) return l;
    throw Refusal("the hierarchy has no level " + name);
  }
  bool has(const std::string& name) const {
    for (const auto& v : levels_)
      if (v.name == name) return true;
    return false;
  }
  std::size_t unit(std::size_t l, const std::string& name) const {
    const auto& u = levels_.at(l).units;
    for (std::size_t i = 0; i < u.size(); i++)
      if (u[i] == name) return i;
    throw Refusal("the level " + levels_.at(l).name + " has no unit " + name);
  }

  // The units of level l − 1 beneath unit u of level l (l ≥ 1).
  const std::vector<std::size_t>& children(std::size_t l, std::size_t u) const {
    if (l == 0) throw Refusal("the finest level has no units beneath it");
    return children_[l].at(u);
  }
  // The unit of level l + 1 above unit u of level l.
  std::size_t parent(std::size_t l, std::size_t u) const {
    if (top(l)) throw Refusal("the top level " + levels_[l].name + " has no level above it");
    return levels_[l].parent.at(u);
  }
  // The unit of level `to` (≥ from) that holds unit u of level `from`.
  std::size_t ancestor(std::size_t from, std::size_t u, std::size_t to) const {
    if (to < from) throw Refusal("an ancestor is in a level above");
    for (std::size_t l = from; l < to; l++) u = parent(l, u);
    return u;
  }
  // The groups of siblings of level l: the children of each unit above, or all units at the top.
  const std::vector<std::vector<std::size_t>>& groups(std::size_t l) const { return groups_.at(l); }
  // The group of siblings unit u of level l belongs to.
  std::size_t group(std::size_t l, std::size_t u) const { return top(l) ? 0 : levels_[l].parent.at(u); }
  const std::vector<std::pair<std::size_t, double>>& neighbours(std::size_t l, std::size_t u) const { return neighbours_.at(l).at(u); }

  // A unit's measure, or 1 when no measure is named.
  double measure(std::size_t l, const std::string& name, std::size_t u) const {
    if (name.empty()) return 1.0;
    const auto& m = levels_.at(l).measures;
    const auto it = m.find(name);
    if (it == m.end()) throw Refusal("the level " + levels_[l].name + " has no measure " + name);
    return it->second.at(u);
  }
  const std::vector<double>& costs(std::size_t l, const std::string& name) const {
    const auto& c = levels_.at(l).costs;
    const auto it = c.find(name);
    if (it == c.end()) throw Refusal("the level " + levels_[l].name + " has no costs named " + name);
    return it->second;
  }

 private:
  std::vector<Level> levels_;
  std::vector<std::vector<std::vector<std::size_t>>> children_;
  std::vector<std::vector<std::vector<std::pair<std::size_t, double>>>> neighbours_;
  std::vector<std::vector<std::vector<std::size_t>>> groups_;
};

}  // namespace kinesis::scales
