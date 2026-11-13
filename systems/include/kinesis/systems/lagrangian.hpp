// Systems given by an action: the derivation `vary`, from an action to the equations of its stationary paths.
//
// A declaration `action S = integral(L)` with `histories stationary S` gives a Lagrangian L(q, q', t) on the
// coordinates q_1 … q_n it declares as states. A path makes S stationary exactly when it satisfies the
// Euler–Lagrange equations (Arnold 1989, Mathematical Methods of Classical Mechanics, 2nd ed., §12–13),
//
//     d/dt ∂L/∂q_i' − ∂L/∂q_i = 0.
//
// With the momenta p_i = ∂L/∂q_i' and D_t the total derivative, D_t p_i = ∂p_i/∂t + Σ_j ∂p_i/∂q_j q_j' + Σ_j
// ∂p_i/∂q_j' q_j'', so the equations are linear in the accelerations,
//
//     M q'' = F,   M_ij = ∂²L/∂q_i'∂q_j',   F_i = ∂L/∂q_i − ∂p_i/∂t − Σ_j ∂p_i/∂q_j q_j'.
//
// Every one of these is a term derived from L by the terms' own derivative, so the derivation is exact. vary()
// returns them, and for n ≤ 3 solves for q'' as terms (Cramer's rule), which gives a declaration of a flow again:
// the derived system is a word of the language, and runs as any flow does. A Lagrangian whose M is singular as a
// term (L linear in a velocity) is refused; more than three coupled coordinates are not solved as terms yet.
//
// The energy E = Σ p_i q_i' − L is conserved along the paths when L does not hold the time (Noether; Arnold §20),
// and is declared so in the derived system, where a check can hold a run to it.
#pragma once

#include <kinesis/systems/dynamics.hpp>

namespace kinesis::systems {

struct Lagrangian {
  std::string name, time;
  std::vector<std::string> coordinates;
  TermPtr lagrangian;                       // L
  std::vector<TermPtr> momenta;             // p_i
  std::vector<std::vector<TermPtr>> mass;   // M
  std::vector<TermPtr> force;               // F
  std::vector<TermPtr> acceleration;        // q_i'' as terms, when solved (n ≤ 3)
  TermPtr energy;                           // E
  bool autonomous = true;                   // L does not hold the time
};

namespace detail {
// The integrand of an action: integral(L), or a sum of such, or L itself.
inline TermPtr integrand(const TermPtr& action) {
  return terms::rewrite(action, [](const TermPtr& x) -> TermPtr { return x->kind == terms::Kind::call && x->name == "integral" && x->args.size() == 1 ? x->args[0] : nullptr; });
}
inline TermPtr determinant(const std::vector<std::vector<TermPtr>>& m) {
  const std::size_t n = m.size();
  if (n == 1) return m[0][0];
  if (n == 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}
}  // namespace detail

inline Lagrangian vary(const System& s) {
  if (form_of(s) != Form::lagrangian) throw kernel::Refusal(s.name + " is a " + std::string(name(form_of(s))) + " system; `vary` needs an action with `histories stationary`");
  const Scope sc = scope(s);
  const language::Named* action = nullptr;
  for (const auto& a : s.change.actions)
    if (a.name == s.structure.stationary) action = &a;
  if (!action) throw ScriptError(s.at, "histories are stationary for '" + s.structure.stationary + "', which is not an action of " + s.name);
  Lagrangian L;
  L.name = s.name, L.time = sc.time, L.coordinates = sc.states;
  L.lagrangian = sc.resolve(detail::integrand(action->value));
  const std::size_t n = L.coordinates.size();
  if (n == 0) throw ScriptError(s.at, s.name + " declares no coordinate");
  for (const auto& [nm, primes] : terms::symbols(L.lagrangian))
    if (sc.is_state(nm) && primes > 1) throw ScriptError(action->at, "the action holds " + nm + std::string(primes, '\'') + "; a Lagrangian of the first order is needed");
  L.autonomous = !terms::depends_on(L.lagrangian, L.time, 0);
  const std::set<std::string> dynamical(L.coordinates.begin(), L.coordinates.end());
  try {
    for (const auto& q : L.coordinates) L.momenta.push_back(terms::derivative(L.lagrangian, q, 1));
    L.mass.assign(n, std::vector<TermPtr>(n));
    for (std::size_t i = 0; i < n; i++) {
      std::vector<TermPtr> f{terms::derivative(L.lagrangian, L.coordinates[i], 0), -terms::derivative(L.momenta[i], L.time, 0)};
      for (std::size_t j = 0; j < n; j++) {
        L.mass[i][j] = terms::derivative(L.momenta[i], L.coordinates[j], 1);
        f.push_back(-(terms::derivative(L.momenta[i], L.coordinates[j], 0) * terms::symbol(L.coordinates[j], 1)));
      }
      L.force.push_back(terms::sum(std::move(f)));
    }
  } catch (const kernel::Refusal& e) {
    throw ScriptError(action->at, std::string("the action cannot be varied: ") + e.what());
  }
  std::vector<TermPtr> e{-L.lagrangian};
  for (std::size_t i = 0; i < n; i++) e.push_back(L.momenta[i] * terms::symbol(L.coordinates[i], 1));
  L.energy = terms::sum(std::move(e));
  if (n <= 3) {
    const TermPtr det = detail::determinant(L.mass);
    if (terms::is_number(det, 0)) throw ScriptError(action->at, "the action is degenerate: its second derivative in the velocities is singular, so it does not determine the accelerations");
    for (std::size_t i = 0; i < n; i++) {
      auto replaced = L.mass;
      for (std::size_t r = 0; r < n; r++) replaced[r][i] = L.force[r];
      L.acceleration.push_back(detail::determinant(replaced) / det);
    }
  }
  return L;
}

// The flow of an action's stationary paths, as a declaration: the coordinates with equations of the second order,
// the parameters as they were, and the energy declared conserved when the Lagrangian does not hold the time.
inline System varied(const System& s, const std::string& as = "") {
  const Lagrangian L = vary(s);
  if (L.acceleration.empty()) throw kernel::Refusal(s.name + " has " + std::to_string(L.coordinates.size()) + " coupled coordinates; their accelerations are not solved as terms beyond three");
  System out;
  out.name = as.empty() ? s.name + "_paths" : as;
  out.at = s.at;
  out.states = s.states;
  out.params = s.params;
  out.inputs = s.inputs;
  out.structure.timed = true, out.structure.time = L.time;
  out.structure.time_domain = s.structure.timed ? s.structure.time_domain : language::Domain{};
  out.difference = s.difference;
  out.constraint = s.constraint;
  for (std::size_t i = 0; i < L.coordinates.size(); i++) out.change.equations.push_back({L.coordinates[i], 2, {}, L.acceleration[i], s.at});
  if (L.autonomous) {
    bool declared = false;
    for (const auto& c : out.constraint.conserved) declared = declared || terms::same(c, L.energy);
    if (!declared) out.constraint.conserved.push_back(L.energy);
  }
  return out;
}

}  // namespace kinesis::systems
