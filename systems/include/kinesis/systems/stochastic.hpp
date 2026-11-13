// Stochastic systems: equations with declared sources of noise, as a deterministic map of state and noise.
//
// In continuous time a system is a stochastic differential
//
//     dX_i = a_i(t, X) dt + Σ_j b_ij(t, X) dW_j,
//
// written as it is on paper: the right side of `dX = …` is linear in the differentials dt and dW_j, and its
// coefficients are read off as terms (a_i = ∂/∂(dt), b_ij = ∂/∂(dW_j); a right side that is not linear in them is
// refused). The calculus is stated in the declaration. A Stratonovich equation is converted once, exactly, to the
// Itô equation with the same solutions (Kloeden and Platen 1992, Numerical Solution of Stochastic Differential
// Equations, Springer, Eq. 4.9.15):
//
//     a_i^Itô = a_i^Strat + ½ Σ_j Σ_k b_kj ∂b_ij/∂x_k,
//
// the derivative taken from the terms. Without a stated calculus an equation is read in Itô's.
//
// Schemes (Kloeden and Platen 1992, §10.2 and §10.3), each a map of the state and the increments ΔW_j ~ N(0, h):
//   Euler–Maruyama   X += a h + Σ_j b_j ΔW_j                              strong order ½, weak order 1
//   Milstein         … + ½ b_i ∂b_i/∂x_i (ΔW_i² − h)                       strong order 1, for diagonal noise
// Milstein's scheme is given for diagonal noise (each source drives one component, and that component's
// coefficient depends on it alone among the driven ones); otherwise it is refused, since the general scheme needs
// iterated integrals that are not sampled here.
//
// In discrete time a system is a map x ↦ f(n, x, ξ) with a draw ξ_j of each source at each step.
//
// Noise as a declared source. Each source has a stream of its own, keyed by (seed, path, source), so that a
// path depends on its key alone and not on the order in which paths are run (calx/random/stream.hpp). A run may
// instead be given its increments, which is how two step sizes are compared on one Brownian path.
#pragma once

#include <kinesis/systems/dynamics.hpp>
#include <kinesis/systems/simulate.hpp>

#include <calx/random/stream.hpp>

#include <cmath>
#include <numbers>

namespace kinesis::systems {

class Stochastic {
 public:
  std::string name;
  bool discrete = false;
  std::string time = "t";
  std::string calculus = "ito";            // as declared; the terms below are Itô's
  std::vector<std::string> states, params;
  std::vector<double> param_values;
  std::vector<Source> sources;
  std::vector<TermPtr> drift_terms;                       // a_i (Itô); for a map, f_i
  std::vector<std::vector<TermPtr>> diffusion_terms;      // b_ij
  std::vector<InitialLaw> initial;
  bool diagonal = false;                                  // Milstein's scheme applies

  std::size_t dimension() const { return states.size(); }
  int state_index(const std::string& n) const {
    const auto it = std::find(states.begin(), states.end(), n);
    return it == states.end() ? -1 : static_cast<int>(it - states.begin());
  }
  void set_param(const std::string& n, double v) {
    const auto it = std::find(params.begin(), params.end(), n);
    if (it == params.end()) throw kernel::Refusal("'" + n + "' is not a parameter of " + name);
    param_values[static_cast<std::size_t>(it - params.begin())] = v;
  }

  // a(t, x), and for a map f(n, x, ξ) with the draws in `noise`.
  void drift(double t, const State& x, State& a, const double* noise = nullptr) const {
    load(t, x, noise);
    for (std::size_t i = 0; i < drift_.size(); i++) a[i] = drift_[i].run(values_.data(), stack_.data());
  }
  // B(t, x), row by row: B[i * m + j] = b_ij.
  void diffusion(double t, const State& x, std::vector<double>& B) const {
    load(t, x, nullptr);
    const std::size_t m = sources.size();
    B.assign(dimension() * m, 0.0);
    for (const auto& e : diffusion_) B[e.row * m + e.column] = e.program.run(values_.data(), stack_.data());
  }
  // One step of length h with the increments dW (one for each source).
  void step(double t, State& x, double h, const double* dW, bool milstein) const {
    const std::size_t n = dimension(), m = sources.size();
    load(t, x, nullptr);
    scratch_.assign(n, 0.0);
    for (std::size_t i = 0; i < n; i++) scratch_[i] = drift_[i].run(values_.data(), stack_.data()) * h;
    for (const auto& e : diffusion_) scratch_[e.row] += e.program.run(values_.data(), stack_.data()) * dW[e.column];
    if (milstein)
      for (const auto& e : correction_) scratch_[e.row] += 0.5 * e.program.run(values_.data(), stack_.data()) * (dW[e.column] * dW[e.column] - h);
    for (std::size_t i = 0; i < n; i++) x[i] += scratch_[i];
    (void)m;
  }

 private:
  friend Stochastic stochastic(const System&);
  void load(double t, const State& x, const double* noise) const {
    const std::size_t n = states.size(), p = params.size(), m = sources.size();
    if (values_.size() != n + p + 1 + m) values_.assign(n + p + 1 + m, 0.0);
    std::copy(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(n), values_.begin());
    std::copy(param_values.begin(), param_values.end(), values_.begin() + static_cast<std::ptrdiff_t>(n));
    values_[n + p] = t;
    for (std::size_t j = 0; j < m; j++) values_[n + p + 1 + j] = noise ? noise[j] : 0.0;
  }
  struct Compiled {
    std::size_t row, column;
    terms::Program program;
  };
  std::vector<terms::Program> drift_;
  std::vector<Compiled> diffusion_, correction_;  // b_ij, and b_i ∂b_i/∂x_i for Milstein
  mutable std::vector<double> values_, stack_, scratch_;
};

inline Stochastic stochastic(const System& s) {
  if (form_of(s) != Form::stochastic) throw kernel::Refusal(s.name + " is a " + std::string(name(form_of(s))) + " system, not a stochastic one");
  const Scope sc = scope(s);
  Stochastic d;
  d.name = s.name, d.discrete = sc.discrete, d.time = sc.time;
  d.calculus = s.structure.calculus.empty() ? "ito" : s.structure.calculus;
  d.states = sc.states, d.params = sc.params, d.param_values = sc.param_values;
  d.sources = sc.sources, d.initial = sc.initial;
  const std::size_t n = d.states.size(), m = d.sources.size();
  terms::Frame frame;
  for (const auto& st : d.states) frame.add(st);
  for (const auto& p : d.params) frame.add(p);
  frame.add(d.time);
  for (const auto& src : d.sources) frame.add(src.name);
  std::size_t deepest = 1;
  const auto compiled = [&](const TermPtr& t, language::Place at, const std::string& what) {
    try {
      terms::Program p = terms::compile(t, frame.resolver());
      deepest = std::max(deepest, p.depth());
      return p;
    } catch (const kernel::Refusal& e) {
      throw ScriptError(at, what + ": " + e.what());
    }
  };
  d.drift_terms.assign(n, nullptr);
  d.diffusion_terms.assign(n, std::vector<TermPtr>(m, terms::number(0)));
  const auto index_of = [&](const std::string& st, language::Place at) {
    const int i = d.state_index(st);
    if (i < 0) throw ScriptError(at, "'" + st + "' is not a state of " + s.name);
    return static_cast<std::size_t>(i);
  };

  if (d.discrete) {
    // A map with draws: x' = f(n, x, ξ).
    for (const auto& e : s.change.equations) {
      const auto v = sc.vectors.find(e.name);
      const std::size_t count = v == sc.vectors.end() ? 1 : v->second.size();
      for (std::size_t k = 0; k < count; k++) {
        const std::size_t i = index_of(v == sc.vectors.end() ? e.name : v->second[k], e.at);
        if (e.primes != 1) throw ScriptError(e.at, "in discrete time write " + e.name + "' for the next value");
        d.drift_terms[i] = sc.resolve(v == sc.vectors.end() ? e.value : component(sc, e.value, k));
      }
    }
  } else {
    const std::string dt = "d" + d.time;
    for (const auto& e : s.change.equations) {
      if (e.primes != 0 || e.name.size() < 2 || e.name[0] != 'd') throw ScriptError(e.at, "a stochastic differential is written d" + (e.name.empty() ? std::string("X") : e.name) + " = a * " + dt + " + b * dW");
      const std::string target = e.name.substr(1);
      const auto v = sc.vectors.find(target);
      const std::size_t count = v == sc.vectors.end() ? 1 : v->second.size();
      for (std::size_t k = 0; k < count; k++) {
        const std::size_t i = index_of(v == sc.vectors.end() ? target : v->second[k], e.at);
        const TermPtr rhs = sc.resolve(v == sc.vectors.end() ? e.value : component(sc, e.value, k, "d"));
        TermPtr rest = rhs;
        const auto coefficient = [&](const std::string& differential) {
          const TermPtr c = terms::derivative(rhs, differential);
          rest = rest - c * terms::symbol(differential);
          return c;
        };
        d.drift_terms[i] = coefficient(dt);
        for (std::size_t j = 0; j < m; j++) d.diffusion_terms[i][j] = coefficient("d" + d.sources[j].name);
        // What is left must vanish, and no coefficient may hold a differential: the right side is linear in them.
        bool linear = true;
        std::set<std::pair<std::string, unsigned>> used = terms::symbols(d.drift_terms[i]);
        for (const auto& b : d.diffusion_terms[i]) terms::symbols(b, used);
        for (const auto& [nm, pr] : used)
          if (nm == dt || std::any_of(d.sources.begin(), d.sources.end(), [&](const Source& src) { return nm == "d" + src.name; })) linear = false;
        if (linear && !terms::is_number(rest, 0)) {
          // The normal form may not have cancelled; decide by values at a point.
          std::size_t salt = 0;
          const double left = terms::evaluate(rest, [&](const std::string&, unsigned) { return 0.37 + 0.11 * static_cast<double>(++salt % 7); });
          linear = std::fabs(left) < 1e-9;
        }
        if (!linear) throw ScriptError(e.at, "the right side of d" + target + " is not linear in the differentials " + dt + " and dW");
      }
    }
  }
  for (std::size_t i = 0; i < n; i++)
    if (!d.drift_terms[i]) throw ScriptError(s.at, "no change is given for the state '" + d.states[i] + "' of " + s.name);

  if (!d.discrete && d.calculus == "stratonovich") {
    for (std::size_t i = 0; i < n; i++) {
      std::vector<TermPtr> parts{d.drift_terms[i]};
      for (std::size_t j = 0; j < m; j++)
        for (std::size_t k = 0; k < n; k++) parts.push_back(terms::number(0.5) * d.diffusion_terms[k][j] * terms::derivative(d.diffusion_terms[i][j], d.states[k]));
      d.drift_terms[i] = terms::sum(std::move(parts));
    }
  }
  // Diagonal noise: each source drives one component, each component is driven by at most one source, and no
  // coefficient depends on another driven component (b_j ∂b_i/∂x_j = 0 for j ≠ i, the condition under which the
  // noise commutes and Milstein's scheme needs no iterated integral; Kloeden and Platen 1992, §10.3). Which
  // components are driven is found first, for every source, and the dependence is asked after.
  d.diagonal = !d.discrete;
  std::vector<int> source_of(n, -1);
  for (std::size_t j = 0; j < m && d.diagonal; j++) {
    int driven = -1;
    for (std::size_t i = 0; i < n; i++)
      if (!terms::is_number(d.diffusion_terms[i][j], 0)) {
        if (driven >= 0 || source_of[i] >= 0) d.diagonal = false;
        driven = static_cast<int>(i), source_of[i] = static_cast<int>(j);
      }
  }
  for (std::size_t i = 0; i < n && d.diagonal; i++) {
    if (source_of[i] < 0) continue;
    for (std::size_t k = 0; k < n; k++)
      if (k != i && source_of[k] >= 0 && terms::depends_on(d.diffusion_terms[i][static_cast<std::size_t>(source_of[i])], d.states[k], 0)) d.diagonal = false;
  }
  for (std::size_t i = 0; i < n; i++) {
    d.drift_.push_back(compiled(d.drift_terms[i], s.at, "the change of " + d.states[i]));
    for (std::size_t j = 0; j < m; j++) {
      const TermPtr& b = d.diffusion_terms[i][j];
      if (terms::is_number(b, 0)) continue;
      d.diffusion_.push_back({i, j, compiled(b, s.at, "the noise of " + d.states[i])});
      if (d.diagonal) {
        const TermPtr c = b * terms::derivative(b, d.states[i]);
        if (!terms::is_number(c, 0)) d.correction_.push_back({i, j, compiled(c, s.at, "Milstein's term of " + d.states[i])});
      }
    }
  }
  d.stack_.assign(deepest, 0.0);
  return d;
}

// ---- draws ----

// A draw of a source's law for a map. Symmetric α-stable draws are Chambers, Mallows and Stuck's (1976, J. Amer.
// Statist. Assoc. 71, 340–344, the symmetric case): with V uniform on (−π/2, π/2) and E exponential,
// X = sin(αV) / cos(V)^{1/α} · (cos((1 − α)V) / E)^{(1 − α)/α}; α = 2 gives a normal of variance 2, α = 1 Cauchy.
inline double draw(const Source& s, calx::random::Stream& r) {
  const auto arg = [&](std::size_t i, double fallback) { return i < s.numbers.size() ? s.numbers[i] : fallback; };
  if (s.kind == "normal" || s.kind == "wiener") return r.normal(arg(0, 0), arg(1, 1));
  if (s.kind == "uniform") return r.uniform(arg(0, 0), arg(1, 1));
  if (s.kind == "exponential") return r.exponential() / arg(0, 1);
  if (s.kind == "poisson") return static_cast<double>(r.poisson(arg(0, 1)));
  if (s.kind == "cauchy") return arg(0, 0) + arg(1, 1) * std::tan(std::numbers::pi * (r.uniform() - 0.5));
  if (s.kind == "levy") {
    const double alpha = arg(0, 1.5);
    if (!(alpha > 0 && alpha <= 2)) throw kernel::Refusal("the index of a stable law lies in (0, 2]");
    const double V = std::numbers::pi * (r.uniform() - 0.5), E = r.exponential();
    if (alpha == 1) return std::tan(V);
    return std::sin(alpha * V) / std::pow(std::cos(V), 1 / alpha) * std::pow(std::cos((1 - alpha) * V) / E, (1 - alpha) / alpha);
  }
  throw kernel::Refusal("the law '" + s.kind + "' of the source " + s.name + " is not known");
}

// Of the options of a run, a stochastic one uses the times, the step, `every` and `transient`; a system with
// noise has no events here, so `event_tolerance` and `max_events` are not read.
struct StochasticOptions : Options {
  std::uint64_t seed = 1;
  std::uint32_t path = 0;             // which path of the ensemble: its draws depend on (seed, path) alone
  bool milstein = false;
};

struct StochasticTrajectory : Trajectory {
  std::vector<State> w;  // the driving paths W_j at the recorded times (continuous time)
};

// A path. `increments`, when given, holds ΔW for each step and source (steps × sources) and replaces the streams.
inline StochasticTrajectory simulate(const Stochastic& d, State x, const StochasticOptions& o = {}, const std::vector<std::vector<double>>* increments = nullptr) {
  if (x.size() != d.dimension()) throw kernel::Refusal(d.name + " has " + std::to_string(d.dimension()) + " components; " + std::to_string(x.size()) + " were given");
  if (o.milstein && !d.diagonal) throw kernel::Refusal("Milstein's scheme is given for diagonal noise; the noise of " + d.name + " is not diagonal");
  const std::size_t m = d.sources.size();
  std::vector<calx::random::Stream> streams;
  for (std::size_t j = 0; j < m; j++) streams.emplace_back(o.seed, 0x6b696e65u, o.path, static_cast<std::uint32_t>(j));  // 'kine'
  StochasticTrajectory out;
  out.names = d.states;
  const std::size_t every = std::max<std::size_t>(o.every, 1);
  std::vector<double> dW(std::max<std::size_t>(m, 1), 0.0);
  State W(m, 0.0);
  double t = o.t0;
  // The transient: the system is run, with draws of its own streams, before recording begins; the driving paths
  // are counted from where recording begins. A run that is given its increments is given them for what is recorded.
  if (o.transient > 0) {
    if (increments) throw kernel::Refusal("a run that is given its increments has no transient: begin it where the increments begin");
    if (d.discrete) {
      State next(x.size());
      for (double k = 0; k < o.transient; k += 1) {
        for (std::size_t j = 0; j < m; j++) dW[j] = draw(d.sources[j], streams[j]);
        d.drift(t, x, next, dW.data());
        x.swap(next), t += 1;
      }
    } else {
      if (!(o.step > 0)) throw kernel::Refusal("the step of a stochastic differential is positive");
      for (double s = 0; s < o.transient - 0.5 * o.step; s += o.step) {
        for (std::size_t j = 0; j < m; j++) dW[j] = std::sqrt(o.step) * streams[j].normal();
        d.step(t, x, o.step, dW.data(), o.milstein);
        t += o.step;
      }
    }
    t = o.t0;
  }
  out.t.push_back(t), out.x.push_back(x), out.w.push_back(W);
  if (d.discrete) {
    const std::size_t count = static_cast<std::size_t>(std::max(0.0, o.t1 - o.t0));
    State next(x.size());
    for (std::size_t k = 1; k <= count; k++) {
      for (std::size_t j = 0; j < m; j++) dW[j] = increments ? (*increments)[k - 1][j] : draw(d.sources[j], streams[j]);
      d.drift(t, x, next, dW.data());
      x.swap(next), t += 1;
      if (k % every == 0 || k == count) out.t.push_back(t), out.x.push_back(x), out.w.push_back(W);
    }
    out.steps = count;
    return out;
  }
  if (!(o.step > 0)) throw kernel::Refusal("the step of a stochastic differential is positive");
  const std::size_t count = static_cast<std::size_t>(std::ceil((o.t1 - o.t0) / o.step - 1e-9));
  for (std::size_t k = 1; k <= count; k++) {
    const double h = std::min(o.step, o.t1 - t);
    for (std::size_t j = 0; j < m; j++) dW[j] = increments ? (*increments)[k - 1][j] : std::sqrt(h) * streams[j].normal(), W[j] += dW[j];
    d.step(t, x, h, dW.data(), o.milstein);
    t += h;
    if (k % every == 0 || k == count) out.t.push_back(t), out.x.push_back(x), out.w.push_back(W);
    for (double v : x)
      if (!std::isfinite(v)) throw kernel::Refusal(d.name + " left the finite numbers at t = " + terms::number_text(t));
  }
  out.steps = count;
  return out;
}

}  // namespace kinesis::systems
