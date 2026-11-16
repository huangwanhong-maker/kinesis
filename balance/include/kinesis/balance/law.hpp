// A balance law on a graph as a sum of terms, each with its coefficient, and its run (decision 0005).
//
// A law is du/dt = Σ_k θ_k T_k[u], every term T_k an operator on the field and every coefficient θ_k a number. The terms:
//
//   constant      1            a density made everywhere at one rate
//   linear        u            growth at a rate (Malthus)
//   quadratic     u²           with linear, logistic growth: θ₁ u + θ₂ u² = r u (1 − u/K), r = θ₁, K = −θ₁/θ₂ (Verhulst 1838)
//   gompertz      u ln u       with linear, Gompertz's growth: θ₁ u + θ₂ u ln u = r u ln(K/u), r = −θ₂, K = e^{−θ₁/θ₂}
//   diffusion     Δu           spread down the gradient (Fick); with logistic growth, Fisher and Kolmogorov, Petrovsky and
//                              Piskunov's equation (Fisher 1937, Ann. Eugen. 7, 355–369), whose fronts run at 2√(rD)
//   crowding      Δ(u²/2)      a diffusivity that grows with the density (the porous-medium equation)
//   drift         −∇·(u∇Φ_k)   carried up the potential Φ_k, one of the covariates given with the law
//   source        Φ_k          made in proportion to a covariate
//
// Every law is linear in its coefficients, so that learning them from a record is least squares (learn.hpp), and every
// flux term is conservative (graph.hpp). A run integrates the law by the classical fourth-order Runge–Kutta step on the
// whole field, with a step held within the stability of the explicit scheme: the step is at most 2/ρ, where ρ bounds
// the spectral radius of the law's Jacobian by Gershgorin's theorem at the field's current values (the region of the
// fourth-order method on the negative axis reaches −2.78; Hairer and Wanner 1996, Solving ODEs II, §IV.2). A density
// that a step would make negative is set to zero, and the quantity so added is reported (`clipped`): the law then left
// what it can describe.
#pragma once

#include <kinesis/balance/graph.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace kinesis::balance {

enum class Kind { constant, linear, quadratic, gompertz, diffusion, crowding, drift, source };

struct Term {
  Kind kind = Kind::linear;
  std::size_t covariate = 0;  // for drift and source: which covariate
};

inline std::string symbol(const Term& t, const std::vector<std::string>& covariates = {}) {
  const std::string phi = t.covariate < covariates.size() ? covariates[t.covariate] : "Φ" + std::to_string(t.covariate + 1);
  switch (t.kind) {
    case Kind::constant: return "1";
    case Kind::linear: return "u";
    case Kind::quadratic: return "u²";
    case Kind::gompertz: return "u ln u";
    case Kind::diffusion: return "Δu";
    case Kind::crowding: return "Δ(u²/2)";
    case Kind::drift: return "−∇·(u ∇" + phi + ")";
    default: return phi;
  }
}

inline const char* name(Kind k) {
  switch (k) {
    case Kind::constant: return "constant";
    case Kind::linear: return "linear";
    case Kind::quadratic: return "quadratic";
    case Kind::gompertz: return "gompertz";
    case Kind::diffusion: return "diffusion";
    case Kind::crowding: return "crowding";
    case Kind::drift: return "drift";
    default: return "source";
  }
}

inline Kind kind_of(const std::string& s) {
  for (Kind k : {Kind::constant, Kind::linear, Kind::quadratic, Kind::gompertz, Kind::diffusion, Kind::crowding, Kind::drift, Kind::source})
    if (s == name(k)) return k;
  throw Refusal("there is no term '" + s + "' (there are constant, linear, quadratic, gompertz, diffusion, crowding, drift and source)");
}

inline bool moves(Kind k) { return k == Kind::diffusion || k == Kind::crowding || k == Kind::drift; }

struct Law {
  std::vector<Term> terms;
  std::vector<double> coefficients;  // one per term
};

using Covariates = std::vector<std::vector<double>>;  // each a value per cell

inline void require_law(const Graph& g, const Law& law, const Covariates& cov) {
  if (law.terms.size() != law.coefficients.size()) throw Refusal("a law needs one coefficient per term");
  for (const Term& t : law.terms)
    if (t.kind == Kind::drift || t.kind == Kind::source) {
      if (t.covariate >= cov.size()) throw Refusal("a term asks for covariate " + std::to_string(t.covariate + 1) + " and " + std::to_string(cov.size()) + " were given");
      require_size(g, cov[t.covariate], "a covariate");
    }
}

// T[u] for one term.
inline std::vector<double> evaluate(const Graph& g, const Term& t, const std::vector<double>& u, const Covariates& cov) {
  require_size(g, u, "the field");
  std::vector<double> out(u.size());
  switch (t.kind) {
    case Kind::constant: std::fill(out.begin(), out.end(), 1.0); return out;
    case Kind::linear: return u;
    case Kind::quadratic: for (std::size_t i = 0; i < u.size(); i++) out[i] = u[i] * u[i]; return out;
    case Kind::gompertz: for (std::size_t i = 0; i < u.size(); i++) out[i] = u[i] > 0 ? u[i] * std::log(u[i]) : 0.0; return out;
    case Kind::diffusion: return laplacian(g, u);
    case Kind::crowding: return crowding(g, u);
    case Kind::drift:
      if (t.covariate >= cov.size()) throw Refusal("the drift asks for a covariate that was not given");
      return drift(g, u, cov[t.covariate]);
    default:
      if (t.covariate >= cov.size()) throw Refusal("the source asks for a covariate that was not given");
      require_size(g, cov[t.covariate], "a covariate");
      return cov[t.covariate];
  }
}

// du/dt under the law, with an extra rate per cell (a planned source, or nothing).
inline std::vector<double> rate(const Graph& g, const Law& law, const std::vector<double>& u, const Covariates& cov,
                                const std::vector<double>& extra = {}) {
  std::vector<double> out = extra.empty() ? std::vector<double>(u.size(), 0.0) : extra;
  if (out.size() != u.size()) throw Refusal("the extra rate is not one value per cell");
  for (std::size_t k = 0; k < law.terms.size(); k++) {
    if (law.coefficients[k] == 0) continue;
    const std::vector<double> t = evaluate(g, law.terms[k], u, cov);
    for (std::size_t i = 0; i < u.size(); i++) out[i] += law.coefficients[k] * t[i];
  }
  return out;
}

// The flux of each edge, from a to b, in units of the quantity per unit time: what the law's moving terms carry.
inline std::vector<double> flux(const Graph& g, const Law& law, const std::vector<double>& u, const Covariates& cov) {
  require_law(g, law, cov);
  require_size(g, u, "the field");
  std::vector<double> f(g.edges.size(), 0.0);
  for (std::size_t k = 0; k < law.terms.size(); k++) {
    const Term& t = law.terms[k];
    const double th = law.coefficients[k];
    if (th == 0 || !moves(t.kind)) continue;
    for (std::size_t e = 0; e < g.edges.size(); e++) {
      const Edge& ed = g.edges[e];
      const double ua = u[ed.a], ub = u[ed.b];
      if (t.kind == Kind::diffusion) f[e] += -th * ed.conductance * (ub - ua);
      else if (t.kind == Kind::crowding) f[e] += -th * ed.conductance * 0.5 * (ub * ub - ua * ua);
      else f[e] += th * ed.conductance * (cov[t.covariate][ed.b] - cov[t.covariate][ed.a]) * 0.5 * (ua + ub);
    }
  }
  return f;
}

// A bound on the spectral radius of the law's Jacobian at u, by Gershgorin's discs.
inline double stiffness(const Graph& g, const Law& law, const std::vector<double>& u, const Covariates& cov) {
  double umax = 0;
  for (double v : u) umax = std::max(umax, std::abs(v));
  std::vector<double> row(g.size(), 0.0);  // Σ_e c_e × the edge's coefficient of exchange, per cell
  double react = 0;
  for (std::size_t k = 0; k < law.terms.size(); k++) {
    const double th = std::abs(law.coefficients[k]);
    if (th == 0) continue;
    const Term& t = law.terms[k];
    switch (t.kind) {
      case Kind::linear: react += th; break;
      case Kind::quadratic: react += 2 * th * umax; break;
      case Kind::gompertz: react += th * (std::abs(std::log(std::max(umax, 1e-12))) + 1 + 30); break;  // ln u + 1, with room for small u
      case Kind::diffusion: for (const Edge& e : g.edges) row[e.a] += th * e.conductance, row[e.b] += th * e.conductance; break;
      case Kind::crowding: for (const Edge& e : g.edges) {
          const double m = th * e.conductance * std::max(std::abs(u[e.a]), std::abs(u[e.b]));
          row[e.a] += m, row[e.b] += m;
        } break;
      case Kind::drift: for (const Edge& e : g.edges) {
          const double m = th * e.conductance * std::abs(cov[t.covariate][e.b] - cov[t.covariate][e.a]);
          row[e.a] += m, row[e.b] += m;
        } break;
      default: break;
    }
  }
  double rho = 0;
  for (std::size_t i = 0; i < row.size(); i++) rho = std::max(rho, 2 * row[i] / g.volume[i]);
  return rho + react;
}

struct RunOptions {
  double most_step = std::numeric_limits<double>::infinity();  // the largest step allowed, besides stability's
  bool clip = true;                                             // set a density made negative to zero, and count it
  std::size_t most_steps = 2'000'000;
};

struct Run {
  std::vector<double> t;
  std::vector<std::vector<double>> u;  // the field at each time asked for
  double clipped = 0;                  // the content added by setting negative densities to zero
  std::size_t steps = 0;
};

// The field at each of `times` (non-decreasing, none before t0), from u0 at t0.
inline Run simulate(const Graph& g, const Law& law, const std::vector<double>& u0, double t0, const std::vector<double>& times,
                    const Covariates& cov = {}, const std::vector<double>& extra = {}, const RunOptions& o = {}) {
  require_law(g, law, cov);
  require_size(g, u0, "the initial field");
  if (!extra.empty()) require_size(g, extra, "the extra rate");
  Run run;
  std::vector<double> u = u0, k1, k2, k3, k4, w(u.size());
  double t = t0;
  for (double target : times) {
    if (target < t - 1e-12 * std::max(1.0, std::abs(t))) throw Refusal("the times asked for must not decrease, nor come before the start");
    while (t < target) {
      const double rho = stiffness(g, law, u, cov);
      double h = std::min(target - t, o.most_step);
      if (rho > 0) h = std::min(h, 2.0 / rho);
      if (!(h > 0)) break;
      if (++run.steps > o.most_steps) throw Refusal("the run needed more than " + std::to_string(o.most_steps) + " steps: the law is too stiff for an explicit step");
      k1 = rate(g, law, u, cov, extra);
      for (std::size_t i = 0; i < u.size(); i++) w[i] = u[i] + 0.5 * h * k1[i];
      k2 = rate(g, law, w, cov, extra);
      for (std::size_t i = 0; i < u.size(); i++) w[i] = u[i] + 0.5 * h * k2[i];
      k3 = rate(g, law, w, cov, extra);
      for (std::size_t i = 0; i < u.size(); i++) w[i] = u[i] + h * k3[i];
      k4 = rate(g, law, w, cov, extra);
      for (std::size_t i = 0; i < u.size(); i++) {
        u[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        if (!std::isfinite(u[i])) throw Refusal("the field became infinite near t = " + std::to_string(t) + ": the law blows up from this start");
        if (o.clip && u[i] < 0) run.clipped += -u[i] * g.volume[i], u[i] = 0;
      }
      t = (target - t - h <= 1e-12 * std::max(1.0, std::abs(target))) ? target : t + h;
    }
    run.t.push_back(target);
    run.u.push_back(u);
  }
  return run;
}

}  // namespace kinesis::balance
