// Unknowns of a system across scales filled from records at any of its levels, and the fit judged on what it was not
// shown (decision 0007).
//
// A parameter marked to be fitted is an unknown. Its value is found by Levenberg–Marquardt (Levenberg 1944, Q. Appl. Math.
// 2, 164–168; Marquardt 1963, J. SIAM 11, 431–441) on the residuals (model − record)/σ of every observation before a time,
// each observation a quantity at one unit and one time, of any level: a city's total and a block's value can be fitted
// together. Each run is a shooting from the initial values, so a fit is a run fitted to the record, not a law fitted to
// differences. Bounds are kept by a change of variable (a logistic between two bounds, an exponential above one).
//
// What the record says of the unknowns. Standard errors come from the Jacobian at the optimum, scaled by the residual
// variance, s² (JᵀJ)⁻¹ with s² = RSS/(n − k): a local, linearised statement. An unknown whose standard error is more than
// half its value, or that moves with another (|correlation| > 0.95), is reported as not identified by this record: a
// quantity to be declared, or a record to be found, rather than a number to be believed.
//
// Judgement. Observations after a time are held aside. The fitted model is judged on them against persistence (the last
// value seen before that time, carried forward, unit by unit): skill = 1 − SSE_model/SSE_persistence, positive when the
// model forecasts better than standing still.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "calx/linalg/dense.hpp"
#include "calx/linalg/eigen.hpp"
#include "calx/linalg/rectangular.hpp"
#include "kinesis/scales/run.hpp"

namespace kinesis::scales {

struct Observation {
  std::string quantity;  // a quantity or derived quantity of the model
  std::string unit;      // a unit of its level
  double time = 0;
  double value = 0;
  double sigma = 1;
};

struct FitOptions {
  double held_after = std::numeric_limits<double>::infinity();  // observations after this time are held aside
  std::size_t most = 40;                                         // iterations
};

struct Fitted {
  std::vector<std::string> names;
  std::vector<double> values, standard_errors;
  std::vector<std::vector<double>> correlation;
  double condition = 0;                    // of the correlation matrix of the estimates
  std::vector<std::string> not_identified;  // each with why
  std::vector<double> predicted;           // for every observation, fitted or held aside
  std::vector<char> held;                  // 1 where the observation was held aside
  double rms = 0;                          // of the weighted residuals fitted
  double held_rms = std::numeric_limits<double>::quiet_NaN(), persistence_rms = std::numeric_limits<double>::quiet_NaN();
  double skill = std::numeric_limits<double>::quiet_NaN();
  std::size_t fitted = 0, held_count = 0, iterations = 0;
  bool converged = false;
  std::string method;
  Inputs inputs;  // the inputs with the unknowns set as estimated, with their standard errors
};

namespace detail {

struct Transform {
  double lower, upper;
  double to_x(double th) const {
    if (std::isfinite(lower) && std::isfinite(upper)) return lower + (upper - lower) / (1 + std::exp(-th));
    if (std::isfinite(lower)) return lower + std::exp(th);
    if (std::isfinite(upper)) return upper - std::exp(th);
    return th;
  }
  double to_theta(double x) const {
    if (std::isfinite(lower) && std::isfinite(upper)) {
      const double f = std::clamp((x - lower) / (upper - lower), 1e-9, 1 - 1e-9);
      return std::log(f / (1 - f));
    }
    if (std::isfinite(lower)) return std::log(std::max(x - lower, 1e-12));
    if (std::isfinite(upper)) return std::log(std::max(upper - x, 1e-12));
    return x;
  }
  double dx_dtheta(double th) const {
    if (std::isfinite(lower) && std::isfinite(upper)) {
      const double s = 1 / (1 + std::exp(-th));
      return (upper - lower) * s * (1 - s);
    }
    if (std::isfinite(lower)) return std::exp(th);
    if (std::isfinite(upper)) return -std::exp(th);
    return 1;
  }
};

}  // namespace detail

inline Fitted fit(const Compiled& c, const Inputs& in, Options o, const std::vector<Observation>& obs, const FitOptions& fo = {}) {
  Fitted out;
  std::vector<detail::Transform> tr;
  std::vector<double> theta;
  for (const auto& p : c.model.parameters) {
    if (!p.fit) continue;
    out.names.push_back(p.name);
    tr.push_back({p.lower, p.upper});
    double start = p.value;
    const auto it = in.parameters.find(p.name);
    if (it != in.parameters.end() && std::isfinite(centre(it->second))) start = centre(it->second);
    if (!std::isfinite(start)) {
      if (std::isfinite(p.lower) && std::isfinite(p.upper)) start = 0.5 * (p.lower + p.upper);
      else if (std::isfinite(p.lower)) start = p.lower + 1;
      else if (std::isfinite(p.upper)) start = p.upper - 1;
      else start = 0;
    }
    theta.push_back(tr.back().to_theta(start));
  }
  const std::size_t k = theta.size();
  if (k == 0) throw Refusal("the model has no unknown to fit: mark a parameter to be fitted");
  if (obs.empty()) throw Refusal("a fit needs observations");

  // Where each observation is read from a run.
  std::set<double> times;
  std::set<std::string> names;
  std::vector<std::size_t> unit_of(obs.size());
  for (std::size_t i = 0; i < obs.size(); i++) {
    const auto lv = c.level_of.find(obs[i].quantity);
    if (lv == c.level_of.end()) throw Refusal("an observation is of " + obs[i].quantity + ", which the model does not have");
    unit_of[i] = c.hierarchy.unit(lv->second, obs[i].unit);
    if (!(obs[i].sigma > 0)) throw Refusal("an observation's σ must be positive");
    times.insert(obs[i].time);
    names.insert(obs[i].quantity);
  }
  o.report.assign(times.begin(), times.end());
  o.record.assign(names.begin(), names.end());
  out.held.assign(obs.size(), 0);
  for (std::size_t i = 0; i < obs.size(); i++) {
    out.held[i] = obs[i].time > fo.held_after ? 1 : 0;
    (out.held[i] ? out.held_count : out.fitted)++;
  }
  if (out.fitted <= k) throw Refusal("the fit has " + std::to_string(out.fitted) + " observations for " + std::to_string(k) + " unknowns");

  const auto with = [&](const std::vector<double>& th) {
    Inputs x = in;
    for (std::size_t j = 0; j < k; j++) x.parameters[out.names[j]] = Value{tr[j].to_x(th[j]), Source::estimated, 0, missing_value, missing_value};
    return x;
  };
  // Predictions of every observation, or nothing when the model refuses this point.
  const auto predict = [&](const std::vector<double>& th) -> std::optional<std::vector<double>> {
    try {
      const Run r = simulate(c, with(th), o);
      if (r.counts.unconverged) return std::nullopt;
      std::vector<double> p(obs.size());
      for (std::size_t i = 0; i < obs.size(); i++) {
        const auto ti = std::find(r.times.begin(), r.times.end(), obs[i].time) - r.times.begin();
        p[i] = r.values.at(obs[i].quantity)[static_cast<std::size_t>(ti)][unit_of[i]];
        if (!std::isfinite(p[i])) return std::nullopt;
      }
      return p;
    } catch (const Refusal&) {
      return std::nullopt;
    }
  };
  const auto residuals = [&](const std::vector<double>& p) {
    std::vector<double> r;
    for (std::size_t i = 0; i < obs.size(); i++)
      if (!out.held[i]) r.push_back((p[i] - obs[i].value) / obs[i].sigma);
    return r;
  };
  const auto cost = [](const std::vector<double>& r) {
    double s = 0;
    for (double v : r) s += v * v;
    return s;
  };

  auto p0 = predict(theta);
  if (!p0) throw Refusal("the model does not run at the unknowns' starting values; give better starts or bounds");
  std::vector<double> r = residuals(*p0);
  double C = cost(r);
  double lambda = 1e-3;
  const std::size_t m = r.size();
  calx::linalg::Rectangular J(m, k);
  for (std::size_t it = 0; it < fo.most; it++) {
    out.iterations = it + 1;
    // The Jacobian by forward differences in the unbounded variables.
    for (std::size_t j = 0; j < k; j++) {
      std::vector<double> th = theta;
      const double h = 1e-6 * (1 + std::abs(th[j]));
      th[j] += h;
      auto pj = predict(th);
      if (!pj) {
        th[j] = theta[j] - h;
        pj = predict(th);
        if (!pj) throw Refusal("the model does not run near the unknowns' current values");
        const std::vector<double> rj = residuals(*pj);
        for (std::size_t i = 0; i < m; i++) J(i, j) = (r[i] - rj[i]) / h;
      } else {
        const std::vector<double> rj = residuals(*pj);
        for (std::size_t i = 0; i < m; i++) J(i, j) = (rj[i] - r[i]) / h;
      }
    }
    calx::linalg::Matrix A(k);
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t b = 0; b < k; b++) {
        double s = 0;
        for (std::size_t i = 0; i < m; i++) s += J(i, a) * J(i, b);
        A(a, b) = s;
      }
    std::vector<double> grad(k, 0.0);
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t i = 0; i < m; i++) grad[a] += J(i, a) * r[i];
    bool stepped = false;
    for (int tries = 0; tries < 12 && !stepped; tries++) {
      calx::linalg::Matrix M = A, rhs(k);
      for (std::size_t a = 0; a < k; a++) {
        M(a, a) += lambda * std::max(A(a, a), 1e-12);
        rhs(a, 0) = -grad[a];
      }
      const calx::linalg::Matrix d = calx::linalg::solve(M, rhs);
      std::vector<double> th = theta;
      double size = 0;
      for (std::size_t a = 0; a < k; a++) th[a] += d(a, 0), size += d(a, 0) * d(a, 0);
      const auto pn = predict(th);
      if (pn) {
        const std::vector<double> rn = residuals(*pn);
        const double Cn = cost(rn);
        if (Cn < C) {
          const double gain = (C - Cn) / std::max(C, 1e-300);
          theta = th, r = rn, C = Cn, p0 = pn;
          lambda = std::max(lambda / 3, 1e-12);
          stepped = true;
          if (gain < 1e-12 || std::sqrt(size) < 1e-10) out.converged = true;
          continue;
        }
      }
      lambda *= 4;
    }
    if (!stepped) out.converged = true;  // no step lowers the cost: a minimum, to the precision of the differences
    if (out.converged) break;
  }
  out.method = "Levenberg–Marquardt on runs from the initial values (shooting), Jacobian by forward differences";
  out.rms = std::sqrt(C / static_cast<double>(m));
  out.values.resize(k);
  for (std::size_t j = 0; j < k; j++) out.values[j] = tr[j].to_x(theta[j]);

  // Standard errors in the original variables.
  calx::linalg::Matrix A(k);
  for (std::size_t a = 0; a < k; a++)
    for (std::size_t b = 0; b < k; b++) {
      double s = 0;
      for (std::size_t i = 0; i < m; i++) s += J(i, a) * J(i, b);
      A(a, b) = s / (tr[a].dx_dtheta(theta[a]) * tr[b].dx_dtheta(theta[b]));
    }
  const double s2 = m > k ? C / static_cast<double>(m - k) : 1.0;
  calx::linalg::Matrix D(k);
  for (std::size_t a = 0; a < k; a++) D(a, a) = A(a, a) > 0 ? 1 / std::sqrt(A(a, a)) : 0;
  calx::linalg::Matrix S = D * A * D;  // the scaled normal matrix, unit diagonal where defined
  const auto eig = calx::linalg::symmetric_eigen(S);
  const double lo = eig.values.front(), hi = eig.values.back();
  out.condition = lo > 0 ? hi / lo : std::numeric_limits<double>::infinity();
  out.standard_errors.assign(k, std::numeric_limits<double>::infinity());
  out.correlation.assign(k, std::vector<double>(k, std::numeric_limits<double>::quiet_NaN()));
  if (lo > 1e-12 * std::max(hi, 1.0)) {
    const calx::linalg::Matrix Si = calx::linalg::inverse(S);
    for (std::size_t a = 0; a < k; a++) {
      out.standard_errors[a] = std::sqrt(s2 * Si(a, a)) * D(a, a);
      for (std::size_t b = 0; b < k; b++) out.correlation[a][b] = Si(a, b) / std::sqrt(Si(a, a) * Si(b, b));
    }
    for (std::size_t a = 0; a < k; a++) {
      if (out.standard_errors[a] > 0.5 * std::abs(out.values[a]))
        out.not_identified.push_back(out.names[a] + ": its standard error is more than half its value");
      for (std::size_t b = a + 1; b < k; b++)
        if (std::abs(out.correlation[a][b]) > 0.95)
          out.not_identified.push_back(out.names[a] + " and " + out.names[b] + " move together (correlation " +
                                       std::to_string(out.correlation[a][b]).substr(0, 5) + "): the record does not tell them apart");
    }
  } else {
    for (const auto& n : out.names) out.not_identified.push_back(n + ": the record does not determine the unknowns together");
  }

  // Predictions of every observation, and the judgement on those held aside.
  out.predicted = *p0;
  double sse = 0, ssp = 0;
  std::size_t judged = 0;
  for (std::size_t i = 0; i < obs.size(); i++) {
    if (!out.held[i]) continue;
    // persistence: the last observation of the same quantity and unit before the hold
    double last_t = -std::numeric_limits<double>::infinity(), last_v = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t j = 0; j < obs.size(); j++)
      if (!out.held[j] && obs[j].quantity == obs[i].quantity && obs[j].unit == obs[i].unit && obs[j].time > last_t)
        last_t = obs[j].time, last_v = obs[j].value;
    if (!std::isfinite(last_v)) continue;
    const double e = (out.predicted[i] - obs[i].value) / obs[i].sigma;
    const double q = (last_v - obs[i].value) / obs[i].sigma;
    sse += e * e, ssp += q * q, judged++;
  }
  if (judged) {
    out.held_rms = std::sqrt(sse / static_cast<double>(judged));
    out.persistence_rms = std::sqrt(ssp / static_cast<double>(judged));
    out.skill = ssp > 0 ? 1 - sse / ssp : std::numeric_limits<double>::quiet_NaN();
  }
  out.inputs = in;
  for (std::size_t j = 0; j < k; j++)
    out.inputs.parameters[out.names[j]] = Value{out.values[j], Source::estimated,
                                                std::isfinite(out.standard_errors[j]) ? out.standard_errors[j] : 0, missing_value, missing_value};
  return out;
}

}  // namespace kinesis::scales
