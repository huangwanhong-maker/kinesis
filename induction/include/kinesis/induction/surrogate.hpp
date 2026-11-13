// Models with no law in them: a forecaster fitted to a record, for a system known only as data or as an oracle.
//
// When no template is proposed, a record can still be continued. The forecaster here is the next-generation
// reservoir computer of Gauthier, Bollt, Griffith and Barbosa (2021, Nat. Commun. 12, 5564): its features are a
// constant, the state at the last k samples, and every product of two of those; the increment of the state over
// one sample is a linear function of the features, fitted by ridge regression,
//
//     W = Y Φᵀ (Φ Φᵀ + αN I)⁻¹,
//
// on features each scaled to a unit root-mean-square over the N samples, so that α is a number without units.
// It is a nonlinear vector autoregression, and has no random part: the same record gives the same model. Run on
// its own output it forecasts; on a chaotic system the forecast leaves the truth at the rate of the largest
// Lyapunov exponent, however good one step is, so that what is reported is the horizon over which it stays
// within a tolerance, in units of the record's time.
// *Caveat.* The ridge is not only against noise. A forecaster with delays is a step of several levels, and has
// modes of its own beside the system's; with too small a ridge one step is fitted better and those modes grow in
// a free run (the check shows both). The ridge is to be chosen by the horizon on a part kept aside, not by the
// error of one step.
// *Caveat.* Its quadratic features are an assumption about the law. A system whose law is not near a quadratic
// polynomial of a few delays is forecast badly, with no warning from the fit (Zhang and Cornelius 2023, Phys.
// Rev. Research 5, 033213).
//
// The model is a black box by construction: it is not a system of the language, and has no terms to analyse.
#pragma once

#include <kinesis/induction/sparse.hpp>

#include <calx/linalg/qr.hpp>
#include <calx/linalg/rectangular.hpp>

#include <cmath>
#include <vector>

namespace kinesis::induction {

class Forecaster {
 public:
  std::size_t dimension = 0, delays = 2, skip = 1;
  double ridge = 1e-6;
  double training_error = 0;   // the root-mean-square of one step on the record, against that of the increments
  std::size_t features = 0;

  // The increment predicted from the last `delays` states, the newest last.
  State increment(const std::vector<State>& recent) const {
    const std::vector<double> phi = features_of(recent);
    State out(dimension, 0.0);
    for (std::size_t i = 0; i < dimension; i++)
      for (std::size_t j = 0; j < phi.size(); j++) out[i] += weights_[i * phi.size() + j] * phi[j];
    return out;
  }
  // The forecast of `steps` samples that follows `history` (its last (delays − 1) · skip + 1 samples are used, and at least as many are asked).
  std::vector<State> forecast(const std::vector<State>& history, std::size_t steps) const {
    if (history.size() < (delays - 1) * skip + 1) throw kernel::Refusal("the forecaster needs the last " + std::to_string((delays - 1) * skip + 1) + " samples");
    std::vector<State> all(history.end() - static_cast<std::ptrdiff_t>((delays - 1) * skip + 1), history.end());
    std::vector<State> out;
    for (std::size_t s = 0; s < steps; s++) {
      std::vector<State> recent;
      for (std::size_t d = 0; d < delays; d++) recent.push_back(all[all.size() - 1 - (delays - 1 - d) * skip]);
      State next = all.back();
      const State dx = increment(recent);
      for (std::size_t i = 0; i < dimension; i++) next[i] += dx[i];
      all.push_back(next), out.push_back(next);
    }
    return out;
  }

 private:
  friend Forecaster forecaster(const Trajectory&, std::size_t, std::size_t, double);
  std::vector<double> features_of(const std::vector<State>& recent) const {
    std::vector<double> linear;
    for (const auto& x : recent) linear.insert(linear.end(), x.begin(), x.end());
    std::vector<double> phi{1.0};
    phi.insert(phi.end(), linear.begin(), linear.end());
    for (std::size_t a = 0; a < linear.size(); a++)
      for (std::size_t b = a; b < linear.size(); b++) phi.push_back(linear[a] * linear[b]);
    return phi;
  }
  std::vector<double> weights_;  // dimension × features
};

inline Forecaster forecaster(const Trajectory& record, std::size_t delays = 2, std::size_t skip = 1, double ridge = 1e-6) {
  Forecaster f;
  f.dimension = record.x.at(0).size(), f.delays = delays, f.skip = skip, f.ridge = ridge;
  const std::size_t first = (delays - 1) * skip, n = record.x.size();
  if (n < first + 10) throw kernel::Refusal("the record is too short for a forecaster of " + std::to_string(delays) + " delays");
  std::vector<std::vector<double>> rows;
  for (std::size_t k = first; k + 1 < n; k++) {
    std::vector<State> recent;
    for (std::size_t d = 0; d < delays; d++) recent.push_back(record.x[k - (delays - 1 - d) * skip]);
    rows.push_back(f.features_of(recent));
  }
  const std::size_t F = rows[0].size(), N = rows.size();
  f.features = F;
  // Ridge regression as least squares on the features, each scaled to a unit root-mean-square, with √(αN) times
  // the identity appended (Φ Φᵀ + αN I in the scaled features): solved by an orthogonal factorisation, since the
  // features at neighbouring delays are nearly dependent and the normal equations would square their condition.
  std::vector<double> scale(F, 0.0);
  for (const auto& r : rows)
    for (std::size_t a = 0; a < F; a++) scale[a] += r[a] * r[a];
  for (double& v : scale) v = v > 0 ? std::sqrt(v / static_cast<double>(N)) : 1.0;
  Rectangular design(N + F, F);
  for (std::size_t k = 0; k < N; k++)
    for (std::size_t a = 0; a < F; a++) design(k, a) = rows[k][a] / scale[a];
  for (std::size_t a = 0; a < F; a++) design(N + a, a) = std::sqrt(ridge * static_cast<double>(N));
  const calx::linalg::Qr factor = calx::linalg::qr(design);
  f.weights_.assign(f.dimension * F, 0.0);
  std::vector<double> y(N + F, 0.0);
  for (std::size_t i = 0; i < f.dimension; i++) {
    for (std::size_t k = 0; k < N; k++) y[k] = record.x[first + k + 1][i] - record.x[first + k][i];
    const std::vector<double> w = calx::linalg::least_squares(factor, y);
    for (std::size_t a = 0; a < F; a++) f.weights_[i * F + a] = w[a] / scale[a];
  }
  double error = 0, size = 0;
  for (std::size_t k = 0; k < N; k++)
    for (std::size_t i = 0; i < f.dimension; i++) {
      double model = 0;
      for (std::size_t a = 0; a < F; a++) model += f.weights_[i * F + a] * rows[k][a];
      const double y = record.x[first + k + 1][i] - record.x[first + k][i];
      error += (model - y) * (model - y), size += y * y;
    }
  f.training_error = std::sqrt(error / std::max(size, 1e-300));
  return f;
}

// The number of samples over which a forecast stays within `tolerance` of the truth, the distance measured
// against the root-mean-square size of the truth's departures from its mean.
inline std::size_t horizon(const std::vector<State>& forecast, const std::vector<State>& truth, double tolerance = 0.3) {
  const std::size_t n = std::min(forecast.size(), truth.size()), dim = truth.at(0).size();
  std::vector<double> mean(dim, 0.0);
  for (const auto& x : truth)
    for (std::size_t i = 0; i < dim; i++) mean[i] += x[i] / static_cast<double>(truth.size());
  double spread = 0;
  for (const auto& x : truth)
    for (std::size_t i = 0; i < dim; i++) spread += (x[i] - mean[i]) * (x[i] - mean[i]) / static_cast<double>(truth.size());
  spread = std::sqrt(spread);
  for (std::size_t k = 0; k < n; k++) {
    double gap = 0;
    for (std::size_t i = 0; i < dim; i++) gap += (forecast[k][i] - truth[k][i]) * (forecast[k][i] - truth[k][i]);
    if (std::sqrt(gap) > tolerance * spread) return k;
  }
  return n;
}

}  // namespace kinesis::induction
