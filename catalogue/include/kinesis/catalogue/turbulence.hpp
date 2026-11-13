// A shell model of turbulence, as a declaration: the Sabra model, written out shell by shell.
//
// Turbulence in a fluid is a cascade: energy put in at large scales is passed, by the nonlinear term of the
// equations of motion, through smaller and smaller scales to those at which viscosity removes it. A shell model
// keeps that and little else. Wavenumbers are cut into shells k_n = k_0 λⁿ⁻¹, each shell holds one complex
// amplitude u_n, and each amplitude is coupled to its two nearest neighbours on either side by a term quadratic
// in them, as the nonlinear term of the Navier–Stokes equations is, with viscosity and a force:
//
//     du_n/dt = i ( a k_{n+1} u_{n+2} u*_{n+1} + b k_n u_{n+1} u*_{n−1} − c k_{n−1} u_{n−1} u_{n−2} ) − ν k_n² u_n + f_n
//
// (L'vov, Podivilov, Pomyalov, Procaccia and Vandembroucq 1998, "Improved shell model of turbulence", Phys. Rev. E
// 58, 1811–1822; the review of Biferale 2003, Annu. Rev. Fluid Mech. 35, 441–468). It is not a discretisation of
// the Navier–Stokes equations: it has no space, and one amplitude stands for all the motion at a scale.
//
// What is exact, and is checked. With a + b + c = 0, and no viscosity and no force, the nonlinear term keeps the
// energy E = Σ |u_n|² and a second quadratic quantity Σ (a/c)ⁿ |u_n|² (with a = 1, b = c = −½ it alternates in
// sign, as helicity does). With viscosity and a force, in a statistically steady state the mean power put in
// equals the mean power dissipated, 2 Re⟨f* u_1⟩ = 2ν Σ k_n² ⟨|u_n|²⟩ = ε, and the mean flux of energy through any
// shell between the forced one and the dissipative ones equals ε: the model's form of Kolmogorov's four-fifths
// law. What is measured, not exact: the scaling ⟨|u_n|²⟩ ∝ k_n^{−ζ₂}, for which Kolmogorov's argument of 1941
// gives ζ₂ = 2/3, and the growth of the flatness ⟨|u_n|⁴⟩/⟨|u_n|²⟩² with n, which is intermittency.
//
// The declaration. The language holds real components, so u_n = x_n + i y_n is two of them. shell_model() writes
// the 2N equations, and what is observed: the energy e_n of each shell, the whole energy and the second quantity,
// the power put in and dissipated, and the flux of energy through two shells. A path started from zero stays on
// the forced shell alone, which is an invariant subspace: start_of() gives a start with every shell moved.
#pragma once

#include <kinesis/terms/term.hpp>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace kinesis::catalogue {

struct ShellModel {
  std::string name = "sabra";
  std::size_t shells = 14;
  double nu = 3e-5;                    // the viscosity: with the force 1 the cascade ends near the twelfth shell
  double forcing = 1;                  // f_1 = forcing (1 + i); the other shells are not forced
  double k0 = 1, ratio = 2;            // k_n = k0 ratio^(n − 1)
  double a = 1, b = -0.5, c = -0.5;    // a + b + c = 0 keeps the energy
  std::vector<std::size_t> flux = {3, 6};  // the shells through which the flux of energy is observed
};

inline double wavenumber(const ShellModel& m, std::size_t n) { return m.k0 * std::pow(m.ratio, static_cast<double>(n) - 1); }

// The declaration of the model, a script of one system.
inline std::string shell_model(const ShellModel& m = {}) {
  const std::size_t N = m.shells;
  if (N < 4) throw kernel::Refusal("a shell model has four shells or more");
  const auto num = [](double v) { return terms::number_text(v); };
  const auto x = [](std::size_t n) { return "x" + std::to_string(n); };
  const auto y = [](std::size_t n) { return "y" + std::to_string(n); };
  const auto k = [&](std::size_t n) { return wavenumber(m, n); };
  // The nonlinear term of shell n, its real and imaginary parts as sums of signed products. With Z the bracket of
  // the equation, the term is iZ: its real part is −Im Z and its imaginary part Re Z.
  std::vector<std::string> re(N + 1), im(N + 1);
  const auto add = [](std::string& to, double coefficient, const std::string& product) {
    if (coefficient == 0) return;
    to += (coefficient < 0 ? " - " : (to.empty() ? "" : " + ")) + terms::number_text(std::fabs(coefficient)) + " * " + product;
  };
  for (std::size_t n = 1; n <= N; n++) {
    if (n + 2 <= N) {  // a k_{n+1} u_{n+2} u*_{n+1}
      add(re[n], -m.a * k(n + 1), "(" + y(n + 2) + " * " + x(n + 1) + " - " + x(n + 2) + " * " + y(n + 1) + ")");
      add(im[n], m.a * k(n + 1), "(" + x(n + 2) + " * " + x(n + 1) + " + " + y(n + 2) + " * " + y(n + 1) + ")");
    }
    if (n + 1 <= N && n >= 2) {  // b k_n u_{n+1} u*_{n−1}
      add(re[n], -m.b * k(n), "(" + y(n + 1) + " * " + x(n - 1) + " - " + x(n + 1) + " * " + y(n - 1) + ")");
      add(im[n], m.b * k(n), "(" + x(n + 1) + " * " + x(n - 1) + " + " + y(n + 1) + " * " + y(n - 1) + ")");
    }
    if (n >= 3) {  // − c k_{n−1} u_{n−1} u_{n−2}
      add(re[n], m.c * k(n - 1), "(" + x(n - 1) + " * " + y(n - 2) + " + " + y(n - 1) + " * " + x(n - 2) + ")");
      add(im[n], -m.c * k(n - 1), "(" + x(n - 1) + " * " + x(n - 2) + " - " + y(n - 1) + " * " + y(n - 2) + ")");
    }
    if (re[n].empty()) re[n] = "0";
    if (im[n].empty()) im[n] = "0";
    // A sum that begins with a minus is written so: " - 2 * (…)" becomes "-2 * (…)".
    for (std::string* s : {&re[n], &im[n]})
      if (s->rfind(" - ", 0) == 0) *s = "-" + s->substr(3);
  }
  const auto energy = [&](std::size_t n) { return "(" + x(n) + "^2 + " + y(n) + "^2)"; };
  std::string out = "system " + m.name + " {\n  state ";
  for (std::size_t n = 1; n <= N; n++) out += (n > 1 ? ", " : "") + x(n) + ", " + y(n);
  out += " : real\n  param nu = " + num(m.nu) + ", f = " + num(m.forcing) + "\n  change {\n";
  for (std::size_t n = 1; n <= N; n++) {
    const std::string damping = " - nu * " + num(k(n) * k(n)) + " * ";
    out += "    " + x(n) + "' = " + re[n] + damping + x(n) + (n == 1 ? " + f" : "") + "\n";
    out += "    " + y(n) + "' = " + im[n] + damping + y(n) + (n == 1 ? " + f" : "") + "\n";
  }
  out += "  }\n  difference {\n    observe ";
  for (std::size_t n = 1; n <= N; n++) out += "e" + std::to_string(n) + " = " + energy(n).substr(1, energy(n).size() - 2) + ", ";
  out += "energy = ";
  for (std::size_t n = 1; n <= N; n++) out += (n > 1 ? " + " : "") + x(n) + "^2 + " + y(n) + "^2";
  out += ", second = ";
  for (std::size_t n = 1; n <= N; n++) {
    const double weight = std::pow(m.a / m.c, static_cast<double>(n));
    out += (weight < 0 ? (n > 1 ? " - " : "-") : (n > 1 ? " + " : "")) + num(std::fabs(weight)) + " * " + energy(n);
  }
  out += ", input = 2 * f * (x1 + y1), dissipation = 2 * nu * (";
  for (std::size_t n = 1; n <= N; n++) out += (n > 1 ? " + " : "") + num(k(n) * k(n)) + " * " + energy(n);
  out += ")";
  for (std::size_t through : m.flux) {
    if (through < 1 || through >= N) throw kernel::Refusal("the flux of energy is observed through a shell that has shells beyond it");
    out += ", flux" + std::to_string(through) + " = -2 * (";
    for (std::size_t n = 1; n <= through; n++) out += (n > 1 ? " + " : "") + x(n) + " * (" + re[n] + ") + " + y(n) + " * (" + im[n] + ")";
    out += ")";
  }
  out += "\n  }\n}\n";
  return out;
}

// A start with every shell moved: amplitudes that fall as k^(−1/3), with phases that differ from shell to shell.
inline std::vector<std::pair<std::string, double>> start_of(const ShellModel& m = {}) {
  std::vector<std::pair<std::string, double>> out;
  for (std::size_t n = 1; n <= m.shells; n++) {
    const double size = 0.5 * std::pow(wavenumber(m, n) / m.k0, -1.0 / 3), phase = 1.7 * static_cast<double>(n * n);
    out.push_back({"x" + std::to_string(n), size * std::cos(phase)});
    out.push_back({"y" + std::to_string(n), size * std::sin(phase)});
  }
  return out;
}

}  // namespace kinesis::catalogue
