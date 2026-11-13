// Checks of kinesis/terms: the normal form of terms, their derivatives, their values and how they are written.
//
// - Normal form: equal terms are collected (x + x is 2x, x·x is x², (x + y) − (y + x) is 0, x / x is 1), a number is
//   distributed over a sum, and the order of the operands does not matter.
// - Derivatives against central differences at random points, for products, quotients, powers with a variable
//   exponent, every named function of one argument, atan2, min and max; and identities in closed form
//   (d/dx tanh x = 1 − tanh² x as terms; the second derivative of sin is −sin).
// - The total derivative in time: D_t (q'² / 2 − cos q) = q' q'' + q' sin q.
// - Programs: the value of a compiled term is the value of the term, to rounding, at random points; a program of a
//   term with no symbol reads no slot; integer powers are exact products.
// - Refusals: a symbol that is not declared, an application that is not a function, an unknown.
// - Writing: text and LaTeX of terms chosen to exercise signs, quotients, powers and names.

#include <kinesis/terms/program.hpp>
#include <kinesis/terms/term.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace kinesis::terms;

int main() {
  calx::verification::Suite s;
  const TermPtr x = symbol("x"), y = symbol("y"), z = symbol("z"), t = symbol("t");
  const auto n = [](double v) { return number(v); };

  // ---- normal form ----
  s.require("x + x is 2 x", same(x + x, n(2) * x), text(x + x));
  s.require("x x is x^2", same(x * x, power(x, n(2))), text(x * x));
  s.require("(x + y) - (y + x) is 0", is_number((x + y) - (y + x), 0), text((x + y) - (y + x)));
  s.require("x / x is 1", is_number(x / x, 1), text(x / x));
  s.require("x y is y x", same(x * y, y * x));
  s.require("x + y + z in any order", same((x + y) + z, z + (y + x)));
  s.require("2 (x + y) - 2 x - 2 y is 0", is_number(n(2) * (x + y) - n(2) * x - n(2) * y, 0), text(n(2) * (x + y) - n(2) * x - n(2) * y));
  s.require("x^0 is 1 and x^1 is x", is_number(power(x, n(0)), 1) && same(power(x, n(1)), x));
  s.require("(x^2)^3 is x^6", same(power(power(x, n(2)), n(3)), power(x, n(6))));
  s.require("(x y)^2 is x^2 y^2", same(power(x * y, n(2)), power(x, n(2)) * power(y, n(2))));
  s.require("0 x is 0 and numbers fold", is_number(n(0) * x, 0) && is_number(n(2) * n(3) + n(1), 7));
  s.require("x - x + 3 is 3", is_number(x - x + n(3), 3));
  s.require("x^2 x^-2 is 1", is_number(power(x, n(2)) * power(x, n(-2)), 1));

  // ---- derivatives against differences ----
  std::mt19937_64 rng(20261112);
  std::uniform_real_distribution<double> u(0.3, 1.7);
  const auto at = [&](const TermPtr& f, double a, double b) {
    return evaluate(f, [&](const std::string& name, unsigned) { return name == "x" ? a : b; });
  };
  const auto against_differences = [&](const std::string& what, const TermPtr& f) {
    const TermPtr dx = derivative(f, "x"), dy = derivative(f, "y");
    double worst = 0;
    for (int k = 0; k < 40; k++) {
      const double a = u(rng), b = u(rng), h = 1e-5;
      const double nx = (at(f, a + h, b) - at(f, a - h, b)) / (2 * h), ny = (at(f, a, b + h) - at(f, a, b - h)) / (2 * h);
      worst = std::max({worst, std::fabs(at(dx, a, b) - nx) / (1 + std::fabs(nx)), std::fabs(at(dy, a, b) - ny) / (1 + std::fabs(ny))});
    }
    s.check("the derivative of " + what + " against central differences", worst, 2e-9);
  };
  against_differences("x^3 sin(x y)", power(x, n(3)) * call("sin", {x * y}));
  against_differences("(x + 2 y) / (1 + x^2)", (x + n(2) * y) / (n(1) + power(x, n(2))));
  against_differences("x^y", power(x, y));
  against_differences("exp(-x^2 / y)", call("exp", {-power(x, n(2)) / y}));
  for (const char* f : {"sin", "cos", "exp", "log", "sqrt", "sinh", "cosh", "tanh", "atan"})
    against_differences(std::string(f) + "(x y)", call(f, {x * y}));
  against_differences("tan(x y / 4)", call("tan", {x * y / n(4)}));  // away from its poles
  against_differences("asin(x y / 4)", call("asin", {x * y / n(4)}));
  against_differences("acos(x y / 4)", call("acos", {x * y / n(4)}));
  against_differences("atan2(y, x)", call("atan2", {y, x}));
  against_differences("abs(x - y) x", call("abs", {x - y}) * x);
  against_differences("min(x^2, y) + max(x, y^2)", call("min", {power(x, n(2)), y}) + call("max", {x, power(y, n(2))}));

  // ---- derivatives in closed form ----
  s.require("d/dx tanh x is 1 - tanh^2 x", same(derivative(call("tanh", {x}), "x"), n(1) - power(call("tanh", {x}), n(2))), text(derivative(call("tanh", {x}), "x")));
  s.require("the second derivative of sin x is -sin x", same(derivative(derivative(call("sin", {x}), "x"), "x"), -call("sin", {x})));
  s.require("d/dx (x^2 y) is 2 x y", same(derivative(power(x, n(2)) * y, "x"), n(2) * x * y));
  s.require("a derivative in a symbol that does not occur is 0", is_number(derivative(call("exp", {y}), "x"), 0));
  {
    const TermPtr q = symbol("q"), qd = symbol("q", 1), qdd = symbol("q", 2);
    const TermPtr L = power(qd, n(2)) / n(2) - call("cos", {q});
    const TermPtr D = total_derivative(L, "t", {"q"});
    s.require("D_t (q'^2 / 2 - cos q) is q' q'' + q' sin q", same(D, qd * qdd + qd * call("sin", {q})), text(D));
    s.require("d/dq' sees q' and not q", same(derivative(L, "q", 1), qd) && same(derivative(L, "q", 0), call("sin", {q})));
  }

  // ---- programs ----
  {
    Frame frame;
    frame.add("x"), frame.add("y"), frame.add("z");
    const std::vector<TermPtr> terms = {
        power(x, n(3)) * call("sin", {x * y}) - z / (n(1) + power(y, n(2))),
        power(x, y) + call("atan2", {y, x}) * call("exp", {-z}),
        power(x + y, n(-3)) * power(z, n(7)),
        call("sqrt", {power(x, n(2)) + power(y, n(2)) + power(z, n(2))}),
        derivative(power(x, n(3)) * call("tanh", {x * y * z}), "x"),
    };
    double worst = 0;
    for (const auto& f : terms) {
      const Program p = compile(f, frame.resolver());
      for (int k = 0; k < 50; k++) {
        const double v[3] = {u(rng), u(rng), u(rng)};
        const double a = p(v), b = evaluate(f, [&](const std::string& name, unsigned) { return v[frame.find(name)]; });
        worst = std::max(worst, std::fabs(a - b) / (1 + std::fabs(b)));
      }
    }
    s.check("a program's value against its term's, 5 terms at 50 points", worst, 1e-14);
    const Program c = compile(n(2) * n(3) + call("cos", {n(0)}), frame.resolver());
    s.require("a program of a term with no symbol reads no slot, and is 7", c.constant() && c(std::vector<double>{}) == 7);
    const Program p7 = compile(power(x, n(7)), frame.resolver());
    const double v[3] = {1.5, 0, 0};
    s.check("x^7 by multiplication at 1.5", p7(v) - 17.0859375, 0);
    s.require("x^-2 at 1.5 is 4/9", std::fabs(compile(power(x, n(-2)), frame.resolver())(v) - 4.0 / 9) < 1e-16);
  }

  // ---- refusals ----
  const auto refused = [&](const std::string& what, const TermPtr& f, const std::string& names) {
    Frame frame;
    frame.add("x");
    try {
      compile(f, frame.resolver());
      s.require(what, false, "compiled");
    } catch (const kinesis::kernel::Refusal& e) {
      s.require(what, std::string(e.what()).find(names) != std::string::npos, e.what());
    }
  };
  refused("a symbol that is not declared is refused by name", x + symbol("w"), "'w'");
  refused("a primed symbol that is not declared is refused", symbol("x", 1), "'x''");
  refused("an application that is not a function is refused", call("beta", {x}), "'beta'");
  refused("an unknown has no value", hole("f", {x}), "?f");
  try {
    derivative(hole("f", {x}), "x");
    s.require("the derivative of an unknown is refused", false);
  } catch (const kinesis::kernel::Refusal&) {
    s.require("the derivative of an unknown is refused", true);
  }

  // ---- writing ----
  const auto written = [&](const TermPtr& f, const std::string& want) { s.require("text: " + want, text(f) == want, text(f)); };
  written(x - y, "x - y");
  written(-x, "-x");
  written(n(2) * x * y, "2 * x * y");
  written(x / y, "x / y");
  written(x / (y * z), "x / y / z");
  written((x + y) * z, "z * (x + y)");
  written(power(x + y, n(2)), "(x + y)^2");
  written(power(x, n(-1)), "1 / x");
  written(-power(x, n(2)), "-x^2");
  written(n(1) - power(call("tanh", {x}), n(2)), "1 - tanh(x)^2");
  written(symbol("q", 2) + call("sin", {symbol("q")}), "q'' + sin(q)");
  written(hole("V", {symbol("q")}), "?V(q)");
  written(n(0.1) * symbol("a.v"), "0.1 * a.v");
  written(n(2.5e-7) * x, "2.5e-07 * x");
  const auto typeset = [&](const TermPtr& f, const std::string& want) { s.require("LaTeX: " + want, latex(f) == want, latex(f)); };
  typeset(symbol("sigma") * (y - x), "\\sigma\\,\\left(-x + y\\right)");
  typeset(x / (n(1) + power(x, n(2))), "\\frac{x}{1 + x^{2}}");
  typeset(symbol("q", 2), "\\ddot{q}");
  typeset(symbol("x1") + symbol("k_on"), "k_{\\mathrm{on}} + x_{1}");
  typeset(call("sqrt", {x}) * call("exp", {-x}), "e^{-x}\\,\\sqrt{x}");
  // Numbers in mathematics: few digits as they are, a fraction of small parts as that fraction, others to six digits;
  // the text keeps every digit.
  typeset(n(0.45) * x, "0.45\\,x");
  typeset(n(8.0 / 3) * x, "\\frac{8}{3}\\,x");
  typeset(n(1.0 / 24) + x, "\\frac{1}{24} + x");
  typeset(n(0.123456789123) * x, "0.123457\\,x");
  written(n(8.0 / 3) * x, "2.6666666666666665 * x");

  return s.done();
}
