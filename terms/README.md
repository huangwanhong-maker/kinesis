# Terms

Status: **draft** (2026-11-13). Module `kinesis::terms`.

## 1. Concepts

A term is an expression of the language of dynamical systems, held as a tree: a number, a symbol (a name with the order of its derivative in time: q, q', q''), a sum, a product, a power, an application f(a, …) or an unknown ?f(a, …). There is no difference and no quotient: a − b is a + (−1)·b and a / b is a·b⁻¹, so that one form stands for each. Terms are built in a normal form by their constructors: sums and products are flat and sorted by a fixed total order, equal terms and equal bases are collected, numbers are folded, and a number times a sum is distributed inside a sum.

Derivatives are symbolic and exact, and the total derivative in time is $D_t f = \partial f/\partial t + \sum \partial f/\partial s^{(k)} \cdot s^{(k+1)}$ over the dynamical symbols s, which is how the Euler–Lagrange equations of an action are derived as terms. Written out, a term is a word of the language again. For evaluation many times a term is compiled once into a straight-line program for a stack machine, each symbol resolved to a slot of an array of values, so that evaluation touches no name and allocates nothing. Theory: [docs/theory.md](../docs/theory.md).

## 2. Content

- `term.hpp`:
  - `Kind{number, symbol, sum, product, power, call, hole}`, `Term`, `TermPtr`; the order `compare`, `same`, `TermLess`.
  - Constructors: `number(v)`, `symbol(name, primes)`, `call(name, args)`, `hole(name, args)`, `sum(terms)`, `product(factors)`, `power(base, exponent)`, and the operators `+`, `-`, `*`, `/`; `is_number`, `is_integer`, `split_coefficient` (a term as coefficient × core).
  - The named functions, `functions()` and `function_index(name, arity)`: sin, cos, tan, exp, log, sqrt, sinh, cosh, tanh, asin, acos, atan, abs, sign, step, floor, and of two arguments atan2, min, max, mod.
  - Reading: `symbols(t)`, `depends_on(t, name, primes)`, `holes(t, out)`, `applications(t, out)` (applications of names that are not functions).
  - `rewrite(t, replace)` (the result rebuilt through the constructors, so in normal form), `substitute(t, name, primes, by)`.
  - `derivative(t, name, primes)`, `total_derivative(t, time, dynamical)`; `evaluate(t, lookup)`, a slow reading for small uses.
  - Writing: `number_text(v)`, `text(t)` (as the language writes it), `latex(t)` (for a page).
- `program.hpp`: `Op`, `Instruction`, `Resolver` (the slot of a symbol, or −1); `Program` with `run(values, stack)`, `operator()(values)`, `depth()`, `size()`, `constant()` (it reads no slot); `compile(term, slot)`; `Frame`, slots by name, with `add`, `find`, `size`, `names`, `resolver`.

*Caveat.* The normal form is what the constructors make, not a decision procedure for equality of expressions: two terms that differ after it may still be equal as functions (sin² + cos² is not folded to 1).

*Caveat.* The derivative of an unknown, and of an application whose derivative is not known, is refused where the term depends on the variable. A symbol that the resolver does not know, an application of a name that is not a function, and an unknown are refused at compilation with the name: a program never evaluates what has no meaning.

## 3. Checks

71: the normal form (x + x is 2x, x·x is x², (x + y) − (y + x) is 0, x / x is 1, a number distributed over a sum, the order of the operands not mattering); derivatives against central differences at random points, for products, quotients, powers with a variable exponent, every named function of one argument, atan2, min and max, and identities in closed form (d/dx tanh x = 1 − tanh² x as terms; the second derivative of sin is −sin); the total derivative D_t(q'²/2 − cos q) = q' q'' + q' sin q; programs, the value of a compiled term being the value of the term to rounding at random points, a program of a term with no symbol reading no slot, integer powers as exact products; the refusals of a symbol that is not declared, an application that is not a function and an unknown; and the text and LaTeX of terms chosen to exercise signs, quotients, powers and names. Numbers in LaTeX are set as mathematics writes them (a number of few digits as it is, 8/3 and 1/24 as fractions, others to six digits), while the text keeps every digit.

## 4. Dependencies

`kinesis::kernel`, for `Refusal`. The checks use calx's `verification`.

## 5. References

The headers cite no source.
