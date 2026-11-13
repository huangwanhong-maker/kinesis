// Checks of kinesis/kernel: what every tool answers with.
//
// - Access is ordered: a presentation gives an oracle, an oracle gives records; a tool asked for more than a
//   system gives refuses, naming what it needs and what it was given.
// - A set of identified values is kept distinct and sorted; it is certified with one value, ambiguous with several,
//   inconsistent with none; a set that is known to be a part certifies nothing, and is not inconsistent when empty.
// - The value of a set that is not certified is refused with the reason.
// - A consistent refinement can only contract the set (Huang 2026, Proposition 7.1): `refines`.
// - The dispositions of a certificate against the one before it, in the cases of the draft's Table 2: issued,
//   retained by refinement, withdrawn, issued again, ambiguous, inconsistent. The seventh, transported, is a
//   caller's to give: the kernel compares within one representation.

#include <kinesis/kernel/answer.hpp>

#include <calx/verification/check.hpp>

#include <string>

using namespace kinesis::kernel;

int main() {
  calx::verification::Suite s;

  // ---- access ----
  {
    bool passed = true;
    try {
      require_access(Access::presentation, Access::oracle, "a tool");
      require_access(Access::oracle, Access::record, "a tool");
      require_access(Access::record, Access::record, "a tool");
    } catch (const Refusal&) {
      passed = false;
    }
    s.require("a system known by more than a tool needs is given to it", passed);
    std::string said;
    try {
      require_access(Access::record, Access::presentation, "Lyapunov's exponents by the tangent map");
    } catch (const Refusal& r) {
      said = r.what();
    }
    s.require("a tool asked for more than a system gives refuses, naming what it needs and what it was given",
              said == "Lyapunov's exponents by the tangent map needs a presentation; this system is known only as a record", said);
    s.require("the levels are named", std::string(name(Access::oracle)) == "oracle" && std::string(name(Disposition::retained)) == "retained by refinement");
  }

  // ---- identified values ----
  Identified<int> none, one, two, part;
  {
    one.add(3);
    for (int v : {5, 3, 5, 4, 3}) two.add(v);
    part.add(1), part.complete = false;
    s.require("values are kept distinct and sorted", two.values == std::vector<int>({3, 4, 5}));
    s.require("one value is certified, several are ambiguous, none is inconsistent",
              one.certified() && !one.ambiguous() && !one.inconsistent() && two.ambiguous() && !two.certified() && none.inconsistent() && !none.ambiguous());
    s.require("a set known to be a part certifies nothing", !part.certified() && part.ambiguous() && !part.inconsistent());
    Identified<int> cut;
    cut.complete = false;
    s.require("and an empty part is not inconsistent: nothing was found, which is not that nothing is there", !cut.inconsistent() && cut.ambiguous());
    s.require("the value of a certified set is its value", one.value() == 3);
    two.reason = "three windings are admitted";
    std::string said;
    try {
      two.value();
    } catch (const Refusal& r) {
      said = r.what();
    }
    s.require("the value of a set that is not certified is refused with the reason", said == "the feature is not identified: three windings are admitted", said);
    Identified<int> wider;
    for (int v : {2, 3, 4, 5, 6}) wider.add(v);
    s.require("a refinement contracts: {3} within {3, 4, 5} within {2, ..., 6}, and not the other way", one.refines(two) && two.refines(wider) && !wider.refines(two) && none.refines(one));
  }

  // ---- dispositions ----
  {
    Identified<int> other;
    other.add(4);
    s.require("a first certificate is issued", disposition(none, one) == Disposition::issued);
    s.require("the same value certified again is retained", disposition(one, one) == Disposition::retained);
    s.require("a certificate that becomes several values is withdrawn", disposition(one, two) == Disposition::withdrawn);
    s.require("several values that become one are issued again, and so is another value", disposition(two, one) == Disposition::reissued && disposition(one, other) == Disposition::reissued);
    s.require("several values that stay several are ambiguous", disposition(two, two) == Disposition::ambiguous && disposition(none, two) == Disposition::ambiguous);
    s.require("no history left is inconsistent, whatever stood before", disposition(one, none) == Disposition::inconsistent && disposition(two, none) == Disposition::inconsistent);
    s.require("a part that becomes one value is issued again, not issued: something stood before", disposition(part, one) == Disposition::reissued);
  }

  return s.done();
}
