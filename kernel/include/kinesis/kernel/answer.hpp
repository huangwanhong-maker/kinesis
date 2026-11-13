// What every tool of kinesis answers with: the levels at which a system is known, a refusal with its reason, and
// the values a record identifies.
//
// Access. A system is known at a level (the theory note, §4): by a presentation, a law one can read; by an oracle,
// something one can drive and read; or by a record, finite data seen through an observation. A tool states the
// access it needs. One that is asked for more than a system gives refuses with the reason, and never answers with
// a made-up value.
//
// Identified values. An answer drawn from a record is not one value but the set of values over every history the
// record admits under stated assumptions,
//
//     I_Φ(D; A) = { Φ(Γ) : Γ ∈ R(D; A) }
//
// (Huang 2026, "Rethinking trajectory-based credibility under symbol inflation", Eqs. 19–20). It is certified when
// the set is a single value, ambiguous when it has several, and inconsistent when no history satisfies the
// assumptions. A certificate keeps the assumptions it was issued under and how it stands to the one before it:
// retained by refinement, transported across a change of representation, withdrawn, or issued again (its §7).
#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace kinesis::kernel {

// The level at which a system is known, ordered: a presentation gives an oracle, an oracle gives records.
enum class Access { record = 0, oracle = 1, presentation = 2 };

inline const char* name(Access a) {
  switch (a) {
    case Access::record: return "record";
    case Access::oracle: return "oracle";
    default: return "presentation";
  }
}

// A request that cannot be answered, with the reason: thrown by tools, reported by their callers.
struct Refusal : std::runtime_error {
  explicit Refusal(const std::string& reason) : std::runtime_error(reason) {}
};

inline void require_access(Access have, Access need, const std::string& tool) {
  if (static_cast<int>(have) < static_cast<int>(need))
    throw Refusal(tool + " needs a " + name(need) + "; this system is known only as a " + name(have));
}

// How a certificate stands to the one before it.
enum class Disposition { issued, retained, transported, withdrawn, reissued, ambiguous, inconsistent };

inline const char* name(Disposition d) {
  switch (d) {
    case Disposition::issued: return "issued";
    case Disposition::retained: return "retained by refinement";
    case Disposition::transported: return "transported";
    case Disposition::withdrawn: return "withdrawn";
    case Disposition::reissued: return "issued again";
    case Disposition::ambiguous: return "ambiguous";
    default: return "inconsistent";
  }
}

// The values of a feature over every admissible history. `complete` is false when the list is known to be a part
// of the set (an enumeration that was cut short): a set that is not complete certifies nothing.
template <class T>
struct Identified {
  std::vector<T> values;   // distinct, sorted
  bool complete = true;
  std::vector<std::string> assumptions;
  std::string reason;      // why it is certified, or what leaves it open

  bool certified() const { return complete && values.size() == 1; }
  bool inconsistent() const { return complete && values.empty(); }
  bool ambiguous() const { return !certified() && !inconsistent(); }
  const T& value() const {
    if (!certified()) throw Refusal("the feature is not identified: " + reason);
    return values.front();
  }
  void add(const T& v) {
    auto it = std::lower_bound(values.begin(), values.end(), v);
    if (it == values.end() || *it != v) values.insert(it, v);
  }
  // Every value of this set is a value of `earlier` (Proposition 7.1: a consistent refinement can only contract).
  bool refines(const Identified& earlier) const {
    return std::includes(earlier.values.begin(), earlier.values.end(), values.begin(), values.end());
  }
};

// How a new set of identified values stands to the one issued before, under the same feature and representation.
template <class T>
Disposition disposition(const Identified<T>& before, const Identified<T>& now) {
  if (now.inconsistent()) return Disposition::inconsistent;
  if (now.certified()) {
    if (before.certified()) return before.values == now.values ? Disposition::retained : Disposition::reissued;
    return before.values.empty() && before.complete ? Disposition::issued : Disposition::reissued;
  }
  return before.certified() ? Disposition::withdrawn : Disposition::ambiguous;
}

}  // namespace kinesis::kernel
