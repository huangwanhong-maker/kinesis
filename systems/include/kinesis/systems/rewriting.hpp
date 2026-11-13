// Rewriting systems: rules on words, as a dynamical system whose histories are derivations.
//
// The state is a word over the alphabet; a change is the application of one rule at one position (the changes are
// typed and partial: a rule applies where its left side occurs); the structure is the free composition of such
// applications. A generative grammar is this system read from its start word (the user's remark that each
// generative grammar is a dynamical system; Mjolsness and Yosiphon 2006, "Stochastic process semantics for
// dynamical grammars", Ann. Math. Artif. Intell. 47, 329–395, give rules a process semantics in the same spirit).
//
// Weights (the idea of 2026-11-12, §6). With Boolean weights a history is possible or not, and what the system
// says is its reachable set: the language is the reachable words that no rule applies to and that hold only
// letters no rule rewrites. With probabilities a rule's weight, normalised over the rules with its left side, is
// the chance that it is the one applied; the weight of a history is the product over its steps.
//
// Histories. explore() builds the graph of forms reachable from the start, within a bound on the length of a
// form. For a context-free system the steps are leftmost (a step rewrites the first position at which any rule
// applies): each derivation tree has exactly one leftmost derivation, so the number of paths from the start to a
// word is the number of its derivations, and their summed weight its probability. With rules that have context
// a leftmost order would lose words (a later step may be needed first), so every step is followed, and a history
// is a sequence of steps. A word reached by several paths is one
// whose history the word does not identify: its parse is ambiguous (Difference, §3.1 of the idea: distinct
// histories that the observation stops telling apart).
//
// Bounds. The graph is finite when the rules do not shorten a form (monotone rules) and forms are bounded in
// length; a system with a shortening rule is explored within a bound on the number of forms, and says when the
// bound cut it short. Nothing is reported as the whole language unless the exploration was complete.
//
// The grammar is morphe's (morphe/formal/grammar.hpp), and membership of a word in a context-free language is
// decided by its Earley parser.
#pragma once

#include <kinesis/systems/model.hpp>

#include <calx/random/stream.hpp>
#include <morphe/formal/earley.hpp>
#include <morphe/formal/grammar.hpp>

#include <deque>
#include <map>

namespace kinesis::systems {

using Word = std::vector<std::string>;

class Rewriting {
 public:
  std::string name;
  std::vector<std::string> alphabet;
  Word start;
  struct Rule {
    std::string name;
    Word left, right;
    double weight = 1;       // as declared
    double probability = 1;  // normalised over the rules with the same left side
  };
  std::vector<Rule> rules;
  std::string weights = "boolean";
  morphe::formal::Grammar grammar;   // the same rules for morphe: terminals are the letters no rule rewrites
  bool context_free = true, monotone = true;

  bool terminal(const std::string& letter) const { return rewritten_.count(letter) == 0; }
  bool terminal(const Word& w) const {
    for (const auto& l : w)
      if (!terminal(l)) return false;
    return true;
  }

  struct Step {
    std::size_t rule, position;
    Word to;
  };
  // Every application of every rule to w; with `leftmost`, only those at the first position where a rule applies.
  std::vector<Step> steps(const Word& w, bool leftmost = false) const {
    std::vector<Step> out;
    for (std::size_t i = 0; i < w.size(); i++) {
      for (std::size_t r = 0; r < rules.size(); r++) {
        const Word& l = rules[r].left;
        if (i + l.size() > w.size() || !std::equal(l.begin(), l.end(), w.begin() + static_cast<std::ptrdiff_t>(i))) continue;
        Word v(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(i));
        v.insert(v.end(), rules[r].right.begin(), rules[r].right.end());
        v.insert(v.end(), w.begin() + static_cast<std::ptrdiff_t>(i + l.size()), w.end());
        out.push_back({r, i, std::move(v)});
      }
      if (leftmost && !out.empty()) break;
    }
    return out;
  }

  // Whether the word is in the language. Decided by Earley's parser for a context-free system.
  bool accepts(const Word& w) const {
    if (!context_free) throw kernel::Refusal("membership is decided here for context-free rules; " + name + " has a rule with a longer left side");
    std::vector<morphe::formal::Symbol> symbols;
    for (const auto& l : w) {
      const auto it = symbol_.find(l);
      if (it == symbol_.end() || !grammar.terminal(it->second)) return false;
      symbols.push_back(it->second);
    }
    return morphe::formal::Earley(grammar).accepts(symbols);
  }

 private:
  friend Rewriting rewriting(const System&);
  std::map<std::string, char> rewritten_;                    // the letters some rule rewrites
  std::map<std::string, morphe::formal::Symbol> symbol_;
};

inline Rewriting rewriting(const System& s) {
  if (form_of(s) != Form::rules) throw kernel::Refusal(s.name + " is a " + std::string(name(form_of(s))) + " system, not one of rules");
  const Scope sc = scope(s);
  Rewriting r;
  r.name = s.name;
  r.weights = s.structure.weights.empty() ? "boolean" : s.structure.weights;
  for (const auto& st : s.states)
    if (st.space.words) r.alphabet = st.space.alphabet, r.start = st.start;
  if (r.alphabet.empty()) throw ScriptError(s.at, s.name + " has rules and no state that is a word: declare `state w : words over { … }, start …`");
  if (r.start.empty()) throw ScriptError(s.at, s.name + " has no start: declare `start …` with its state");
  const auto known = [&](const Word& w, language::Place at) {
    for (const auto& l : w)
      if (std::find(r.alphabet.begin(), r.alphabet.end(), l) == r.alphabet.end()) throw ScriptError(at, "'" + l + "' is not a letter of the alphabet of " + s.name);
  };
  known(r.start, s.at);
  std::map<Word, double> total;
  for (const auto& rule : s.change.rules) {
    known(rule.left, rule.at), known(rule.right, rule.at);
    const double w = rule.weight ? sc.constant(rule.weight, rule.at, "the weight of the rule " + rule.name) : 1.0;
    if (w < 0) throw ScriptError(rule.at, "the weight of the rule " + rule.name + " is negative");
    r.rules.push_back({rule.name, rule.left, rule.right, w, w});
    total[rule.left] += w;
    // The letters a rule rewrites are those of its left side that are not context: what is left of the left side
    // after the prefix and the suffix it shares with the right side (in `a B -> a b` the letter a is context).
    std::size_t prefix = 0, suffix = 0;
    const std::size_t shorter = std::min(rule.left.size(), rule.right.size());
    while (prefix < shorter && rule.left[prefix] == rule.right[prefix]) prefix++;
    while (suffix < shorter - prefix && rule.left[rule.left.size() - 1 - suffix] == rule.right[rule.right.size() - 1 - suffix]) suffix++;
    for (std::size_t i = prefix; i + suffix < rule.left.size(); i++) r.rewritten_[rule.left[i]] = 1;
    if (rule.left.size() != 1) r.context_free = false;
    if (rule.right.size() < rule.left.size()) r.monotone = false;
  }
  for (auto& rule : r.rules) rule.probability = total[rule.left] > 0 ? rule.weight / total[rule.left] : 0;
  // morphe's grammar: terminals first, then the letters that are rewritten, the start letter first of those.
  std::vector<std::string> terminals, nonterminals;
  for (const auto& l : r.alphabet) (r.terminal(l) ? terminals : nonterminals).push_back(l);
  const bool fresh = r.start.size() != 1 || r.terminal(r.start[0]);
  if (!fresh) {
    const auto it = std::find(nonterminals.begin(), nonterminals.end(), r.start[0]);
    std::rotate(nonterminals.begin(), it, it + 1);
  }
  r.grammar.terminals = terminals.size();
  for (std::size_t i = 0; i < terminals.size(); i++) r.symbol_[terminals[i]] = static_cast<morphe::formal::Symbol>(i);
  morphe::formal::Symbol start = 0;
  if (fresh) start = r.grammar.add_nonterminal("#start");
  for (const auto& l : nonterminals) r.symbol_[l] = r.grammar.add_nonterminal(l);
  const auto form = [&](const Word& w) {
    morphe::formal::Form f;
    for (const auto& l : w) f.push_back(r.symbol_.at(l));
    return f;
  };
  if (fresh) r.grammar.add(start, form(r.start));
  for (const auto& rule : r.rules) r.grammar.add(form(rule.left), form(rule.right));
  return r;
}

// The graph of forms reachable from the start: by leftmost steps when every rule rewrites one letter (the
// histories are then the leftmost derivations), by every step when a rule has context.
struct Histories {
  std::vector<Word> forms;                  // forms[0] is the start
  struct Edge {
    std::size_t from, to, rule, position;
  };
  std::vector<Edge> edges;
  std::vector<double> paths;                // the number of histories from the start to each form
  std::vector<double> probability;          // their summed weight (rules' probabilities multiplied along a history)
  std::vector<char> final;                  // no rule applies
  bool complete = true;                     // nothing was cut by a bound
  std::string cut;                          // what cut it, when it is not complete
};

struct Bounds {
  std::size_t max_length = 12;     // of a form
  std::size_t max_forms = 200000;
};

inline Histories explore(const Rewriting& r, const Bounds& b = {}) {
  Histories h;
  std::map<Word, std::size_t> index;
  h.forms.push_back(r.start), index[r.start] = 0;
  // Breadth first; then the counts in an order in which every form follows the forms that lead to it. The leftmost
  // graph of a monotone system is acyclic in length-then-discovery order only when no rule keeps the length and
  // returns, so the counts are accumulated by a topological pass, and a cycle is reported.
  for (std::size_t at = 0; at < h.forms.size(); at++) {
    const Word w = h.forms[at];
    for (auto& step : r.steps(w, r.context_free)) {  // with context, a history is a sequence of steps at any position
      if (step.to.size() > b.max_length) {
        h.complete = false, h.cut = "forms longer than " + std::to_string(b.max_length) + " letters were not followed";
        continue;
      }
      auto it = index.find(step.to);
      if (it == index.end()) {
        if (h.forms.size() >= b.max_forms) {
          h.complete = false, h.cut = "more than " + std::to_string(b.max_forms) + " forms";
          continue;
        }
        it = index.emplace(step.to, h.forms.size()).first;
        h.forms.push_back(std::move(step.to));
      }
      h.edges.push_back({at, it->second, step.rule, step.position});
    }
  }
  const std::size_t n = h.forms.size();
  h.final.assign(n, 1);
  std::vector<std::size_t> indegree(n, 0);
  std::vector<std::vector<std::size_t>> out(n);
  for (std::size_t e = 0; e < h.edges.size(); e++) h.final[h.edges[e].from] = 0, indegree[h.edges[e].to]++, out[h.edges[e].from].push_back(e);
  h.paths.assign(n, 0.0), h.probability.assign(n, 0.0);
  h.paths[0] = 1, h.probability[0] = 1;
  std::deque<std::size_t> ready;
  for (std::size_t i = 0; i < n; i++)
    if (indegree[i] == 0) ready.push_back(i);
  std::size_t done = 0;
  while (!ready.empty()) {
    const std::size_t u = ready.front();
    ready.pop_front(), done++;
    for (std::size_t e : out[u]) {
      const auto& edge = h.edges[e];
      h.paths[edge.to] += h.paths[u];
      h.probability[edge.to] += h.probability[u] * r.rules[edge.rule].probability;
      if (--indegree[edge.to] == 0) ready.push_back(edge.to);
    }
  }
  if (done != n) h.complete = false, h.cut = "the forms return to themselves: histories are without number";
  return h;
}

// A word of the language with the histories that lead to it.
struct Sentence {
  Word word;
  double histories = 0;    // the number of derivations: more than one is an ambiguous word
  double probability = 0;
};

// The words of the language within the bounds, shortest first: final forms of letters that no rule rewrites.
inline std::vector<Sentence> language(const Rewriting& r, const Histories& h) {
  std::vector<Sentence> out;
  for (std::size_t i = 0; i < h.forms.size(); i++)
    if (h.final[i] && r.terminal(h.forms[i])) out.push_back({h.forms[i], h.paths[i], h.probability[i]});
  std::sort(out.begin(), out.end(), [](const Sentence& a, const Sentence& b) { return a.word.size() != b.word.size() ? a.word.size() < b.word.size() : a.word < b.word; });
  return out;
}

// One history drawn by the rules' probabilities, leftmost: the forms it passes through.
inline std::vector<Word> sample(const Rewriting& r, calx::random::Stream& stream, std::size_t max_steps = 10000) {
  std::vector<Word> history{r.start};
  for (std::size_t k = 0; k < max_steps; k++) {
    const auto steps = r.steps(history.back(), true);
    if (steps.empty()) break;
    double total = 0;
    for (const auto& st : steps) total += r.rules[st.rule].probability;
    if (total <= 0) break;
    double u = stream.uniform() * total;
    std::size_t chosen = steps.size() - 1;
    for (std::size_t i = 0; i < steps.size(); i++) {
      u -= r.rules[steps[i].rule].probability;
      if (u < 0) {
        chosen = i;
        break;
      }
    }
    history.push_back(steps[chosen].to);
  }
  return history;
}

}  // namespace kinesis::systems
