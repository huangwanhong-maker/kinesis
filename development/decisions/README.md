# Decision Log

Each record captures one decision that shapes the structure, interfaces or tooling of kinesis.

- Records are numbered and never rewritten.
- When a decision changes, a new record supersedes the old one, and the old record's status is updated to point to it.

## 1. Index

| No. | Decision | Status |
|---|---|---|
| [0001](0001-record-architecture-decisions.md) | Record architecture decisions | Accepted |
| [0002](0002-scope-boundaries-and-dependencies.md) | Scope, boundaries and dependencies: what a generator does over time, and generators recovered from records | Accepted |
| [0003](0003-the-first-modules.md) | The first modules, and the language | Accepted, provisionally |
| [0004](0004-averages-coupling-by-an-action-and-a-cascade.md) | Averages of what a system observes, coupling by an interaction action, and a shell model of turbulence | Accepted, provisionally |
| [0005](0005-balance-laws-on-graphs.md) | Balance laws of a density on a graph, run and learned from a few snapshots (the integral form, shooting, well-posedness), judged as a whole and as a pattern, and the modes of a field from a few snapshots | Accepted, provisionally |
| [0006](0006-prisms.md) | Prisms: what a record of positions identifies of a path under a bound on speed (Hägerstrand's space-time prism); the deepest point of discs; consistency with a witness or the pair that breaks it; reach and presence in a region, an alibi certified | Accepted, provisionally |
| [0007](0007-systems-across-scales.md) | Systems across scales: quantities on a hierarchy of units, carried between levels by named operators, advanced each at its own pace; provenance on every input and what is missing refused; ensembles and screening over declarations; fits at any level with identification said; perfect aggregation | Accepted, provisionally |
| [0008](0008-critical-transitions-in-structured-states.md) | Critical transitions in structured states: a chart for every kind of state (vectors, compositions, histograms, covariances), a neighbourhood and a hierarchy; indicators that say what they see and are blind to; the change tested against the fitted null with error rates against a declared alternative; verdicts of signal, no signal or undecidable with reasons | Accepted, provisionally |
| [0009](0009-densities-and-learned-laws-as-evaluators.md) | Densities moved by a generator (drift and diffusion varying by cell, in switching types; implicit flow, stationary density, interventions as maps), a value equation whose policy drives them, and learned laws as evaluators (held to baselines, the record's support, an ensemble's standing on an order, bagged sparse laws) | Accepted, provisionally |

## 2. Template

```markdown
# NNNN. Title

- **Status:** Proposed | Accepted | Superseded by NNNN
- **Date:** YYYY-MM-DD
- **Decided by:** …

## 1. Context
## 2. Decision
## 3. Alternatives Considered
## 4. Consequences
```
