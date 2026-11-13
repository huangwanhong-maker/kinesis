# 0002. Scope, Boundaries and Dependencies

- **Status:** Accepted
- **Date:** 2026-11-13
- **Decided by:** the user (an eleventh foundation: "yes"; "braid group is sehr important"; "we even need create a 'domain specific langauge' that can define dynamical systems"; the name kinesis, whose repository they create) and Claude; the infrastructure's ADR 0103

## 1. Context

The user opened a project on dynamical systems: to analyse and draw systems classical and open quantum, deterministic and stochastic, continuous and not; to induce a structured representation from data; to treat a system known only as a black box; and to study trajectories by topology and algebra, braid groups first. The discussion is recorded in the infrastructure's idea of 2026-11-12, *A Language of Dynamical Systems*, and a survey of the literature in its `research/literature/reports/Dynamical systems models and braids.md`.

The foundations already held parts of the subject: calx its integrators, physis its spectra of Lyapunov exponents, its continuation and the laws of physics as generators, morphe its grammars, tyche its time series. What none held was a way to say what a system is that serves every kind, and the study of a system as such.

Three things had to be settled: whether the subject is a foundation; where it stands against physis; and what the library promises of an answer drawn from data.

## 2. Decision

1. **A foundation.** The study of what a law does over time, and the recovery of a law from a record, name no object in nature, and every domain whose objects change may use them. The library is a foundation, kinesis. The application, the Dynamics Lab, is a work of the infrastructure; the drawings are the visualization layer's.
2. **On calx, morphe, physis and tyche.** kinesis uses calx's linear algebra and streams, morphe's BNF reader, grammars and braid groups, physis's Runge–Kutta step, spectra of exponents and Lindblad generator, and tyche's statistics. It depends on nothing else.
3. **The rule against physis.** physis writes generators: the laws of mechanics, fields and quanta. kinesis studies a generator of any kind (what it does over time, what it keeps, where it changes its behaviour) and recovers generators from records. A law of physics is never written in kinesis; a method that would serve a grammar as well as a pendulum is.
4. **Braids in two places.** The braid groups, their reduction, their coordinates and their representations are algebra, and are morphe's (its decision 0014). Reading a braid from a record of trajectories, and what such a record certifies, is kinesis's, on the user's draft *Rethinking trajectory-based credibility under symbol inflation*.
5. **What an answer from data is.** A system is known by a presentation, an oracle or a record, and a tool states which it needs. An answer drawn from a record is the set of values over the histories the record admits under stated assumptions, certified only when it is one value. This is the form of every such answer in the library, not only of braids.
6. **Randomness.** The user agreed that samplers of laws belong in tyche. Steppers of stochastic differentials are part of running a declared system and are here; they use calx's keyed streams.
7. **Out of scope:** the laws of physics (physis); numerical methods as such (calx); distributions and estimation (tyche); search (heuresis); what a system is of (the domains); drawing (the infrastructure's `visualization/`).

## 3. Alternatives Considered

- **A module of physis.** physis's `dynamics` holds flows and their exponents. Rejected with the user: the subject includes grammars, stochastic models, data and inference, which are not physics; and a language with its reader is a different kind of thing from a module of laws. The agreed move of physis's general analysis up into kinesis follows from the same rule and is deferred (§4).
- **Several libraries** (a language, an analysis, an inference). The declaration is what joins them: a fitted model must be a system that every tool reads. One library keeps that true.
- **A domain.** A domain names objects in nature. Lorenz's system in the catalogue is a form of law, and the atmosphere is meteorology's.
- **Answers from data as plain values.** Rejected: a braid read from samples is one of many the samples admit, and a library that reports it without its assumptions will sometimes report what is false.

## 4. Consequences

- A new kind of system is added as a form of the language, and each tool states whether it serves it.
- **Deferred, and said to be.** The move of physis's general analysis (spectra of exponents, Floquet multipliers, continuation of equilibria) into kinesis, as panels moved from nomos to tyche, was agreed and is not done: kinesis calls physis's. It is to be proposed to the user as a change of its own, since it moves code and checks between two repositories.
- **Deferred.** Samplers in tyche are not written; nothing in kinesis needs one yet.
- Domains that hold dynamics of their own (neuroscience's spiking networks, urban science's flows, physical computing's substrates) keep them. Where one is a general method, it is proposed here.
