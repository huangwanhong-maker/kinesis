# 0008. Critical Transitions in Structured States: Detectors That Say What They Can See, Read in Any Chart

- **Status:** Accepted, provisionally
- **Date:** 2026-11-19
- **Decided by:** the user ("this is an important infrastructure, no only in climate change detection, but also city goverance, like population critical transition, population/economics structural critical change / brain-dynamics disorder detection, etc. So it should support various dynamical system. So state variable is highly structural, right? and even multi-dimensional. Anyways, survey and proceed.") and Claude; the infrastructure's ADR 0135

## 1. Context

kinesis held signs of criticality for one coordinate (`analysis/criticality.hpp`): rolling variance, lag-1 autocorrelation and skewness with Kendall's τ against red-noise surrogates, detrended fluctuation analysis, avalanches and the branching ratio. The infrastructure's survey *Critical transitions in structured systems* (2026, five parts written from the literature) found three things that module did not meet:

- **Kinds of transition.** Slowing down is a theorem only for local bifurcations of codimension one under slow forcing. Transitions driven by noise, by the rate of a parameter, through global bifurcations, crises or thresholds give no such warning, and a detector must say so before it looks.
- **Error rates.** The empirical record splits. Controlled and physically guaranteed transitions are often anticipated; series chosen because a transition happened mostly are not (in field ecology 9–13% true positives against 47–53% false; no slowing before financial crashes; conflicting results before seizures; tipping times of the Atlantic overturning spanning millennia). So every warning must carry its false-alarm and detection rates at the record's length, and an absence of warning means something only where the test could have seen one.
- **Structured states.** The states that matter (a population by age and place, an economy by sector, a city by block, a network, a multichannel recording) are not scalars. Every multivariate indicator estimates the spectrum of the lag-1 propagator, the leading eigenpair of the covariance, or a projection on the critical direction, and each goes silent in ways the structure predicts. Nothing published defines slowing down in Aitchison's geometry, in Wasserstein's, or on covariance matrices.

## 2. Decision

- **A new module, `kinesis::transitions`,** after `scales` in the order, using `analysis` (the smoother) and `scales` (the hierarchy).
- **A state space declares a chart** into linear coordinates: identity; isometric log-ratios for compositions, with the log total as a coordinate when asked; the quantile function for histograms (W₂ on the line); the whitened matrix logarithm for covariances (an isometry of the affine-invariant metric at its reference). It may also declare a neighbourhood of its units and a hierarchy of them. A detector reads the chart's coordinates, so one battery serves every kind of state.
- **Every indicator declares** what it estimates, the transitions it sees, those it is blind to and what it assumes. Every assessment lists what none of its detectors can see.
- **The test is of the change** between the first window and the last, against the stationary linear process fitted to the residuals. The same pipeline is run on that null and on a declared alternative (the leading mode slowing to a tenth of its rate), giving the false-alarm rate, detection rate and ROC at the record's length. Further nulls and a grid of windows and bandwidths check robustness. Kendall's τ is reported but not tested, since it has little power with long windows (the theory note, §6.1).
- **Verdicts are three:**
  - *signal*;
  - *no signal*, only when the test was powered against the alternative;
  - *undecidable*, with every reason. The conditions that force it include a lag-1 autocorrelation already near one, and a variance that rises while a powered estimate of the restoring rate does not fall (a change of noise).
- **The restoring rate is estimated directly** by four estimators: AR(1), increments with AR(1) noise, red noise, and the drift's slope. The red-noise estimator reports the unordered pair of rates it identifies.
- **Further detectors:**
  - a likelihood test of a drifting against a constant rate, which gives a tipping time only when significant and bounded;
  - the dominant eigenvalue's path (towards +1, −1 or the circle);
  - wells and switching for noise-induced transitions;
  - on a hierarchy, the totals and the contrasts within groups at every level.
- **Test models with known answers** sit in `models.hpp`, positive and negative.

## 3. Alternatives Considered

- **Extending `analysis/criticality.hpp`.** Rejected: the structured charts, the nulls and operating characteristics, and the verdicts are a subject of their own. The older header stays for avalanches, branching and DFA, and its `early_warnings` remains for compatibility.
- **Testing Kendall's τ, as the literature does.** Rejected after the checks. A fold whose variance rises sevenfold was not significant by τ, because a stationary record's windowed series is monotone too often. τ is reported beside the test.
- **Classifiers learnt from normal forms (Bury et al. 2021 and successors).** Deferred. Their failures out of distribution would need a training manifest and a refusal, and their error rates the same nulls.
- **Detecting rate-induced tipping from a record.** Not done: the survey found no validated method from records alone. It is said to be detectable only with a model of the forcing.

## 4. Consequences

- **Checks:** 77, against closed forms and positive and negative cases. Among them:
  - false alarms of a stationary process held to the stated level;
  - rate-induced tipping at constant r showing no slowing;
  - reddening noise moving the AR(1) rate and not the red-noise rate;
  - a mode the noise does not excite giving no signal;
  - a hierarchy's totals carrying a slowing its contrasts within do not.
- **What the infrastructure can now compose:** detectors for any of its domains, whose objects are named only there. The climate domain's energy balance and overturning are the first test models with a fold in nature (ADR 0135).
- **Provisional:** the module's name; the default window (half the record) and bandwidth (a tenth); the declared alternative (a tenth of the rate); the required power (0.8); the rule of half the further nulls and half the grid.
