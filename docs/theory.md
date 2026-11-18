# Dynamical Systems Declared, Run and Recovered: What kinesis Computes, and on What It Rests

Status: **draft** (2026-11-13).

This note states what a dynamical system is taken to be in kinesis, how one is declared, what each module computes from a declaration or from a record, and why its answers are right. Statements proved here are short; the others are cited. The headers of each module carry the same sources beside the code.

## 1. Four Primitives

The design starts from four primitives, given by the user for the generative relations of their own papers. Let $\mathcal M$ be a monoid of changes and $\mathcal S$ a set of states, with an equivalence $\sim$ on $\mathcal S$.

$$
\begin{aligned}
\textbf{(Change)}\quad & \alpha : \mathcal{M} \times \mathcal{S} \to \mathcal{S} \\
\textbf{(Structure)}\quad & \alpha(m_1 m_2,\, x) = \alpha(m_1,\, \alpha(m_2, x)) \\
\textbf{(Constraint)}\quad & x \sim y \;\Longrightarrow\; \alpha(m, x) \sim \alpha(m, y) \\
\textbf{(Difference)}\quad & \Delta(m) = \{\, (x, y) : x \not\sim y \ \wedge\ \alpha(m, x) \sim \alpha(m, y) \,\}
\end{aligned}
$$

Change and Structure say that $\mathcal M$ acts on $\mathcal S$: a dynamical system is a representation of a monoid of changes. Time is one such monoid ($\mathbb N$, $\mathbb R_{\ge 0}$); the words over a set of rules are another. Constraint says that $\sim$ is a congruence of the action, which then descends to $\bar\alpha(m) : \mathcal S/{\sim} \to \mathcal S/{\sim}$. Difference collects what a change stops telling apart: under Constraint, $\Delta(m) = \emptyset$ exactly when $\bar\alpha(m)$ is injective.

The equations are philosophical, as their author says, and three generalisations make them cover the systems of the library.

1. **A change carries a weight.** For a semiring $K$, a change sends a state to a $K$-weighted family of states, and Structure becomes the sum over intermediate states, $\alpha(m_1 m_2)(x, z) = \sum_y \alpha(m_2)(x, y)\,\alpha(m_1)(y, z)$. Booleans give possibility, probabilities a stochastic process, non-negative reals histories weighted by $e^{-S}$, complex numbers amplitudes.
2. **Changes may be typed.** A rule applies only where it matches, and an event may happen only where its guard holds. $\mathcal M$ is then a category and $\alpha$ a functor on it.
3. **A constraint is declared, and may be checked.** A conserved quantity or an invariance in a declaration is a claim about the law.

The four primitives are the four parts of a declaration (§4).

*Caveat.* Constraint and Difference are stated for unweighted changes. Their weighted forms are not written here. For probabilities the weighted Constraint should be lumpability, which physis holds; that is to be shown.

## 2. The Lift to Operators

The space is replaced by its algebra of observables, and a change becomes an operator on that algebra. Let $A$ be the algebra of all functions $\mathcal S \to k$ into a field, $A_\sim \subseteq A$ the functions constant on the classes of $\sim$, and $\alpha^*(m) f = f \circ \alpha(m, \cdot)$.

**Proposition 2.1.** For an unweighted system:

1. each $\alpha^*(m)$ is an endomorphism of the algebra $A$;
2. $\alpha^*(m_1 m_2) = \alpha^*(m_2)\,\alpha^*(m_1)$;
3. $\sim$ is respected by $m$ if and only if $\alpha^*(m) A_\sim \subseteq A_\sim$;
4. under Constraint, $\Delta(m) = \emptyset$ if and only if $\alpha^*(m)$ maps $A_\sim$ onto $A_\sim$.

*Proof.* (1) and (2) are immediate: $\alpha^*(m_1 m_2) f = f \circ \alpha(m_1) \circ \alpha(m_2) = \alpha^*(m_2)(\alpha^*(m_1) f)$. (3) If $\sim$ is respected and $f \in A_\sim$, then $x \sim y$ gives $f(\alpha(m, x)) = f(\alpha(m, y))$. Conversely, if $x \sim y$ but $\alpha(m, x) \not\sim \alpha(m, y)$, the indicator of the class of $\alpha(m, x)$ lies in $A_\sim$ and its image does not. (4) $A_\sim$ is the algebra of functions on $Q = \mathcal S/{\sim}$, on which $\alpha^*(m)$ is composition with $\varphi = \bar\alpha(m)$. If $\varphi$ is injective, any $g$ on $Q$ is $f \circ \varphi$ for $f(\varphi(q)) = g(q)$, extended by zero. If $\varphi(q_1) = \varphi(q_2)$ with $q_1 \neq q_2$, every $f \circ \varphi$ agrees at $q_1$ and $q_2$, and the indicator of $q_1$ is not reached. $\square$

*Caveat.* The proposition is for the algebra of all functions. For continuous or measurable observables, (4) needs hypotheses on the spaces and maps that are not stated here.

For a classical system $A$ is commutative and $\alpha^*$ is Koopman's operator; for a weighted one $\alpha^*(m) f(x) = \sum_y \alpha(m)(x, y) f(y)$; for a quantum one $A$ is not commutative and $\alpha^*(m)$ is a completely positive map that keeps the unit. The library computes three lifts: Ulam's matrix of a map (`analysis/operators.hpp`), in which Difference is read from the operator, since a change is invertible on what a partition tells apart exactly when the matrix is a permutation; the generator of an open quantum system (`systems/quantum.hpp`, on physis); and, from a record, a finite section of Koopman's operator on a dictionary (`induction/operators.hpp`).

## 3. Levels of Access, and What a Record Identifies

A system is known at a level: by a *presentation*, a law one can read; by an *oracle*, something one can drive and read; or by a *record*, finite data seen through an observation. A presentation gives an oracle and an oracle gives records, not the reverse. A tool states the access it needs, and one asked for more than a system gives refuses with the reason (`kernel::require_access`).

**Definition 3.1 (Huang 2026, Eqs. 19–20).** For a feature $\Phi$ of histories, a record $D$ and assumptions $A$, with $\mathcal R(D; A)$ the histories the record admits under the assumptions, the *identified set* is

$$\mathcal I_\Phi(D; A) = \{\, \Phi(\Gamma) : \Gamma \in \mathcal R(D; A) \,\}.$$

It is *certified* when it is a single value, *ambiguous* when it has several, and *inconsistent* when no history satisfies the assumptions. `kernel::Identified` holds the set with its assumptions; a set known to be only a part of the identified set certifies nothing.

**Proposition 3.2 (Huang 2026, Proposition 7.1).** If $\mathcal R(D'; A') \subseteq \mathcal R(D; A)$, then $\mathcal I_\Phi(D'; A') \subseteq \mathcal I_\Phi(D; A)$: a consistent refinement can only contract the set. $\square$

A certificate therefore has a lineage against the one before it: issued, retained by refinement, withdrawn, issued again, or left ambiguous or inconsistent (`kernel::disposition`). The same object appears wherever a history is reconstructed: the derivations of a word are the histories that could have produced it, and an ambiguous word is one whose history the word does not identify (§5.5).

## 4. A System as a Declaration

A declaration names its carrier and gives the four primitives.

| Part | What it declares |
|---|---|
| `state`, `param`, `noise`, `unknown`, `input` | the carrier, its parameters, its sources of noise, what is not yet known, what another system may drive |
| `change` | the generators and how each acts: equations, an action, a Hamiltonian and jumps, rules, events |
| `structure` | how changes compose into histories and what a history weighs: the time, derivations, stationary or weighted histories, the semiring of weights, the stochastic calculus |
| `constraint` | what every change respects: a background geometry, an invariance, a conserved quantity, an identity |
| `difference` | what is told apart: the observations |
| `part`, `connect` | a system made of systems, wired at their inputs |

The background is a constraint: `x on lattice(L, L)` says that an argument ranges over a fixed geometry. A system with no such clause declares no background.

The grammar is a text in morphe's dialect of BNF over bytes (`language/grammar.hpp`), read by morphe's Earley parser, so that a script that is wrong is answered with its line, its column and what may follow there. Every expression is a *term* (`terms`): a tree in a normal form with exact symbolic derivatives, compiled to a small stack program for evaluation. A term may hold an unknown, written `?name`; a system with unknowns is a *template*.

The *form* of a system is read from its declaration (`systems::form_of`): rules make a rewriting system; a Hamiltonian or jumps an open quantum system; an action with stationary histories a Lagrangian system; an action with weighted histories a statistical field; equations with noise a stochastic system; equations in discrete time a map; equations otherwise a flow; parts a composite; `system B = derivation A` a derived system.

*Caveat.* The language reads more than the library runs. Of the derivations `vary`, `legendre`, `expand`, `reverse`, `lift` and `unravel`, only `vary` is carried out. A statistical field (`histories weight exp(-S)`), with or without a background, is read and described, and not run. An unknown function is read and not filled. Terms carry no physical dimensions. Each of these is answered with a note that says so, not with an error and never with a made-up value.

## 5. What a Declared Law Does over Time

### 5.1 Flows, Maps and Events

A flow $x' = f(t, x)$ is stepped by the classical Runge–Kutta method of physis at a fixed step; a map is iterated. The Jacobian is the derivative of the terms, exact. An equation of order $k$ makes the components $q, q', \dots, q^{(k-1)}$.

An event `when g >= 0 : resets` of a flow is a crossing: it happens when $g$ passes from negative to non-negative within a step, its time is found by bisection to a tolerance, and the state jumps by the resets. An event of a map is a level, tested after each iteration.

### 5.2 Stochastic Differentials

A system with noise is a deterministic map of state and noise, the noise a declared source with a stream of its own keyed by (seed, path, source): a path depends on its key alone. A stochastic differential is written as on paper, $dX_i = a_i\,dt + \sum_j b_{ij}\,dW_j$; its right side must be linear in the differentials, and the coefficients are read off as terms. A Stratonovich equation is converted once, exactly, to Itô's with the same solutions, $a_i^{\text{Itô}} = a_i^{\text{Strat}} + \tfrac12 \sum_{j,k} b_{kj}\,\partial b_{ij}/\partial x_k$ (Kloeden and Platen 1992, Eq. 4.9.15).

Paths are stepped by Euler–Maruyama (strong order ½) or, for diagonal noise, by Milstein's scheme (strong order 1). Noise is diagonal when each source drives one component, each component is driven by at most one source, and no coefficient depends on another driven component: then $b_j\,\partial b_i/\partial x_j = 0$ for $j \neq i$, the noise commutes, and the scheme needs no iterated integral (Kloeden and Platen 1992, §10.3). Otherwise Milstein's scheme is refused.

### 5.3 Systems Given by an Action

For `action S = integral(L)` with stationary histories, `vary` forms the Euler–Lagrange equations from the terms, $M\,q'' = F$ with $M_{ij} = \partial^2 L/\partial q_i'\partial q_j'$, and solves them for the accelerations as terms by Cramer's rule, for at most three coordinates. The result is a declaration of a flow again, with the energy $\sum_i q_i'\,\partial L/\partial q_i' - L$ declared conserved when $L$ does not hold the time. A Lagrangian whose $M$ is singular as a term is refused.

### 5.4 Open Quantum Systems

A Hamiltonian $H$ and jump operators $L_k$ give the generator of Gorini, Kossakowski, Sudarshan and Lindblad,

$$\dot\rho = -i[H, \rho] + \sum_k \Big( L_k \rho L_k^\dagger - \tfrac12 \{ L_k^\dagger L_k, \rho \} \Big),$$

which physis holds as a real matrix on the space of density matrices. A state is evolved by the exact propagator of one step; the stationary states are the kernel of the generator, and the gap is the least rate of decay among its other eigenvalues. Operators are terms of named operators (`sx`, `sm`, `a`, `adag`, `num`, …) on the factors of the state space; since terms commute, a product of two operators on one factor is written `ordered(…)`.

### 5.5 Rules

A rewriting system is a dynamical system whose state is a word, whose changes are rules applied at places, and whose histories are derivations. `explore` builds the graph of forms within bounds and counts the histories to each; a word with more than one history is one whose history is not identified in the sense of §3. For rules that rewrite one letter the histories counted are the leftmost derivations, which correspond one to one with parse trees; with context, every step is followed. Nothing is reported as the whole language unless the exploration was complete.

## 6. What a System Does over a Long Time

**Lyapunov's exponents.** Tangent vectors are carried by the exact Jacobian and orthonormalised (Benettin, Galgani, Giorgilli and Strelcyn 1980).

**Proposition 6.1.** The exponents of a flow sum to the time average of the divergence of its field along the path; those of a map to the average of $\log|\det Df|$. For a field of constant divergence the sum is therefore known exactly, whatever the path.

*Proof.* The volume of the parallelepiped of the tangent vectors obeys $\dot V = (\operatorname{div} f)\,V$ (Liouville), and the sum of the exponents is the growth rate of $\log V$. $\square$

This is the check of the method: Lorenz's exponents sum to $-(\sigma + 1 + \beta)$, Hénon's to $\log b$, the standard map's to zero.

An exponent of a chaotic path settles as one over the square root of the time. With a running record, the error of a finite time is reported: the time is cut into blocks, the mean rate over each is read from the record, and the standard error of their mean is given, the first block left out. *Caveat.* The blocks are taken as independent, which they are when each is long against the time the path remembers; otherwise the error is understated.

**Equilibria** are found by Newton's method from seeds spread over a box, and classified by the eigenvalues of the exact Jacobian. What is found is every root the seeds converged to; the answer never says there is no other.

**Basins.** Each state of a grid is run, and the end states are clustered. An end state is run on for four times as long before it is counted as an attractor of its own, and each end is named: at rest, by its linearisation (a stable node, a saddle, …), or still moving. An end at rest that is not stable is reached only from the boundary between basins. The share of a box that leads to an attractor is its basin stability (Menck, Heitzig, Marwan and Kurths 2013), with the standard error $\sqrt{S(1-S)/N}$ from $N$ independent draws.

**Landscapes.** The word covers several functions, which agree only in special cases (Zhou, Aliyu, Aurell and Huang 2012), and each answer says which it is. A flow of one dimension has the potential $U = -\int f\,dx$; in more dimensions a flow is a gradient exactly when its Jacobian is symmetric on a simply connected domain, which is tested at points, and is otherwise refused a potential. For a stochastic differential of one dimension the stationary landscape $-\ln p$ is exact, $p \propto b^{-2}\exp\!\int 2a/b^2\,dx$ (Gardiner 2009, §5.3.3). For additive noise it is $2U/\sigma^2$, and Kramers's time over a barrier is $\tau = 2\pi\,(U''(x_a)\,|U''(x_b)|)^{-1/2}\exp(2\Delta U/\sigma^2)$ in the limit $\Delta U \gg \sigma^2/2$. In more dimensions the landscape is estimated from long paths as minus the logarithm of a histogram, and is not a Lyapunov function of the deterministic part in general.

**Averages.** A system that neither settles nor repeats is described by its statistics. `averages` follows the observations a declaration names along one long path: each mean comes with the standard error of its block means (Flyvbjerg and Petersen 1989), and with its second moment over its mean squared, which is 1 for what does not fluctuate and, for the squared size of a complex amplitude, the flatness of that amplitude (2 when it is Gaussian). Means that follow a power law of a scale are fitted on their logarithms, with the error of the slope. *Caveat.* The error of a mean is right when a block is long against what the observation remembers; the error of an exponent is of the line through the means, over the octaves fitted and no further.

**Signs of criticality** are read from records, and none is a certificate. Slowing down before a transition is measured by the trends of the variance and of the autocorrelation in a rolling window, by Kendall's $\tau$, with a significance against red-noise surrogates (Dakos et al. 2012). Correlations without a scale are measured by detrended fluctuation analysis; avalanches by their sizes, with the exponent of a power law fitted by maximum likelihood (Clauset, Shalizi and Newman 2009); and the branching ratio from the decay of the regression slopes at increasing lags, which is unbiased when only a share of the units is observed (Wilting and Priesemann 2018). *Caveat.* A rising indicator is not a transition, and a power law over a decade or two does not show criticality; the headers cite what each is known to mistake.

## 7. Braids of Trajectories

The setting is that of Huang (2026), cited as the draft. A record holds the positions of $n$ labelled points of the plane at finitely many times. A history is a continuous, collision-free motion of the points through the recorded positions, and its class is a morphism of the fundamental groupoid of the configuration space. A record does not determine it: between two samples a point may turn about another any number of times (the draft, Proposition 6.6).

**The reading.** The braid of the piecewise-linear interpolant is read by projection on a line: the points are ordered by their projection, and each time two neighbours exchange places a letter $\sigma_i^{\pm 1}$ is written, $i$ the position of the left one, the sign positive when the left one passes below (the counterclockwise half-exchange). Configurations with distinct projections form contractible cells (the top cells of Fox and Neuwirth 1962), so a path between two such configurations has a braid word that is well defined once the line is chosen. The order of two points at a sample is the sign of the difference of their projections, zero counted as positive at both ends of every interval, so that an exchange that falls on a sample is read once and a touch that does not pass is read as two letters that cancel.

**Proposition 7.1 (the certificate).** Let each point $i$ move no faster than $V_i$ and be observed within $\eta$. Let $r_k, r_{k+1}$ be the observed positions of $j$ relative to $i$ at the two ends of an interval of length $h$. If

$$|r_k| + |r_{k+1}| - (V_i + V_j)\,h - 8\eta > 0 \tag{$*$}$$

for every interval and every pair, then every history the assumptions admit is homotopic, in the configuration space, to the interpolant.

*Proof.* The true relative position stays in the ellipse with foci at the true $r_k, r_{k+1}$ and string $(V_i + V_j)h$, which is convex. With $|\tilde r - r| \le 2\eta$ for each observed relative position, $(*)$ puts the origin outside that ellipse, so the straight-line homotopy from the true path to the true interpolant never meets the origin, for all pairs at once. The observed and the true interpolants differ by at most $2\eta$ in each relative position, and $(*)$ keeps the observed chord farther than $2\eta$ from the origin, so the straight-line homotopy between them is collision-free too. $\square$

With $\eta = 0$ and one bound $V$, $(*)$ is implied by the draft's Corollary 6.4 and is weaker; for two points it is sharp, in that where it fails strictly another class is admissible. The ends of the two interpolants lie in one cell when the projections at the first and last samples are farther apart than $2\eta$, which is asked as well. The certificate was derived for this library, and the draft's coarser margin is reported beside it.

**The windings of a pair.** For two points with exact samples and a least separation $\rho$, the classes of the relative path on an interval are its total angles $\Theta = \theta_0 + 2\pi m$. The shortest path of total angle $\Theta$ that keeps $|r| \ge \rho$ is the chord when it has that angle and clears the disc, and otherwise two tangents and an arc, of length

$$\sqrt{a^2 - \rho^2} + \sqrt{b^2 - \rho^2} + \rho\,\big(|\Theta| - \arccos(\rho/a) - \arccos(\rho/b)\big), \qquad a = |r_k|,\ b = |r_{k+1}|,$$

and a class is admissible exactly when this length is at most $(V_i + V_j)h$. The winding of a closed relative path is the sum over the intervals, so its identified set is a sum of finite sets of integers, computed exactly. With $\rho = 0$ one admissible turn admits every number of turns, and the set is reported as not complete.

**Entropy.** A braid is a mapping class of the punctured disc, and its topological entropy is the logarithm of its dilatation when the class is pseudo-Anosov. It is estimated by the growth of Dynnikov's coordinates of the generating loops under repeated application (Moussafir 2006; Thiffeault 2010), in exact integers, and bounded from below by the spectral radius of the Burau matrix on the unit circle (Fried 1986; Kolev 1989; Band and Boyland 2007). The braid groups themselves, Dehornoy's handle reduction, the Dynnikov action and the Burau representation are morphe's.

*Remark.* A braid lifted to an operator is its Burau matrix, and the bound above is a statement about that operator's spectrum: the one case of §2 built for braids.

## 8. Laws from Records

**One mechanism.** A template is a declaration with unknowns, and a model is matched, discovered or fitted by filling them (`induction::fit`). The equations a record must satisfy are those of the declaration: for a flow each equation is held at every sample, with differences of the record on the left; for a map the left side is the next sample; for a system given by an action they are its Euler–Lagrange equations, so that an unknown inside an action is fitted through the equations it implies. When every equation is affine in the unknowns they solve a linear least-squares problem; otherwise the Levenberg–Marquardt iteration is used, with the Jacobian derived from the terms. The result is a system of the language again.

**Discovery.** With no law proposed, the template is a sparse combination of a library of monomials, thresholded sequentially (Brunton, Proctor and Kutz 2016); in the weak form the equations are tested against windows and the regression takes no derivative of the record (Schaeffer and McCalla 2017; Messenger and Bortz 2021).

**Stochastic templates.** The drift is the first Kramers–Moyal coefficient of the increments and the noise the second (Risken 1989); both are of the first order in the step of the record.

**Matching.** Candidates are each fitted to the first part of a record and run on the rest, and ordered by Schwarz's criterion on the one-step errors, $N \ln(\text{mean square error}) + k \ln N$, so that a model with more freedom must earn it. A model of the catalogue becomes a candidate with its parameters freed and its states named after the record's, in each order of the components (§10).

**No law at all.** A record can be continued by a forecaster with no law in it (the next-generation reservoir computer of Gauthier, Bollt, Griffith and Barbosa 2021), and a linear operator can be fitted to it (dynamic mode decomposition; Schmid 2010).

*Caveat.* A residual measures how well equations hold on a record, not whether a template is right. A ranking is among the candidates tried. On a chaotic record no model stays near for long, whatever its merit.

### 8.1 Balance Laws on Graphs, from a Few Snapshots

A density $u$ on the cells of a graph, cell $i$ of area $V_i$, cells sharing a face joined by an edge $e$ of conductance $c_e$ (the face's width over the distance between the centres), changes by what crosses the faces and what is made within. With the flux density $q = -D\nabla u - D_2\, u \nabla u + \chi\, u \nabla\Phi$, the finite-volume flux from $a$ to $b$ is

$$F_e = c_e\Big[-D\,(u_b - u_a) - D_2\,\tfrac{1}{2}(u_b^2 - u_a^2) + \chi\,(\Phi_b - \Phi_a)\,\tfrac{1}{2}(u_a + u_b)\Big],$$

and $\dot u_a = -V_a^{-1}\sum_{e \ni a} F_e + r u_a + s\, u_a^2 + g\, u_a\ln u_a + \dots$ Every flux is conservative by construction, so that $\sum_i V_i \dot u_i$ is what the reaction makes. The law is linear in its coefficients $\theta = (D, D_2, \chi, r, s, g, \dots)$: $\dot u = \sum_k \theta_k T_k[u]$.

**Proposition 8.1** (well-posedness). *With $D < 0$ and no other term the law is the backward heat equation, whose Fourier mode of wavenumber $k$ grows as $e^{|D|k^2 t}$; on the grid the largest growth rate is $8|D|/h^2$, and a run amplifies the smallest wiggle of the record by $e^{8|D|t/h^2}$.* So a negative diffusion found by least squares is set aside, and the other terms fitted again (an active set). *Proof.* The five-point Laplacian has eigenvalues $-(4/h^2)(\sin^2(k_xh/2) + \sin^2(k_yh/2)) \in [-8/h^2, 0]$. $\square$

**The integral form.** Between snapshots $t_k$ and $t_{k+1}$, $u(t_{k+1}) - u(t_k) = \sum_j \theta_j \int T_j[u]\,dt \approx \sum_j \theta_j\, \tfrac{\Delta t}{2}\big(T_j[u(t_k)] + T_j[u(t_{k+1})]\big)$: one linear equation per cell and interval, with no derivative of the record. Its bias is the trapezoid's, $O(\Delta t^3\, \dddot u)$ per interval, small when the field changes little between snapshots.

**Shooting.** The law is run from $u(t_0)$ and $\theta$ chosen to minimise $\sum_{k \ge 1}\sum_i w_{ik}^2\,(u_i^\theta(t_k) - u_i(t_k))^2$ by Levenberg–Marquardt, from the integral form's estimate. It has no bias of the step: its error is the record's.

**Judging a forecast.** On snapshots kept aside, skill is $1 - \text{SSE}/\text{SSE}_0$ against carrying the last snapshot forward. The *pattern's* skill rescales each forecast, and the carried snapshot, to the record's total at each time before comparing: it measures where the quantity went apart from how much there was. Uniform growth has pattern skill zero by construction, since its forecast is the carried snapshot rescaled.

**Modes from a few snapshots.** With $n$ cells far more than the $m$ snapshots, exact DMD fits the operator on the snapshots' span: $X^\top X = V\Sigma^2 V^\top$ (Sirovich's method of snapshots), $\tilde A = \Sigma^{-1}V^\top (X^\top Y)V\Sigma^{-1}$, its eigenpairs $(\lambda_k, w_k)$, the modes $\phi_k = Y V \Sigma^{-1} w_k/\lambda_k$ (Tu et al. 2014). For a field moved by a linear law the modes and eigenvalues are the law's.

*Caveat.* What a law with few shared coefficients can describe is what it holds. A city's density is pinned by its buildings: a diffusion strong enough to move people a kilometre a year would flatten its centre within a year, and is not found in its record. The standard errors treat neighbouring cells as independent, and so understate.

## 9. Systems of Systems

A composite declares parts and wires an input of one to a term of the others' states.

**Proposition 9.1.** When every part is a flow, the composite is the flow whose field is the parts' fields with each input replaced by the term wired to it; `merged` builds its declaration, and no error is added by the composition. $\square$

**Joined by an action.** When every part is given by its action, the composite may give the coupling as an action of its own. The whole is then the system whose Lagrangian is $L = \sum_{\text{parts}} L_{\text{part}} + L_{\text{int}}$, and `joined_action` builds it as one declaration, which is varied as any action is.

**Proposition 9.2.** A coupling given by an interaction action is reciprocal, and the whole keeps what its Lagrangian is invariant under: its energy when $L$ does not hold the time; its total momentum $\sum_a \partial L/\partial q_a'$ when $L_{\text{int}}$ depends on differences of coordinates alone and the parts are free. A coupling given by a wire need not keep either.

*Proof.* The Euler–Lagrange equation of a coordinate $q_a$ of one part holds the force $\partial L_{\text{int}}/\partial q_a$, and that of $q_b$ the force $\partial L_{\text{int}}/\partial q_b$: both come from the one term, and when it depends on $q_a - q_b$ alone they are equal and opposite. The conservation laws are Noether's theorem for the whole Lagrangian (Arnold 1989, §20). For the converse, the wire $b.u = k(a.q - b.q)$ on two free masses moves $b$ and leaves $a$ in uniform motion, so the total momentum changes: a check of the library. $\square$

*Remark (the question of a field theory of couplings).* An interaction action is therefore the form of couplings that are forces. Couplings between living things are mostly not: what eats is not eaten back. They are wires. Every wired composite of flows is nonetheless the stationarity of one action in doubled variables, the response-field action $S[\tilde x, x] = \int \sum_a \tilde x_a\,(\dot x_a - f_a(x_a) - g_a(x))\,dt$ of Martin, Siggia and Rose (1973), stationary in $\tilde x$: there a wire is the term $-\int \tilde x_a\,g_a(x)\,dt$ of an interaction action, and noise adds a term quadratic in $\tilde x$. That form treats reciprocal and one-way couplings alike. It is recorded here and not built.

Parts of different forms have no common presentation and are advanced side by side: at each communication step $H$ the wires are read on the states of all parts, the inputs are held, and each part is advanced over $H$ by its own method (the explicit scheme of Jacobi; Kübler and Schiehlen 2000). Its error is of the first order in $H$ where a wired signal varies within a step, and it is exact where the signals are constant between steps, as when a map drives a flow. The order of the parts does not matter.

**Systems across scales (decision 0007).** A system across scales lives on a hierarchy $\mathcal H$: levels $L_0, \dots, L_m$, the finest first, with maps $\pi_\ell : L_\ell \to L_{\ell+1}$ onto the level above, so that the fibres $\pi_\ell^{-1}(v)$, the children of $v$, partition $L_\ell$. A quantity of level $\ell$ is a function $x : L_\ell \to \mathbb R$. Between levels there are two kinds of linear operator, and a model is consistent across its levels when it uses them for every crossing:

- *carrying up*, $(\Sigma x)(v) = \sum_{u \in \pi^{-1}(v)} x(u)$, and the weighted mean $(M_w x)(v) = \Sigma(wx)(v)/\Sigma w(v)$;
- *allocating down*, $(A_s X)(u) = s(u)\,X(\pi(u))$ with shares $s \ge 0$ summing to one on every fibre (proportional, $s = y/(\Sigma y)\circ\pi$; logit, $s = e^{y}/(\Sigma e^{y})\circ\pi$).

**Proposition 9.3.** $\Sigma A_s = \mathrm{id}$ on quantities of the level above: what is allocated sums to what was there. For extensive $a, b$ the index $\Sigma a / \Sigma b$ equals the weighted mean $M_b(a/b)$ of the units' ratios, which in general differs from their plain mean.

*Proof.* $\Sigma A_s X(v) = X(v)\sum_{u\in\pi^{-1}(v)} s(u) = X(v)$; and $M_b(a/b) = \Sigma(b\cdot a/b)/\Sigma b = \Sigma a/\Sigma b$. $\square$

A non-linear index of the units beneath (Gini's coefficient, the dissimilarity, a concentration) is a reduction $\rho(x|_{\pi^{-1}(v)})$, recomputed at each level from the level below and never averaged; a decomposable one (Theil's) has within and between parts.

**Paces.** Each law has its step $h_k = n_k h$. In the base step $[t_n, t_{n+1}]$ every law whose step begins at $t_n$ reads the states of all others as they are at $t_n$ and advances over its own step; its result is applied when its step ends. This is the explicit (Jacobi) scheme of co-simulation above; its error is of the first order in a slow law's step wherever a fast state varies within it (a check on a linear two-rate system finds the error halving with the step), and when all laws share one step it is Runge–Kutta's classical step on the whole.

**Proposition 9.4 (Tikhonov; Fenichel).** For $\dot x = f(x, y)$, $\varepsilon\,\dot y = g(x, y)$ with $\partial g/\partial y$ uniformly negative definite on a compact piece of the slow manifold $g = 0$, the solution from a point near that manifold stays within $O(\varepsilon)$ of the reduced system $\dot x = f(x, h(x))$, $g(x, h(x)) = 0$, on finite intervals (Tikhonov 1952; Fenichel 1979). A law in *equilibrium* form solves its fixed point at every instant and is that reduced system; a check finds the error of the fast flow falling in proportion to $\varepsilon$.

**Proposition 9.5 (perfect aggregation).** A linear law $\dot x = Mx$ and an aggregation $y = Ax$ ($A$ of full row rank) give a closed law $\dot y = \tilde M y$ for every initial $x$ if and only if $AM = \tilde M A$ for some $\tilde M$, that is, iff the row space of $A$ is invariant under $M^{\mathsf T}$; then $\tilde M = AMA^{\mathsf T}(AA^{\mathsf T})^{-1}$ (Iwasa, Andreasen and Levin 1987). For a partition summed, the condition is that for every two groups $I, J$ the column sums $\sum_{i\in I} M_{ij}$ agree for all $j \in J$.

*Proof.* If $AM = \tilde MA$ then $\dot y = A\dot x = AMx = \tilde M Ax = \tilde M y$. Conversely, if $\dot y$ is a function of $y$ alone for every $x$, it is linear and $AMx = \tilde M A x$ for all $x$. Multiplying $AM = \tilde MA$ on the right by $A^{\mathsf T}(AA^{\mathsf T})^{-1}$ gives $\tilde M$. $\square$

**How each input is known, and what is not.** Every value carries its provenance (observed, estimated with a standard error, declared with a range, modelled). Before a run the model lists what it needs and lacks; a run with anything lacking is refused. This is invariant 4 applied to inputs: no tool answers with a made-up value, and a quantity no record holds is the person's to declare.

## 10. The Catalogue

A model of the catalogue, of which there are thirty, is a form of law from the literature, written as a declaration with the values of its parameters in common use and its source (`catalogue/models.hpp`). It is a starting point for a study and a hypothesis for a record. **A cascade of turbulence.** The catalogue holds the Sabra shell model (L'vov, Podivilov, Pomyalov, Procaccia and Vandembroucq 1998): wavenumbers $k_n = k_0\lambda^{n-1}$, one complex amplitude $u_n$ a shell, and

$$\dot u_n = i\,\big(a\,k_{n+1}\,u_{n+2}u_{n+1}^* + b\,k_n\,u_{n+1}u_{n-1}^* - c\,k_{n-1}\,u_{n-1}u_{n-2}\big) - \nu k_n^2 u_n + f_n, \qquad a + b + c = 0.$$

It is a model of the cascade and not a discretisation of the equations of a fluid: it has no space. What is exact of it is checked. With no viscosity and no force the nonlinear term keeps $\sum |u_n|^2$ and $\sum (a/c)^n |u_n|^2$, as an identity of the field at every state. In a statistically steady state the mean power put in equals the mean power dissipated, and the mean flux of energy through a shell is the power put in less what the shells up to it dissipate, so that through the shells of the cascade it is one number: the model's form of Kolmogorov's four-fifths law. What is measured is given with its error: the exponent of $\langle |u_n|^2\rangle \propto k_n^{-\zeta_2}$, for which Kolmogorov's argument gives $2/3$, and the growth of the flatness down the cascade, which is intermittency.

Since a record's components have no names a catalogue knows, a plain model (a flow or a map given by first-order equations alone) is offered as a candidate with its states named after the record's in every order, for three components or fewer.

## 11. What Is Not Built

- The derivations `expand` (a field theory into the rules that glue its vertices), `reverse` (the time reversal of a stochastic differential, with its score), `lift`, `legendre` and `unravel`; statistical fields as systems that run.
- Unknown functions, and learning them; the score of a diffusion model.
- Oracles as declarations: a physical or simulated system driven as a black box.
- The response-field action of a wired composite (§9); families of components indexed by a number in the grammar (fourteen shells as one line; systems across scales are declared as structures, decision 0007); higher-order multirate schemes, the profile likelihood and history matching for systems across scales; a discretised fluid beyond balance laws of a density on a graph (§8.1, decision 0005): momentum, pressure and the equations of a fluid itself.
- The weighted forms of Constraint and Difference; a check that a declared invariance holds.
- Samplers of laws known up to normalisation, agreed to belong in tyche.
- Garside's normal form and conjugacy of braids; records in which points arrive and leave, beyond cutting them into segments; representations of braids beyond Burau's.
- Physical dimensions of terms.
- The move of physis's general analysis (its spectra of exponents, multipliers and continuation) into this library, agreed and deferred: kinesis uses physis's.

## 12. References

- Band, G. and Boyland, P. (2007). The Burau estimate for the entropy of a braid. *Algebr. Geom. Topol.* 7.
- Benettin, G., Galgani, L., Giorgilli, A. and Strelcyn, J.-M. (1980). Lyapunov characteristic exponents for smooth dynamical systems and for Hamiltonian systems. *Meccanica* 15.
- Brunton, S. L., Proctor, J. L. and Kutz, J. N. (2016). Discovering governing equations from data by sparse identification of nonlinear dynamical systems. *Proc. Natl. Acad. Sci.* 113.
- Clauset, A., Shalizi, C. R. and Newman, M. E. J. (2009). Power-law distributions in empirical data. *SIAM Rev.* 51.
- Dakos, V. et al. (2012). Methods for detecting early warnings of critical transitions in time series. *PLoS ONE* 7.
- Fox, R. and Neuwirth, L. (1962). The braid groups. *Math. Scand.* 10.
- Flyvbjerg, H. and Petersen, H. G. (1989). Error estimates on averages of correlated data. *J. Chem. Phys.* 91.
- Gardiner, C. (2009). *Stochastic Methods*, 4th ed. Springer.
- Gauthier, D. J., Bollt, E., Griffith, A. and Barbosa, W. A. S. (2021). Next generation reservoir computing. *Nat. Commun.* 12.
- Huang, W. (2026). Rethinking trajectory-based credibility under symbol inflation: possibilities for building epistemic infrastructure through topological methods with braid group theory as an example. Draft.
- Kloeden, P. E. and Platen, E. (1992). *Numerical Solution of Stochastic Differential Equations*. Springer.
- Kübler, R. and Schiehlen, W. (2000). Two methods of simulator coupling. *Math. Comput. Model. Dyn. Syst.* 6.
- L'vov, V. S., Podivilov, E., Pomyalov, A., Procaccia, I. and Vandembroucq, D. (1998). Improved shell model of turbulence. *Phys. Rev. E* 58.
- Martin, P. C., Siggia, E. D. and Rose, H. A. (1973). Statistical dynamics of classical systems. *Phys. Rev. A* 8.
- Menck, P. J., Heitzig, J., Marwan, N. and Kurths, J. (2013). How basin stability complements the linear-stability paradigm. *Nature Phys.* 9.
- Risken, H. (1989). *The Fokker–Planck Equation*, 2nd ed. Springer.
- Schmid, P. J. (2010). Dynamic mode decomposition of numerical and experimental data. *J. Fluid Mech.* 656.
- Thiffeault, J.-L. (2010). Braids of entangled particle trajectories. *Chaos* 20.
- Wilting, J. and Priesemann, V. (2018). Inferring collective dynamical states from widely unobserved systems. *Nat. Commun.* 9.
- Zhou, J. X., Aliyu, M. D. S., Aurell, E. and Huang, S. (2012). Quasi-potential landscape in complex multi-stable systems. *J. R. Soc. Interface* 9.

The headers cite the further sources of each method.
