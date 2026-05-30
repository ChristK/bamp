# bamp 2.2.0 — Detailed guide to the changes since 2.1.3

*A technical companion to `NEWS.md`, written to be read end to end. It states the
model and notation once, then explains each new feature and bug fix: the problem
it addresses, the method (with the relevant mathematics), how it is exposed to
the user, and the evidence that it is correct. Sections are largely
self-contained; read §1–§2 first.*

---

## 1. Background, notation, and the problem this release solves

`bamp` fits the Bayesian age–period–cohort (APC) model of Knorr-Held and Rainer
(2001) and Schmid and Held (2007). The data are counts on a Lexis diagram: for
each age group *i* (i = 1, …, I) and period *j* (j = 1, …, J) we observe a number
of cases `y_ij` out of a population at risk `n_ij`. The first-stage model is
binomial,

  y_ij ~ Binomial(n_ij, p_ij),  logit(p_ij) = η_ij,

with the linear predictor decomposed into an intercept and smooth age, period and
cohort effects,

  η_ij = μ + θ_i + φ_j + ψ_k,  k = cohort index = (I − i)·M + j,

where M = `periods_per_agegroup` is the number of periods spanned by one age
group. (Cohort *k* runs from 1 to K = M·(I − 1) + J.) Optional extensions add a
per-cell heterogeneity/overdispersion term and covariates; these are treated in
§4 and §5.

The effects are given **intrinsic Gaussian Markov random-field (GMRF) priors** —
random walks of first or second order (RW1 / RW2). For an effect *x* of length L
with difference operator D of order d (d = 1 for RW1, d = 2 for RW2), the prior
density is

  p(x | κ) ∝ κ^((L − d)/2) · exp( −(κ/2) · xᵀ K x ),  K = DᵀD,

where κ is a precision (smoothing) parameter with a Gamma(a, b) hyper-prior. K is
rank-deficient (rank L − d): it does not penalise a constant (RW1) or a constant
and a linear trend (RW2), which is what makes the prior "intrinsic" and is the
source of the identifiability issues discussed in §7.

**The problem.** The APC likelihood is only weakly informative about each effect
separately, and for incidence/mortality data the cells are highly heterogeneous:
some have large populations and many events, others (young ages, edge cohorts,
rare causes) have huge populations but almost no events. The legacy sampler
(`method = "iwls"`) updates the blocks one at a time with an
iteratively-reweighted-least-squares Metropolis–Hastings proposal. On such data
it mixes poorly: for RW2 models it routinely failed the convergence check, and on
sparse / zero-cell data its automatic safeguards could prune every chain and
return nothing. The diagnosis that motivated this release found that *part* of
the apparent non-convergence was a reporting artefact (§7) but a substantial part
was genuine, and the cure is a better sampler (§2–§3).

Throughout, write β = (μ, θ, φ, ψ) for the stacked parameter vector of length
P = 1 + I + J + K (extended in §4–§5), and z_ij = (y_ij − n_ij/2) for the
"centred" response that recurs in the Polya-Gamma algebra.

---

## 2. The new default engine: Polya-Gamma Gibbs with a Laplace–Newton refinement

### 2.1 Why Polya-Gamma

Polson, Scott and Windle (2013) showed that a logistic likelihood becomes
conditionally Gaussian under a Polya-Gamma (PG) latent variable. If
ω_ij ~ PG(n_ij, η_ij), then conditional on ω the contribution of cell (i, j) to
the log-posterior is exactly a Gaussian in η_ij with "working response"
z_ij / ω_ij and weight ω_ij. Because η is linear in β, the **entire field β
becomes conditionally Gaussian**:

  β | ω, κ ~ N( Q⁻¹ b, Q⁻¹ ),  Q = Xᵀ diag(ω) X + prior precision,  b = Xᵀ z,

where X is the design matrix mapping β to the I·J cell predictors. This lets us
draw the intercept and all three effects **jointly in one exact step** (the
"one-block" GMRF sampler of Rue and Held, 2005), rather than one block at a time.
Joint sampling removes the cross-block autocorrelation that cripples
one-at-a-time updates when the data are informative — exactly the
large-population regime of incidence/mortality data.

This package draws the PG weights with a fast Gaussian (mean/variance-matched,
central-limit) approximation to PG(n, η), which is essentially exact for the
large counts `n` here and analytically continuous as η → 0.

### 2.2 Why a Newton refinement is also needed

Pure PG-Gibbs has a known weakness: in cells where the data are very
informative, the augmented full-conditional of β is far *tighter* than its
marginal, so the sampler takes tiny steps and mixes slowly. On the package's own
example this left ~25% of cells with a fitted-value Gelman–Rubin statistic near
2 even after hundreds of thousands of iterations. The remedy is to interleave the
Gibbs sweep with a **joint Laplace (Newton) Metropolis–Hastings proposal** built
from the *true* binomial Fisher information W_ij = n_ij p_ij (1 − p_ij) (which is
appropriately *wide* where the data are sparse), accepted against the exact
binomial likelihood. Proposing in the constraint null-space coordinates (§7)
makes the move respect the sum-to-zero constraints without a Jacobian
correction. This step is what converges the rare-event cells, in a few thousand
iterations rather than never.

### 2.3 The sweep, step by step

Each iteration performs, in this fixed order:

1. **Augment.** Compute η, draw ω_ij ~ PG(n_ij, η_ij) for every cell.
2. **Joint Gibbs draw** of β from N(Q⁻¹ b, Q⁻¹) subject to the linear constraints
   Aβ = 0, by *conditioning by Kriging* (Rue and Held, 2005, Alg. 2.6): draw the
   unconstrained β, then correct it by Q⁻¹Aᵀ(AQ⁻¹Aᵀ)⁻¹Aβ.
3. **Precision updates.** Draw each κ from its conjugate Gamma full conditional.
4. **ASIS interweaving** (Yu and Meng, 2011). Re-draw each precision in the
   *non-centred* parameterisation (rescaling the effect and its precision
   together), which breaks the precision–effect coupling that otherwise slows
   mixing. This is a Metropolis step on log κ.
5. **Laplace–Newton Metropolis–Hastings** refinement (§2.2), proposing the whole
   of β jointly in the constraint-free coordinates.
6. **Accumulate** the fitted predictor and deviance for the stored draws.

There is no Metropolis tuning, no acceptance-rate-triggered restart, and no chain
pruning: the sampler is parameter-free from the user's point of view.

### 2.4 Identifiability handling inside the sampler

Each effect is constrained to sum to zero (a level constraint), imposed jointly
in step 2. When the **period** effect is a second-order random walk in a *full*
APC model (i.e. all three effects present and `period = "rw2"`), the shared
linear-trend (drift) direction is improper — the RW2 prior does not penalise it
and the likelihood cannot separate it (see §7) — so the sampler pins it with a
single zero-slope constraint on the period effect. The pin is applied only in
that case: a full model with an RW1 period (for example the headline
`age = "rw2", period = "rw1", cohort = "rw1"`) gets no such constraint, because
its RW1 period prior already (weakly) identifies the trend. This is an internal
reporting convention; the fitted rates, the curvatures and the net drift are
unaffected.

### 2.5 What the user sees

```r
m <- bamp(cases, population, age = "rw2", period = "rw1", cohort = "rw1",
          periods_per_agegroup = 5)            # method = "pg" by default
```

The returned object, the stored samples, `print`, `effects`, `plot`,
`predict_apc` and the DIC are all exactly as before — the new engine is a drop-in
replacement. To use the legacy sampler, pass `method = "iwls"`.

### 2.6 Evidence

On the package example and on simulated RW1/RW2 data the new engine converges
(maximum fitted-value Gelman–Rubin ≈ 1.01, zero restarts, zero pruned chains)
where the legacy sampler did not, and it recovers the simulated effects with the
same accuracy. On RW1 models, where the legacy sampler does converge, the two
agree on the fitted values and the DIC to Monte-Carlo error. On real England
mortality data (six cause × sex series, including zero-heavy coronary and stroke
data) the new engine converged on all of them, whereas the legacy sampler failed
or pruned all chains on the majority.

---

## 3. The compiled engine (`pg_engine`)

The sampler exists in two implementations selected by `pg_engine`:

* `pg_engine = "R"` — a readable reference implementation in `R/pg_engine.R`.
* `pg_engine = "C"` — a compiled port of the inner loop in `src/pg_engine.c`
  (the default).

**Design.** The C engine implements the §2.3 sweep using only the R C API and the
LAPACK/BLAS routines the package already links (`dpotrf` for the Cholesky,
`dtrsv`/`dposv`/`dgemv`/`dgemm` for the solves and projections). There is **no new
dependency** — in particular no Rcpp and no RcppEigen. (A proof-of-concept
established that a sparse Cholesky, the usual reason to bring in RcppEigen, is
actually *slower* than dense LAPACK at the problem sizes here, P of order 100–350,
because the dense factorisation is small and already optimal.) The R driver
performs all the one-time setup — building and scaling the structure matrices,
the constraint matrix A and its null-space basis, the empirical starting values,
and the index maps — and passes them to the C routine, which returns the same
per-chain object the R engine returns. The R engine is retained as the
correctness oracle.

**Reproducibility.** The C code consumes R's random-number stream in exactly the
same order and count as the R engine (the draw order is documented at the top of
`src/pg_engine.c`: the per-cell PG weights in column-major order, then the joint
draw's P normals, then the overdispersion and precision draws, then the ASIS and
Metropolis uniforms). Consequently `set.seed(s)` followed by either engine
produces the same chain. Verified to about 1e-10 — the residual is
last-bit floating-point reassociation between the BLAS kernels and R's own
linear algebra, propagated through the (chaotic) chain, far inside any
statistically meaningful tolerance — across every supported model
(plain RW1/RW2, overdispersion, each heterogeneity effect, period and cohort
covariates, and their combinations), end to end through `bamp()` on real data.

**Speed.** Roughly 1.8× faster on the coronary/stroke series (P = 148) and up to
about 3.2× on the smaller "non-modelled deaths" series (P = 112), measured as
full four-chain `bamp()` wall-clock at a fixed seed. The speedup grows as the
model shrinks, because the dense Cholesky — the one cost the C code cannot beat,
since R already runs it through LAPACK — is then a smaller fraction of the work.

```r
m_fast <- bamp(cases, population, age = "rw2", period = "rw1", cohort = "rw1",
               periods_per_agegroup = 5)                    # pg_engine = "C"
m_ref  <- bamp(cases, population, age = "rw2", period = "rw1", cohort = "rw1",
               periods_per_agegroup = 5, pg_engine = "R")   # identical draws
```

---

## 4. Overdispersion and heterogeneity in the new engine

### 4.1 Overdispersion (`overdisp = TRUE`)

The model adds an independent per-cell effect on the logit scale,

  η_ij = μ + θ_i + φ_j + ψ_k + δ_ij,  δ_ij ~ N(0, 1/ζ),  ζ ~ Gamma,

which absorbs extra-binomial variation. Under PG augmentation both δ and ζ have
closed-form conditionals (Gaussian for δ, Gamma for ζ), so no Metropolis tuning is
needed.

The subtlety is that δ_ij and the smooth effects compete to explain the same cell
residual; updating them in alternation mixes very slowly (a structured /
unstructured confounding familiar from BYM models). The engine therefore uses a
**collapsed Gibbs** step: it *marginalises* δ out of the smooth-block update.
Given ω and ζ, integrating δ_ij gives the smooth predictor a per-cell working
precision ω_ij·ζ / (ω_ij + ζ) and a correspondingly reweighted response; the
smooth block is drawn from this marginal, and δ | (smooth, rest) is then drawn
from its Gaussian conditional. This is an exact joint draw of (smooth, δ) and
removes the confounding. Only the scalar precision ζ is stored (as
`samples$overdispersion`); the cell effects are re-sampled at prediction time,
matching the legacy output contract.

*Evidence.* On data simulated with a known ζ the engine recovers it; on simulated
and real data it agrees with the legacy engine on ζ, the DIC and the smooth
effect curves to Monte-Carlo error; and on the strongly overdispersed real
series (where adding overdispersion drops the DIC by hundreds to thousands) it
converges, whereas the legacy engine with overdispersion often does not.

### 4.2 Heterogeneity (`"rw1+het"`, `"rw2+het"`)

A heterogeneity model adds, to a chosen effect, an i.i.d. Gaussian component
sharing that effect's index — e.g. for age, θ*_i = θ_i + θ2_i with
θ2_i ~ N(0, 1/κ2). Because the smooth and i.i.d. components enter the likelihood
through the *same* index, a separate-block update would confound them (as with
overdispersion). The engine therefore draws the i.i.d. component **jointly with
the smooth effects in the one-block Gaussian step**: the precision matrix Q is
extended so that the smooth and heterogeneity indices of an effect share the
likelihood coupling and differ only in their prior block (the GMRF precision κ·K
for the smooth part, the identity precision κ2·I for the i.i.d. part). The
heterogeneity precisions get conjugate Gamma conditionals.

*Display.* `effects()` and `plot()` gain a `combined` argument. By default
(`combined = FALSE`) they show the smooth component θ_i, as before; with
`combined = TRUE` they return the full effect θ_i + θ2_i.

*Evidence.* On data simulated from a known heterogeneity precision the posterior
covers the truth and the chain converges; the engine agrees with the legacy
engine on the *combined* effect curve and on the DIC. (The split of the combined
effect into smooth and i.i.d. parts can differ from the legacy engine, which
mixes that direction poorly with separate blocks; the joint draw used here is the
one validated against ground truth.)

### 4.3 Composition

Overdispersion and heterogeneity compose: δ stays marginalised out of the smooth
block while the heterogeneity component is carried inside it, so a model with
both is drawn correctly.

---

## 5. Covariates in the new engine

A period covariate `period_covariate` (a known, positive vector x of length J;
cohort covariates are symmetric with a length-K vector) enters the model by
**scaling that effect's contribution** to the linear predictor:

  η_ij = μ + θ_i + φ_j · x_j + ψ_k.

So φ_j is a *relative* coefficient and φ_j·x_j is the *absolute* contribution of
period j. The covariate enters the joint Gaussian draw as a **design-column
scaling**: in Q the period block's likelihood couplings pick up one factor of x
per scaled effect (the cross-blocks ×x, the diagonal ×x²); the response vector b
picks up one factor of x; and the Laplace–Newton gradient and Hessian scale to
match (∂η_ij/∂φ_j = x_j). The random-walk prior, the precision draws and the
sum-to-zero / zero-slope constraints all act on the *relative* coefficient φ.

**Storage convention (important for interpreting output).** The stored effect
(`samples$period`) is the *absolute* contribution φ_j·x_j, because
`predict_apc` adds it directly to the linear predictor and `plot` divides the
stored effect by the raw covariate to display the relative coefficient ("raw
period covariate effect"). The covariate is internally normalised to mean 1,
which does not change the absolute contribution and hence is output-invariant.

**Scope.** A per-effect covariate cannot co-occur with heterogeneity on the *same*
effect (the model has no term for it), so a covariate only ever scales the smooth
block. The display gauge of §7 is **not** applied to covariate models (the
linear-trend transform assumes an additive effect), exactly as for non-full-APC
models.

**Validation oracle.** The legacy IWLS covariate path is unreliable, so it is
*not* used to validate covariates. Instead the engine is checked against
`apcSimulate` ground truth: with a known covariate it recovers the relative
effect (the fitted log-odds match to about 0.01 RMSE; the relative coefficient
to correlation ≈ 0.997 for RW1), and with a constant covariate x ≡ 1 it
reproduces the no-covariate fit bit for bit (a regression guard that the
covariate code path does not perturb the base model).

---

## 6. Usability tools

### 6.1 `selectModel()` — automatic model selection

A new exported function answers the practical question *which model specification
do my data support?* — is a first- or second-order random walk more appropriate
for each effect, are the data overdispersed, is extra heterogeneity warranted —
without the user fitting every combination by hand.

**Strategy: greedy forward selection by complexity, on the DIC.** Start from the
simplest model (RW1 for every present effect, no overdispersion). In each round,
fit every candidate that is exactly *one step more complex* than the current best
— an effect upgraded RW1 → RW2, overdispersion switched on, or (optionally)
heterogeneity added to an effect — and adopt the lowest-DIC candidate, **but only
if** it (a) converged and (b) improves the DIC by at least a parsimony margin
(`dic_margin`, default 4, a conventional "clearly better" difference). The
parsimony margin in (b) applies once the current best model has itself converged;
while the running best has *not* yet converged, the first converged candidate is
adopted regardless of the margin (any converged model is preferable to an
unconverged one). Repeat until no admissible improvement remains. This costs a
handful of fits rather than the full grid of ≈ 250 models, follows an
interpretable path, and prefers the simpler model unless the data clearly favour
more structure. Each distinct specification is fitted at most once.

**Two design decisions worth dwelling on.**

* *Convergence gating.* The selection criterion is the DIC, but **a model whose
  chains have not mixed can never be selected**, because a low DIC computed from a
  chain that has not converged is not trustworthy. Convergence is judged by the
  same identified-quantity Gelman–Rubin criterion as `checkConvergence` (§7):
  the maximum fitted-value statistic must be at or below `psrf_tol` (default 1.1).
  This is not academic. On strongly overdispersed but well-populated data the
  overdispersion model both fits far better and converges, so it is selected. On
  zero-heavy data the overdispersion model may have the *lowest* DIC yet fail to
  converge at a moderate run length; `selectModel` then correctly declines it and
  reports the converged binomial model, surfacing the mixing difficulty rather
  than recommending an untrustworthy fit. If *nothing* converges at the screening
  length, the function warns explicitly and falls back to a DIC-only ranking
  flagged as not convergence-backed.

* *Screen short, refit long.* All candidates are fitted with the same moderate
  `screen` MCMC settings for a fast, fair comparison; the selected model is then
  refitted with the longer `final` settings (the data-adaptive default of §6.2)
  before being returned.

```r
sel <- selectModel(cases, population, periods_per_agegroup = 5)
sel                 # ranked comparison table + the chosen specification
plot(sel$model)     # the refitted best model
```

Axes can be pinned to exclude them from the search: `age = "rw2"` fixes the age
effect, `age = " "` removes it, `overdispersion = FALSE` forbids overdispersion;
`try_heterogeneity = TRUE` opens the heterogeneity axis. The result is an
`"apcselect"` object with the ranked `table`, the selected `best` specification,
the refitted `model`, and the adopted `path`.

### 6.2 Data-adaptive MCMC length

Each of `number_of_iterations`, `burn_in` and `step` in `mcmc.options` may now be
a number **or the string `"auto"`** (the new default). With `"auto"`:

* the iteration count is set from a **rarity score** — the fraction of zero cells,
  or half the fraction of cells with ≤ 5 events, whichever is larger — because the
  rare-event, zero-heavy cells are exactly the ones whose PG augmentation mixes
  slowly and so need more iterations. A rarity of 0 (well-populated data) maps to
  40000 iterations, a rarity of 1 (almost every cell empty/rare) to 120000,
  rounded to the nearest thousand;
* `burn_in` defaults to half the iterations;
* `step` is set to keep about 1000 stored samples per chain, so the stored-sample
  count is roughly constant across data sets regardless of the chosen length.

Any element supplied as a number is used exactly as given, so an explicit
`mcmc.options` reproduces the previous behaviour; partial lists work (supply only
`number_of_iterations` and the rest are derived). `verbose = TRUE` reports the
chosen settings.

### 6.3 The `prior_scale` documentation

`prior_scale = TRUE` applies the Sorbye and Rue (2014) scaling: the structure
matrix K is rescaled so that the geometric mean of the (generalised) marginal
variances of the effect equals one. After scaling, 1/√κ is — to a good
approximation — the marginal standard deviation of a typical effect element *on
the logit scale*, **independently of the random-walk order, the number of time
points and the grid spacing**.

Why it matters: with an *unscaled* K the meaning of a fixed κ depends on the order
and size of the random walk (the eigenvalues of K grow with both), so a single
Gamma hyper-prior implies very different prior smoothness across models — a
hyper-prior tuned on one model silently means something else on another. Scaling
gives **portable, interpretable hyper-priors** and **fairer model comparison**
(e.g. RW1 vs RW2 by DIC). The `?bamp` help now contains a dedicated section
explaining this, with a runnable example that tabulates the implied prior effect
standard deviation for a fixed κ across RW orders and grid sizes: it ranges from
about 1.2 to 14.6 without scaling and is 1.0 throughout with scaling.

The default is `prior_scale = FALSE`, which keeps the parameterisation and default
hyper-parameters of the legacy engine; if you turn scaling on you should choose
hyper-parameters appropriate for the scaled prior (where κ ≈ 1/variance).

---

## 7. Identifiability, the convergence check, and the display gauge

The APC model has a well-known non-identifiability (Clayton and Schifflers,
1987): a linear trend can be moved between the age, period and cohort effects
without changing the fitted rates,

  θ_i → θ_i + c·M·i,  φ_j → φ_j − c·j,  ψ_k → ψ_k + c·k,

for any slope c (with a compensating change in the intercept). Only the *full*
APC model (all three effects present) has this aliasing; AP and AC models do not.

**Consequence for convergence diagnostics.** When the raw effect samples are
summarised directly, this non-identified trend drifts between chains and between
runs even when the model has fully converged. The previous `checkConvergence`
applied Gelman's R to the raw effects and therefore reported spurious
non-convergence. In this release **`checkConvergence` assesses the *identified*
quantities** — the smoothing precisions, the intercept, and the fitted log-odds
η_ij in each Lexis cell (these are invariant to the trend re-gauging) — which
removes the false alarm. The raw per-effect diagnostic is still available with
`info = TRUE`.

**Consequence for display.** Because the split of the trend among the three
curves is arbitrary, the plotted effect curves are not reproducible run to run.
`effects()` and `plot()` gain a `convention` argument that fixes this single
degree of freedom *for display only*, applied per MCMC draw:

* `"age"` (default) — remove the linear slope of the age effect (age is shown as
  curvature about a zero trend; the drift is carried by period and cohort);
* `"period"`, `"cohort"` — pin the corresponding effect's slope to zero instead;
* `"none"` — the previous, un-gauged behaviour.

The gauge is display-only: the stored samples, the fitted rates, the `predict_apc`
projections and the DIC are all invariant to it, and it is ignored for models that
are not full APC (and for covariate models, where the multiplicative covariate
breaks the additive transform). Note that because the gauge is applied per draw
and quantiles are non-linear, the *summarised* median curve has only approximately
(not exactly) zero slope.

---

## 8. Bug fixes

* **`predict_apc(object, periods = 0)` crashed.** The retrospective,
  no-forecast case (used for model checking) errored, and with
  `overdisp = TRUE` aborted with "missing values and NaN's not allowed if
  'na.rm' is FALSE" (preceded by many non-fatal "NaNs produced" warnings). The
  cause was a classic R sequence trap: the random-walk extrapolation helper
  looped over `(n1+1):n2`, and when `periods = 0` the future horizon `n2` equals
  the observed `n1`, so `(n1+1):n2` is `(n1+1):n1`, which in R counts *downwards*
  to `c(n1+1, n1)` instead of being empty. This appended a spurious period (and
  cohort) step and shifted every field in the packed per-draw parameter vector,
  so the per-cell predictor read the wrong cohort effect and the
  overdispersion-precision slot landed on a (frequently negative) effect, whence
  √(negative) only *warns* "NaNs produced" — the fatal abort came one step later,
  when the downstream `quantile()` call met those NaNs. The three loops are now
  guarded with `if (n2 > n1)`. The bug was pre-existing and engine-independent
  (it reproduced on `method = "iwls"`); `periods ≥ 1` was
  unaffected.

* **`plot()` failed for some quantile counts.** A single quantile errored with
  "incorrect number of dimensions" (the summary was dropped to a vector), and
  vectors of length 2 or 4 errored with "invalid line type" (the line-type rule
  was defined only for 1, 3 or 5 quantiles). The summary now keeps its matrix
  shape and the line types follow a symmetric distance-from-median rule valid for
  any count (reproducing the historical appearance for 3 and 5 quantiles). The
  covariate panels also now tolerate a zero in the covariate without aborting on
  a non-finite axis range.

* **`cohort = "rw2+het"` cohort effect never updated (legacy engine).** A typo,
  `if (cohort_block == 3 || cohort_block == 3)` in `src/bamp.cc`, should have
  tested `== 4`, so the RW2-plus-heterogeneity cohort block was never sampled
  under the default (no-overdispersion) sampler. Corrected to `== 4`.

* **`method = "pg"` core allocation.** The new engine now honours a numeric
  `parallel` value as the requested number of cores (it was previously capped at
  two, so `parallel = 4` ran four chains on only two cores, roughly doubling the
  wall-clock). The legacy IWLS automatic chain-removal notice is now issued as a
  catchable `warning()` rather than printed to the console, so callers can detect
  it programmatically.

---

## 9. Compatibility notes

* The **default sampler changed** from `"iwls"` to `"pg"`. Code that did not set
  `method` will now use the new engine. To reproduce pre-2.2.0 results exactly,
  pass `method = "iwls"`.
* The **default MCMC length is now data-adaptive** (`"auto"`). To reproduce a
  fixed previous setting, pass an explicit numeric `mcmc.options`.
* `effects()` and `plot()` now **gauge the trend by default** (`convention =
  "age"`); pass `convention = "none"` for the previous un-gauged curves. This is a
  display change only and does not affect fitted rates, predictions or the DIC.
* No new package dependencies were introduced. The compiled engine uses the
  existing LAPACK/BLAS linkage; `RhpcBLASctl` is not required.

---

## References

* Clayton, D. and Schifflers, E. (1987). Models for temporal variation in cancer
  rates. *Statistics in Medicine* 6, 449–481. DOI:10.1002/sim.4780060406.
* Knorr-Held, L. and Rainer, E. (2001). Projections of lung cancer mortality in
  West Germany: a case study in Bayesian prediction. *Biostatistics* 2, 109–129.
  DOI:10.1093/biostatistics/2.1.109.
* Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian inference for
  logistic models using Polya-Gamma latent variables. *JASA* 108, 1339–1349.
* Rue, H. and Held, L. (2005). *Gaussian Markov Random Fields: Theory and
  Applications*. Chapman & Hall/CRC.
* Schmid, V. and Held, L. (2007). Bayesian age-period-cohort modeling and
  prediction — BAMP. *Journal of Statistical Software* 21(8).
  DOI:10.18637/jss.v021.i08.
* Sorbye, S. H. and Rue, H. (2014). Scaling intrinsic Gaussian Markov random
  field priors in spatial modelling. *JASA* (Spatial Statistics).
  DOI:10.1080/01621459.2013.866549.
* Yu, Y. and Meng, X.-L. (2011). To center or not to center: that is not the
  question — an ancillarity-sufficiency interweaving strategy (ASIS).
  *Journal of Computational and Graphical Statistics* 20, 531–570.
