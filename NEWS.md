# bamp 2.2.0
* The `prior_scale` argument is now documented in detail: a new
  "Scaling the random-walk priors" section in `?bamp` explains what the
  Sorbye-Rue unit-variance scaling does, why it makes the precision (smoothing)
  hyper-prior portable and interpretable across random-walk orders and grid
  sizes, when to use it, and the trade-off with the default hyper-parameters. A
  short, instantly-running example shows that for a fixed precision the implied
  prior effect standard deviation ranges from ~1.2 to ~14.6 across RW orders and
  grid sizes without scaling, but is 1.0 throughout with scaling.
* The MCMC length is now chosen automatically from the data by default. Each of
  `number_of_iterations`, `burn_in` and `step` in `mcmc.options` may be a number
  or the string `"auto"` (the new default). With `"auto"`, rare or zero-heavy
  data -- whose rare-event cells mix more slowly -- is given more iterations
  (from 40000 for well-populated data up to 120000 when almost every cell is
  empty or has very few events), `burn_in` is set to half the iterations, and
  `step` is set to keep about 1000 stored samples per chain. The rarity score is
  the fraction of zero cells (or half the fraction of cells with <= 5 events,
  whichever is larger). Any element supplied as a number is used exactly as
  given, so passing an explicit `mcmc.options` reproduces the previous behaviour.
* **`method = "pg"` is now the default sampler** (previously `"iwls"`). The
  Polya-Gamma sampler is feature-complete (RW1/RW2, heterogeneity, overdispersion,
  period/cohort covariates), markedly more robust on the sparse, zero-cell,
  rare-event data typical of incidence/mortality modelling -- where the legacy
  IWLS sampler can fail to converge or prune all of its chains -- and, with the
  compiled engine below, comparable in speed. The original IWLS block
  Metropolis-Hastings sampler remains fully available via `method = "iwls"` and
  can still be the faster choice on well-behaved (non rare-event) data. Existing
  calls that did not set `method` will now use `"pg"`; pass `method = "iwls"`
  explicitly to keep the previous behaviour.
* The `method = "pg"` sampler now has a compiled C implementation of its inner
  loop, selected by the new `pg_engine` argument (`"C"`, the default, or `"R"`).
  Both engines run the identical algorithm and, for a given seed, produce the
  same draws to floating-point tolerance (verified to ~1e-10 across the full
  feature matrix -- plain RW1/RW2, overdispersion, heterogeneity, period/cohort
  covariates and their combinations -- end to end through `bamp()`); the C engine
  is roughly twice as fast on real incidence/mortality data. It uses only the
  R C API and the LAPACK/BLAS the package already links -- no new dependency
  (no Rcpp/RcppEigen). The readable `pg_engine = "R"` reference implementation is
  retained as the verification oracle. The C and R engines share the random-number
  stream draw-for-draw, so results are reproducible and comparable across engines.
* `predict_apc(object, periods = 0)` (the retrospective, no-forecast case used for
  model checking) no longer errors. The random-walk extrapolation helper looped
  over `(n1+1):n2`, which counts *downward* when `n2 == n1` (i.e. `periods = 0`),
  appending a spurious extra period/cohort step; that shifted the packed parameter
  vector so the per-cell predictor read the wrong cohort effect, and with
  `overdisp = TRUE` it took the square root of a (negative) effect in place of the
  overdispersion precision, aborting with "NaNs produced". The loops are now
  guarded with `n2 > n1`. Affects both engines (`"iwls"` and `"pg"`); `periods >= 1`
  was unaffected.
* The `method = "pg"` engine now supports period and cohort covariates
  (`period_covariate` / `cohort_covariate`) natively, completing its model
  coverage -- it now handles the full `bamp` model space with no fallback to
  `method = "iwls"`. A covariate scales its effect's contribution to the linear
  predictor (`phi_j * x_j`); it enters the joint Gaussian draw as a column
  scaling of the design (the affected block's likelihood couplings scale by the
  covariate, its diagonal by the covariate squared). The stored effect is the
  absolute (scaled) contribution, so `plot()` recovers the relative coefficient
  by dividing by the covariate. On data simulated from `apcSimulate` with a known
  covariate the engine recovers the relative effect; with a constant covariate it
  reproduces the no-covariate fit exactly. (The legacy `method = "iwls"`
  covariate path is unreliable, so `apcSimulate` -- not `iwls` -- is the
  validation oracle; `"pg"` is the supported covariate engine.) The display gauge
  `convention` in `effects()` / `plot()` is not applied to covariate models (the
  linear-trend transform assumes an additive effect), as for non-full-APC models.
* The `method = "pg"` engine now supports heterogeneity models (`"rw1+het"`,
  `"rw2+het"`) natively instead of falling back to `method = "iwls"`. The iid
  heterogeneity component of each effect is drawn jointly with the intercept and
  the smooth age/period/cohort effects in the one-block Gaussian step (a separate
  block would mix badly because the smooth and heterogeneity components share the
  same index). On data simulated from a known heterogeneity precision the engine
  recovers it (the posterior covers the truth) and converges (Gelman-R ~1.00);
  it agrees with the iwls engine on the combined effect curve and the DIC. By
  default `effects()` and `plot()` still show the smooth component only (as
  before); the new `combined = TRUE` argument returns the full effect (smooth +
  heterogeneity).
* The `method = "pg"` engine now supports overdispersion (`overdisp = TRUE`)
  natively instead of falling back to `method = "iwls"`. The cell-level random
  effect is marginalised out of the joint draw of the intercept and the age,
  period and cohort effects (a collapsed Gibbs step, so the smooth effects and
  the cell effects do not have to be untangled by slow alternating updates),
  then drawn from its closed-form Gaussian conditional, with a Gamma
  full-conditional for its precision; there is still no Metropolis tuning. It
  reproduces the iwls overdispersion estimates (the overdispersion precision,
  the DIC and the smooth effect curves all agree to within Monte-Carlo error on
  simulated and real data) and converges on the strongly overdispersed
  incidence/mortality data where the binomial-only model is badly misspecified.
* `plot()` of an apc object now works for any number of plotted `quantiles`:
  a single quantile no longer errors with "incorrect number of dimensions",
  and vectors of length 2 or 4 no longer error with "invalid line type" (the
  line types were previously only defined for 1, 3 or 5 quantiles). The
  established appearance for 3 and 5 quantiles is unchanged. The covariate
  panels also tolerate a zero in the covariate without aborting on a
  non-finite axis range.
* `plot()` and `effects()` now take a `convention` argument that fixes the
  non-identified linear trend (drift) of a full age-period-cohort model for
  display. The three effects are identifiable only up to a shared linear trend,
  so summarising the raw samples lets the curves drift between runs. The default
  `convention = "age"` removes the linear slope of the age effect (showing age
  as curvature about a zero trend and putting the drift in the period and cohort
  effects), which removes the dominant source of run-to-run non-reproducibility;
  `"period"` and `"cohort"` pin the corresponding effect's slope instead, and
  `convention = "none"` reproduces the previous (un-gauged) behaviour. The gauge
  is display-only and applied per MCMC draw: the stored samples, the fitted
  rates, the `predict()` projections and the DIC are all unchanged, and it is
  ignored for models that are not full APC.
* New `method = "pg"` engine: a joint sampler that combines Polya-Gamma data
  augmentation with a Laplace (Newton) Metropolis-Hastings refinement. Each
  sweep draws the intercept and the age, period and cohort effects jointly in
  one exact Gibbs step (so there is no Metropolis tuning, no acceptance-rate
  restart and no chain pruning), then a joint Newton proposal against the true
  binomial likelihood refines them. The Gibbs step gives robustness (it fixes
  the RW2 models, where the old engine routinely failed to converge); the
  Newton step gives fast mixing in the cells that pure Polya-Gamma moves only in
  tiny steps -- in particular the high-population, rare-event cells typical of
  incidence/mortality data, which the Gibbs-only sampler could not converge even
  in hundreds of thousands of iterations but the hybrid converges in a few
  thousand. It matches the original engine's estimates on RW1 models. Optional
  Sorbye & Rue (2014) scaling of the random-walk priors is available via
  `prior_scale = TRUE`. The engine supports plain RW1/RW2 models; heterogeneity,
  overdispersion and covariate models fall back to `method = "iwls"`.
* `checkConvergence()` now assesses the *identified* quantities (smoothing
  precisions, intercept and the fitted log-odds in each Lexis cell) instead of
  the raw age/period/cohort effects. Because the three effects share a
  non-identified linear-trend (drift) direction, the raw effects can drift
  apart between chains even when the model has fully converged, which made the
  previous check report spurious non-convergence. The raw per-effect diagnostic
  is still shown with `info = TRUE`.
* Fixed a typo that prevented the cohort effect from being updated for the
  `cohort = "rw2+het"` model under the default (no-overdispersion) sampler.

# bamp 2.1.3
* invalid UTF-8 in comment removed

# bamp 2.1.2
* Adapted to R 4.2

# bamp 2.1.1
* USE_FC_LEN_T for Fortran code

# bamp 2.1.0
* Better default settings (burn in times ten, more informative prior for age).
* Add warnings for failed convergence checks and removed chains in bamp(), including suggestions.
* Add warnings for failed convergence checks in print.apc().
* Fixed unwanted doubling of MCMC when verbose=2.

# bamp 2.0.8
* Minor bug fixes, fix "additional issues".

# bamp 2.0.7
* Better initial setting for restarting iterations; helps with RW2 priors.

# bamp 2.0.6
* Introductory vignette renamed (double vignette name warning from CRAN).

# bamp 2.0.5
* Removed ambiguities (mail Brian Ripley) and clean up in C code.

# bamp 2.0.4
* Added examples to all functions.

# bamp 2.0.3
* Added more details to help pages.

# bamp 2.0.2
* Reference in description changed.

# bamp 2.0.1
* Smaller vignettes.

# bamp 2.0.0
* R package.