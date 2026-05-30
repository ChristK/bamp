# bamp 2.2.0

This is a major release. Its central theme is a new, more robust MCMC engine
that fixes the long-standing convergence failures on the sparse, rare-event data
typical of incidence and mortality modelling, together with tools that make the
package easier to use well. A detailed, self-contained account of every change
(with the statistical reasoning and validation evidence) is in
`CHANGES-2.2.0.md`. The highlights:

## New default sampler: Polya-Gamma + Laplace (`method = "pg"`)

* **A new MCMC engine, `method = "pg"`, is now the default** (previously
  `"iwls"`). Each sweep draws the intercept and the age, period and cohort
  effects *jointly* in one exact Gibbs step using Polya-Gamma data augmentation
  (Polson, Scott & Windle 2013), then refines them with a joint Laplace (Newton)
  Metropolis-Hastings step against the exact binomial likelihood. There is no
  Metropolis tuning, no acceptance-rate restart and no chain pruning. The joint
  Gibbs step fixes the second-order random-walk (RW2) models, where the legacy
  sampler routinely failed to converge; the Newton step gives fast mixing in the
  high-population, rare-event cells that pure data augmentation moves only in
  tiny steps (cells the old sampler could not converge in hundreds of thousands
  of iterations now converge in a few thousand). It is feature-complete --
  RW1/RW2 priors, heterogeneity, overdispersion and period/cohort covariates,
  with no fallback to another engine -- and reproduces the legacy estimates where
  both converge. The original IWLS block Metropolis-Hastings sampler remains
  fully available via `method = "iwls"` and can be faster on well-behaved (non
  rare-event) data. **Existing calls that did not set `method` will now use
  `"pg"`; pass `method = "iwls"` to retain the previous behaviour.**

* **Compiled engine (`pg_engine = "C"`, the default).** The Polya-Gamma sampler
  has both a readable reference implementation in R (`pg_engine = "R"`) and a
  compiled C implementation (`pg_engine = "C"`). The two run the identical
  algorithm and share the random-number stream draw for draw, so for a given
  seed they agree to floating-point tolerance (verified to about 1e-10 end to
  end through `bamp()` across every supported model); the C engine is roughly
  1.8 to 3.2 times faster on real data. The C code uses only the R C API and the
  LAPACK/BLAS the package already links -- there is **no new package
  dependency** (in particular no Rcpp or RcppEigen).

* **Full model coverage of the new engine.** `method = "pg"` supports, natively
  and without falling back:
  - *overdispersion* (`overdisp = TRUE`): the cell-level random effect is
    marginalised out of the joint draw (a collapsed Gibbs step), then drawn from
    its closed-form Gaussian conditional with a Gamma conditional for its
    precision. It reproduces the IWLS overdispersion estimates and converges on
    the strongly overdispersed data where the binomial-only model is misspecified.
  - *heterogeneity* (`"rw1+het"`, `"rw2+het"`): the i.i.d. component of each
    effect is drawn jointly with the smooth effects in the one-block Gaussian
    step. On data simulated from a known heterogeneity precision the engine
    recovers it and converges. `effects()` and `plot()` gain a `combined`
    argument: `FALSE` (default) shows the smooth component, `TRUE` the full
    effect (smooth + heterogeneity).
  - *period and cohort covariates* (`period_covariate`, `cohort_covariate`): a
    covariate scales its effect's contribution to the linear predictor and
    enters the joint draw as a design-column scaling. The stored effect is the
    covariate-scaled contribution, so `plot()` recovers the relative coefficient
    by dividing by the covariate; a constant covariate reproduces the
    no-covariate fit exactly. (Validated against `apcSimulate` ground truth; the
    legacy IWLS covariate path is unreliable and is not used for validation.)

* Optional Sorbye & Rue (2014) unit-variance scaling of the random-walk priors
  via `prior_scale = TRUE`, so a single precision hyper-prior is comparable and
  interpretable across random-walk orders and grid sizes. The `?bamp` help now
  has a dedicated "Scaling the random-walk priors" section explaining the
  rationale, benefits and trade-off, with a short runnable demonstration.

## Usability

* **`selectModel()` (new)**: automatic model selection. It searches over the
  random-walk order of each effect, whether to include overdispersion, and
  optionally heterogeneity, by greedy forward selection on the Deviance
  Information Criterion (DIC), and returns a ranked comparison table plus the
  refitted best model. The search is convergence-gated -- a model whose chains
  have not mixed can never be selected -- and uses a parsimony margin so a more
  complex model is adopted only when clearly better. Axes can be pinned to
  restrict the search (e.g. `age = "rw2"`).

* **Parallel chains on Windows (`method = "pg"`).** The chains of the default
  engine now run in parallel on Windows too. On Unix and macOS the chains are run
  in forked workers (`parallel::mclapply`) as before; on Windows, where forking
  is unavailable, the engine now starts a PSOCK cluster
  (`parallel::makePSOCKcluster` + `parLapply`) instead of falling back to running
  the chains serially. This is reproducible and mechanism-independent: the
  per-chain seeds are drawn in the main process, so a given `set.seed()` produces
  bit-identical chains whether they run forked, on a socket cluster, or serially
  (verified). If a cluster cannot be created the engine falls back to serial
  rather than failing. (The legacy `method = "iwls"` engine still runs serially on
  Windows.)

* **Data-adaptive MCMC length (new default).** Each of `number_of_iterations`,
  `burn_in` and `step` in `mcmc.options` may be a number or the string `"auto"`
  (the new default). With `"auto"`, the iteration count is chosen from the
  rarity of the data -- more iterations for zero-heavy, rare-event data whose
  cells mix more slowly (about 40000 for well-populated data up to 120000 when
  almost every cell is empty) -- `burn_in` defaults to half the iterations, and
  `step` is set to keep about 1000 stored samples per chain. Any element given
  as a number is used exactly as supplied, so an explicit `mcmc.options`
  reproduces the previous behaviour.

* **Interpretable, reproducible effect displays.** `effects()` and `plot()` gain
  a `convention` argument that fixes the non-identified shared linear trend
  (drift) of a full age-period-cohort model for display. In a full APC model the
  three effects are identifiable only up to a common trend, so summarising the
  raw samples lets the curves drift between runs; `convention = "age"` (default)
  removes the age effect's slope (showing age as curvature about a zero trend and
  putting the drift in period and cohort), `"period"`/`"cohort"` pin those
  instead, and `"none"` reproduces the previous behaviour. The gauge is
  display-only and applied per draw: the stored samples, fitted rates, `predict()`
  projections and DIC are all unchanged, and it is ignored for non-full-APC
  models.

* **`checkConvergence()` now assesses the *identified* quantities** (the
  smoothing precisions, the intercept and the fitted log-odds in each Lexis
  cell) instead of the raw age/period/cohort effects. Because the three effects
  share a non-identified linear-trend direction, the raw effects can differ
  between chains even when the model has fully converged, which made the previous
  check report spurious non-convergence. The raw per-effect diagnostic is still
  shown with `info = TRUE`.

## Bug fixes

* `predict_apc(object, periods = 0)` (the retrospective, no-forecast case used
  for model checking) no longer errors. The random-walk extrapolation helper
  looped over `(n1+1):n2`, which counts *downwards* when `n2 == n1`, appending a
  spurious extra step and shifting the packed parameter vector; with
  `overdisp = TRUE` this took the square root of a negative value and aborted
  with "NaNs produced". The loops are now guarded with `n2 > n1`. The bug was
  pre-existing and affected both engines; `periods >= 1` was unaffected.

* `plot()` of an `apc` object now works for any number of plotted `quantiles`: a
  single quantile no longer errors with "incorrect number of dimensions", and
  vectors of length 2 or 4 no longer error with "invalid line type" (line types
  were previously defined only for 1, 3 or 5 quantiles). The established
  appearance for 3 and 5 quantiles is unchanged, and the covariate panels
  tolerate a zero in the covariate without aborting on a non-finite axis range.

* Fixed a typo (`cohort_block == 3 || cohort_block == 3`) that prevented the
  cohort effect from being updated for the `cohort = "rw2+het"` model under the
  legacy (no-overdispersion) sampler.

* `method = "pg"` now honours a numeric `parallel` value as the requested number
  of cores (it was previously capped at two, so `parallel = 4` ran four chains on
  two cores). The IWLS automatic chain-removal notice is now a catchable
  `warning()` rather than a printed message.

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