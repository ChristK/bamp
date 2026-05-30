# bamp 2.2.0
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
  heterogeneity). Covariate models still fall back to `method = "iwls"`.
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
  Heterogeneity and covariate models still fall back to `method = "iwls"`.
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