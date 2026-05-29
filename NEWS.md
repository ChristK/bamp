# bamp 2.2.0
* New Polya-Gamma Gibbs engine, selected with `method = "pg"` in `bamp()`. It
  draws the intercept and the age, period and cohort effects jointly in one
  exact Gibbs step, so it has no Metropolis tuning, never restarts on low
  acceptance and does not prune chains. It is much more robust for RW2 priors
  (where the old engine often failed to converge) and matches the original
  engine's estimates on RW1 models. Optional Sorbye & Rue (2014) scaling of the
  random-walk priors is available via `prior_scale = TRUE`.
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