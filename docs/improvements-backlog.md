# bamp — improvements backlog (coherent / competing-risk forecasting)

Tick `[x]` the items you want implemented, then tell me (or just say the IDs, e.g. "A1, C1, C4, D1").

**Effort:** S ≈ hours · M ≈ a day or two · L ≈ multi-day / large.
**Value** is for *epidemiological projection quality* unless noted.

---

## Already shipped (context — not for ticking)
- Sex / general S-strata coherence (shared common factor + mean-reverting AR1 deviation, sampled ρ) — `bamp_coherent`
- Competing causes (stick-breaking multinomial-logit APC), cross-cause coupling via Wishart Ω (period + cohort) + optional low-rank factor — `bamp_multicause`
- Disease-taxonomy cascade + full [age × period × sex × group × disease] tensor — `bamp_cascade`, `bamp_sex_cascade`
- Hazard emission for microsim, internal-total reconciliation, energy/variogram scoring, sex-only out-of-sample backtest
- Sparse-Cholesky engine foundation (`engine="sparse"`, multicause)

---

## A. Trend / forecast realism  *(what the projection actually does at the horizon)*
- [x] **A1** Optional damped-drift / mean-reverting RW2 on the shared trend (vs flat RW1 or un-damped linear RW2) — S/M — **high** at 20–30y
- [x] **A2** Per-effect RW-order selection with guidance / score-based auto-choice — S
- [x] **A3** Long-horizon drift-uncertainty shrinkage (stop bands exploding) — S
- [x] **A4** Age-specific period drift (relax the single shared period so ages can diverge) — M
- [x] **A5** Changepoint / structural-break period trend — L
- [ ] **A6** Oldest-old / open-ended top age-group handling — M
- [ ] **A7** Partially-observed / edge-cohort handling (fragile recent & oldest cohorts) — M
- [ ] **A8** Pandemic-shock (COVID 2020–21) outlier / intervention handling — S/M

## B. Drivers / scenarios  *(covariates — the policy-model capability)*
- [ ] **B1** Declared risk-factor covariates linking exposure trends → period/cohort drift — L — **the big substantive gap**
- [ ] **B2** Scenario / counterfactual projection API ("set smoking to X% by 2035") — M *(needs B1)*
- [ ] **B3** Lagged / latency exposure effects — M
- [ ] **B4** Cross-group coupling via *shared* RF covariates (IMPACTncd-relevant) — M
- [ ] **B5** Time-varying covariate effects — M

## C. Validation / calibration  *(is the projection actually right?)*
- [x] **C1** Interval coverage / PIT calibration metric (are the 95% bands 95%?) — S — **high** — ✅ shipped (calibration)
- [x] **C2** Extend the backtest harness to multicause + cascade (currently sex-only) — M
- [x] **C3** Rolling-origin (multiple origins), not a single holdout — S
- [x] **C4** Run the backtest on real England mortality data + report — M — **high**
- [x] **C5** Baseline comparisons (Lee–Carter, naïve, independent APC) — M
- [x] **C6** Marginal CRPS / log-score complements to the multivariate scores — S — ✅ shipped (crps_sample/logs_sample)
- [x] **C7** Per-disease automated calibration report — M

## D. Sampler diagnostics / robustness  *(running 30 diseases × strata)*
- [x] **D1** R-hat / ESS reported in coherent / multicause / cascade (none today) — S — ✅ shipped (bamp_diagnostics)
- [x] **D2** Multi-chain support + parallel chains — M
- [x] **D3** Automated non-convergence flagging across many fits — S
- [x] **D4** Trace / diagnostic plots for the new samplers — S — ✅ shipped (bamp_traceplot)
- [x] **D5** Adaptive tuning of the ρ Metropolis step size — S
- [x] **D6** Rare-disease / zero-cell leaf-share stress test — S

## E. External integration / outputs  *(handoff to IMPACTncd)*
- [x] **E1** Optional reconcile to an *external* official total (ONS / GBD), not just the internal fit — M
- [x] **E2** External population-projection denominators + their uncertainty → counts/burden — M
- [x] **E3** Life-table / life-expectancy emission from the coherent hazards — M
- [x] **E4** Burden helpers (YLL / DALY) — M
- [x] **E5** Documented end-to-end IMPACTncd workflow vignette — S/M

## F. Computational / scaling
- [x] **F1** Native-C port (`cport`) of the coherent + multicause samplers — L
- [x] **F2** Sparse engine for `bamp_coherent` (currently multicause-only) — M
- [x] **F3** Blocked / Kronecker Gibbs for the cross-cause coupling (the *real* scaling win) — L
- [x] **F4** Cached-factor sparse Cholesky (the ~12% + min-Matrix-version pin) — S
- [x] **F5** Parallelise across diseases / strata — S

## G. Statistical extensions
- [x] **G1** Spatial / regional coherence (regions of England as exhaustive strata) — M
- [x] **G2** Small-area / subnational projection — L
- [x] **G3** Joint incidence + mortality + survival modelling — L
- [x] **G4** Mechanical competing-risk hooks (shared-survivor coupling) for the microsim — M

## H. Usability / packaging
- [x] **H1** Plot methods for the coherent / cascade objects — S
- [x] **H2** One-call convergence + validation summary report — S
- [x] **H3** CRAN-clean `R CMD check` (Matrix, vignettes, examples) — M
- [x] **H4** Consolidated user guide tying the pieces together — S

---

### My suggested high-value, low-cost starter set
**C1 + C4** (does it forecast well & are the bands calibrated, on real data) · **A1** (damped drift — biggest single lever on the horizon) · **D1/D3** (you'll need these to run 30 diseases) · then **B1** (the big substantive build) when you're ready for scenarios.
