# Graph Report - /home/ckyprid/GH_projects/bamp  (2026-05-29)

## Corpus Check
- Corpus is ~24,845 words - fits in a single context window. You may not need a graph.

## Summary
- 150 nodes · 343 edges · 15 communities (13 shown, 2 thin omitted)
- Extraction: 88% EXTRACTED · 12% INFERRED · 0% AMBIGUOUS · INFERRED: 41 edges (avg confidence: 0.78)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Linear Algebra & Cholesky|Linear Algebra & Cholesky]]
- [[_COMMUNITY_Precision & RNG Sampling|Precision & RNG Sampling]]
- [[_COMMUNITY_MCMC Engine & Fit Diagnostics|MCMC Engine & Fit Diagnostics]]
- [[_COMMUNITY_Init & Utility Helpers|Init & Utility Helpers]]
- [[_COMMUNITY_Covariate Block Updates|Covariate Block Updates]]
- [[_COMMUNITY_APC Model Concepts|APC Model Concepts]]
- [[_COMMUNITY_Taylor-Approx MH Updates|Taylor-Approx MH Updates]]
- [[_COMMUNITY_Core Gibbs Block Update|Core Gibbs Block Update]]
- [[_COMMUNITY_Literature & Model Foundations|Literature & Model Foundations]]
- [[_COMMUNITY_R apc Class & Convergence|R apc Class & Convergence]]
- [[_COMMUNITY_Overdispersion Block Update|Overdispersion Block Update]]
- [[_COMMUNITY_R Prediction & Simulation|R Prediction & Simulation]]
- [[_COMMUNITY_Package Startup Hooks|Package Startup Hooks]]
- [[_COMMUNITY_R Centering Helper|R Centering Helper]]

## God Nodes (most connected - your core abstractions)
1. `bamp()` - 43 edges
2. `coh()` - 23 edges
3. `blocken()` - 19 edges
4. `blocken2()` - 14 edges
5. `blockupdate()` - 11 edges
6. `blockupdate_a2()` - 11 edges
7. `loese()` - 11 edges
8. `loese2()` - 11 edges
9. `blockupdateplus()` - 10 edges
10. `blockupdate_S()` - 10 edges

## Surprising Connections (you probably didn't know these)
- `bamp()` --cites--> `Clayton and Schifflers (1987)`  [EXTRACTED]
  src/bamp.cc → README.md
- `bamp()` --cites--> `Knorr-Held and Rainer (2001)`  [EXTRACTED]
  src/bamp.cc → README.md
- `bamp()` --implements--> `Bayesian Age-Period-Cohort model`  [INFERRED]
  src/bamp.cc → README.md
- `bamp()` --references--> `Hyperparameters / Smoothing Parameters`  [EXTRACTED]
  src/bamp.cc → vignettes/modeling.html
- `apcSimulate (simulate APC data)` --implements--> `Bayesian Age-Period-Cohort model`  [INFERRED]
  R/simulate_apc.R → README.md

## Hyperedges (group relationships)
- **apc object lifecycle: fit, check, summarize, plot, predict** — src_bamp_bamp, check_apc_checkconvergence, print_apc_print_apc, plot_apc_plot_apc, predict_apc_predict_apc, effects_apc_effects_apc, apc_class_apc [INFERRED 0.85]
- **coh cohort-index used across simulation and prediction** — coh_coh, simulate_apc_apcsimulate, predict_apc_ksi_prognose [EXTRACTED 0.90]
- **MCMC Gibbs/MH sweep (mode 0, no overdispersion)** — src_bamp_bamp, src_block_blocken, src_block_blocken2, src_update_my_update_my_mh, src_praez_hyper, src_l_coh [INFERRED 0.85]
- **MCMC sweep with overdispersion (mode 1)** — src_bamp_bamp, src_block_blockupdate, src_zz_ksi_zz_aus_fc_von_ksi0, src_update_my_update_my_1, src_praez_delta_berechnen, src_l_z_aus_ksi_berechnen [INFERRED 0.80]
- **Banded GMRF sampling kernel (Cholesky + triangular solves + Gaussian draw)** — src_mxs_cholesky, src_mxs_loese, src_mxs_loese2, src_block_gausssample, src_block_berechneq [INFERRED 0.85]
- **APC Model Components (Age, Period, Cohort)** — bamp_age_effect, bamp_period_effect, bamp_cohort_effect, bamp_apc_model [EXTRACTED 1.00]
- **Modeling to Prediction to Simulation Workflow** — src_bamp_bamp, predict_apc_predict_apc, simulate_apc_apcsimulate [INFERRED 0.80]
- **MCMC Convergence Assessment** — bamp_mcmc, bamp_gelman_rubin_convergence_diagnostic, check_apc_checkconvergence [EXTRACTED 1.00]

## Communities (15 total, 2 thin omitted)

### Community 0 - "Linear Algebra & Cholesky"
Cohesion: 0.19
Nodes (19): bedinge(), bedinge_lik(), bedinge_lik2(), blockupdate_a2(), lik_bedingt(), ABS(), min(), invers() (+11 more)

### Community 1 - "Precision & RNG Sampling"
Cohesion: 0.17
Nodes (10): delta_berechnen(), delta_berechnen_S(), hyper(), hyper2(), hyper_a(), tau_berechnen(), normal(), nulleins() (+2 more)

### Community 2 - "MCMC Engine & Fit Diagnostics"
Cohesion: 0.14
Nodes (16): APC Data Example (apc dataset), bamp (compiled C MCMC backend), Gelman and Rubin's Convergence Diagnostic, Markov Chain Monte Carlo (MCMC), singlerun (single MCMC chain runner), Schmid (2004, in German), Deviance Information Criterion (DIC), Overdispersion / heterogeneity parameter (+8 more)

### Community 3 - "Init & Utility Helpers"
Cohesion: 0.15
Nodes (6): ksi_berechnen(), logit(), sort(), sortieren(), start(), z_aus_ksi_berechnen()

### Community 4 - "Covariate Block Updates"
Cohesion: 0.20
Nodes (7): berechneBcohortplus(), berechneBplus(), berechneQcohort2(), berechneQcohortplus(), berechneQplus(), blockupdateplus(), gausssample()

### Community 5 - "APC Model Concepts"
Cohesion: 0.29
Nodes (12): Age Effect, Age-Period-Cohort (APC) Model, Cohort Effect, Period Effect, effects.apc (extract APC effects), Cohort Covariate, Hyperparameters / Smoothing Parameters, Period Covariate (+4 more)

### Community 6 - "Taylor-Approx MH Updates"
Cohesion: 0.38
Nodes (12): berechneBtaylor(), berechneBtaylorcohort(), blocken(), loglikelihood(), MausQphi(), MausQpsi(), MausQtheta(), taylor1() (+4 more)

### Community 7 - "Core Gibbs Block Update"
Cohesion: 0.25
Nodes (9): berechneB(), berechneB_S(), berechneBcohort(), berechneBcohort_S(), berechneQ(), berechneQcohort(), berechneQspace(), blockupdate() (+1 more)

### Community 8 - "Literature & Model Foundations"
Cohesion: 0.29
Nodes (7): Berzuini and Clayton (1994), Besag, Green, Higdon and Mengersen (1995), Knorr-Held and Rainer (2001), Bayesian Age-Period-Cohort model, Lexis diagram, Random walk smoothing priors (RW1/RW2), predict_rw (random walk forecast)

### Community 9 - "R apc Class & Convergence"
Cohesion: 0.40
Nodes (6): apc S3 class constructor, checkConvergence (Gelman-Rubin R convergence check), Clayton and Schifflers (1987), APC likelihood non-identifiability, plot.apc (plot APC effects), print.apc (print APC model summary)

### Community 10 - "Overdispersion Block Update"
Cohesion: 0.33
Nodes (6): berechneQ2(), blocken2(), detQ(), loglikelihood2(), loglikelihood2o(), machQ2()

### Community 11 - "R Prediction & Simulation"
Cohesion: 0.60
Nodes (5): coh (cohort index from age and period), pr_aus_ksi_berechnen (logit-inverse rate from ksi), ksi_prognose (predicted linear predictor), predict_apc (APC prediction), apcSimulate (simulate APC data)

## Knowledge Gaps
- **14 isolated node(s):** `zentrieren (centering helper)`, `bamp (compiled C MCMC backend)`, `apc.data (example APC dataset)`, `.onAttach (startup message, onattach.R)`, `plot.apc (plot APC effects)` (+9 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `bamp()` connect `MCMC Engine & Fit Diagnostics` to `Linear Algebra & Cholesky`, `Precision & RNG Sampling`, `Init & Utility Helpers`, `Covariate Block Updates`, `APC Model Concepts`, `Taylor-Approx MH Updates`, `Core Gibbs Block Update`, `Literature & Model Foundations`, `R apc Class & Convergence`, `Overdispersion Block Update`, `R Prediction & Simulation`?**
  _High betweenness centrality (0.541) - this node is a cross-community bridge._
- **Why does `coh()` connect `Taylor-Approx MH Updates` to `Precision & RNG Sampling`, `MCMC Engine & Fit Diagnostics`, `Init & Utility Helpers`, `Covariate Block Updates`, `Core Gibbs Block Update`, `Overdispersion Block Update`?**
  _High betweenness centrality (0.120) - this node is a cross-community bridge._
- **Why does `blocken()` connect `Taylor-Approx MH Updates` to `Linear Algebra & Cholesky`, `Precision & RNG Sampling`, `MCMC Engine & Fit Diagnostics`, `Covariate Block Updates`, `Core Gibbs Block Update`, `Overdispersion Block Update`?**
  _High betweenness centrality (0.081) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `bamp()` (e.g. with `apc.data (example APC dataset)` and `apcSimulate (simulate APC data)`) actually correct?**
  _`bamp()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `zentrieren (centering helper)`, `bamp (compiled C MCMC backend)`, `apc.data (example APC dataset)` to the rest of the system?**
  _14 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `MCMC Engine & Fit Diagnostics` be split into smaller, more focused modules?**
  _Cohesion score 0.13970588235294118 - nodes in this community are weakly interconnected._