# Coherent & competing-risk forecasting in `bamp` — design note

Status: design / not yet implemented. Audience: `bamp` maintainer and IMPACTncd lead.
This is a decision document, not code. File:line references are to the current tree.

Two extensions are scoped here:

- **Problem A — stratum coherence.** Project mortality/incidence by a stratifier (sex,
  education, region, …) so the strata (i) do not diverge implausibly and (ii) aggregate
  consistently to the unstratified projection.
- **Problem B — dependent multi-disease / competing risk.** Project several diseases
  without assuming independence, respecting the shared population at risk and coherence
  with all-cause.

Both reduce, in `bamp` terms, to the **same** small set of changes: a shared common
factor + mean-reverting/stationary stratum deviations in the in-sample field, and — the
load-bearing half — a rewrite of the projection step so deviations are extrapolated as
stationary processes rather than free random walks.

---

## 0. The governing principle

In `bamp`, **coherence is won or lost in `predict_apc()` / `predict_rw()`
(`R/predict_apc.R:58-87`), not in the model fit.** The Pólya-Gamma field can fit several
strata or diseases beautifully in-sample and they will still re-diverge in projection,
because each period/cohort effect is extrapolated as a *free* RW1/RW2 carrying only that
draw's own scalar precision (`period` λ at `R/predict_apc.R:103`, `cohort` ν at `:116`).
There is no slot for a stationary mean, an AR coefficient, or a cross-stratum innovation
covariance.

Treat the in-sample prior change as the easy half and the **projection change as the
load-bearing half**. Under RW2 in particular, the terminal slope of each stratum's period
effect is extrapolated indefinitely, so any small in-sample slope difference becomes an
unbounded gap — this is the failure that naïve stratification produces.

### Two distinct meanings of "coherent with the total"

Exhaustiveness buys only the first of these for free:

- **Internal / aggregation coherence** — the total *is* `Σ_s pop_s · rate_s`. Automatic for
  any exhaustive partition (you sum the parts). This is what the shared-factor model gives.
- **External coherence** — the aggregate exactly equals a *separately-fitted* all-population
  model. **Not** automatic, because the logit link is nonlinear: sharing the period/cohort
  effect on the log-odds scale does not make the population-weighted aggregate rate match a
  fitted total on the rate scale. If you need this, use the **total-plus-share / cascade**
  construction (fit the total, model strata as shares of it).

The extra condition beyond "exhaustive" is therefore: **you need correct population weights
(denominators) per stratum, in-sample and projected.**

---

## 1. Problem A — stratum coherence (sex as the canonical case)

### Recommended: shared common factor + mean-reverting deviation (native PG)

This is Riebler & Held's multivariate APC expressed in `bamp`'s existing one-block PG
Gaussian field — the Bayesian-APC form of Li-Lee's augmented common factor and Hyndman's
product-ratio. For sex `s ∈ {F, M}`:

```
logit p_{sij} = mu_s + theta_{s,i} + phi_j + psi_k + d^phi_{s,j} + d^psi_{s,k}
```

- `phi_j, psi_k` — **shared** period/cohort effects; existing intrinsic RW1/RW2 priors,
  unchanged.
- `d^phi_{s,·}, d^psi_{s,·}` — **stratum deviations** with a *proper, mean-reverting* prior
  (AR1 `d_t = ρ d_{t-1} + ε`, `|ρ| < 1`; or, as a first cut, iid ridge = `ρ = 0`),
  sum-to-zero over strata so `d_M = −d_F`.

Coherence is structural: both strata share the same `phi/psi` draws, so the
population-weighted aggregate is internally consistent; the log stratum-ratio (`= 2d` for
two strata) is stationary, so the gap cannot drift even under RW2. It also borrows strength
across strata — the only option that does — which directly helps the rare-event England
cells where the `iwls` method collapses.

### In-sample changes

- **Input** `R/bamp.R:300-318`: accept a 3-D `[S × J × I]` array or a list of matrices;
  flatten to a stacked `Ymat/Nmat`.
- **Index layout** `R/pg_engine.R:124-136`: extend `beta` to
  `(mu_s, theta_s, phi, psi, d^phi_s, d^psi_s)` — `phi/psi` stay a single shared block; add
  per-stratum `theta` blocks and one deviation block per stratified effect.
- **Prior** `R/pg_engine.R:217-260` (`assemble_prec`): `phi/psi` design columns are shared
  (all strata's rows load them); deviation columns load only that stratum's cells. Add a new
  `.pg_ARmat(L, ρ)` tridiagonal prior beside `.pg_Kmat` (`R/pg_engine.R:25`) for the
  deviations.
- **Constraints** `R/pg_engine.R:138-160`: add `Σ_s d_{s,·} = 0` rows; replicate per-effect
  sum-to-zero per stratum for the `theta` blocks; add a stratum-level contrast so `mu_s` and
  `theta_s` do not alias across strata.
- **Precisions** `R/pg_engine.R:367-369`: add a Gamma draw for `λ_d`; if AR1, a small
  griddy-Gibbs/MH draw for `ρ`.

### The decisive `predict_apc` change

- Project shared `phi/psi` with the existing free-RW `predict_rw` (the total is allowed to
  trend; this is the behaviour already trusted).
- Project each deviation with a **new `predict_ar(ρ, λ_d)`**: `d_t = ρ d_{t-1} + N(0, 1/λ_d)`
  — mean-reverts to 0, so the stratum gap stops widening. Carry `ρ` and `λ_d` per draw the
  way `λ`/`ν` are carried at `:103`/`:116` — this requires **adding those slots to the
  per-draw output**, which currently holds a single scalar precision.
- Reconstruct `rate_s = expit(mu_s + theta_{s,i} + phi_j + psi_k + d_{s,j} + d_{s,k})`;
  the aggregate is the population-weighted sum.

### Traps (all verified against the code)

1. **Cost is cubic, not "2-3×".** `assemble_prec` builds a *dense* P×P precision; the
   bottleneck is `chol(Q)` (`R/pg_engine.R:57`), O(P³). Doubling the field is ~6-10× per
   sweep. The C port becomes mandatory and is *not* a mirror of the het-block port (het adds
   a diagonal; this adds correlated cross-blocks).
2. **Silent-bias trap.** The Laplace-MH refinement (`R/pg_engine.R:396-458`) and ASIS step
   (`:375-394`) are hard-wired to the current effect set. If copied unchanged the chain runs
   and looks plausible but the MH ratio is silently wrong → biased posterior, no crash.
   Rederive `state()` and `nc_step` for the new parameterisation.
3. **Do not Sørbye-Rue-scale the AR1 block.** `.pg_scale` (`R/pg_engine.R:32-38`) is for
   *intrinsic* K only; the proper AR1 precision is full-rank.
4. **Constraint completeness.** An incomplete `A` does not error; it surfaces as a
   `1e-6`-ridged near-singular draw that mixes terribly. Validate
   `ncol(Zbasis) = P − rank(A)` on a tiny synthetic set first.
5. **Jallbjørn caution.** Forcing `d` stationary assumes the long-run gap reverts to the
   common trend — for two-sex models this mathematically imposes a unimodal life-expectancy
   gap. Expose `ρ` as a documented dial and choose it from data.

### Lighter-weight fallback: total-plus-share

Fit one ordinary `bamp()` to the summed total (the coherent total), then a second `bamp()`
on the stratum share by passing `cases_tot` as the `population` argument — the share model
is then a binomial-logit fit with `N := cases_tot`. Exact **count** coherence every draw, no
engine change. Correctness condition: pick one formula — if `q` is the share of *events*,
then `cases_s = cases_tot · q` and the stratum **rate** `= cases_s / pop_s`, **never**
`pr_tot · pr_share` (that mixes a share identity with a rate identity and is wrong by
`pop_tot / pop_s` at every non-equal split). Best role: a fast first deliverable and an
oracle to validate the structural model.

---

## 2. Problem B — competing risk / multi-disease

### 2.1 Split the two couplings (this is the whole game)

| Coupling | Nature | Owner | Object exchanged |
|---|---|---|---|
| **Mechanical competing risk** ("who is alive", shared survivors) | within-survivor, hazard/exposure | **IMPACTncd life-table** | per-person-year **additive cause-specific hazards** |
| **Statistical trend coupling** (correlated cause trends, all-cause coherence) | between-trend, prior + projection | **`bamp`** | coherently-projected cause-specific hazards/shares |

Cause-specific *hazards* are separable — a full competing-risks process can be simulated
from them alone. The "shared survivors" coupling enters only when hazards become
risks/counts, through the shared survival function. So the microsimulation already handles
the mechanical coupling correctly **provided `bamp` emits cause-specific hazards, not
independent whole-population probabilities** (which can sum > 1 — the double-counting trap).

Concretely: the current `rbinom`-on-`expit(eta)` path (`R/predict_apc.R:154-155`) is the
literal double-counting form. Add an `emit='hazard'` branch at `R/predict_apc.R:148`
returning `h = -log(1 - p)`, **and** specify the bin-width/exposure correction so emitted
hazards are **per-person-year and additive** across causes (otherwise period-binned hazards,
`periods_per_agegroup > 1`, are not additive at the sim's granularity).

### 2.2 Recommended statistical model: stick-breaking multinomial-logit APC via PG

Keep the all-cause sub-model unchanged (`y_{ij+} ~ Binomial(N_ij, q_ij)`); split causes by
**stick-breaking**: for `c = 1..C-1`, `y_{ijc} ~ Binomial(R_{ijc}, π_{ijc})` where `R` is
the running remainder and each `logit(π^c)` is its own APC field. Cross-cause trend
correlation enters via a separable prior `κ (K ⊗ Ω)` with a Wishart hyperprior on `Ω`. Each
sub-model is binomial-logit + PG + RW-GMRF, so the existing block draw, constraint geometry,
ASIS, and C engine all carry over.

Honest limits: "coherent by construction" is exact only in-sample; coherent *trends* require
carrying the cross-cause innovation covariance through `predict_rw`. Replace the PG
Normal-CLT draw (`.pg_rpg`, `R/pg_engine.R:43-53`) with an exact PG sampler for small
remainders; order causes by prevalence; replicate sum-to-zero **and** RW2 zero-slope
constraints on every cause field and drop deviation intercepts. `K ⊗ Ω` makes Q C²-dense →
O((C·P)³); C port mandatory.

### 2.3 Non-invasive fallback: forecast reconciliation

Run `bamp` once per cause + once for all-cause, then reconcile so cause-specific **hazards**
sum to the all-cause hazard. Zero engine/C changes. Conditions: reconcile on the **hazard**
scale (`emit='hazard'`), not `expit` risks; the fits have independent RNG streams, so
reconcile **marginally** (don't pretend chains share a draw index).

---

## 3. Arbitrary exhaustive strata + nested cascade

Problem A generalises to **any** mutually-exclusive, collectively-exhaustive (MECE)
partition of the population. Sex is just the binary case `S = 2`. Riebler & Held's
multivariate APC was built for arbitrary strata (their examples: 3 regions, 3 countries).
Education, region, deprivation quintile, ethnicity — all fit the same machinery, and
hierarchical stratifiers (e.g. **education within sex**) fit the **nested cascade** form
below, which is the population-side analogue of the disease cascade in §2.

### 3.1 Flat `S`-stratum form

For a single exhaustive partition into `S` strata, replace the binary `d_M = −d_F` with a
**sum-to-zero over the `S` strata**:

```
logit p_{s i j} = mu_s + theta_{s,i} + phi_j + psi_k + d^phi_{s,j} + d^psi_{s,k}
constraint:  Σ_s w_s · d^·_{s,·} = 0          (population-weighted; see 3.3)
prior:       d^·_{s,·} ~ proper mean-reverting GMRF (AR1 / iid ridge), precision λ_d
```

Both target properties survive for any `S`:

1. **Non-divergence** — between-stratum log-ratios are stationary, so no two strata drift
   apart without bound in projection.
2. **Aggregation coherence** — the total is reconstructed as the population-weighted sum of
   strata.

### 3.2 Nested cascade form (e.g. education within sex)

When the stratifier is hierarchical, do **not** flatten to a single `S = 2 × K` partition and
lose the structure. Cascade the deviations so coherence holds at **every** level
simultaneously. With sex `s` and education `e`:

```
logit p_{s e i j} = mu_{s e} + theta_{s e, i} + phi_j + psi_k
                    + d^sex_{s, j}            (level 1: sex deviation)
                    + d^edu_{s e, j}          (level 2: education-within-sex deviation)
                    + d^sex_{s, k} + d^edu_{s e, k}        (cohort analogues)

nested constraints (population-weighted):
   Σ_s  w_s    · d^sex_{s, ·}     = 0                  → sexes cohere to the grand total
   Σ_e  w_{s e} · d^edu_{s e, ·}  = 0   for each s     → education coheres to each sex total
```

This is structurally the GBD/Foreman cascade (sum-to-parent at every level), applied to a
population hierarchy instead of a disease hierarchy. The hazard/rate tensor then sums
correctly on **both** margins by construction of the two nested partitions:
`Σ_e pop · rate_{s e} = pop · rate_s` and `Σ_s pop · rate_s = pop · rate_total`.

Each deviation level has its own precision `λ_d^level` and its own mean-reversion strength
`ρ^level`, projected independently in `predict_apc` with `predict_ar`.

### 3.3 Population-weighted sum-to-zero (necessary for skewed strata)

The constraint must be **population-weighted** (`Σ_e w_{s e} · d^edu = 0`), not unweighted,
where `w_{s e} = pop_{s e} / Σ_e pop_{s e}`. With a roughly balanced stratifier (sex ≈ 50/50)
the weighted and unweighted constraints nearly coincide and the distinction is cosmetic. With
a skewed stratifier (education, deprivation) an unweighted constraint treats a small
high-education stratum and a large low-education one symmetrically, and the implied total
drifts off the population mean. Implementation: scale the constraint rows in `A`
(`R/pg_engine.R:138-160`) by the population weights; `Zbasis` follows from `svd(A)` as now.

Note the weights `w_{s e}` are themselves time-varying and must be projected (§3.4).

### 3.4 Education-specific caveats that do **not** apply to sex

These are the substantive differences — education is not a free drop-in:

1. **Do not default to strong reversion — the Jallbjørn sign flips.** The sex gap has been
   *narrowing* (mean-reversion plausible), but **educational mortality gaps have been
   widening** for decades in most high-income countries. A strong-coherence prior (`ρ` small,
   fast reversion) imposes the opposite of the data. For education you likely want `ρ` near 1
   (weak, slow reversion) or a shared-innovation structure. **Let the data choose `ρ` per
   level** — the validation harness (§5) is mandatory here, not optional.

2. **Moving denominators (educational expansion).** The share of the population in each
   education stratum changes strongly across cohorts/periods, so `pop_{s e}` cannot be treated
   as fixed the way `pop_F/pop_M` is. You need population projections **by education**;
   wrong denominators break aggregation coherence even with a perfect rate model. The
   population-weighted constraints (§3.3) depend on these projected weights.

3. **Age-range validity.** Attained education is undefined for children, so the partition is
   exhaustive only within the adult population (≈ 25+). Handle young ages outside the
   education stratification (e.g. stratify by sex only below the education-onset age).

4. **Selective survival is the microsim's job, not `bamp`'s.** Lower-educated people die
   earlier, so survivors at old ages are selectively higher-educated — the same "who is
   alive" mechanical coupling assigned to the life-table in §2.1. In a microsimulation with
   education as an individual attribute this is handled naturally; an aggregate model would
   bias old-age stratum rates. The architecture already places this on the correct side of
   the boundary.

5. **Parameter/cost blow-up.** `S` strata multiply the field; `chol(Q)` is O(P³). Sex
   doubles; sex × 3-level education (`S = 6`) is 6× the strata. Beyond `S ≈ 3-4`, use a
   **low-rank deviation** — e.g. a shared education-gradient *shape* scaled per sex/age, or a
   CBD-style 1-2 coefficient age deviation — instead of a full per-stratum `theta`, and move
   to the sparse-GMRF solver sooner (§4).

### 3.5 Code mapping for the general / nested case

- **Input** `R/bamp.R:300-318`: accept an N-D array `[S₁ × … × S_L × J × I]` or a nested
  list; record the level structure and the per-cell population weights (including their
  projected trajectories).
- **Index layout** `R/pg_engine.R:124-136`: one shared `phi/psi` block at the root; one
  deviation block per (level, stratified-effect) pair. For the nested cascade, level-2 blocks
  are indexed within their level-1 parent.
- **Constraints** `R/pg_engine.R:138-160`: one population-weighted sum-to-zero row set per
  level (nested within parent for level ≥ 2); replicate per-stratum sum-to-zero for the
  `theta` blocks; stratum-level contrasts so the `mu` aliasing is broken at every level.
  **Validate `ncol(Zbasis) = P − rank(A)`** — the APC drift identification problem recurs
  once per added axis, so the constraint count must be checked, not assumed.
- **Precisions** `R/pg_engine.R:367-369`: a Gamma draw for each `λ_d^level` and a `ρ^level`
  draw (or fixed hyperparameter) per level.
- **Projection** `R/predict_apc.R:58-118`: shared `phi/psi` via free-RW `predict_rw`; each
  deviation level via `predict_ar(ρ^level, λ_d^level)`; reconstruct each leaf rate and
  aggregate up the cascade using the projected population weights.

### 3.6 Worked example — sex × education (3 levels)

`S = 2 sexes × 3 education levels = 6 leaf strata`. Root shares `phi/psi`. Level-1: two sex
deviations (`Σ_s w_s d^sex = 0`). Level-2: within each sex, three education deviations
(`Σ_e w_{s e} d^edu_{s e} = 0`). Projection: free-RW for `phi/psi`; `ρ^sex` moderate (gap
narrowing is plausible), `ρ^edu` near 1 (gaps widening — do not force reversion). Aggregation:
education leaves sum (population-weighted) to each sex total; sexes sum to the grand total;
the grand total matches the unstratified free-RW projection internally. External coherence to
a separately-fitted national total, if required, comes from the total-plus-share/cascade
construction rather than the shared-factor model.

---

## 4. Cross-cutting: layer, do not build a monolithic high-D model

The end-state is a tensor (age × period × cohort × **sex × education × cause × …**), and
IMPACTncd needs a hazard tensor that sums correctly on **all** stratum margins and the cause
margin simultaneously. Whether the recommended per-axis strategies *compose* — stacked
constraint matrix full row-rank, augmentation order, feasibility of the multiply-Kronecker
`chol(Q)` — is unproven and likely forces a **sparse/banded GMRF solver** rewrite
(`R/pg_engine.R:57` and `dpotrf` in `src/pg_engine.c` are dense, no sparse path).

Recommendation: **layer.** Shared `phi/psi` at the root; population stratifiers nested by the
cascade of §3 (binary/stable stratifiers like sex innermost; skewed/trending stratifiers like
education as outer levels); cause stick-breaking layered on a stratum-coherent total; nested
disease groups (CVD > CHD/stroke) as a 2-level cause cascade. This lets each axis ship and be
validated independently, reuses the same `predict_ar` / `emit='hazard'` plumbing, and defers
the hardest linear algebra until evidence (a held-out **joint** score showing real
cross-margin interaction) shows the fully-joint field is needed.

**IMPLEMENTED — the disease cascade** (`bamp_cascade()` / `predict_cascade()`,
`R/coherent_cascade.R`): the 2-level taxonomy cascade as nested layered fits — a group-level
`bamp_multicause` (groups partition all-cause) plus a leaf-level `bamp_multicause` within each
group — combined by `leaf = all-cause × group-share × leaf-share`. Coherent at **every** level by
construction (validated machine-exact ~1e-17: leaves→group→all-cause for rates *and* additive
hazards). Because each fit is small it scales where a flat multinomial cannot: **30 diseases
(5 groups × 6) fit + predict end-to-end in ~60 s** on the dense reference engine. It couples trends
*within* a group; the cross-group, risk-factor-driven coupling (smoking across CVD/cancer/COPD) is
the tree-orthogonal piece — see the cross-cutting note below.

**IMPLEMENTED — latent factor cross-cause coupling** (`bamp_multicause(..., factor = R)`): a LOW-RANK
factor model for the cross-cause period covariance (`Σ = ΛΛ' + Ψ`, `R` latent factors) replacing the
full Wishart `Ω`. The factors are shared, **cross-cutting** latent drivers (risk-factor proxies) that
the cascade's tree cannot express, and they make the coupling identifiable for many causes
(`R×C` loadings, not `C²/2` correlations). Implemented as a one-sweep Bayesian factor analysis of the
period *increments* (`δ̃ = √s_p·DΦ`, so `δ̃'δ̃ = Φ'K_pΦ`), run *standardised* so the priors are
scale-appropriate, then `Σ` rescaled; the field draw and projection are unchanged (they just use the
factor-implied `Ω`). Rotation-invariant (only `Σ` is used), so loadings need no identification.
Validated: recovers a block (within-factor `0.6-0.9` vs cross-factor `~0.13`) correlation structure —
including a factor loading a *cross-group* subset — with far fewer parameters; `predict_multicause`
returns the posterior-mean `loadings`. This is the **latent-factor fallback** for cross-cutting
coupling; the alternative — **declared risk-factor covariates** (IMPACTncd-natural, scenario-able) —
remains a planned extension (you supply RF trajectories; correlation comes from a named, intervenable
cause; mechanically just extra design columns + Gaussian priors on the coefficients).

---

## 5. Validation harness (the biggest gap — budget real time)

Coherence gains are invisible to the wrong metric:

- **Multivariate scoring, not marginal CRPS.** Hold out the last *h* periods, refit, score the
  **joint** forecast with an energy / Dawid-Sebastiani variogram score. Marginal RMSE can
  worsen slightly while the joint score and the gap/share calibration improve — coherence is a
  dependence property.
- **No-coupling baseline arm.** Independent fits with no coherence. Jallbjørn's point is that
  the independent model is sometimes right; without this arm there is nothing to measure the
  gain against.
- **Hold-out-a-stratum.** Drop one stratum (a sex, a rare education × cause cell) entirely;
  check the joint model reconstructs it from shared field + deviation. Directly tests the
  borrowing-strength claim in the rare-event regime where `iwls` collapses.
- **Data-driven coherence strength, per level.** Cross-validate / WAIC / posterior-predictive
  check on the realised gap trajectory to choose each `ρ^level` (especially `ρ^edu`). A tunable
  knob with no selection procedure just relocates the modelling decision.
- **Constraint sanity.** Assert `ncol(Zbasis) = P − rank(A)` on a tiny synthetic dataset.

---

## 6. Phased implementation plan

- **Phase 0 (days, low risk) — IMPLEMENTED** (branch `phase0-coherent-forecasting`). Gives
  IMPACTncd coherent additive hazards now and becomes the validation oracle. Delivered:
  - `predict_apc(..., hazard=TRUE, period_length=)` — returns the cumulative cause-specific
    hazard `-log(1-p)/period_length` as `$hazard` (quantiles) and `$samples$hazard`, additive
    across competing causes. Backward-compatible (off by default). `R/predict_apc.R`.
  - `bamp_strata()` + `predict_strata()` — total-plus-share construction for any exhaustive
    partition; exact count coherence (`Σ_s cases_s = cases_total` every draw) with each stratum
    rate computed as `cases_s / pop_s` (denominator bug avoided). `R/coherent_strata.R`.
  - `reconcile_apc()` — OLS forecast reconciliation of cause-specific **hazards** to the
    all-cause hazard; marginal (per mean and per quantile), non-negative and coherent by
    construction. `R/reconcile_apc.R`.
  - `disaggregate_hazard()` — expand a `[period, agegroup]` (or `[period, agegroup, draw]`)
    hazard to single year of age / single calendar year by piecewise-constant replication, for
    the IMPACTncd handoff. `R/disaggregate_hazard.R`.
  - Rare cells (`cases_total = 0`) are now passed to `bamp` as zero-trial `N = 0` cells under the
    default `method = "pg"` (no information; the smooth prior interpolates), not floored; `iwls`
    still floors with a warning.
  - Tests: `tinytest` suite `inst/tinytest/test_coherent.R` (24 checks, all pass) + runner
    `tests/tinytest.R`; `tinytest` added to `Suggests`. (Project convention: tinytest, never
    testthat.)
  - Vignette: `vignettes/microsimulation-hazards.Rmd` — runnable end-to-end example (per-cause
    hazards → reconcile → disaggregate → competing-risk microsimulation, incl. the
    "halve one cause, others rise via shared survivors" demonstration, plus sex-coherent hazards).
- **Phase 1 (R engine, S = 2 sex) — PROTOTYPE IMPLEMENTED** (branch `phase0-coherent-forecasting`).
  `bamp_coherent()` / `predict_coherent()` (`R/coherent_joint.R`): one joint posterior for two sexes,
  `logit p_{s,i,j} = mu0 + a_s + theta_{s,i} + phi_j + psi_k + sgn(s)·delta_j`, with shared
  period/cohort effects (free-RW projection) and a sex-specific period **deviation** `delta_j` under a
  proper mean-reverting prior (`iid` ridge, or `ar1` with a `rho` coherence dial). Implemented as a
  clean, auditable **dense one-block Polya-Gamma Gibbs** reusing the engine helpers
  (`.pg_rpg`, `.pg_draw_block`, `.pg_Kmat`, `.pg_scale`) and **deliberately omitting the ASIS / Laplace-MH
  refinement** — sidestepping the silent-bias trap rather than re-deriving those hard-wired steps.
  Validated: in-sample fitted rates recover empirical (corr ~0.99); projected sex-gap variance is
  **flat across horizon (×1.02 over 10 periods)** while two independent `bamp` fits diverge (×5.4,
  3.2× wider at h=10); the population-weighted total is coherent by construction. tinytest-covered
  (31/31 incl. Phase 1). **HARDENED** (commits on branch, design+adversarially-verified in
  `docs/hardening-plan.md`): `rho` is now **sampled** (determinant-correct logit-MH); an optional
  **cohort-axis deviation** (`deviation_cohort=`, default off = bit-exact); and **`S > 2`** strata via a
  contr.sum period-deviation general sampler (S=2 legacy path unchanged). **Still deferred:** sparse-GMRF
  solver + native-C port (the two large performance efforts — current sampler is dense O(P³), fine as a
  reference, not production); the Laplace-MH/ASIS re-derivation for the production engine.
- **Phase 1 (R engine, competing CAUSES) — PROTOTYPE IMPLEMENTED** (branch
  `phase0-coherent-forecasting`). `bamp_multicause()` / `predict_multicause()` (`R/coherent_cause.R`):
  the Strategy-0 statistical cause model. The deaths are stick-broken into `C-1` conditional
  binomial-logit APC shares (multinomial-PG); cause-specific age/cohort effects; the cause **period
  trends share a multivariate random walk with cross-cause innovation precision `Omega` (Wishart
  prior, `K_p ⊗ Omega`)**. Cause rates/hazards are `share × total`, so **coherent with all-cause by
  construction** (validated machine-exact, ~1e-19; cause hazards sum to the all-cause hazard).
  Crucially `Omega` admits **negative** cross-cause correlation — *cause replacement* — which the
  Phase 0 `bamp_strata` fallback and independent fits cannot represent: on simulated data with a true
  innovation correlation of −0.8, the model recovers −0.71 (90% CrI [−0.88, −0.45], excludes 0). The
  correlation is carried into projection (shares forecast with `Omega`-correlated innovations). Same
  dense-Gibbs reference design (no ASIS/Laplace-MH). **HARDENED:** the **cohort** trends are now also
  cross-cause coupled (second Wishart `Omega_psi`, `K_c ⊗ Omega_psi`; `cor_omega_psi` reported); and a
  **prevalence ordering** default + `order_sensitivity()` diagnostic (predict returns causes in the
  user's original order). tinytest-covered. **Still deferred:** age cross-cause coupling (off the
  forecast axis); sparse-GMRF + native-C port.
- **Phase 2 (validation) — IMPLEMENTED** (branch `phase0-coherent-forecasting`). The §5 harness:
  `energy_score()` / `variogram_score()` (proper *multivariate* scoring rules, verified against
  analytic values; `R/scoring.R`); `coherence_backtest()` (hold out the last *h* periods, refit each
  model, score the **joint** sex-by-age forecast against held-out truth — coherent vs the no-coupling
  **independent** baseline vs the Phase 0 `totalshare`; reports joint energy/variogram, a *marginal*
  energy for contrast, and a `gap_growth` non-divergence diagnostic; `R/coherence_backtest.R`);
  `select_rho()` (data-driven AR1 coherence strength by held-out score — answers the Jallbjorn "how
  strong should coherence be" question from data). tinytest-covered (41/41 total).
  **Findings:** (a) *no bug* — on the rate scale coherent and independent forecast equally well
  (in-sample/held-out RMSE within ~1e-5); the large *lograte* gap on `data(apc)` was a rare-cell clip
  artifact (rates 0–0.001 with zeros), so `scale = "rate"` is the default. (b) **Accuracy is a tie**
  (joint energy within ~0.4% on controlled data) — exactly the literature's "at least as accurate"
  (Hyndman 2013). (c) The coherent model's value is **non-divergence**: across a held-out block the
  projected sex-gap variance grows ~1.5× vs ~3.9× for independent fits (2.6× slower divergence),
  consistent with the Phase 1 long-horizon result. So coherence buys plausibility/non-divergence, not
  short-horizon error — and the harness measures both honestly. **Still open in §5:** hold-out-a-stratum
  is degenerate for S = 2 (a fully-missing sex has an unidentified level) — meaningful only for the
  multi-stratum / rare-cell regime, deferred with it.
- **Phase 3 (high, dominant cost, gated on need).** Port to C after the R reference validates
  to ~1e-10; scope a sparse GMRF solver here, not as an afterthought.

---

## 7. Code-seam index

- `R/pg_engine.R:25` `.pg_Kmat` — add `.pg_ARmat` (proper mean-reverting block).
- `R/pg_engine.R:32-38` `.pg_scale` — do **not** apply to AR/proper blocks.
- `R/pg_engine.R:57-67` dense `chol(Q)` — will not scale; sparse solver needed for stacked models.
- `R/pg_engine.R:124-136` parameter index layout — shared vs per-level deviation blocks.
- `R/pg_engine.R:138-160` constraints / `Zbasis` — population-weighted, nested sum-to-zero.
- `R/pg_engine.R:152` RW2 drift pin — recurs once per stacked axis.
- `R/pg_engine.R:217-260` `assemble_prec` — block layout, `K ⊗ Ω`.
- `R/pg_engine.R:367-369` precision draws — add `λ_d^level`, `ρ^level`, Wishart `Ω`.
- `R/pg_engine.R:375-394`, `:396-458` ASIS / Laplace-MH — must be rederived.
- `R/predict_apc.R:58-87` `predict_rw` — the decisive rewrite; add `predict_ar`.
- `R/predict_apc.R:148` `emit='hazard'` branch.
- `R/predict_apc.R:154-155` `rbinom`-on-risk — the literal double-counting form.
- `R/bamp.R:300-318` input validation — N-D array / nested list.
- `R/bamp.R:884-893` output assembly — deviation samples, `ρ`, `λ_d`, `Ω`.
- `src/pg_engine.c:538` `dpotrf` — dense, no sparse path.

---

## References

Coherence across strata: Hyndman, Booth & Yasmeen 2013 (product-ratio, *Demography*);
Li & Lee 2005 (augmented common factor, *Demography*); Riebler & Held 2011 (multivariate APC,
*Ann. Appl. Stat.*); Cairns et al. 2011 (Bayesian two-population APC, *ASTIN*);
Chen et al. 2018 (2-tier ACF, sex × country, *Eur. Actuar. J.*);
Bergeron-Boucher et al. 2017 (compositional, *Demographic Research*);
Jallbjørn et al. 2022 (sex-gap unimodality caution, *Forecasting*);
Vanchev et al. 2017 (common-age-effect, *Scand. Actuar. J.*).

Competing risk / multi-cause: Beyersmann et al. 2009 (simulating competing risks,
*Stat. Med.*); Andersen et al. 2012 (competing risks in epi, *IJE*); Putter et al. 2007
(competing risks & multistate tutorial, *Stat. Med.*); Nigri et al. 2026
(Dirichlet-Multinomial-Poisson); Dong et al. 2025 (α-transformation for zeros,
*Ann. Actuar. Sci.*); Han Li et al. 2019 (forecast reconciliation, *IME*);
Arnold-Gaille et al. 2013 / Arnold et al. 2021 (VECM / cointegration, *NAAJ*);
Hong Li et al. 2018 (hierarchical Archimedean copula, *Scand. Actuar. J.*);
Foreman et al. 2017 (Bayesian hierarchical cause forecasts, *JRSS-C*);
Vollset et al. 2024 (GBD 2021 cascading cause forecasts, *The Lancet*).
