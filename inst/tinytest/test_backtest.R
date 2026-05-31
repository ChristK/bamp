## Batch 2 -- validation harness: baselines, rolling-origin backtest, reports, rare-cell robustness.
suppressMessages(library(bamp)); options(mc.cores = 1L)

set.seed(21); I <- 7; J <- 15; pop <- matrix(2e5, J, I)
ax <- seq(-1.2, 1.2, length.out = I); kt <- cumsum(rnorm(J, -0.03, 0.02))
base <- plogis(-3 + outer(rep(1, J), ax) + kt)
tot <- matrix(rbinom(I * J, as.integer(pop), base), J, I)
l1 <- round(tot * .35); l2 <- round(tot * .25); l3 <- round(tot * .25); l4 <- tot - l1 - l2 - l3
mc <- list(iterations = 600, burn_in = 200, thin = 2)
mb <- list(number_of_iterations = 600, burn_in = 200, step = 2, tuning = 120)

## C5: baseline forecasters return draws x agegroups per horizon step
fn <- forecast_naive(l1, pop, h = 2, draws = 200); fl <- forecast_leecarter(l1, pop, h = 3, draws = 200)
expect_equal(length(fn), 2L); expect_equal(dim(fn[[1]]), c(200L, I))
expect_equal(length(fl), 3L); expect_equal(dim(fl[[2]]), c(200L, I))
expect_true(all(fl[[1]] > 0 & fl[[1]] < 1))                      # rates in (0,1)

## C2 + C3 + C7: multicause backtest, 1 origin, coupled vs cheap baselines
cases <- list(cvd = l1, cancer = l2, other = tot - l1 - l2)
bt <- multicause_backtest(cases, pop, holdout = 2, n_origins = 1, periods_per_agegroup = 1,
        order = 1:3, models = c("coupled", "naive", "leecarter"), draws_baseline = 200,
        mcmc = mc, mcmc_bamp = mb)
expect_inherits(bt, "bamp_backtest")
expect_true(all(c("coupled", "naive", "leecarter") %in% bt$overall$model))
expect_true(all(is.finite(bt$overall$energy)))
expect_true(all(bt$overall$cov90 >= 0 & bt$overall$cov90 <= 1))
expect_equal(nrow(bt$per_cause), 3L)                            # C7 per-cause breakdown
expect_true(all(bt$per_cause$cause == c("cvd", "cancer", "other")))

## C3: rolling-origin (naive only = no MCMC, fast) yields one per_origin row per origin
bt2 <- multicause_backtest(cases, pop, holdout = 2, n_origins = 3, periods_per_agegroup = 1,
        order = 1:3, models = "naive", draws_baseline = 200)
expect_true(bt2$n_origins >= 2L)
expect_equal(nrow(bt2$per_origin), bt2$n_origins)

## C2 cascade: per-leaf breakdown over a 2-group taxonomy
tax <- list(circ = c("ihd", "stroke"), neo = c("lung", "colo"))
cl <- list(ihd = l1, stroke = l2, lung = l3, colo = l4)
ctb <- cascade_backtest(tax, cl, pop, holdout = 2, n_origins = 1, periods_per_agegroup = 1,
        models = c("cascade", "naive"), draws_baseline = 200, mcmc = mc, mcmc_bamp = mb)
expect_inherits(ctb, "bamp_backtest")
expect_equal(nrow(ctb$per_cause), 4L)                          # per-leaf
expect_true("cascade" %in% ctb$overall$model)

## D3: convergence report across several fits
f1 <- suppressMessages(bamp_multicause(cases, pop, periods_per_agegroup = 1, order = 1:3, mcmc = mc))
f2 <- suppressMessages(bamp_coherent(list(m = l1, w = l2), list(m = pop, w = pop),
        periods_per_agegroup = 1, mcmc = mc))
cr <- convergence_report(list(multicause = f1, coherent = f2))
expect_inherits(cr, "bamp_convergence_report")
expect_equal(nrow(cr$table), 2L)
expect_true(all(c("fit", "max_rhat", "min_ess", "n_flagged", "ok") %in% names(cr$table)))
expect_true(is.logical(cr$all_ok))

## H2: validation report bundles diagnostics + verdict
vr <- validation_report(f1, backtest = bt)
expect_inherits(vr, "bamp_validation")
expect_true(is.logical(vr$verdict))
expect_false(is.null(vr$backtest))

## D6: rare-disease / zero-cell robustness -- no NaN, coherence preserved
set.seed(5); rare <- matrix(rbinom(I * J, as.integer(pop), 1e-5), J, I)
big <- tot - rare; big[big < 0] <- 0
fr <- suppressMessages(bamp_multicause(list(rare = rare, big = big), pop,
        periods_per_agegroup = 1, order = 1:2, mcmc = mc))
prr <- predict_multicause(fr, periods = 1)
expect_true(sum(rare == 0) > 0L)                               # genuinely has zero cells
expect_false(any(!is.finite(prr$rare$samples$rate)))
expect_true(prr$coherence_maxerr < 1e-7)
