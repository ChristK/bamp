## Batch 1 -- measurement backbone: marginal scores, calibration, convergence diagnostics.
suppressMessages(library(bamp)); options(mc.cores = 1L)

## ---- C6: CRPS fair estimator equals the brute-force O(m^2) form ----
set.seed(1); y <- 0.3; x <- rnorm(400, 0, 1)
bf <- mean(abs(x - y)) - sum(abs(outer(x, x, "-"))) / (2 * length(x)^2)
expect_equal(crps_sample(y, matrix(x, ncol = 1)), bf, tolerance = 1e-12)
expect_true(logs_sample(0, matrix(rnorm(500), ncol = 1)) > 0)         # neg log density of N(0,1) at 0 ~ 0.92

## ---- C1: calibration separates a calibrated from an over-confident ensemble ----
mk <- function(sd_fac, K = 500, m = 300) {
  obs <- as.list(rnorm(K)); ens <- lapply(seq_len(K), function(i) matrix(rnorm(m, 0, sd_fac), ncol = 1))
  list(obs, ens)
}
set.seed(7); g <- mk(1.0); cg <- calibration(g[[1]], g[[2]])
set.seed(7); b <- mk(0.5); cb <- calibration(b[[1]], b[[2]])
cov90 <- function(c) c$coverage$empirical[c$coverage$level == 0.9]
expect_true(cov90(cg) > 0.84 && cov90(cg) < 0.95)                     # calibrated ~ nominal 0.90
expect_true(cov90(cb) < 0.72)                                        # over-confident under-covers
expect_true(cg$pit_uniformity$p_value > 0.05)                        # PIT consistent with uniform
expect_true(cb$pit_uniformity$p_value < 0.01)                        # PIT rejects uniform
expect_equal(length(pit_values(c(0, 1), matrix(rnorm(20), 10, 2))), 2L)
expect_true(all(interval_coverage(c(0, 99), matrix(rnorm(20), 10, 2), 0.9) %in% c(0, 1)))

## ---- D1: convergence diagnostics on multicause + coherent fits ----
I <- 8; J <- 14; pop <- matrix(2e5, J, I); set.seed(2)
base <- plogis(-3 + rep(seq(-1, 1, length.out = I), each = J) + rnorm(I * J, 0, .1))
tot <- matrix(rbinom(I * J, as.integer(pop), base), J, I)
ca <- round(tot * .5); cbb <- round(tot * .3); cc <- tot - ca - cbb
mh <- list(iterations = 1200, burn_in = 400, thin = 2)
fm <- suppressMessages(bamp_multicause(list(a = ca, b = cbb, c = cc), pop,
        periods_per_agegroup = 1, order = 1:3, mcmc = mh))
dm <- bamp_diagnostics(fm)
expect_inherits(dm, "bamp_diagnostics")
expect_true(all(c("by_param", "worst", "summary", "ok") %in% names(dm)))
expect_true(all(c("block", "parameter", "rhat", "ess", "flagged") %in% names(dm$by_param)))
expect_true(nrow(dm$by_param) > 20L)
expect_true(is.finite(dm$summary$max_rhat) && dm$summary$max_rhat < 1.3)   # short chain, but not diverging
expect_true(is.logical(dm$ok))
expect_equal(dm$summary$n_flagged, sum(dm$by_param$flagged))

fc <- suppressMessages(bamp_coherent(list(men = ca, women = cbb),
        list(men = pop, women = pop), periods_per_agegroup = 1, mcmc = mh))
dc <- bamp_diagnostics(fc)
expect_true(dc$summary$n_params > 0L)
expect_true(any(grepl("mu0|^a", dc$by_param$parameter)))             # coherent-specific blocks present

## D4: traceplot returns the names it drew (no on-screen device in tests)
pdf(NULL); labs <- bamp_traceplot(fm, n = 2, diag = dm); grDevices::dev.off()
expect_equal(length(labs), 2L)

## error path: not a fitted object
expect_error(bamp_diagnostics(list(a = 1)))
