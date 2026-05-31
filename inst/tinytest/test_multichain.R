## Batch 3 -- sampler robustness: multi-chain between-chain Rhat (D2), adaptive rho (D5).
suppressMessages(library(bamp)); options(mc.cores = 1L)

set.seed(31); I <- 8; J <- 13; pop <- matrix(2e5, J, I)
base <- plogis(-3 + rep(seq(-1, 1, length.out = I), each = J) + cumsum(rnorm(J, -.02, .02)))
tot <- matrix(rbinom(I * J, as.integer(pop), base), J, I)
ca <- round(tot * .5); cb <- round(tot * .3); cc <- tot - ca - cb
mc <- list(iterations = 1200, burn_in = 400, thin = 2)

## D2: run several chains (serial in tests), between-chain Rhat, pooled prediction
mch <- run_chains(bamp_multicause,
        args = list(cases = list(a = ca, b = cb, c = cc), population = pop,
                    periods_per_agegroup = 1, order = 1:3, mcmc = mc),
        chains = 3, parallel = FALSE)
expect_inherits(mch, "bamp_multichain")
expect_equal(mch$n_chains, 3L)
expect_equal(length(mch$chains), 3L)

dm <- bamp_diagnostics(mch)
expect_inherits(dm, "bamp_diagnostics")
expect_true(is.finite(dm$summary$max_rhat) && dm$summary$max_rhat < 1.3)
expect_true(dm$summary$n_params > 20L)

cf <- combine_chains(mch)
nk <- nrow(mch$chains[[1]]$samples$phi)
expect_equal(dim(cf$samples$phi)[1], 3L * nk)             # draws pooled across chains
pp <- predict_multicause(cf, periods = 1)
expect_true(pp$coherence_maxerr < 1e-7)                   # pooled fit still predicts coherently

## seeds default and custom
expect_equal(mch$seeds, 1:3)
mch2 <- run_chains(bamp_coherent,
         args = list(cases = list(m = ca, w = cb), population = list(m = pop, w = pop),
                     periods_per_agegroup = 1, mcmc = mc),
         chains = 2, seeds = c(11, 22), parallel = FALSE)
expect_equal(mch2$seeds, c(11, 22))
expect_true(bamp_diagnostics(mch2)$summary$n_params > 0L)

## D5: adaptive rho tunes the proposal toward ~0.234 acceptance during burn-in
fa <- suppressMessages(bamp_coherent(list(m = ca, w = cb), list(m = pop, w = pop),
        periods_per_agegroup = 1, deviation = "ar1", rho = 0.5,
        mcmc = list(iterations = 3000, burn_in = 1200, thin = 2), adapt_rho = TRUE))
fo <- suppressMessages(bamp_coherent(list(m = ca, w = cb), list(m = pop, w = pop),
        periods_per_agegroup = 1, deviation = "ar1", rho = 0.5,
        mcmc = list(iterations = 3000, burn_in = 1200, thin = 2),
        adapt_rho = FALSE, mh_sd_rho = 0.3))
expect_true(fa$model$mh_sd_rho_adapted)
expect_false(fo$model$mh_sd_rho_adapted)
expect_true(fa$model$rho_accept > 0.10 && fa$model$rho_accept < 0.5)   # healthy acceptance band
## adaptive lands closer to the 0.234 target than the (badly) fixed proposal here
expect_true(abs(fa$model$rho_accept - 0.234) < abs(fo$model$rho_accept - 0.234))
expect_true(length(fa$samples$rho) > 0L)                  # rho draws stored
