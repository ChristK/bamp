## Batch 4 (part 1) -- trend realism: damped drift (A1), variance shrinkage (A3), order selection (A2).
suppressMessages(library(bamp)); options(mc.cores = 1L)

set.seed(41); I <- 8; J <- 15; pop <- matrix(2e5, J, I)
base <- plogis(-3 + rep(seq(-1, 1, length.out = I), each = J) + cumsum(rnorm(J, -.03, .02)))
tot <- matrix(rbinom(I * J, as.integer(pop), base), J, I)
ca <- round(tot * .5); cb <- round(tot * .3); cc <- tot - ca - cb
mc <- list(iterations = 1200, burn_in = 400, thin = 2)
fco <- suppressMessages(bamp_coherent(list(m = ca, w = cb), list(m = pop, w = pop),
        periods_per_agegroup = 1, age = "rw1", period = "rw2", cohort = "rw1", mcmc = mc))
H <- 12; last <- J + H

## A1/A3: defaults reproduce the free random walk bit-for-bit
set.seed(7); a <- predict_coherent(fco, periods = H)
set.seed(7); b <- predict_coherent(fco, periods = H, damping = 1, var_damping = 1)
expect_identical(a$m$samples$rate, b$m$samples$rate)

## A1: damping reduces the final-horizon point-forecast drift
set.seed(7); pu <- predict_coherent(fco, periods = H, damping = 1.0, var_damping = 1.0)
set.seed(7); pd <- predict_coherent(fco, periods = H, damping = 0.5, var_damping = 1.0)
drift <- function(p) abs(mean(p$m$rate["50%", last, ]) - mean(p$m$rate["50%", J, ]))
expect_true(drift(pd) < drift(pu))

## A3: var_damping narrows the predictive band at long horizon
set.seed(7); pv <- predict_coherent(fco, periods = H, damping = 1.0, var_damping = 0.85)
bw <- function(p) mean(log(pmax(p$m$rate["95%", last, ], 1e-9)) - log(pmax(p$m$rate["5%", last, ], 1e-9)))
expect_true(bw(pv) < bw(pu))

## argument validation
expect_error(predict_coherent(fco, periods = 1, damping = 1.5))
expect_error(predict_coherent(fco, periods = 1, var_damping = 0))

## multicause damped path stays coherent
fmc <- suppressMessages(bamp_multicause(list(a = ca, b = cb, c = cc), pop,
        periods_per_agegroup = 1, order = 1:3, period = "rw2", mcmc = mc))
pm <- predict_multicause(fmc, periods = H, damping = 0.7, var_damping = 0.95)
expect_true(pm$coherence_maxerr < 1e-7)
## predict_apc damping default reproduces and validates
set.seed(3); q0 <- predict_apc(fmc$total, periods = 4)
set.seed(3); q1 <- predict_apc(fmc$total, periods = 4, damping = 1, var_damping = 1)
expect_identical(q0$samples$pr, q1$samples$pr)
expect_error(predict_apc(fmc$total, periods = 1, damping = -0.1))

## A2: select_rw_order returns a ranked table and a best in the candidate set
s <- select_rw_order(ca, pop, periods_per_agegroup = 1, holdout = 3, n_origins = 1,
       period_orders = c("rw1", "rw2"), cohort_orders = "rw1",
       mcmc_bamp = list(number_of_iterations = 1000, burn_in = 300, step = 2, tuning = 120))
expect_equal(nrow(s$table), 2L)
expect_true(all(c("period", "cohort", "energy", "crps") %in% names(s$table)))
expect_true(s$best$period %in% c("rw1", "rw2"))
expect_true(s$table$energy[1] <= s$table$energy[2])   # sorted by score (default energy)
