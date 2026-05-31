# Phase 0 coherent / competing-risk forecasting: hazards, total-plus-share strata,
# and forecast reconciliation. Run with tinytest::test_package("bamp").

options(mc.cores = 1L)          # make predict_apc's mclapply == lapply (reproducible draws)
set.seed(1)

e <- new.env()
data("apc", package = "bamp", envir = e)
cases <- e$cases; population <- e$population

mco <- list(number_of_iterations = 800, burn_in = 200, step = 2, tuning = 100)
ppa <- 5

m <- bamp(cases, population, age = "rw1", period = "rw1", cohort = "rw1",
          periods_per_agegroup = ppa, mcmc.options = mco, parallel = FALSE, verbose = FALSE)

## ---- Deliverable 1: hazard emission ----
p0 <- predict_apc(m, periods = 2)
ph <- predict_apc(m, periods = 2, hazard = TRUE, period_length = ppa)

expect_true(is.null(p0$hazard))                                     # off by default
expect_true(!is.null(ph$hazard))
expect_identical(dim(ph$hazard), dim(ph$pr))                        # hazard shaped like pr
expect_equal(ph$samples$hazard, -log1p(-ph$samples$pr) / ppa, tolerance = 1e-12)  # identity

set.seed(7); h1 <- predict_apc(m, periods = 2, hazard = TRUE, period_length = 1)
set.seed(7); h5 <- predict_apc(m, periods = 2, hazard = TRUE, period_length = ppa)
expect_equal(h1$samples$pr, h5$samples$pr, tolerance = 1e-12)       # same draws under same seed
expect_equal(h1$samples$hazard / ppa, h5$samples$hazard, tolerance = 1e-12)  # h(L) = h(1)/L

# competing-risk additivity: -log(1-p_all) == sum_c -log(1-p_c)
p1 <- 0.03; p2 <- 0.05; p_all <- 1 - (1 - p1) * (1 - p2)
expect_equal(-log1p(-p_all), (-log1p(-p1)) + (-log1p(-p2)), tolerance = 1e-12)

expect_error(predict_apc(m, periods = 1, hazard = TRUE, period_length = 0),
             pattern = "period_length")

## ---- Deliverable 2: total-plus-share strata (exact count coherence) ----
ca <- round(cases * 0.6); cb <- cases - ca
pa <- round(population * 0.55); pb <- population - pa
n1 <- nrow(cases); pers <- 2
fut <- function(P) rbind(P, P[rep(n1, pers), , drop = FALSE])

fit_s <- bamp_strata(list(F = ca, M = cb), list(F = pa, M = pb),
                     age = "rw1", period = "rw1", cohort = "rw1",
                     periods_per_agegroup = ppa, mcmc.options = mco,
                     parallel = FALSE, verbose = FALSE)
ps <- predict_strata(fit_s, periods = pers, population = list(F = fut(pa), M = fut(pb)),
                     hazard = TRUE, period_length = ppa)

expect_equal(ps$coherence_maxerr, 0)                               # internal coherence flag
expect_equal(max(abs(ps$F$samples$cases + ps$M$samples$cases - ps$total$samples$cases)), 0)
expect_equal(ps$F$samples$hazard, -log1p(-ps$F$samples$rate) / ppa, tolerance = 1e-12)

## ---- Gap 1: rare cells (N = 0) handled (no floor) under default method = "pg" ----
caz <- ca; cbz <- cb; paz <- pa; pbz <- pb
caz[1, 1] <- 0L; cbz[1, 1] <- 0L; paz[1, 1] <- 0L; pbz[1, 1] <- 0L   # an empty cell
fit_z <- bamp_strata(list(F = caz, M = cbz), list(F = paz, M = pbz),
                     age = "rw1", period = "rw1", cohort = "rw1",
                     periods_per_agegroup = ppa, mcmc.options = mco,
                     parallel = FALSE, verbose = FALSE)
expect_inherits(fit_z, "apc_strata")
expect_false(any(is.nan(fit_z$shares[[1]]$samples$period[[1]])))

## ---- Gap 3: disaggregation to single year of age / calendar year ----
hd <- disaggregate_hazard(ph$samples$hazard, agegroup_width = 5, period_width = ppa)
expect_equal(dim(hd)[1], dim(ph$samples$hazard)[1] * ppa)          # periods -> years
expect_equal(dim(hd)[2], dim(ph$samples$hazard)[2] * 5)            # groups  -> single ages
expect_equal(dim(hd)[3], dim(ph$samples$hazard)[3])                # draws preserved
# piecewise-constant: the first age group / first period block all equal the source cell
expect_equal(hd[1, 1, 1], ph$samples$hazard[1, 1, 1], tolerance = 1e-12)
expect_equal(hd[ppa, 5, 1], ph$samples$hazard[1, 1, 1], tolerance = 1e-12)
expect_error(disaggregate_hazard(ph$samples$hazard, agegroup_width = 5, period_width = 0))

## ---- Deliverable 3: reconcile cause hazards to all-cause ----
m_c1 <- bamp(round(cases * 0.6), population, age = "rw1", period = "rw1", cohort = "rw1",
             periods_per_agegroup = ppa, mcmc.options = mco, parallel = FALSE, verbose = FALSE)
m_c2 <- bamp(round(cases * 0.55), population, age = "rw1", period = "rw1", cohort = "rw1",
             periods_per_agegroup = ppa, mcmc.options = mco, parallel = FALSE, verbose = FALSE)
pt  <- predict_apc(m,    periods = 2, hazard = TRUE, period_length = ppa)
pc1 <- predict_apc(m_c1, periods = 2, hazard = TRUE, period_length = ppa)
pc2 <- predict_apc(m_c2, periods = 2, hazard = TRUE, period_length = ppa)

rec <- reconcile_apc(pt, list(IHD = pc1, stroke = pc2))
expect_true(rec$coherence_maxerr < 1e-8)                           # mean coherence
expect_equal(rec$causes$IHD$hazard_mean + rec$causes$stroke$hazard_mean,
             rec$total$hazard_mean, tolerance = 1e-8)
expect_true(max(abs(rec$causes$IHD$hazard + rec$causes$stroke$hazard - rec$total$hazard)) < 1e-8)
expect_true(min(rec$causes$IHD$hazard_mean, rec$causes$stroke$hazard_mean) >= 0)
expect_error(reconcile_apc(predict_apc(m, periods = 1), list(a = pc1)),
             pattern = "hazard")                                   # needs hazard=TRUE inputs

## ---- Phase 1: joint sex-coherent model (bamp_coherent / predict_coherent) ----
fut8 <- function(P) rbind(P, P[rep(n1, 8), , drop = FALSE])
fc <- bamp_coherent(list(female = ca, male = cb), list(female = pa, male = pb),
                    age = "rw1", period = "rw1", cohort = "rw1", periods_per_agegroup = ppa,
                    deviation = "iid", mcmc = list(iterations = 3000, burn_in = 1000, thin = 2))
expect_inherits(fc, "apc_coherent")

# in-sample recovery: fitted median tracks the empirical rate
fit_in <- predict_coherent(fc, periods = 0)
expect_true(cor(as.vector(fit_in$female$rate["50%", , ]), as.vector(ca / pa)) > 0.9)

pcoh <- predict_coherent(fc, periods = 8,
                         population = list(female = fut8(pa), male = fut8(pb)),
                         hazard = TRUE, period_length = ppa)

# NON-divergence: the period-component of the sex gap (= 2*delta) has ~flat variance
# across the projection horizon (independent fits would grow several-fold).
gv <- apply(2 * pcoh$deviation$samples, 1, var)
expect_true(gv[n1 + 8] / gv[n1 + 1] < 1.5)

# aggregation coherence: total rate is the population-weighted mean of the two sexes
np <- dim(pcoh$total$samples$rate)[1]
wF <- fut8(pa)[seq_len(np), ]; wM <- fut8(pb)[seq_len(np), ]
wmean <- array(0, dim(pcoh$total$samples$rate))
for (d in seq_len(dim(wmean)[3]))
  wmean[, , d] <- (wF * pcoh$female$samples$rate[seq_len(np), , d] +
                   wM * pcoh$male$samples$rate[seq_len(np), , d]) / (wF + wM)
expect_equal(max(abs(wmean - pcoh$total$samples$rate)), 0, tolerance = 1e-9)

# hazards on the coherent output behave like predict_apc's
expect_equal(pcoh$female$samples$hazard,
             -log1p(-pcoh$female$samples$rate) / ppa, tolerance = 1e-12)

# AR1 deviation variant also runs and is a valid object
fc_ar <- bamp_coherent(list(female = ca, male = cb), list(female = pa, male = pb),
                       age = "rw1", period = "rw1", cohort = "rw1", periods_per_agegroup = ppa,
                       deviation = "ar1", rho = 0.7,
                       mcmc = list(iterations = 1500, burn_in = 500, thin = 2))
expect_inherits(fc_ar, "apc_coherent")
expect_error(bamp_coherent(list(female = ca, male = cb), list(female = pa, male = pb),
                           age = "rw1", period = "rw1", cohort = "rw1",
                           periods_per_agegroup = ppa, deviation = "ar1", rho = 1.2))

## ---- Phase 2: multivariate scoring rules + backtest harness ----
# proper scoring rules against hand-computed values
expect_equal(energy_score(c(0, 0), matrix(c(3, 4), 1, 2)), 5, tolerance = 1e-9)
expect_equal(energy_score(c(1, 0), rbind(c(0, 0), c(2, 0))), 0.5, tolerance = 1e-9)
expect_equal(variogram_score(c(0, 2), rbind(c(0, 0)), p = 1), 4, tolerance = 1e-9)
expect_equal(energy_score(c(2, 5, 1), rbind(c(2, 5, 1), c(2, 5, 1))), 0)        # perfect
expect_error(energy_score(c(0, 0), matrix(1, 1, 3)))                            # dim mismatch

# backtest harness runs and returns a valid comparison (tools, not the science)
bt <- suppressMessages(coherence_backtest(
  list(female = ca, male = cb), list(female = pa, male = pb),
  holdout = 2, periods_per_agegroup = ppa, models = c("coherent", "independent"),
  scale = "rate", mcmc_coherent = list(iterations = 1500, burn_in = 500, thin = 2),
  mcmc_bamp = list(number_of_iterations = 1500, burn_in = 500, step = 2, tuning = 200)))
expect_true(all(c("model", "energy", "variogram", "gap_growth") %in% names(bt)))
expect_equal(nrow(bt), 2L)
expect_true(all(is.finite(bt$energy)) && all(bt$energy >= 0))

# data-driven rho selection returns a value from the grid
sr <- suppressMessages(select_rho(
  list(female = ca, male = cb), list(female = pa, male = pb),
  holdout = 2, periods_per_agegroup = ppa, rho_grid = c(0, 0.6),
  mcmc_coherent = list(iterations = 1200, burn_in = 400, thin = 2)))
expect_true(sr$best_rho %in% c(0, 0.6))
expect_equal(nrow(sr$table), 2L)
