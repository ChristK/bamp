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

## ---- Phase 1 cause: joint multinomial competing-cause model ----
mc <- suppressMessages(bamp_multicause(
  list(ihd = round(cases * 0.5), stroke = round(cases * 0.3),
       other = cases - round(cases * 0.5) - round(cases * 0.3)), population,
  age = "rw1", period = "rw1", cohort = "rw1", periods_per_agegroup = ppa,
  mcmc = list(iterations = 2000, burn_in = 600, thin = 2)))
expect_inherits(mc, "apc_multicause")
pmc <- predict_multicause(mc, periods = 2, hazard = TRUE, period_length = ppa)
# coherence by construction: cause rates sum to the all-cause total (machine exact)
expect_true(pmc$coherence_maxerr < 1e-9)
# cause-specific hazards sum to the all-cause hazard
sumhz <- Reduce(`+`, lapply(pmc$causes, function(nm) pmc[[nm]]$samples$hazard))
expect_true(max(abs(sumhz - pmc$total$samples$hazard)) < 1e-9)
expect_equal(dim(pmc$cor_omega), c(2L, 2L))

# cross-cause correlation recovery: simulated cause-replacement (negative) is identified
set.seed(3); Ic <- 8; Jc <- 14; Nc <- 40000
Lc <- chol(0.25^2 * matrix(c(1, -0.8, -0.8, 1), 2, 2))
ph <- matrix(0, Jc, 2); for (j in 2:Jc) ph[j, ] <- ph[j - 1, ] + as.numeric(rnorm(2) %*% Lc)
muc <- c(0.2, -0.3); thc <- cbind(seq(-.5, .5, length.out = Ic), seq(.3, -.3, length.out = Ic))
thA <- seq(-1.2, 1.2, length.out = Ic); phiT <- -0.04 * (1:Jc); phiT <- phiT - mean(phiT)
a1 <- a2 <- a3 <- matrix(0, Jc, Ic); popc <- matrix(Nc, Jc, Ic)
for (j in 1:Jc) for (i in 1:Ic) {
  tr <- plogis(log(0.03 / 0.97) + thA[i] + phiT[j])
  p1 <- plogis(muc[1] + thc[i, 1] + ph[j, 1]); p2 <- plogis(muc[2] + thc[i, 2] + ph[j, 2])
  sp <- as.numeric(rmultinom(1, rbinom(1, Nc, tr), c(p1, (1 - p1) * p2, (1 - p1) * (1 - p2))))
  a1[j, i] <- sp[1]; a2[j, i] <- sp[2]; a3[j, i] <- sp[3]
}
fr <- suppressMessages(bamp_multicause(list(a = a1, b = a2, c = a3), popc, periods_per_agegroup = 1,
        age = "rw1", period = "rw1", cohort = "rw1", order = 1:3,   # keep simulated stick-break order
        mcmc = list(iterations = 3500, burn_in = 1000, thin = 2)))
# reference by dimname (robust to ordering); a<->b were simulated with corr -0.8
expect_true(predict_multicause(fr, periods = 0)$cor_omega["a", "b"] < -0.2)  # cause replacement recovered

## ---- Hardening: sampled rho, cohort deviation, cohort coupling, ordering ----
mh <- list(iterations = 800, burn_in = 200, thin = 2)
# coh-rho: iid keeps rho at 0 (backward-compatible RNG path); ar1 samples & stores rho
f_iid <- suppressMessages(bamp_coherent(list(f = ca, m = cb), list(f = pa, m = pb),
   periods_per_agegroup = ppa, deviation = "iid", mcmc = mh))
expect_true(all(f_iid$samples$rho == 0))
f_ar <- suppressMessages(bamp_coherent(list(f = ca, m = cb), list(f = pa, m = pb),
   periods_per_agegroup = ppa, deviation = "ar1", mcmc = mh))
expect_equal(length(f_ar$samples$rho), length(f_ar$samples$mu0))
expect_true(all(f_ar$samples$rho >= 0 & f_ar$samples$rho < 1))
expect_true(is.finite(f_ar$model$rho_accept) && f_ar$model$rho_accept > 0)
# coh-rho determinant correctness: constrained rho log-posterior peaks near truth at long J
ar1p <- bamp:::.ar1_prec
lp <- function(de, lam, r) { J <- length(de); Z <- svd(matrix(1, 1, J), nu = 0, nv = J)$v[, 2:J, drop = FALSE]
  sum(log(diag(chol(crossprod(Z, ar1p(J, r) %*% Z))))) - 0.5 * lam * sum(de * (ar1p(J, r) %*% de)) }
g <- seq(0.02, 0.96, 0.02)
set.seed(2); deH <- as.numeric(arima.sim(list(ar = 0.9), 40, sd = 0.25)); deH <- deH - mean(deH)
expect_true(g[which.max(sapply(g, function(r) lp(deH, 1 / var(diff(deH)), r)))] > 0.6)
set.seed(3); deL <- as.numeric(arima.sim(list(ar = 0.0), 40, sd = 0.25)); deL <- deL - mean(deL)
expect_true(g[which.max(sapply(g, function(r) lp(deL, 1 / var(diff(deL)), r)))] < 0.3)

# coh-cohdev: default off has no dpsi; turning it on runs and stays coherent
expect_true(is.null(f_iid$samples$dpsi))
f_cd <- suppressMessages(bamp_coherent(list(f = ca, m = cb), list(f = pa, m = pb),
   periods_per_agegroup = ppa, deviation_cohort = "iid", mcmc = mh))
expect_true(isTRUE(f_cd$model$use_dpsi) && !is.null(f_cd$samples$dpsi))
pcd <- predict_coherent(f_cd, periods = 2)
npc <- dim(pcd$total$samples$rate)[1]                       # total covers observed periods (pop=NULL)
rf <- pcd$f$samples$rate[seq_len(npc), , ]; rm_ <- pcd$m$samples$rate[seq_len(npc), , ]
expect_true(all(pcd$total$samples$rate >= pmin(rf, rm_) - 1e-9 &
                pcd$total$samples$rate <= pmax(rf, rm_) + 1e-9))

# mc-cohcoupling: cohort cross-cause correlation reported; Cm=1 (C=2) path works
mc3 <- suppressMessages(bamp_multicause(list(x = round(cases * .5), y = round(cases * .3),
   z = cases - round(cases * .5) - round(cases * .3)), population, periods_per_agegroup = ppa, mcmc = mh))
expect_equal(dim(predict_multicause(mc3, periods = 2)$cor_omega_psi), c(2L, 2L))
mc2 <- suppressMessages(bamp_multicause(list(x = round(cases * .5), y = cases - round(cases * .5)),
   population, periods_per_agegroup = ppa, mcmc = mh))
expect_equal(dim(predict_multicause(mc2, periods = 1)$cor_omega), c(1L, 1L))

# mc-order: prevalence reorders for fitting but predict reports original order; access by name correct
mco <- suppressMessages(bamp_multicause(list(rare = round(cases * .1), big = round(cases * .6),
   mid = cases - round(cases * .1) - round(cases * .6)), population, periods_per_agegroup = ppa,
   order = "prevalence", mcmc = mh))
expect_identical(mco$data$causes[1], "big")               # most prevalent fitted first
ppo <- predict_multicause(mco, periods = 1)
expect_identical(ppo$causes, c("rare", "big", "mid"))     # original order reported
expect_false(is.null(ppo$rare$rate))                      # access by name correct

# coh-Sgt2: general S>2 strata; S=2 stays on the legacy path
expect_true(is.null(f_iid$model$S))                       # S=2 legacy: no $S flag
caL <- list(low = round(cases * .5), mid = round(cases * .3),
            high = cases - round(cases * .5) - round(cases * .3))
paL <- list(low = round(population * .5), mid = round(population * .3),
            high = population - round(population * .5) - round(population * .3))
f3 <- suppressMessages(bamp_coherent(caL, paL, periods_per_agegroup = ppa,
                                     deviation = "ar1", rho = 0.6, mcmc = mh))
expect_equal(f3$model$S, 3L)
expect_equal(dim(f3$samples$d), c(length(f3$samples$mu0), 3L, nrow(cases)))
expect_true(max(abs(apply(f3$samples$d, c(1, 3), sum))) < 1e-8)   # per-period sum_s d = 0
fut3 <- function(P) rbind(P, P[rep(nrow(cases), 2), , drop = FALSE])
p3 <- predict_coherent(f3, periods = 2, population = lapply(paL, fut3))
expect_true(all(c("low", "mid", "high", "total") %in% names(p3)))
np3 <- dim(p3$total$samples$rate)[1]
rL <- p3$low$samples$rate[seq_len(np3), , ]; rM <- p3$mid$samples$rate[seq_len(np3), , ]
rH <- p3$high$samples$rate[seq_len(np3), , ]
expect_true(all(p3$total$samples$rate >= pmin(rL, rM, rH) - 1e-9 &
                p3$total$samples$rate <= pmax(rL, rM, rH) + 1e-9))

## ---- cascade: nested disease taxonomy, coherent at every level ----
d1 <- round(cases * .3); d2 <- round(cases * .25); d3 <- round(cases * .25)
d4 <- cases - d1 - d2 - d3; d4[d4 < 0] <- 0
casc <- suppressMessages(bamp_cascade(list(A = c("d1", "d2"), B = c("d3", "d4")),
   list(d1 = d1, d2 = d2, d3 = d3, d4 = d4), population, periods_per_agegroup = ppa, mcmc = mh))
expect_inherits(casc, "apc_cascade")
pcasc <- predict_cascade(casc, periods = 2, hazard = TRUE, period_length = ppa)
expect_identical(sort(pcasc$leaves), c("d1", "d2", "d3", "d4"))
expect_true(pcasc$coherence_maxerr < 1e-9)               # leaves sum to all-cause
gA <- pcasc$d1$samples$rate + pcasc$d2$samples$rate      # within-group: d1+d2 == group A
Dc <- min(dim(gA)[3], dim(pcasc$groups$A$samples$rate)[3])
expect_true(max(abs(gA[, , seq_len(Dc)] - pcasc$groups$A$samples$rate[, , seq_len(Dc)])) < 1e-9)
shz <- Reduce(`+`, lapply(pcasc$leaves, function(l) pcasc[[l]]$samples$hazard))
expect_true(max(abs(shz - pcasc$total$samples$hazard)) < 1e-9)   # leaf hazards sum to all-cause
# a singleton group (leaf == group) works and stays coherent
casc2 <- suppressMessages(bamp_cascade(list(A = c("d1", "d2"), C = "d34"),
   list(d1 = d1, d2 = d2, d34 = d3 + d4), population, periods_per_agegroup = ppa, mcmc = mh))
p2 <- predict_cascade(casc2, periods = 1)
expect_true(p2$coherence_maxerr < 1e-9)
expect_false(is.null(p2$d34$rate))
