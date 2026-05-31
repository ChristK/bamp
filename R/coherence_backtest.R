## ===========================================================================
## Phase 2 validation harness: out-of-sample, MULTIVARIATE-scored comparison of
## coherent vs no-coupling (independent) vs total-plus-share sex forecasts, and
## data-driven selection of the AR1 coherence strength rho.
##
## The point (design note, section 5): coherence is a joint property, so it must
## be scored with a proper MULTIVARIATE rule (energy / variogram) against a
## no-coupling baseline. Marginal scores can even prefer the wrong model.
## ===========================================================================

## per-model ensemble of joint [draws x 2I] sex rate forecasts at the held-out
## periods. Returns a list (length = holdout) of [m x 2I] matrices.
.cb_ens <- function(model, ctr, ptr, popF_full, popM_full, h, J, I, strata,
                    age, period, cohort, ppa, deviation, rho,
                    mcmc_coherent, mcmc_bamp, seed) {
  tidx <- (J - h) + seq_len(h)                       # global period indices of the test block
  joint <- function(F_ID, M_ID) {                    # F_ID/M_ID: [age x draw] at one period
    md <- min(ncol(F_ID), ncol(M_ID))
    cbind(t(F_ID[, seq_len(md), drop = FALSE]), t(M_ID[, seq_len(md), drop = FALSE]))
  }
  if (model == "coherent") {
    fit <- bamp_coherent(stats::setNames(ctr, strata), stats::setNames(ptr, strata),
                         age = age, period = period, cohort = cohort,
                         periods_per_agegroup = ppa, deviation = deviation, rho = rho,
                         mcmc = mcmc_coherent, seed = seed)
    pc <- predict_coherent(fit, periods = h,
                           population = stats::setNames(list(popF_full, popM_full), strata))
    lapply(tidx, function(t) joint(pc[[strata[1]]]$samples$rate[t, , ],
                                   pc[[strata[2]]]$samples$rate[t, , ]))
  } else if (model == "independent") {
    fF <- bamp(ctr[[1]], ptr[[1]], age = age, period = period, cohort = cohort,
               periods_per_agegroup = ppa, mcmc.options = mcmc_bamp, parallel = FALSE, verbose = FALSE)
    fM <- bamp(ctr[[2]], ptr[[2]], age = age, period = period, cohort = cohort,
               periods_per_agegroup = ppa, mcmc.options = mcmc_bamp, parallel = FALSE, verbose = FALSE)
    pF <- predict_apc(fF, periods = h, population = popF_full)
    pM <- predict_apc(fM, periods = h, population = popM_full)
    lapply(tidx, function(t) joint(pF$samples$pr[t, , ], pM$samples$pr[t, , ]))
  } else if (model == "totalshare") {
    fs <- bamp_strata(stats::setNames(ctr, strata), stats::setNames(ptr, strata),
                      age = age, period = period, cohort = cohort,
                      periods_per_agegroup = ppa, mcmc.options = mcmc_bamp,
                      parallel = FALSE, verbose = FALSE)
    ps <- predict_strata(fs, periods = h,
                         population = stats::setNames(list(popF_full, popM_full), strata))
    lapply(tidx, function(t) joint(ps[[strata[1]]]$samples$rate[t, , ],
                                   ps[[strata[2]]]$samples$rate[t, , ]))
  } else stop("unknown model '", model, "'")
}

#' Out-of-sample multivariate backtest of coherent sex forecasts
#'
#' @description
#' Hold out the last \code{holdout} periods, refit each model on the rest, project, and score the
#' \emph{joint} sex-by-age forecast at the held-out periods with proper multivariate rules
#' (\code{\link{energy_score}}, \code{\link{variogram_score}}) against the observed rates. Compares
#' the coherent joint model against the no-coupling \code{independent} baseline (and optionally the
#' Phase 0 \code{totalshare} wrapper). A marginal energy score is also reported to show that marginal
#' scoring is blind to the coherence (dependence) gain.
#'
#' @param cases,population named lists of two \code{[periods x agegroups]} matrices (the sexes).
#' @param holdout number of final periods to hold out and forecast.
#' @param periods_per_agegroup,age,period,cohort model settings (as in \code{\link{bamp}}).
#' @param models which models to compare: any of \code{"coherent"}, \code{"independent"},
#'   \code{"totalshare"}.
#' @param deviation,rho deviation prior for the coherent model (see \code{\link{bamp_coherent}}).
#' @param scale score on the \code{"rate"} or \code{"lograte"} scale (log rates clipped at 1e-6).
#' @param vs_p order of the variogram score.
#' @param mcmc_coherent,mcmc_bamp MCMC settings for \code{\link{bamp_coherent}} and \code{\link{bamp}}.
#' @param seed RNG seed.
#'
#' @return a data.frame with one row per model: mean held-out \code{energy} and \code{variogram}
#'   scores (joint; lower = better), \code{energy_marginal} (sum of per-sex energy scores), and a
#'   non-divergence diagnostic \code{gap_growth} = across-draw variance of the projected sex gap at
#'   the last held-out period divided by the first (\eqn{\approx 1} for the coherent model, \eqn{>1}
#'   for diverging independent fits) with its components \code{gap_var_first}/\code{gap_var_last}.
#' @seealso \code{\link{energy_score}}, \code{\link{variogram_score}}, \code{\link{select_rho}}
#' @export
coherence_backtest <- function(cases, population, holdout = 2, periods_per_agegroup,
                               age = "rw1", period = "rw1", cohort = "rw1",
                               models = c("coherent", "independent", "totalshare"),
                               deviation = "iid", rho = 0,
                               scale = c("rate", "lograte"), vs_p = 0.5,
                               mcmc_coherent = list(iterations = 4000, burn_in = 1000, thin = 2),
                               mcmc_bamp = list(number_of_iterations = 4000, burn_in = 1000,
                                                step = 2, tuning = 300),
                               seed = 1) {
  scale <- match.arg(scale)
  if (is.null(names(cases))) names(cases) <- c("sex1", "sex2")
  names(population) <- names(cases); strata <- names(cases)
  cases <- lapply(cases, as.matrix); population <- lapply(population, as.matrix)
  J <- nrow(cases[[1]]); I <- ncol(cases[[1]]); h <- holdout
  if (h < 1 || h >= J) stop("'holdout' must be between 1 and the number of periods - 1.")

  ctr <- lapply(cases, function(x) x[seq_len(J - h), , drop = FALSE])
  ptr <- lapply(population, function(x) x[seq_len(J - h), , drop = FALSE])
  tf <- if (scale == "lograte") function(x) log(pmax(x, 1e-6)) else identity

  ## held-out observations: joint c(empF, empM) per test period, transformed
  empF <- cases[[1]] / population[[1]]; empM <- cases[[2]] / population[[2]]
  tidx <- (J - h) + seq_len(h)
  obs <- lapply(tidx, function(t) tf(c(empF[t, ], empM[t, ])))

  res <- lapply(models, function(model) {
    ens <- .cb_ens(model, ctr, ptr, population[[1]], population[[2]], h, J, I, strata,
                   age, period, cohort, periods_per_agegroup, deviation, rho,
                   mcmc_coherent, mcmc_bamp, seed)
    es <- vs <- esm <- gapvar <- numeric(h)
    clip <- function(x) pmin(pmax(x, 1e-9), 1 - 1e-9)
    for (tt in seq_len(h)) {
      E <- tf(ens[[tt]]); o <- obs[[tt]]
      es[tt]  <- energy_score(o, E)
      vs[tt]  <- variogram_score(o, E, p = vs_p)
      esm[tt] <- energy_score(o[1:I], E[, 1:I, drop = FALSE]) +
                 energy_score(o[I + 1:I], E[, I + 1:I, drop = FALSE])
      # non-divergence diagnostic: across-draw variance of the period-component of
      # the sex gap (mean-over-ages logit difference), on the raw rate-scale ensemble.
      R <- ens[[tt]]
      gapvar[tt] <- stats::var(rowMeans(stats::qlogis(clip(R[, 1:I, drop = FALSE]))) -
                               rowMeans(stats::qlogis(clip(R[, I + 1:I, drop = FALSE]))))
    }
    data.frame(model = model, energy = mean(es), variogram = mean(vs),
               energy_marginal = mean(esm),
               gap_var_first = gapvar[1], gap_var_last = gapvar[h],
               gap_growth = gapvar[h] / gapvar[1], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res); rownames(out) <- NULL
  out
}

#' Select the AR1 coherence strength rho by out-of-sample score
#'
#' @description
#' Choose the deviation AR1 coefficient \code{rho} (the coherence dial of \code{\link{bamp_coherent}})
#' from the data, instead of fixing it, by minimising a held-out multivariate score over a grid. This
#' addresses the Jallbjorn caution: rather than asserting how strongly the sex gap reverts, let the
#' out-of-sample fit decide. \code{rho = 0} corresponds to the \code{"iid"} deviation.
#'
#' @param cases,population,holdout,periods_per_agegroup,age,period,cohort as in
#'   \code{\link{coherence_backtest}}.
#' @param rho_grid grid of AR1 coefficients to evaluate (each in \code{[0, 1)}).
#' @param score which score to minimise: \code{"energy"} or \code{"variogram"}.
#' @param scale,vs_p,mcmc_coherent,seed as in \code{\link{coherence_backtest}}.
#'
#' @return list with \code{table} (data.frame of \code{rho} and its held-out scores) and
#'   \code{best_rho} (the minimiser of the chosen score).
#' @seealso \code{\link{coherence_backtest}}, \code{\link{bamp_coherent}}
#' @export
select_rho <- function(cases, population, holdout = 2, periods_per_agegroup,
                       age = "rw1", period = "rw1", cohort = "rw1",
                       rho_grid = c(0, 0.3, 0.6, 0.9),
                       score = c("energy", "variogram"),
                       scale = c("rate", "lograte"), vs_p = 0.5,
                       mcmc_coherent = list(iterations = 4000, burn_in = 1000, thin = 2),
                       seed = 1) {
  score <- match.arg(score); scale <- match.arg(scale)
  rows <- lapply(rho_grid, function(r) {
    bt <- coherence_backtest(cases, population, holdout = holdout,
                             periods_per_agegroup = periods_per_agegroup,
                             age = age, period = period, cohort = cohort,
                             models = "coherent",
                             deviation = if (r == 0) "iid" else "ar1", rho = r,
                             scale = scale, vs_p = vs_p, mcmc_coherent = mcmc_coherent, seed = seed)
    data.frame(rho = r, energy = bt$energy, variogram = bt$variogram)
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  list(table = tab, best_rho = tab$rho[which.min(tab[[score]])])
}
