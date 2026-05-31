## ===========================================================================
## Out-of-sample validation harness for the competing-cause / cascade forecasts,
## with rolling origins (C3), demographic baselines (C5), and a per-cause
## calibration breakdown (C7). Reuses the Batch-1 scoring + calibration layer.
##
## Coherence and good mixing are necessary but not sufficient: the only test that
## a forecast is USABLE is that it forecasts held-out data well AND its intervals
## are calibrated. multicause_backtest() refits each model on a training span,
## projects the held-out periods, and scores the JOINT cause-by-age forecast with
## the proper multivariate rules plus marginal CRPS / coverage / PIT -- against
## independent-APC, naive random-walk, and Lee-Carter baselines.
## ===========================================================================

## rolling-origin train/test period index sets (C3).
## scheme: expanding training window (default) or a fixed sliding 'window' length.
.bt_origins <- function(J, holdout, n_origins, window = NULL) {
  if (holdout < 1 || holdout >= J) stop("'holdout' must be in 1..(periods-1).")
  if (n_origins < 1) stop("'n_origins' must be >= 1.")
  starts <- (J - holdout) - seq_len(n_origins) + 1L                  # last train period per origin
  starts <- starts[starts >= max(4L, holdout)]                      # need enough history to fit
  lapply(starts, function(last_train) {
    train_from <- if (is.null(window)) 1L else max(1L, last_train - window + 1L)
    list(train = train_from:last_train, test = last_train + seq_len(holdout))
  })
}

#' Naive random-walk-with-drift forecast ensemble (per-age baseline)
#'
#' @description
#' Independent random-walk-with-drift on each age's log (or logit) central rate:
#' drift and innovation sd estimated from the training first-differences, forecast
#' \eqn{s} steps ahead as \eqn{x_J + s\,\hat\delta + \sqrt{s}\,\hat\sigma Z}. The
#' simplest honest baseline a real model must beat.
#'
#' @param cases,population training \code{[periods x agegroups]} count / population matrices.
#' @param h forecast horizon (number of future periods).
#' @param draws ensemble size.
#' @param scale \code{"log"} (log rate) or \code{"logit"} (log odds).
#' @param eps small rate floor before the log/logit transform.
#' @return a list of length \code{h}; element \code{s} is a \code{draws x agegroups}
#'   rate ensemble for the \code{s}-th future period.
#' @seealso \code{\link{forecast_leecarter}}, \code{\link{multicause_backtest}}
#' @export
forecast_naive <- function(cases, population, h, draws = 2000,
                           scale = c("log", "logit"), eps = 1e-6) {
  scale <- match.arg(scale)
  g  <- if (scale == "logit") stats::qlogis else log
  gi <- if (scale == "logit") stats::plogis else exp
  m <- pmin(pmax(as.matrix(cases) / as.matrix(population), eps), 1 - eps)
  x <- g(m); J <- nrow(x); I <- ncol(x)
  d <- diff(x)                                                       # [(J-1) x I] first differences
  drift <- colMeans(d); sigma <- apply(d, 2, stats::sd); sigma[!is.finite(sigma)] <- 0
  last <- x[J, ]
  lapply(seq_len(h), function(s) {
    Z <- matrix(stats::rnorm(draws * I), draws, I)
    fx <- sweep(Z * rep(sqrt(s) * sigma, each = draws), 2, last + s * drift, "+")
    gi(fx)
  })
}

#' Lee-Carter forecast ensemble (demographic baseline)
#'
#' @description
#' Classic Lee-Carter (1992) on the log central rates \eqn{\log m_{x,t}=a_x+b_x k_t}
#' (first SVD component, \eqn{\sum_x b_x=1}, \eqn{\sum_t k_t=0}), with a
#' random-walk-with-drift on the mortality index \eqn{k_t}. Forecast uncertainty comes
#' from the \eqn{k_t} drift innovations -- the standard demographic forecasting yardstick.
#'
#' @param cases,population training \code{[periods x agegroups]} count / population matrices.
#' @param h forecast horizon.
#' @param draws ensemble size.
#' @param eps small rate floor before the log transform.
#' @return a list of length \code{h} of \code{draws x agegroups} rate ensembles.
#' @seealso \code{\link{forecast_naive}}, \code{\link{multicause_backtest}}
#' @export
forecast_leecarter <- function(cases, population, h, draws = 2000, eps = 1e-6) {
  m <- pmin(pmax(as.matrix(cases) / as.matrix(population), eps), 1 - eps)
  lm <- log(m); J <- nrow(lm); I <- ncol(lm)
  ax <- colMeans(lm)
  clm <- sweep(lm, 2, ax, "-")
  sv <- svd(clm)
  bx <- sv$v[, 1]; kt <- sv$u[, 1] * sv$d[1]
  sb <- sum(bx)
  if (sb == 0) sb <- 1
  bx <- bx / sb; kt <- kt * sb                                      # impose sum_x b_x = 1
  km <- mean(kt); kt <- kt - km; ax <- ax + bx * km                # impose sum_t k_t = 0
  drift <- (kt[J] - kt[1]) / (J - 1)
  sigma <- stats::sd(diff(kt)); if (!is.finite(sigma)) sigma <- 0
  lapply(seq_len(h), function(s) {
    kfs <- kt[J] + s * drift + sqrt(s) * sigma * stats::rnorm(draws)
    exp(matrix(ax, draws, I, byrow = TRUE) + outer(kfs, bx))        # [draws x I]
  })
}

## per-test-period [draws x (Cn*I)] joint ensembles + obs, for one train/test split.
.bt_mc_block <- function(model, cases_tr, pop_tr, pop_pred, test, J, I, Cn, cause_names,
                         order, ppa, age, period, cohort, mcmc, mcmc_bamp,
                         draws_baseline, seed) {
  Jtr <- nrow(cases_tr[[1]])
  hh <- max(test) - Jtr                                             # periods beyond training end
  off <- test - Jtr                                                 # 1..hh local future indices
  asm <- function(per_cause_list) {                                 # list[cause] -> list[test] [draws x Cn*I]
    lapply(seq_along(off), function(tt) {
      do.call(cbind, lapply(per_cause_list, function(E) E[[off[tt]]]))
    })
  }
  if (model == "coupled") {
    fit <- bamp_multicause(stats::setNames(cases_tr, cause_names), pop_tr,
                           periods_per_agegroup = ppa, order = order,
                           mcmc = mcmc, seed = seed)
    pr <- predict_multicause(fit, periods = hh, population = pop_pred)
    list(ens = lapply(off, function(s)
           do.call(cbind, lapply(cause_names, function(nm)
             t(pr[[nm]]$samples$rate[Jtr + s, , ])))), fit = fit)
  } else if (model == "independent") {
    per <- lapply(cause_names, function(nm) {
      f <- bamp(cases_tr[[nm]], pop_tr, age = age, period = period, cohort = cohort,
                periods_per_agegroup = ppa, mcmc.options = mcmc_bamp,
                parallel = FALSE, verbose = FALSE)
      p <- predict_apc(f, periods = hh, population = pop_pred)
      lapply(off, function(s) t(p$samples$pr[Jtr + s, , ]))
    })
    list(ens = asm(per), fit = NULL)
  } else if (model %in% c("naive", "leecarter")) {
    fc <- if (model == "naive") forecast_naive else forecast_leecarter
    per <- lapply(cause_names, function(nm) fc(cases_tr[[nm]], pop_tr, hh, draws = draws_baseline))
    list(ens = asm(per), fit = NULL)
  } else stop("unknown model '", model, "'")
}

#' Out-of-sample backtest of competing-cause forecasts
#'
#' @description
#' Rolling-origin, multivariate-scored backtest of the coupled competing-cause model
#' (\code{\link{bamp_multicause}}) against independent-APC, naive random-walk and
#' Lee-Carter baselines. For each origin the models are refit on a training span and
#' the held-out periods projected; the \emph{joint} cause-by-age forecast is scored with
#' \code{\link{energy_score}} / \code{\link{variogram_score}} plus marginal
#' \code{\link{crps_sample}} and interval \code{\link{calibration}}. A per-cause
#' breakdown is returned so you can see which diseases are well- or badly-calibrated.
#'
#' @param cases named list of \code{Cn} cause-specific \code{[periods x agegroups]} count matrices.
#' @param population a single \code{[periods x agegroups]} population matrix (shared denominator).
#' @param holdout number of final periods to forecast at each origin.
#' @param n_origins number of rolling origins (each steps the training end back by one period).
#' @param window fixed training-window length for a sliding origin; \code{NULL} = expanding.
#' @param models any of \code{"coupled"}, \code{"independent"}, \code{"naive"}, \code{"leecarter"}.
#' @param order cause ordering passed to \code{\link{bamp_multicause}}.
#' @param scale score on the \code{"rate"} or \code{"lograte"} scale.
#' @param vs_p order of the variogram score.
#' @param levels nominal central-interval levels for coverage.
#' @param draws_baseline ensemble size for the naive / Lee-Carter baselines.
#' @param periods_per_agegroup,age,period,cohort model settings (as in \code{\link{bamp}}).
#' @param mcmc,mcmc_bamp MCMC settings for \code{\link{bamp_multicause}} and \code{\link{bamp}}.
#' @param seed RNG seed.
#'
#' @return a list of class \code{bamp_backtest} with \code{overall} (data.frame: one row per
#'   model, origin-averaged \code{energy}, \code{variogram}, \code{crps}, \code{logs},
#'   \code{cov90}, \code{pit_p}), \code{per_origin} (per-origin scores), and \code{per_cause}
#'   (per-disease \code{energy}/\code{crps}/\code{cov90} for the coupled model; the C7 report).
#' @seealso \code{\link{coherence_backtest}}, \code{\link{calibration}}, \code{\link{validation_report}}
#' @export
multicause_backtest <- function(cases, population, holdout = 2, n_origins = 1, window = NULL,
                                models = c("coupled", "independent", "naive", "leecarter"),
                                order = "prevalence", scale = c("rate", "lograte"), vs_p = 0.5,
                                levels = c(0.5, 0.8, 0.9, 0.95), draws_baseline = 2000,
                                periods_per_agegroup, age = "rw1", period = "rw1", cohort = "rw1",
                                mcmc = list(iterations = 3000, burn_in = 1000, thin = 2),
                                mcmc_bamp = list(number_of_iterations = 3000, burn_in = 1000,
                                                 step = 2, tuning = 300),
                                seed = 1) {
  scale <- match.arg(scale)
  if (is.null(names(cases))) names(cases) <- paste0("cause", seq_along(cases))
  cause_names <- names(cases); Cn <- length(cases)
  cases <- lapply(cases, as.matrix); population <- as.matrix(population)
  J <- nrow(cases[[1]]); I <- ncol(cases[[1]])
  tf <- if (scale == "lograte") function(x) log(pmax(x, 1e-6)) else identity
  emp <- lapply(cases, function(x) x / population)
  origins <- .bt_origins(J, holdout, n_origins, window)
  if (!length(origins)) stop("not enough history for the requested origins/holdout.")

  perm <- list()
  for (oi in seq_along(origins)) {
    tr <- origins[[oi]]$train; te <- origins[[oi]]$test
    cases_tr <- lapply(cases, function(x) x[tr, , drop = FALSE])
    pop_tr <- population[tr, , drop = FALSE]
    pop_pred <- population[seq_len(max(te)), , drop = FALSE]         # train + held-out rows
    obs <- lapply(te, function(t) tf(unlist(lapply(emp, function(e) e[t, ]))))   # length Cn*I
    for (model in models) {
      bl <- .bt_mc_block(model, cases_tr, pop_tr, pop_pred, te, J, I, Cn, cause_names, order,
                         periods_per_agegroup, age, period, cohort, mcmc, mcmc_bamp,
                         draws_baseline, seed)
      ensT <- lapply(bl$ens, tf)
      es <- vapply(seq_along(te), function(tt) energy_score(obs[[tt]], ensT[[tt]]), numeric(1))
      vs <- vapply(seq_along(te), function(tt) variogram_score(obs[[tt]], ensT[[tt]], p = vs_p), numeric(1))
      cal <- calibration(obs, ensT, levels = levels)
      cov90 <- cal$coverage$empirical[cal$coverage$level == 0.9]
      perm[[length(perm) + 1L]] <- data.frame(
        origin = oi, model = model, energy = mean(es), variogram = mean(vs),
        crps = cal$crps, logs = cal$logs,
        cov90 = if (length(cov90)) cov90 else NA_real_,
        pit_p = cal$pit_uniformity$p_value, stringsAsFactors = FALSE)
      ## per-cause breakdown (C7) for the coupled model, last origin
      if (model == "coupled" && oi == length(origins)) {
        pc <- lapply(seq_len(Cn), function(c) {
          idx <- (c - 1L) * I + seq_len(I)
          oC <- lapply(obs, function(o) o[idx]); eC <- lapply(ensT, function(E) E[, idx, drop = FALSE])
          calc <- calibration(oC, eC, levels = levels)
          data.frame(cause = cause_names[c],
                     energy = mean(vapply(seq_along(te), function(tt) energy_score(oC[[tt]], eC[[tt]]), numeric(1))),
                     crps = calc$crps,
                     cov90 = calc$coverage$empirical[calc$coverage$level == 0.9],
                     pit_p = calc$pit_uniformity$p_value, stringsAsFactors = FALSE)
        })
        per_cause <- do.call(rbind, pc)
      }
    }
  }
  per_origin <- do.call(rbind, perm); rownames(per_origin) <- NULL
  agg <- stats::aggregate(cbind(energy, variogram, crps, logs, cov90, pit_p) ~ model,
                          data = per_origin, FUN = function(z) mean(z, na.rm = TRUE),
                          na.action = stats::na.pass)
  agg <- agg[order(agg$energy), ]; rownames(agg) <- NULL
  structure(list(overall = agg, per_origin = per_origin,
                 per_cause = if (exists("per_cause")) per_cause else NULL,
                 per_unit_label = "per-cause (coupled model)",
                 scale = scale, n_origins = length(origins)),
            class = "bamp_backtest")
}

#' @method print bamp_backtest
#' @export
print.bamp_backtest <- function(x, ...) {
  cat(sprintf("bamp out-of-sample backtest (%d origin(s), %s scale; lower energy/crps = better)\n",
              x$n_origins, x$scale))
  print(x$overall, row.names = FALSE)
  lab <- if (is.null(x$per_unit_label)) "per-unit" else x$per_unit_label
  if (!is.null(x$per_cause)) { cat("\n", lab, ":\n", sep = ""); print(x$per_cause, row.names = FALSE) }
  invisible(x)
}

## shared scorer: list of obs vectors + aligned ensembles -> one summary row
.bt_score_one <- function(obs, ensT, vs_p, levels) {
  h <- length(obs)
  es <- vapply(seq_len(h), function(tt) energy_score(obs[[tt]], ensT[[tt]]), numeric(1))
  vs <- vapply(seq_len(h), function(tt) variogram_score(obs[[tt]], ensT[[tt]], p = vs_p), numeric(1))
  cal <- calibration(obs, ensT, levels = levels)
  c90 <- cal$coverage$empirical[cal$coverage$level == 0.9]
  data.frame(energy = mean(es), variogram = mean(vs), crps = cal$crps, logs = cal$logs,
             cov90 = if (length(c90)) c90 else NA_real_,
             pit_p = cal$pit_uniformity$p_value, stringsAsFactors = FALSE)
}

## per-test-period leaf-by-age joint ensembles for one cascade train/test split.
.bt_casc_block <- function(model, taxonomy, cases_tr, pop_tr, pop_pred, off, leaves,
                           age, period, cohort, ppa, mcmc, mcmc_bamp, draws_baseline, seed) {
  Jtr <- nrow(pop_tr)
  if (model == "cascade") {
    fit <- bamp_cascade(taxonomy, stats::setNames(cases_tr, leaves), pop_tr,
                        age = age, period = period, cohort = cohort,
                        periods_per_agegroup = ppa, mcmc = mcmc, seed = seed)
    pr <- predict_cascade(fit, periods = max(off), population = pop_pred)
    Dmin <- min(vapply(leaves, function(l) dim(pr[[l]]$samples$rate)[3], integer(1)))
    list(ens = lapply(off, function(s) do.call(cbind, lapply(leaves, function(l)
           t(pr[[l]]$samples$rate[Jtr + s, , seq_len(Dmin)])))), fit = fit)
  } else if (model == "independent") {
    per <- lapply(leaves, function(l) {
      f <- bamp(cases_tr[[l]], pop_tr, age = age, period = period, cohort = cohort,
                periods_per_agegroup = ppa, mcmc.options = mcmc_bamp, parallel = FALSE, verbose = FALSE)
      p <- predict_apc(f, periods = max(off), population = pop_pred)
      lapply(off, function(s) t(p$samples$pr[Jtr + s, , ]))
    })
    list(ens = lapply(seq_along(off), function(tt) do.call(cbind, lapply(per, function(E) E[[tt]]))), fit = NULL)
  } else {
    fc <- if (model == "naive") forecast_naive else forecast_leecarter
    per <- lapply(leaves, function(l) fc(cases_tr[[l]], pop_tr, max(off), draws = draws_baseline))
    list(ens = lapply(seq_along(off), function(tt) do.call(cbind, lapply(per, function(E) E[[tt]]))), fit = NULL)
  }
}

#' Out-of-sample backtest of a disease-taxonomy cascade forecast
#'
#' @description
#' Rolling-origin, multivariate-scored backtest of the coherent disease cascade
#' (\code{\link{bamp_cascade}}) against independent-APC, naive and Lee-Carter baselines,
#' scoring the \emph{joint} leaf-by-age forecast plus a per-leaf calibration breakdown.
#' The cascade analogue of \code{\link{multicause_backtest}}.
#'
#' @param taxonomy named list of groups, each a vector of leaf-disease names (as in \code{\link{bamp_cascade}}).
#' @param cases named list of leaf-disease count matrices; names must cover the taxonomy leaves.
#' @param population a single \code{[periods x agegroups]} population matrix.
#' @param holdout,n_origins,window rolling-origin controls (see \code{\link{multicause_backtest}}).
#' @param models any of \code{"cascade"}, \code{"independent"}, \code{"naive"}, \code{"leecarter"}.
#' @param scale,vs_p,levels,draws_baseline scoring controls (see \code{\link{multicause_backtest}}).
#' @param periods_per_agegroup,age,period,cohort,mcmc,mcmc_bamp,seed model / MCMC settings.
#'
#' @return a \code{bamp_backtest} object (\code{overall}, \code{per_origin}, and \code{per_cause}
#'   = the per-leaf breakdown for the cascade model).
#' @seealso \code{\link{multicause_backtest}}, \code{\link{validation_report}}
#' @export
cascade_backtest <- function(taxonomy, cases, population, holdout = 2, n_origins = 1, window = NULL,
                             models = c("cascade", "independent", "naive", "leecarter"),
                             scale = c("rate", "lograte"), vs_p = 0.5,
                             levels = c(0.5, 0.8, 0.9, 0.95), draws_baseline = 2000,
                             periods_per_agegroup, age = "rw1", period = "rw1", cohort = "rw1",
                             mcmc = list(iterations = 3000, burn_in = 1000, thin = 2),
                             mcmc_bamp = list(number_of_iterations = 3000, burn_in = 1000,
                                              step = 2, tuning = 300),
                             seed = 1) {
  scale <- match.arg(scale)
  leaves <- unlist(taxonomy, use.names = FALSE)
  cases <- lapply(cases, as.matrix)[leaves]; population <- as.matrix(population)
  J <- nrow(population); I <- ncol(population); nl <- length(leaves)
  tf <- if (scale == "lograte") function(x) log(pmax(x, 1e-6)) else identity
  emp <- lapply(cases, function(x) x / population)
  origins <- .bt_origins(J, holdout, n_origins, window)
  if (!length(origins)) stop("not enough history for the requested origins/holdout.")

  perm <- list(); per_leaf <- NULL
  for (oi in seq_along(origins)) {
    tr <- origins[[oi]]$train; te <- origins[[oi]]$test; off <- te - max(tr)
    cases_tr <- lapply(cases, function(x) x[tr, , drop = FALSE])
    pop_tr <- population[tr, , drop = FALSE]; pop_pred <- population[seq_len(max(te)), , drop = FALSE]
    obs <- lapply(te, function(t) tf(unlist(lapply(emp, function(e) e[t, ]))))
    for (model in models) {
      bl <- .bt_casc_block(model, taxonomy, cases_tr, pop_tr, pop_pred, off, leaves,
                           age, period, cohort, periods_per_agegroup, mcmc, mcmc_bamp, draws_baseline, seed)
      ensT <- lapply(bl$ens, tf)
      row <- .bt_score_one(obs, ensT, vs_p, levels)
      perm[[length(perm) + 1L]] <- cbind(data.frame(origin = oi, model = model), row)
      if (model == "cascade" && oi == length(origins)) {
        per_leaf <- do.call(rbind, lapply(seq_len(nl), function(c) {
          idx <- (c - 1L) * I + seq_len(I)
          r <- .bt_score_one(lapply(obs, function(o) o[idx]),
                             lapply(ensT, function(E) E[, idx, drop = FALSE]), vs_p, levels)
          cbind(data.frame(leaf = leaves[c]), r[, c("energy", "crps", "cov90", "pit_p")])
        }))
      }
    }
  }
  per_origin <- do.call(rbind, perm); rownames(per_origin) <- NULL
  agg <- stats::aggregate(cbind(energy, variogram, crps, logs, cov90, pit_p) ~ model,
                          data = per_origin, FUN = function(z) mean(z, na.rm = TRUE),
                          na.action = stats::na.pass)
  agg <- agg[order(agg$energy), ]; rownames(agg) <- NULL
  structure(list(overall = agg, per_origin = per_origin, per_cause = per_leaf,
                 per_unit_label = "per-leaf (cascade model)",
                 scale = scale, n_origins = length(origins)), class = "bamp_backtest")
}

#' Select random-walk orders by out-of-sample score
#'
#' @description
#' Choose the period and cohort prior order (\code{"rw1"} vs \code{"rw2"}) for a single-population
#' APC model from the data rather than by assumption, by rolling-origin out-of-sample scoring. RW1
#' forecasts a flat trend (constant + noise); RW2 continues the last linear trend (and over-extrapolates
#' at long horizons unless damped, see \code{\link{predict_apc}}'s \code{damping}). This refits each
#' candidate on a training span and scores the held-out age profile, returning the score table and the
#' minimiser -- the principled answer to "RW1 or RW2?".
#'
#' @param cases,population single \code{[periods x agegroups]} count / population matrices.
#' @param periods_per_agegroup period/age grid ratio (as in \code{\link{bamp}}).
#' @param holdout,n_origins,window rolling-origin controls (see \code{\link{multicause_backtest}}).
#' @param age fixed age order (usually \code{"rw1"} or \code{"rw2"}).
#' @param period_orders,cohort_orders candidate orders to compare.
#' @param score \code{"energy"} (joint over ages) or \code{"crps"} (marginal).
#' @param scale score on the \code{"lograte"} (default) or \code{"rate"} scale.
#' @param mcmc_bamp,seed MCMC settings / seed for \code{\link{bamp}}.
#' @return a list with \code{table} (one row per period x cohort combination, mean held-out
#'   \code{energy} and \code{crps}) and \code{best} (the minimising \code{period}/\code{cohort}).
#' @seealso \code{\link{multicause_backtest}}, \code{\link{predict_apc}}
#' @export
select_rw_order <- function(cases, population, periods_per_agegroup, holdout = 3, n_origins = 1,
                            window = NULL, age = "rw1",
                            period_orders = c("rw1", "rw2"), cohort_orders = c("rw1", "rw2"),
                            score = c("energy", "crps"), scale = c("lograte", "rate"),
                            mcmc_bamp = list(number_of_iterations = 3000, burn_in = 1000,
                                             step = 2, tuning = 300), seed = 1) {
  score <- match.arg(score); scale <- match.arg(scale)
  cases <- as.matrix(cases); population <- as.matrix(population)
  J <- nrow(cases); I <- ncol(cases)
  tf <- if (scale == "lograte") function(x) log(pmax(x, 1e-6)) else identity
  emp <- cases / population
  origins <- .bt_origins(J, holdout, n_origins, window)
  combos <- expand.grid(period = period_orders, cohort = cohort_orders, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(combos)), function(ci) {
    es <- cr <- numeric(0)
    for (o in origins) {
      tr <- o$train; te <- o$test; Jtr <- max(tr)
      f <- bamp(cases[tr, , drop = FALSE], population[tr, , drop = FALSE], age = age,
                period = combos$period[ci], cohort = combos$cohort[ci],
                periods_per_agegroup = periods_per_agegroup, mcmc.options = mcmc_bamp,
                parallel = FALSE, verbose = FALSE)
      p <- predict_apc(f, periods = max(te) - Jtr, population = population[seq_len(max(te)), , drop = FALSE])
      for (tt in te) {
        o_v <- tf(emp[tt, ]); E <- tf(t(p$samples$pr[tt, , ]))
        es <- c(es, energy_score(o_v, E)); cr <- c(cr, crps_sample(o_v, E))
      }
    }
    data.frame(period = combos$period[ci], cohort = combos$cohort[ci],
               energy = mean(es), crps = mean(cr), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab[[score]]), ]; rownames(tab) <- NULL
  list(table = tab, best = list(period = tab$period[1], cohort = tab$cohort[1]))
}

#' Convergence pass/fail report across many fits
#'
#' @description
#' Run \code{\link{bamp_diagnostics}} on each of a named list of fits and tabulate the
#' worst split-\eqn{\hat R}, smallest ESS and flag count per fit -- the operational tool
#' for a 30-disease x strata sweep, where you cannot inspect every chain. Sorts the
#' worst fits first and reports an overall pass/fail.
#'
#' @param fits a named list of fitted objects (coherent / multicause / cascade ...).
#' @param rhat_max,ess_min thresholds passed to \code{\link{bamp_diagnostics}}.
#' @return a list of class \code{bamp_convergence_report} with \code{table} (one row per fit:
#'   \code{fit}, \code{n_params}, \code{max_rhat}, \code{min_ess}, \code{n_flagged}, \code{ok})
#'   and \code{all_ok} (logical).
#' @seealso \code{\link{bamp_diagnostics}}, \code{\link{validation_report}}
#' @export
convergence_report <- function(fits, rhat_max = 1.01, ess_min = 400) {
  if (is.null(names(fits))) names(fits) <- paste0("fit", seq_along(fits))
  rows <- lapply(names(fits), function(nm) {
    d <- tryCatch(bamp_diagnostics(fits[[nm]], rhat_max = rhat_max, ess_min = ess_min),
                  error = function(e) NULL)
    if (is.null(d)) return(data.frame(fit = nm, n_params = NA_integer_, max_rhat = NA_real_,
                                      min_ess = NA_real_, n_flagged = NA_integer_, ok = NA))
    cbind(data.frame(fit = nm), d$summary, ok = d$ok)
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  tab <- tab[order(-replace(tab$max_rhat, is.na(tab$max_rhat), -Inf)), ]
  structure(list(table = tab, all_ok = isTRUE(all(tab$ok))),
            class = "bamp_convergence_report")
}

#' @method print bamp_convergence_report
#' @export
print.bamp_convergence_report <- function(x, ...) {
  cat(sprintf("bamp convergence report: %d fit(s), %s\n", nrow(x$table),
              if (x$all_ok) "ALL OK" else "SOME FLAGGED (worst first)"))
  print(x$table, row.names = FALSE); invisible(x)
}

#' One-call convergence + calibration validation report
#'
#' @description
#' Bundle the two questions you must answer before trusting a forecast -- "did the
#' sampler converge?" (\code{\link{bamp_diagnostics}}) and "is it calibrated / does it
#' beat the baselines out-of-sample?" (\code{\link{multicause_backtest}} /
#' \code{\link{cascade_backtest}}) -- into a single object with a one-screen summary.
#'
#' @param fit a fitted coherent / multicause / cascade object.
#' @param backtest optional result of \code{\link{multicause_backtest}} /
#'   \code{\link{cascade_backtest}} for the same data.
#' @param rhat_max,ess_min convergence thresholds.
#' @return a list of class \code{bamp_validation} with \code{diagnostics}, \code{backtest},
#'   and \code{verdict} (logical: converged AND, if a backtest is supplied, the model is not
#'   the worst-calibrated by PIT).
#' @seealso \code{\link{bamp_diagnostics}}, \code{\link{multicause_backtest}}
#' @export
validation_report <- function(fit, backtest = NULL, rhat_max = 1.01, ess_min = 400) {
  diag <- bamp_diagnostics(fit, rhat_max = rhat_max, ess_min = ess_min)
  verdict <- diag$ok
  structure(list(diagnostics = diag, backtest = backtest, verdict = verdict),
            class = "bamp_validation")
}

#' @method print bamp_validation
#' @export
print.bamp_validation <- function(x, ...) {
  cat("== bamp validation report ==\n[1] convergence\n")
  print(x$diagnostics)
  if (!is.null(x$backtest)) { cat("\n[2] out-of-sample calibration / skill\n"); print(x$backtest) }
  cat(sprintf("\nverdict: %s\n", if (x$verdict) "converged" else "NOT converged -- do not trust forecasts yet"))
  invisible(x)
}
