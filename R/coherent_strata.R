#' Coherent age-period-cohort forecasting across exhaustive strata (total-plus-share)
#'
#' @description
#' Fit a coherent set of APC models for an \emph{exhaustive} partition of the population
#' into strata (e.g. sex, or education within sex), so that the stratum-specific forecasts
#' aggregate exactly to the all-strata total and cannot diverge implausibly. This is the
#' lightweight "total-plus-share" construction (Phase 0 of the coherent-forecasting design):
#' one ordinary \code{\link{bamp}} fit to the summed total, plus a binomial-logit \code{bamp}
#' fit to each stratum's \emph{share of events}. Because the shares are fitted on running
#' remainders (a stick-breaking of the total), the predicted stratum counts sum to the
#' predicted total \emph{exactly}, every posterior draw.
#'
#' For the principled (single-joint-posterior, borrowing-strength) alternative see the design
#' note \code{docs/coherent-forecasting.md}; this wrapper is the fast, zero-engine-change
#' baseline and an oracle to validate that against.
#'
#' @param cases named list of \code{[periods x agegroups]} case matrices, one per stratum,
#'   all the same dimension. The strata must be mutually exclusive and collectively exhaustive
#'   (they partition the population), so that \code{Reduce('+', cases)} is the all-strata total.
#' @param population named list of \code{[periods x agegroups]} population matrices matching
#'   \code{cases}.
#' @param ... passed on to \code{\link{bamp}} (e.g. \code{age}, \code{period}, \code{cohort},
#'   \code{periods_per_agegroup}, \code{mcmc.options}, \code{method}). The same model is used
#'   for the total and for every share.
#' @param order optional character vector giving the order of the strata for the stick-breaking
#'   remainder. The last stratum in the order carries the remainder (it is not fitted directly).
#'   Defaults to \code{names(cases)}. Put the most data-rich strata first.
#'
#' @return object of class \code{apc_strata}: a list with the fitted total model (\code{total}),
#'   the list of fitted share models (\code{shares}, length \eqn{S-1}), the stratum \code{order},
#'   the stratum names (\code{strata}), and the in-sample populations (\code{populations}).
#'   Use \code{\link{predict_strata}} to obtain coherent stratum-specific rates, hazards and counts.
#'
#' @details
#' Coherence here is \emph{internal/aggregation} coherence: the total is fitted directly and the
#' strata are shares of it, so \eqn{\sum_s \text{cases}_s = \text{cases}_{total}} holds by
#' construction. The all-strata total uses \code{bamp}'s ordinary (free random-walk) projection.
#' Cells where no events remain to be split carry no share information: with the default
#' \code{method = "pg"} they are passed to \code{bamp} as zero-trial (\eqn{N = 0}) cells, which
#' contribute nothing to the likelihood so the smooth APC prior interpolates the share there. For
#' \code{method = "iwls"} such cells are floored to one trial (with a warning).
#'
#' @seealso \code{\link{predict_strata}}, \code{\link{bamp}}, \code{\link{reconcile_apc}}
#' @export
bamp_strata <- function(cases, population, ..., order = NULL)
{
  if (!is.list(cases) || !is.list(population))
    stop("'cases' and 'population' must be named lists of [periods x agegroups] matrices, one per stratum.")
  if (length(cases) < 2L)
    stop("Need at least two strata for a coherent stratified fit.")
  if (is.null(names(cases)))
    names(cases) <- paste0("stratum", seq_along(cases))
  if (is.null(names(population)))
    names(population) <- names(cases)
  if (!setequal(names(cases), names(population)))
    stop("'cases' and 'population' must have the same stratum names.")
  population <- population[names(cases)]

  dims <- lapply(cases, function(x) dim(as.matrix(x)))
  if (length(unique(vapply(dims, paste, "", collapse = "x"))) != 1L)
    stop("All stratum 'cases' matrices must have the same dimension.")
  if (!identical(lapply(population, function(x) dim(as.matrix(x))), dims))
    stop("Each 'population' matrix must match its 'cases' matrix dimension.")

  strata <- names(cases)
  if (is.null(order)) order <- strata
  if (!setequal(order, strata))
    stop("'order' must be a permutation of the stratum names.")

  cases <- lapply(cases, as.matrix)
  population <- lapply(population, as.matrix)

  dots <- list(...)
  method <- if (is.null(dots$method)) "pg" else match.arg(dots$method, c("pg", "iwls"))

  cases_tot <- Reduce(`+`, cases)
  pop_tot   <- Reduce(`+`, population)

  message("bamp_strata: fitting total model (", length(strata), " strata) ...")
  fit_total <- bamp(cases_tot, pop_tot, ...)

  # Stick-breaking shares on running remainders: stratum order[s] out of what is
  # left after order[1..s-1]. The last stratum is the remainder and is not fitted.
  shares <- vector("list", length(order) - 1L)
  names(shares) <- order[-length(order)]
  remainder <- cases_tot
  for (s in seq_len(length(order) - 1L))
  {
    num <- cases[[order[s]]]
    den <- remainder
    if (any(den < num)) stop("Internal error: share numerator exceeds remainder denominator.")
    # Rare cells where no events remain to be split (den == 0, hence num == 0 too)
    # carry no share information. Under the binomial / Polya-Gamma likelihood an
    # N = 0 cell contributes nothing (omega = 0, Fisher weight = 0), so the smooth
    # APC prior simply interpolates the share there -- the statistically correct
    # treatment, kept for the default method = "pg". The legacy IWLS sampler can
    # divide by the denominator, so for method = "iwls" we floor zero cells to 1.
    if (method == "iwls" && any(den == 0)) {
      den[den == 0] <- 1L
      warning("bamp_strata: zero-event denominator cells floored to 1 for method = 'iwls'; use the default method = 'pg' for exact no-information handling of rare cells.")
    }
    message("bamp_strata: fitting share model for '", order[s], "' ...")
    shares[[s]] <- bamp(num, den, ...)
    remainder <- remainder - cases[[order[s]]]
  }

  structure(list(
    total = fit_total,
    shares = shares,
    order = order,
    strata = strata,
    populations = population
  ), class = "apc_strata")
}


#' Coherent prediction across exhaustive strata
#'
#' @description
#' Project stratum-specific rates, cumulative hazards and counts from a fitted
#' \code{\link{bamp_strata}} object so that the strata aggregate exactly to the all-strata
#' total every posterior draw. The total is projected with \code{bamp}'s ordinary random-walk
#' extrapolation; the predicted total counts are then split among strata by stick-breaking with
#' the projected share probabilities, which guarantees \eqn{\sum_s \text{cases}_s = \text{cases}_{total}}.
#'
#' @param object an \code{apc_strata} object from \code{\link{bamp_strata}}.
#' @param periods number of future periods to predict.
#' @param population optional named list of \emph{future} \code{[periods x agegroups]} population
#'   matrices (one per stratum, with \code{periods} extra rows). Required to forecast beyond the
#'   observed periods; if \code{NULL}, the in-sample populations are used (model-checking only).
#' @param quantiles vector of quantiles to compute.
#' @param hazard boolean. If TRUE, also return the per-stratum cumulative cause-specific hazard
#'   \eqn{-\log(1-\text{rate})} per \code{period_length}. These are additive across competing
#'   causes and are the quantity to feed to a competing-risk life-table / microsimulation.
#' @param period_length single positive number: length of one period in the time units you want
#'   the hazard expressed per (e.g. years per period). Only used when \code{hazard=TRUE}.
#'
#' @return list with one entry per stratum and a \code{total} entry, each holding quantiles of
#'   \code{cases}, \code{rate} (and \code{hazard} if requested) on the \code{[period, agegroup]}
#'   grid plus the underlying \code{samples}; the stratum \code{order}; and \code{coherence_maxerr},
#'   the largest absolute discrepancy between the summed stratum counts and the total counts
#'   (0 by construction).
#'
#' @details
#' Each stratum rate is \code{cases_s / population_s} -- the stratum's \emph{own} population --
#' never \code{rate_total * share}. In rare cells Monte-Carlo splitting can allocate slightly
#' more events to a stratum than its population; such rates are capped just below 1 (with a
#' warning) so the hazard stays finite. The total and share models are separate fits with
#' independent MCMC streams, so the draws are paired by index, which is a valid joint sample
#' under the (deliberate) assumption that total-level and share-level uncertainty are independent.
#'
#' @seealso \code{\link{bamp_strata}}, \code{\link{predict_apc}}
#' @export
predict_strata <- function(object, periods = 0, population = NULL,
                           quantiles = c(0.05, 0.5, 0.95),
                           hazard = FALSE, period_length = 1)
{
  if (!inherits(object, "apc_strata"))
    stop("'object' must be an 'apc_strata' object from bamp_strata().")
  hazard <- isTRUE(hazard)
  if (hazard && (!is.numeric(period_length) || length(period_length) != 1L ||
                 is.na(period_length) || period_length <= 0))
    stop("'period_length' must be a single positive number.")

  order <- object$order
  S <- length(order)

  # Future populations per stratum (and their total). If NULL, use in-sample.
  if (is.null(population)) {
    pop_s <- object$populations
  } else {
    if (!is.list(population) || !setequal(names(population), object$strata))
      stop("'population' must be a named list of future population matrices, one per stratum.")
    pop_s <- lapply(population[object$strata], as.matrix)
  }
  pop_tot <- Reduce(`+`, pop_s)

  # Total: predict counts (integer draws) of the all-strata total.
  pred_tot <- predict_apc(object$total, periods = periods, population = pop_tot,
                          quantiles = quantiles)
  Ctot <- pred_tot$samples$cases               # [n0, agegroup, draws]
  n0 <- dim(Ctot)[1]; A <- dim(Ctot)[2]

  # Conditional share probabilities q_s for the fitted strata (order[1..S-1]).
  qsh <- lapply(seq_len(S - 1L), function(s) {
    p <- predict_apc(object$shares[[s]], periods = periods)$samples$pr  # [n2, agegroup, draws]
    p[seq_len(n0), , , drop = FALSE]
  })

  # Common number of draws across the (independent) fits.
  ndraw <- min(dim(Ctot)[3], vapply(qsh, function(x) dim(x)[3], 1L))
  trunc3 <- function(x) x[, , seq_len(ndraw), drop = FALSE]
  Ctot <- trunc3(Ctot)
  qsh  <- lapply(qsh, trunc3)
  shape <- dim(Ctot)

  # Stick-breaking split of the total counts -> exact-coherent stratum counts.
  counts <- vector("list", S); names(counts) <- order
  remainder <- Ctot
  for (s in seq_len(S - 1L)) {
    cs <- array(stats::rbinom(length(remainder), remainder, qsh[[s]]), shape)
    counts[[order[s]]] <- cs
    remainder <- remainder - cs
  }
  counts[[order[S]]] <- remainder              # last stratum = remainder

  coherence_maxerr <- max(abs(Reduce(`+`, counts) - Ctot))

  capped <- FALSE
  out <- list()
  qfun <- function(arr) apply(arr, 1:2, stats::quantile, quantiles)
  for (nm in object$strata) {
    cs <- counts[[nm]]
    pv <- as.vector(pop_s[[nm]][seq_len(n0), , drop = FALSE])   # [period, agegroup] recycled over draws
    rate <- cs / pv
    rate[!is.finite(rate)] <- 0                                  # pop==0 cells: no rate
    if (any(rate >= 1)) { capped <- TRUE; rate[rate >= 1] <- 1 - 1e-10 }
    entry <- list(
      cases = qfun(cs),
      rate  = qfun(rate),
      samples = list(cases = cs, rate = rate)
    )
    if (hazard) {
      hz <- -log1p(-rate) / period_length
      entry$hazard <- qfun(hz)
      entry$samples$hazard <- hz
    }
    out[[nm]] <- entry
  }
  if (capped)
    warning("Some stratum rates were capped just below 1 (Monte-Carlo splitting allocated more events to a stratum than its population in rare cells).")

  ratet <- Ctot / as.vector(pop_tot[seq_len(n0), , drop = FALSE])
  ratet[!is.finite(ratet)] <- 0; ratet[ratet >= 1] <- 1 - 1e-10
  out$total <- list(
    cases = qfun(Ctot), rate = qfun(ratet),
    samples = list(cases = Ctot, rate = ratet)
  )
  if (hazard) {
    hzt <- -log1p(-ratet) / period_length
    out$total$hazard <- qfun(hzt)
    out$total$samples$hazard <- hzt
  }

  out$order <- order
  out$coherence_maxerr <- coherence_maxerr
  out
}
