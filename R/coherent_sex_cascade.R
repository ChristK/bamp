## ===========================================================================
## Per-stratum disease cascade on a stratum-coherent total -- the full
## [period x age x sex x group x disease] coherent hazard tensor for a
## microsimulation. Composition of the two coherent ideas:
##   * bamp_coherent() sets the LEVEL: the all-cause total per sex, made
##     sex-coherent (non-diverging) by the shared trend + mean-reverting deviation.
##   * bamp_cascade() (per sex) sets the COMPOSITION: the disease shares of that
##     sex's all-cause, coherent down the taxonomy (group then leaf).
## Layered (leaf = sex-coherent total x cascade disease-share), the result is
## coherent on BOTH margins: diseases sum to each sex's total, and the sex totals
## stay plausibly related. (This is the vignette Part-C layering with a cascade.)
## ===========================================================================

#' Per-stratum disease cascade on a stratum-coherent total
#'
#' @description
#' Fit the full sex (or other exhaustive stratum) by disease-taxonomy model coherently: a
#' \code{\link{bamp_coherent}} model for the stratum all-cause totals (so the strata do not diverge)
#' plus a \code{\link{bamp_cascade}} disease taxonomy within each stratum. \code{\link{predict_sex_cascade}}
#' layers them into a hazard tensor that is coherent on both the stratum margin and the disease
#' margin -- the object an NCD microsimulation (e.g. IMPACTncd) consumes.
#'
#' @param taxonomy disease taxonomy, as in \code{\link{bamp_cascade}}.
#' @param cases named list with one element per stratum (e.g. \code{list(men=, women=)}); each element
#'   is a named list of leaf-disease count matrices for that stratum (names = taxonomy leaves).
#' @param population named list with one \code{[periods x agegroups]} population matrix per stratum.
#' @param age,period,cohort,periods_per_agegroup model settings.
#' @param deviation deviation prior for the stratum-coherent total (\code{\link{bamp_coherent}}).
#' @param mcmc_total,mcmc_cascade MCMC settings for the coherent total and each cascade fit.
#' @param factor,order,prior_scale,seed passed to the cascade fits (\code{\link{bamp_cascade}}).
#'
#' @return object of class \code{apc_sex_cascade}: the stratum-coherent total model, the per-stratum
#'   cascades, and metadata; used by \code{\link{predict_sex_cascade}}.
#' @seealso \code{\link{predict_sex_cascade}}, \code{\link{bamp_coherent}}, \code{\link{bamp_cascade}}
#' @export
bamp_sex_cascade <- function(taxonomy, cases, population, age = "rw1", period = "rw1", cohort = "rw1",
                             periods_per_agegroup, deviation = "ar1",
                             mcmc_total = list(iterations = 4000, burn_in = 1000, thin = 2),
                             mcmc_cascade = list(iterations = 4000, burn_in = 1000, thin = 2),
                             factor = NULL, order = "prevalence", prior_scale = TRUE, seed = 1) {
  if (!is.list(cases) || length(cases) < 2L || is.null(names(cases)))
    stop("'cases' must be a NAMED list with one element per stratum (each a named list of leaf counts).")
  sexes <- names(cases)
  if (!is.list(population) || !setequal(names(population), sexes))
    stop("'population' must be a named list of one matrix per stratum, matching 'cases'.")
  leaves <- unlist(taxonomy, use.names = FALSE)

  ## stratum all-cause totals -> one sex-coherent bamp_coherent fit (the LEVEL)
  message("bamp_sex_cascade: fitting the stratum-coherent all-cause total ...")
  totals <- stats::setNames(lapply(sexes, function(sx) Reduce(`+`, lapply(cases[[sx]], as.matrix))), sexes)
  fit_total <- bamp_coherent(totals, population[sexes], age = age, period = period, cohort = cohort,
                             periods_per_agegroup = periods_per_agegroup, deviation = deviation,
                             mcmc = mcmc_total, prior_scale = prior_scale, seed = seed)

  ## a disease cascade within each stratum (the COMPOSITION)
  fit_cascade <- stats::setNames(lapply(sexes, function(sx) {
    message("bamp_sex_cascade: fitting the disease cascade for '", sx, "' ...")
    bamp_cascade(taxonomy, cases[[sx]], population[[sx]], age = age, period = period, cohort = cohort,
                 periods_per_agegroup = periods_per_agegroup, mcmc = mcmc_cascade,
                 order = order, prior_scale = prior_scale, seed = seed, factor = factor)
  }), sexes)

  structure(list(total = fit_total, cascade = fit_cascade, sexes = sexes,
                 taxonomy = taxonomy, leaves = leaves), class = "apc_sex_cascade")
}


#' Coherent [stratum x disease] projection from a per-stratum disease cascade
#'
#' @description
#' Project the full stratum-by-disease hazard tensor from \code{\link{bamp_sex_cascade}}: the
#' stratum-coherent all-cause total (Part B) times each stratum's cascade disease shares (Part A,
#' nested). Diseases sum to each stratum's total (disease coherence) and the stratum totals stay
#' bounded/non-diverging (stratum coherence). Cause-specific hazards likewise sum on both margins.
#'
#' @param object an \code{apc_sex_cascade} object.
#' @param periods number of future periods to project.
#' @param population named list of \emph{future} \code{[periods x agegroups]} population per stratum.
#' @param quantiles quantiles to summarise.
#' @param hazard,period_length if \code{hazard=TRUE} also return additive per-person-year hazards
#'   (stratum-coherent total hazard times the cascade disease share; coherent on both margins).
#'
#' @return list with one element per stratum; each holds an entry per leaf disease and per group, a
#'   stratum \code{total}, each with quantiles of \code{rate} (and \code{hazard}) and \code{samples};
#'   plus \code{sexes}, \code{leaves}, and \code{coherence_maxerr} (largest within-stratum deviation
#'   of the summed diseases from the stratum total; ~0 by construction).
#' @seealso \code{\link{bamp_sex_cascade}}
#' @export
predict_sex_cascade <- function(object, periods = 0, population = NULL,
                                quantiles = c(0.05, 0.5, 0.95), hazard = FALSE, period_length = 1) {
  if (!inherits(object, "apc_sex_cascade")) stop("'object' must come from bamp_sex_cascade().")
  hazard <- isTRUE(hazard); sexes <- object$sexes
  qf <- function(a) apply(a, 1:2, stats::quantile, quantiles)
  pack <- function(rt, hz = NULL) {
    e <- list(rate = qf(rt), samples = list(rate = rt))
    if (hazard && !is.null(hz)) { e$hazard <- qf(hz); e$samples$hazard <- hz }; e
  }
  ## stratum-coherent all-cause total (non-diverging across strata)
  pc <- predict_coherent(object$total, periods = periods, population = population,
                         hazard = hazard, period_length = period_length)
  out <- list(); cmax <- 0
  for (sx in sexes) {
    pcs <- predict_cascade(object$cascade[[sx]], periods = periods,
                           population = if (is.null(population)) NULL else population[[sx]])
    ## align the coherent-total and cascade draws to a common count so the tensor is consistent
    Ds <- min(dim(pc[[sx]]$samples$rate)[3], dim(pcs$total$samples$rate)[3])
    tr <- function(a) a[, , seq_len(Ds), drop = FALSE]
    sex_rate <- tr(pc[[sx]]$samples$rate)
    sex_haz  <- if (hazard) tr(pc[[sx]]$samples$hazard) else NULL
    ctot <- tr(pcs$total$samples$rate)
    shr <- function(num) { s <- tr(num) / ctot; s[!is.finite(s)] <- 0; s }   # cascade share, aligned
    e <- list()
    for (l in object$leaves) {                            # disease share (within stratum) x stratum total
      sh <- shr(pcs[[l]]$samples$rate)
      e[[l]] <- pack(sex_rate * sh, if (hazard) sex_haz * sh else NULL)
    }
    e$groups <- stats::setNames(lapply(names(object$taxonomy), function(g) {
      gr <- sex_rate * shr(pcs$groups[[g]]$samples$rate); list(rate = qf(gr), samples = list(rate = gr))
    }), names(object$taxonomy))
    e$total <- pack(sex_rate, sex_haz)
    summed <- Reduce(`+`, lapply(object$leaves, function(l) e[[l]]$samples$rate))
    cmax <- max(cmax, max(abs(summed - sex_rate)))        # diseases sum to the stratum total (exact)
    out[[sx]] <- e
  }
  out$sexes <- sexes; out$leaves <- object$leaves; out$taxonomy <- object$taxonomy
  out$coherence_maxerr <- cmax
  out
}
