## ===========================================================================
## Nested disease-taxonomy CASCADE (scaffold) -- the immediately-runnable route
## to MANY competing diseases. Generalises the two-level "total x share" layering
## to a tree: all-cause -> group shares -> within-group leaf shares. Coherence
## telescopes (leaf = all-cause x group-share x leaf-share; shares sum to 1 at
## every level), so leaves sum to their group and groups sum to all-cause, by
## construction. Each level/group is a SMALL bamp_multicause() fit, so 30 diseases
## run on the dense reference engine without ever forming one huge joint field.
##
## The cascade couples disease trends WITHIN a group (they share a group trend).
## Cross-group, risk-factor-driven coupling (smoking across cancer + CVD + COPD)
## is a tree-orthogonal dependence -- a planned extension via declared risk-factor
## covariates or a latent factor on the innovations (see docs/coherent-forecasting.md).
## ===========================================================================

#' Coherent competing-disease projection over a disease taxonomy (cascade)
#'
#' @description
#' Fit many competing diseases organised in a 2-level taxonomy (groups of leaf diseases) as a
#' coherent cascade: a group-level \code{\link{bamp_multicause}} (groups partition all-cause) plus,
#' within each group, a leaf-level \code{\link{bamp_multicause}} (leaves partition the group). The
#' projected leaf rates/hazards are coherent at \emph{every} level by construction. Because each fit
#' is small (a handful of groups; a handful of leaves per group), this scales to dozens of diseases
#' on the dense reference engine, where a single flat multinomial over all leaves would not.
#'
#' @param taxonomy named list mapping each group to a character vector of its leaf-disease names
#'   (e.g. \code{list(CVD = c("CHD","stroke"), cancer = c("lung","colorectal"))}). Need >= 2 groups.
#' @param cases named list of leaf-disease count matrices (\code{[periods x agegroups]}); the names
#'   must cover exactly the leaves in \code{taxonomy}. Leaves partition all deaths.
#' @param population a single shared \code{[periods x agegroups]} population matrix.
#' @param age,period,cohort,periods_per_agegroup model settings, as in \code{\link{bamp_multicause}}.
#' @param mcmc,order,hyper,prior_scale,seed passed to each \code{\link{bamp_multicause}} fit.
#' @param factor optional integer: low-rank cross-cause coupling (see \code{\link{bamp_multicause}}),
#'   applied within each fit that has more shares than \code{factor} (else full Wishart).
#'
#' @return object of class \code{apc_cascade}: the fitted group model, the per-group leaf models, the
#'   taxonomy, and metadata; used by \code{\link{predict_cascade}}.
#' @seealso \code{\link{predict_cascade}}, \code{\link{bamp_multicause}}, \code{\link{bamp_coherent}}
#' @export
bamp_cascade <- function(taxonomy, cases, population, age = "rw1", period = "rw1", cohort = "rw1",
                         periods_per_agegroup,
                         mcmc = list(iterations = 4000, burn_in = 1000, thin = 2),
                         order = "prevalence", hyper = NULL, prior_scale = TRUE, seed = 1,
                         factor = NULL) {
  if (!is.list(taxonomy) || length(taxonomy) < 2L)
    stop("'taxonomy' must be a named list of >= 2 groups, each a vector of leaf-disease names.")
  if (is.null(names(taxonomy))) names(taxonomy) <- paste0("group", seq_along(taxonomy))
  leaves_all <- unlist(taxonomy, use.names = FALSE)
  if (anyDuplicated(leaves_all)) stop("each leaf disease may belong to only one group.")
  if (!setequal(leaves_all, names(cases)))
    stop("names(cases) must be exactly the leaf diseases listed in 'taxonomy'.")
  cases <- lapply(cases, as.matrix)
  hy <- if (is.null(hyper)) list(age = c(1, 0.5), omega = c(2, 1e-4), omega_c = c(2, 1e-4)) else hyper
  ## factor=R adds low-rank cross-cause coupling within a fit (the cross-cutting drivers); applied
  ## per fit only where it has enough causes (R < number of shares), else falls back to Wishart.
  fit_mc <- function(cl) {
    R <- if (!is.null(factor) && (length(cl) - 1L) > as.integer(factor)) as.integer(factor) else NULL
    bamp_multicause(cl, population, age = age, period = period, cohort = cohort,
                    periods_per_agegroup = periods_per_agegroup, order = order,
                    mcmc = mcmc, hyper = hy, prior_scale = prior_scale, seed = seed, factor = R)
  }

  ## group level: each group = sum of its leaves; groups partition all-cause
  message("bamp_cascade: fitting group level (", length(taxonomy), " groups) ...")
  group_counts <- lapply(taxonomy, function(lv) Reduce(`+`, cases[lv]))
  group_fit <- fit_mc(group_counts)

  ## within each group: leaves partition the group total (singletons need no split)
  within_fits <- stats::setNames(lapply(names(taxonomy), function(g) {
    lv <- taxonomy[[g]]
    if (length(lv) < 2L) return(NULL)
    message("bamp_cascade: fitting leaves within '", g, "' (", length(lv), " diseases) ...")
    fit_mc(cases[lv])
  }), names(taxonomy))

  structure(list(group = group_fit, within = within_fits, taxonomy = taxonomy,
                 data = list(population = population, leaves = leaves_all,
                             periods_per_agegroup = as.integer(periods_per_agegroup))),
            class = "apc_cascade")
}


#' Coherent projection from a disease-taxonomy cascade
#'
#' @description
#' Project every leaf disease from a \code{\link{bamp_cascade}} object by layering the group share
#' (group / all-cause) and the within-group leaf share (leaf / group). Leaf rates sum to their group
#' rate, group rates sum to all-cause; cause-specific hazards likewise sum to the all-cause hazard at
#' every level.
#'
#' @param object an \code{apc_cascade} object.
#' @param periods number of future periods to project.
#' @param population optional future \code{[periods x agegroups]} population (NULL uses in-sample).
#' @param quantiles quantiles to summarise.
#' @param hazard,period_length if \code{hazard=TRUE} also return additive per-person-year hazards
#'   (\code{leaf-share x group-share x all-cause hazard}; coherent at every level).
#'
#' @return list with one entry per leaf disease, a \code{groups} list (one entry per group), a
#'   \code{total} (all-cause), each holding quantiles of \code{rate} (and \code{hazard}) on the
#'   \code{[period, agegroup]} grid plus \code{samples}; the taxonomy \code{leaves}/\code{taxonomy};
#'   and \code{coherence_maxerr} (largest deviation of summed leaves from all-cause; ~0).
#' @seealso \code{\link{bamp_cascade}}
#' @export
predict_cascade <- function(object, periods = 0, population = NULL,
                            quantiles = c(0.05, 0.5, 0.95), hazard = FALSE, period_length = 1) {
  if (!inherits(object, "apc_cascade")) stop("'object' must come from bamp_cascade().")
  hazard <- isTRUE(hazard)
  qf <- function(a) apply(a, 1:2, stats::quantile, quantiles)
  pack <- function(rt, hz = NULL) {
    e <- list(rate = qf(rt), samples = list(rate = rt))
    if (hazard && !is.null(hz)) { e$hazard <- qf(hz); e$samples$hazard <- hz }
    e
  }
  mul <- function(a, b) {            # elementwise product of two [period, age, draw] arrays, aligned
    D <- min(dim(a)[3], dim(b)[3])
    a[, , seq_len(D), drop = FALSE] * b[, , seq_len(D), drop = FALSE]
  }

  pg <- predict_multicause(object$group, periods = periods, population = population,
                           hazard = hazard, period_length = period_length)
  out <- list()
  for (g in names(object$taxonomy)) {
    lv <- object$taxonomy[[g]]
    grate <- pg[[g]]$samples$rate                      # group rate = group-share x all-cause
    ghaz  <- if (hazard) pg[[g]]$samples$hazard else NULL
    if (length(lv) < 2L) {                             # singleton: leaf == group
      out[[lv]] <- pack(grate, ghaz)
    } else {
      pw <- predict_multicause(object$within[[g]], periods = periods, population = population)
      gtot <- pw$total$samples$rate                    # within-fit's group total
      for (l in lv) {
        share <- pw[[l]]$samples$rate / gtot; share[!is.finite(share)] <- 0   # leaf share within group
        out[[l]] <- pack(mul(grate, share), if (hazard) mul(ghaz, share) else NULL)
      }
    }
  }
  out$groups <- stats::setNames(lapply(names(object$taxonomy), function(g)
    list(rate = pg[[g]]$rate, samples = pg[[g]]$samples)), names(object$taxonomy))
  out$total <- pack(pg$total$samples$rate, if (hazard) pg$total$samples$hazard else NULL)

  summed <- Reduce(`+`, lapply(object$data$leaves, function(l) out[[l]]$samples$rate))
  D <- min(dim(summed)[3], dim(out$total$samples$rate)[3])
  out$coherence_maxerr <- max(abs(summed[, , seq_len(D)] - out$total$samples$rate[, , seq_len(D)]))
  out$leaves <- object$data$leaves; out$taxonomy <- object$taxonomy
  out
}
