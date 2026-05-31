## ===========================================================================
## Proper MULTIVARIATE scoring rules for evaluating coherent forecasts.
##
## Coherence is a dependence property: a model can get each stratum's marginal
## forecast right yet get the joint (e.g. the sex gap) wrong. Marginal CRPS is
## blind to that; the energy score and the variogram score are proper scoring
## rules for the whole forecast vector and reward correct cross-stratum
## dependence. Both are negatively oriented (lower = better).
## ===========================================================================

#' Energy score of a multivariate ensemble forecast
#'
#' @description
#' Proper multivariate scoring rule (Gneiting & Raftery 2007). For an ensemble
#' forecast \eqn{X_1,\dots,X_m} (rows) and observation \eqn{y},
#' \deqn{ES = \frac1m\sum_i \lVert X_i - y\rVert - \frac1{2m^2}\sum_{i,j}\lVert X_i - X_j\rVert.}
#' Lower is better. Sensitive to the dependence structure of the forecast, so it
#' detects coherence gains that marginal scores cannot.
#'
#' @param obs numeric vector, the observed multivariate outcome (length d).
#' @param ens an \code{m x d} matrix of ensemble (e.g. posterior predictive) draws.
#' @param max_n optional cap on the number of ensemble members used (subsampled
#'   without replacement) to bound the \eqn{O(m^2)} pairwise term; \code{NULL} uses all.
#' @return the energy score (scalar; lower is better).
#' @seealso \code{\link{variogram_score}}, \code{\link{coherence_backtest}}
#' @export
energy_score <- function(obs, ens, max_n = 2000) {
  obs <- as.numeric(obs); ens <- as.matrix(ens)
  if (ncol(ens) != length(obs)) stop("ncol(ens) must equal length(obs).")
  m <- nrow(ens)
  if (!is.null(max_n) && m > max_n) { ens <- ens[sample.int(m, max_n), , drop = FALSE]; m <- max_n }
  term1 <- mean(sqrt(rowSums((ens - rep(obs, each = m))^2)))
  term2 <- sum(as.matrix(stats::dist(ens))) / (2 * m^2)
  term1 - term2
}

#' Variogram score of a multivariate ensemble forecast
#'
#' @description
#' Proper multivariate scoring rule of order \code{p} (Scheuerer & Hamill 2015):
#' \deqn{VS_p = \sum_{k<l} w_{kl}\,\bigl(|y_k-y_l|^p - \tfrac1m\sum_i |X_{i,k}-X_{i,l}|^p\bigr)^2.}
#' It compares observed and forecast pairwise differences, so it targets the
#' dependence (here, cross-stratum / cross-age) structure directly. Lower is better.
#'
#' @param obs numeric observation vector (length d).
#' @param ens an \code{m x d} ensemble matrix.
#' @param p order of the variogram (default 0.5, robust to heavy tails).
#' @param weights optional \code{d x d} matrix of non-negative pair weights;
#'   \code{NULL} weights all pairs equally.
#' @return the variogram score (scalar; lower is better).
#' @seealso \code{\link{energy_score}}, \code{\link{coherence_backtest}}
#' @export
variogram_score <- function(obs, ens, p = 0.5, weights = NULL) {
  obs <- as.numeric(obs); ens <- as.matrix(ens); d <- length(obs)
  if (ncol(ens) != d) stop("ncol(ens) must equal length(obs).")
  if (!is.null(weights) && !all(dim(weights) == c(d, d))) stop("'weights' must be d x d.")
  vs <- 0
  for (k in 1:(d - 1)) for (l in (k + 1):d) {
    ekl <- mean(abs(ens[, k] - ens[, l])^p)
    okl <- abs(obs[k] - obs[l])^p
    w <- if (is.null(weights)) 1 else weights[k, l]
    vs <- vs + w * (okl - ekl)^2
  }
  vs
}

## ---------------------------------------------------------------------------
## Marginal scores + calibration. The energy/variogram scores above reward the
## JOINT (dependence) structure; CRPS and the log score are their MARGINAL
## complements (per-component, then averaged), and the calibration tools answer
## the orthogonal question "are the predictive intervals the right width?" --
## a sharp, coherent forecast can still be mis-calibrated.
## ---------------------------------------------------------------------------

#' Continuous ranked probability score of an ensemble forecast (marginal)
#'
#' @description
#' Mean per-component CRPS using the fair (unbiased) ensemble estimator
#' \eqn{CRPS = \frac1m\sum_i|x_i-y| - \frac1{2m^2}\sum_{i,j}|x_i-x_j|}, averaged over
#' the \code{d} components. The marginal complement to \code{\link{energy_score}}: it
#' is blind to cross-component dependence (use it alongside, not instead of, the energy
#' score). Lower is better.
#'
#' @param obs numeric observation vector (length d).
#' @param ens an \code{m x d} ensemble matrix.
#' @return mean per-component CRPS (scalar; lower is better).
#' @seealso \code{\link{energy_score}}, \code{\link{logs_sample}}, \code{\link{calibration}}
#' @export
crps_sample <- function(obs, ens) {
  obs <- as.numeric(obs); ens <- as.matrix(ens); m <- nrow(ens)
  if (ncol(ens) != length(obs)) stop("ncol(ens) must equal length(obs).")
  t1 <- colMeans(abs(ens - rep(obs, each = m)))
  # 0.5 E|X-X'| per column via the sorted-order formula (O(m log m), no m x m matrix)
  t2 <- vapply(seq_len(ncol(ens)), function(k) {
    x <- sort(ens[, k]); sum((2 * seq_len(m) - m - 1) * x) / (m^2)
  }, numeric(1))
  mean(t1 - t2)
}

#' Gaussian log score of an ensemble forecast (marginal)
#'
#' @description
#' Mean per-component negative log predictive density under a Gaussian fit to each
#' ensemble margin (mean / sd of the draws). A strictly proper marginal score; more
#' sensitive than CRPS to tail mis-fit but assumes approximate marginal normality.
#' Lower is better.
#'
#' @param obs numeric observation vector (length d).
#' @param ens an \code{m x d} ensemble matrix.
#' @param min_sd floor on the per-margin sd to avoid \eqn{-\infty} on degenerate columns.
#' @return mean per-component negative log density (scalar; lower is better).
#' @seealso \code{\link{crps_sample}}, \code{\link{calibration}}
#' @export
logs_sample <- function(obs, ens, min_sd = 1e-9) {
  obs <- as.numeric(obs); ens <- as.matrix(ens)
  if (ncol(ens) != length(obs)) stop("ncol(ens) must equal length(obs).")
  mu <- colMeans(ens); sdv <- pmax(apply(ens, 2, stats::sd), min_sd)
  mean(-stats::dnorm(obs, mu, sdv, log = TRUE))
}

#' Probability integral transform values of an ensemble forecast (per component)
#'
#' @description
#' \eqn{PIT_k = } fraction of ensemble draws \eqn{\le} the observation in component \eqn{k}.
#' If the forecast is calibrated, pooled PIT values are Uniform(0,1); systematic
#' departures diagnose bias (PIT mass at the ends) or over/under-dispersion (U- or
#' n-shaped histograms). Pool the returned vectors across many forecast cases before
#' assessing uniformity (see \code{\link{calibration}}).
#'
#' @param obs numeric observation vector (length d).
#' @param ens an \code{m x d} ensemble matrix.
#' @return numeric vector of length \code{d} of PIT values in \[0,1\].
#' @seealso \code{\link{calibration}}
#' @export
pit_values <- function(obs, ens) {
  obs <- as.numeric(obs); ens <- as.matrix(ens)
  if (ncol(ens) != length(obs)) stop("ncol(ens) must equal length(obs).")
  vapply(seq_along(obs), function(k) mean(ens[, k] <= obs[k]), numeric(1))
}

#' Central-interval coverage indicators of an ensemble forecast (per component)
#'
#' @param obs numeric observation vector (length d).
#' @param ens an \code{m x d} ensemble matrix.
#' @param level nominal central-interval probability (e.g. 0.9 for the 5\%--95\% interval).
#' @return numeric 0/1 vector of length \code{d}: 1 if the observation falls in the
#'   central \code{level} predictive interval. Pool and average across cases for coverage.
#' @seealso \code{\link{calibration}}
#' @export
interval_coverage <- function(obs, ens, level = 0.9) {
  obs <- as.numeric(obs); ens <- as.matrix(ens)
  a <- (1 - level) / 2
  lo <- apply(ens, 2, stats::quantile, probs = a, names = FALSE)
  hi <- apply(ens, 2, stats::quantile, probs = 1 - a, names = FALSE)
  as.numeric(obs >= lo & obs <= hi)
}

#' Calibration report pooled over many ensemble forecast cases
#'
#' @description
#' Pool PIT values and interval-coverage indicators across a list of forecast cases
#' (e.g. all held-out cells of a backtest) and summarise calibration: empirical vs
#' nominal central-interval coverage, a PIT histogram with a binned chi-square test of
#' Uniform(0,1)-ness, and the mean marginal CRPS / log score. Complements the joint
#' \code{\link{energy_score}} / \code{\link{variogram_score}} with the "right-width?" check.
#'
#' @param obs_list list of observation vectors (one per forecast case).
#' @param ens_list list of \code{m x d} ensemble matrices, aligned with \code{obs_list}.
#' @param levels nominal central-interval probabilities to report coverage for.
#' @param pit_bins number of equal-width PIT histogram bins for the uniformity test.
#' @return a list with \code{coverage} (data.frame: \code{level}, \code{nominal},
#'   \code{empirical}), \code{crps}, \code{logs}, \code{pit} (pooled vector),
#'   \code{pit_hist} (bin counts), and \code{pit_uniformity}
#'   (\code{chisq}, \code{df}, \code{p_value}; large p = consistent with calibration).
#' @seealso \code{\link{crps_sample}}, \code{\link{pit_values}}, \code{\link{coherence_backtest}}
#' @export
calibration <- function(obs_list, ens_list, levels = c(0.5, 0.8, 0.9, 0.95), pit_bins = 10) {
  if (length(obs_list) != length(ens_list)) stop("obs_list and ens_list must be the same length.")
  pit <- unlist(Map(pit_values, obs_list, ens_list), use.names = FALSE)
  cov <- vapply(levels, function(L)
    mean(unlist(Map(function(o, e) interval_coverage(o, e, L), obs_list, ens_list), use.names = FALSE)),
    numeric(1))
  crps <- mean(vapply(seq_along(obs_list), function(i) crps_sample(obs_list[[i]], ens_list[[i]]), numeric(1)))
  logs <- mean(vapply(seq_along(obs_list), function(i) logs_sample(obs_list[[i]], ens_list[[i]]), numeric(1)))
  br <- seq(0, 1, length.out = pit_bins + 1)
  counts <- tabulate(findInterval(pmin(pmax(pit, 0), 1), br, rightmost.closed = TRUE), nbins = pit_bins)
  expct <- length(pit) / pit_bins
  chisq <- if (expct > 0) sum((counts - expct)^2 / expct) else NA_real_
  pval <- if (is.na(chisq)) NA_real_ else stats::pchisq(chisq, pit_bins - 1, lower.tail = FALSE)
  list(coverage = data.frame(level = levels, nominal = levels, empirical = cov),
       crps = crps, logs = logs, pit = pit, pit_hist = counts,
       pit_uniformity = list(chisq = chisq, df = pit_bins - 1, p_value = pval))
}
