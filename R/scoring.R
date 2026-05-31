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
