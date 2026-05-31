#' Reconcile cause-specific forecasts to be coherent with all-cause (hazard scale)
#'
#' @description
#' Non-invasive post-processor that takes independently-fitted cause-specific forecasts and an
#' all-cause forecast and reconciles them so the cause-specific \emph{cumulative hazards} sum to
#' the all-cause hazard. This is the zero-engine-change fallback of the coherent-forecasting
#' design (Phase 0): you keep using ordinary \code{\link{bamp}}/\code{\link{predict_apc}} once per
#' cause and once for all-cause, then call \code{reconcile_apc()} to make them add up.
#'
#' Reconciliation is done on the \strong{hazard} scale (\eqn{h=-\log(1-p)}), because cause-specific
#' hazards are additive across competing causes whereas risks are not -- reconciling probabilities
#' would be coherence on the wrong scale. It is performed \strong{marginally} (on the posterior
#' mean and on each quantile band), not per draw, because the cause and all-cause fits have
#' independent MCMC streams and pairing them by draw index would fabricate posterior dependence.
#'
#' @param total a prediction from \code{predict_apc(all_cause_fit, ..., hazard=TRUE)} (must contain
#'   \code{samples$hazard}).
#' @param causes a named list of predictions from \code{predict_apc(cause_fit, ..., hazard=TRUE)},
#'   one per cause, all on the same \code{[period, agegroup]} grid as \code{total}.
#' @param weights optional positive numeric vector of length \code{1 + length(causes)} giving the
#'   diagonal of the reconciliation weight matrix \eqn{W} (relative forecast-error variances, in the
#'   order \code{total, causes...}). Default \code{NULL} uses ordinary least squares (\eqn{W=I}).
#' @param quantiles vector of quantiles to recompute for the reconciled hazards.
#'
#' @return list with the reconciled all-cause forecast (\code{total}) and reconciled cause-specific
#'   forecasts (\code{causes}, named), each holding \code{hazard_mean} (\code{[period, agegroup]})
#'   and \code{hazard} quantiles (\code{[quantile, period, agegroup]}); the summing matrix \code{S}
#'   and projection \code{P} used; and \code{coherence_maxerr}, the largest absolute discrepancy
#'   between the summed reconciled cause hazards and the reconciled total (0 by construction).
#'
#' @details
#' With base hazards \eqn{b=(h_0, h_1, \dots, h_C)} per cell (all-cause then causes) and summing
#' matrix \eqn{S=[\mathbf{1}_C^\top;\, I_C]}, the GLS reconciliation projects onto the coherent
#' subspace via \eqn{P = S (S^\top W^{-1} S)^{-1} S^\top W^{-1}}. The map \eqn{P} is the same for
#' every cell, so it is applied to the mean and to each quantile surface. Reconciled cause hazards
#' are floored at 0 (a projection can overshoot slightly negative in rare cells) and the reconciled
#' total is then re-summed from the floored causes, which keeps the result both non-negative and
#' coherent. For the principled single-posterior alternative (a multinomial/stick-breaking APC that
#' is coherent by construction) see the design note \code{docs/coherent-forecasting.md}.
#'
#' @seealso \code{\link{predict_apc}}, \code{\link{bamp_strata}}
#' @export
reconcile_apc <- function(total, causes, weights = NULL, quantiles = c(0.05, 0.5, 0.95))
{
  if (is.null(total$samples$hazard))
    stop("'total' must be a predict_apc() result with hazards: call predict_apc(..., hazard=TRUE).")
  if (!is.list(causes) || length(causes) < 1L)
    stop("'causes' must be a non-empty (named) list of predict_apc() results.")
  if (is.null(names(causes)) || any(names(causes) == ""))
    names(causes) <- paste0("cause", seq_along(causes))
  if (any(vapply(causes, function(x) is.null(x$samples$hazard), logical(1))))
    stop("Each element of 'causes' must contain hazards: call predict_apc(..., hazard=TRUE).")

  C <- length(causes)
  hz_list <- c(list(total), causes)                       # length C+1, total first

  # All on the same [period, agegroup] grid.
  grid <- dim(hz_list[[1]]$samples$hazard)[1:2]
  if (any(vapply(hz_list, function(x) !identical(dim(x$samples$hazard)[1:2], grid), logical(1))))
    stop("'total' and all 'causes' must share the same [period, agegroup] grid.")
  np <- grid[1]; na <- grid[2]; ncell <- np * na

  # Summing matrix and projection onto the coherent subspace.
  S <- rbind(matrix(1, 1, C), diag(C))                    # (C+1) x C
  if (is.null(weights)) {
    Winv <- diag(C + 1L)
  } else {
    if (length(weights) != C + 1L || any(weights <= 0))
      stop("'weights' must be a positive vector of length 1 + number of causes (total first).")
    Winv <- diag(1 / weights)
  }
  G <- solve(t(S) %*% Winv %*% S) %*% t(S) %*% Winv        # C x (C+1)
  P <- S %*% G                                            # (C+1) x (C+1) projection

  # Reconcile a stack of [C+1, ncell] hazards: project, floor causes at 0, re-sum total.
  reconcile_stack <- function(M) {                        # M: (C+1) x ncell, rows = total, causes
    R <- P %*% M
    bottom <- pmax(R[-1, , drop = FALSE], 0)              # non-negative cause hazards
    rbind(colSums(bottom), bottom)                        # coherent total = sum of causes
  }

  # Posterior means.
  meanmat <- vapply(hz_list, function(x) apply(x$samples$hazard, 1:2, mean),
                    matrix(0, np, na))                    # [np, na, C+1]
  meanmat <- matrix(aperm(meanmat, c(3, 1, 2)), C + 1L, ncell)
  rec_mean <- reconcile_stack(meanmat)

  # Quantile surfaces: project each quantile level (marginal approximation).
  qsurf <- lapply(hz_list, function(x) apply(x$samples$hazard, 1:2, stats::quantile, quantiles))
  nq <- length(quantiles)
  rec_q <- array(0, c(C + 1L, nq, np, na))
  for (qi in seq_len(nq)) {
    M <- matrix(vapply(qsurf, function(a) as.vector(a[qi, , ]), numeric(ncell)), ncol = C + 1L)
    # reconcile_stack returns (C+1) x ncell with cells column-major (period fastest),
    # which array() folds straight back into [C+1, period, age] with no transpose.
    rec_q[, qi, , ] <- array(reconcile_stack(t(M)), c(C + 1L, np, na))
  }

  unpack <- function(idx) list(
    hazard_mean = matrix(rec_mean[idx, ], np, na),
    hazard = array(rec_q[idx, , , ], c(nq, np, na),
                   dimnames = list(paste0(quantiles * 100, "%"), NULL, NULL))
  )

  out <- list(total = unpack(1L))
  out$causes <- stats::setNames(lapply(seq_len(C), function(c) unpack(c + 1L)), names(causes))
  out$S <- S; out$P <- P
  out$coherence_maxerr <- max(abs(Reduce(`+`, lapply(out$causes, `[[`, "hazard_mean")) -
                                 out$total$hazard_mean))
  out
}
