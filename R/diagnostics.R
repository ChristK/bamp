## ===========================================================================
## Convergence diagnostics for the coherent / competing-cause / cascade samplers.
##
## The new one-block PG samplers (bamp_coherent, bamp_multicause, bamp_cascade,
## bamp_sex_cascade) store their kept draws under fit$samples as a named list of
## [nkeep x ...] arrays (iteration is always the FIRST dimension). Running 30
## diseases x strata, you cannot eyeball every chain, so bamp_diagnostics()
## flattens every stored parameter to a scalar chain and reports, per parameter:
##   * rank-normalised split-Rhat (Vehtari et al. 2021) -- robust to non-normal,
##     heavy-tailed posteriors; splits the single chain in half so it still
##     detects within-run non-stationarity without needing multiple chains;
##   * bulk effective sample size (coda::effectiveSize).
## It then flags parameters breaching Rhat/ESS thresholds. Multi-chain Rhat
## (true between-chain) arrives with D2; the split-Rhat here is the single-chain
## stand-in and uses the same machinery.
## ===========================================================================

## --- single-chain rank-normalised split-Rhat + ESS ------------------------
.diag_split_halves <- function(x) {                       # two non-overlapping halves
  n <- length(x); h <- n %/% 2L
  list(x[seq_len(h)], x[(n - h + 1L):n])
}
.diag_zscale <- function(x) {                             # rank-normalise (to normal scores)
  r <- rank(x, ties.method = "average")
  stats::qnorm((r - 0.5) / length(r))
}
.diag_rhat_from_chains <- function(chains) {              # Gelman split-Rhat on equal-length chains
  n <- length(chains[[1]])
  if (n < 2L) return(NA_real_)
  means <- vapply(chains, mean, numeric(1)); vars <- vapply(chains, stats::var, numeric(1))
  W <- mean(vars); if (!is.finite(W) || W <= 0) return(NA_real_)
  B <- n * stats::var(means)
  sqrt(((n - 1) / n * W + B / n) / W)
}
## rank-normalised split-Rhat + bulk ESS. `mat` is a vector (one chain, split in
## half) or an [n x C] matrix (C chains -> true between-chain Rhat, each split in half).
.diag_rhat_ess <- function(mat) {
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
  subch <- list()
  for (c in seq_len(ncol(mat))) {                         # split each chain into 2 halves
    x <- mat[, c]; x <- x[is.finite(x)]; h <- length(x) %/% 2L
    if (h < 2L) next
    subch[[length(subch) + 1L]] <- x[seq_len(h)]
    subch[[length(subch) + 1L]] <- x[(length(x) - h + 1L):length(x)]
  }
  if (length(subch) < 2L) return(c(rhat = NA_real_, ess = NA_real_))
  n <- min(vapply(subch, length, integer(1)))
  subch <- lapply(subch, function(s) s[seq_len(n)])
  pooled <- unlist(subch)
  if (stats::var(pooled) <= 0) return(c(rhat = NA_real_, ess = NA_real_))
  z <- .diag_zscale(pooled)                               # rank-normalise across all split chains
  zch <- split(z, rep(seq_along(subch), each = n))
  rhat <- .diag_rhat_from_chains(zch)
  ess <- tryCatch(sum(vapply(seq_len(ncol(mat)),
           function(c) as.numeric(coda::effectiveSize(mat[, c])), numeric(1))),
           error = function(e) NA_real_)
  c(rhat = rhat, ess = ess)
}

## flatten one fit$samples list to a [nkeep x P] matrix of scalar chains.
## iteration is dim 1; everything else becomes a column with an indexed name.
.diag_param_matrix <- function(samples, prefix = "") {
  if (!is.list(samples) || !length(samples)) return(NULL)
  lens <- vapply(samples, function(a) if (is.null(dim(a))) length(a) else dim(a)[1], integer(1))
  nkeep <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])   # modal first-dim = nkeep
  cols <- list()
  for (nm in names(samples)) {
    a <- samples[[nm]]
    d <- dim(a)
    if (is.null(d)) {                                     # plain vector, one chain
      if (length(a) != nkeep) next
      cols[[paste0(prefix, nm)]] <- as.numeric(a)
    } else {
      if (d[1] != nkeep) next
      mat <- matrix(as.numeric(a), nrow = nkeep)          # collapse dims 2..k column-major
      idx <- if (length(d) > 1) do.call(expand.grid, lapply(d[-1], seq_len)) else NULL
      for (j in seq_len(ncol(mat))) {
        tag <- if (is.null(idx)) "" else paste0("[", paste(idx[j, ], collapse = ","), "]")
        cols[[paste0(prefix, nm, tag)]] <- mat[, j]
      }
    }
  }
  if (!length(cols)) return(NULL)
  do.call(cbind, cols)
}

## collect every fit-level $samples in a (possibly nested, e.g. cascade) fit
.diag_collect <- function(obj, path = "") {
  out <- list()
  if (is.list(obj)) {
    s <- obj[["samples"]]
    if (is.list(s) && length(s) &&
        any(vapply(s, function(a) is.numeric(a) && !is.null(dim(a)), logical(1))))
      out[[if (nzchar(path)) path else "model"]] <- s
    for (nm in names(obj)) if (!identical(nm, "samples"))
      out <- c(out, .diag_collect(obj[[nm]], if (nzchar(path)) paste(path, nm, sep = "/") else nm))
  }
  out
}

#' Convergence diagnostics for coherent / competing-cause / cascade fits
#'
#' @description
#' Flatten every stored parameter draw of a \code{bamp_coherent}, \code{bamp_multicause},
#' \code{bamp_cascade} or \code{bamp_sex_cascade} fit to a scalar chain and report
#' rank-normalised split-\eqn{\hat R} (Vehtari et al. 2021) and bulk effective sample
#' size per parameter, then flag any breaching the thresholds. Nested fits (cascade
#' sub-models) are traversed and labelled by path. Designed for batch use across many
#' disease/stratum fits: read \code{$summary} for the worst \eqn{\hat R} / smallest ESS
#' and \code{$ok} for a single pass/fail.
#'
#' @param fit a fitted object with draws under \code{$samples} (possibly nested).
#' @param rhat_max flag parameters with split-\eqn{\hat R} above this (default 1.01).
#' @param ess_min flag parameters with bulk ESS below this (default 400).
#' @param top number of worst parameters (by \eqn{\hat R}, then ESS) to return in \code{$worst}.
#'
#' @return a list of class \code{bamp_diagnostics} with: \code{by_param} (data.frame:
#'   \code{block}, \code{parameter}, \code{rhat}, \code{ess}, \code{flagged}); \code{worst}
#'   (the \code{top} most suspect rows); \code{summary} (\code{n_params}, \code{max_rhat},
#'   \code{min_ess}, \code{n_flagged}); and \code{ok} (logical: nothing flagged).
#' @seealso \code{\link{bamp_traceplot}}, \code{\link{calibration}}
#' @export
bamp_diagnostics <- function(fit, rhat_max = 1.01, ess_min = 400, top = 10) {
  rows <- list()
  if (inherits(fit, "bamp_multichain")) {                 # D2: true between-chain Rhat
    cblocks <- lapply(fit$chains, .diag_collect)
    for (bn in names(cblocks[[1]])) {
      Ms <- lapply(cblocks, function(cb) .diag_param_matrix(cb[[bn]]))
      Ms <- Filter(Negate(is.null), Ms)
      if (length(Ms) < 1L) next
      pars <- colnames(Ms[[1]]); ncol1 <- length(pars)
      re <- t(vapply(seq_len(ncol1), function(j) {
        mat <- vapply(Ms, function(M) M[, j], numeric(nrow(Ms[[1]])))   # nkeep x chains
        .diag_rhat_ess(mat)
      }, numeric(2)))
      rows[[bn]] <- data.frame(block = bn, parameter = pars,
                               rhat = re[, 1], ess = re[, 2], stringsAsFactors = FALSE)
    }
  } else {
    blocks <- .diag_collect(fit)
    if (!length(blocks)) stop("no fit-level $samples found; is this a fitted bamp object?")
    for (bn in names(blocks)) {
      M <- .diag_param_matrix(blocks[[bn]])
      if (is.null(M)) next
      re <- t(apply(M, 2, .diag_rhat_ess))
      rows[[bn]] <- data.frame(block = bn, parameter = colnames(M),
                               rhat = re[, "rhat"], ess = re[, "ess"],
                               stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("no fit-level $samples found; is this a fitted bamp object?")
  by_param <- do.call(rbind, rows); rownames(by_param) <- NULL
  by_param$flagged <- (is.finite(by_param$rhat) & by_param$rhat > rhat_max) |
                      (is.finite(by_param$ess)  & by_param$ess  < ess_min)
  ord <- order(-replace(by_param$rhat, !is.finite(by_param$rhat), -Inf),
               replace(by_param$ess, !is.finite(by_param$ess), Inf))
  worst <- utils::head(by_param[ord, ], top)
  summary <- data.frame(
    n_params  = nrow(by_param),
    max_rhat  = if (any(is.finite(by_param$rhat))) max(by_param$rhat, na.rm = TRUE) else NA_real_,
    min_ess   = if (any(is.finite(by_param$ess)))  min(by_param$ess,  na.rm = TRUE) else NA_real_,
    n_flagged = sum(by_param$flagged))
  structure(list(by_param = by_param, worst = worst, summary = summary,
                 ok = summary$n_flagged == 0L,
                 thresholds = c(rhat_max = rhat_max, ess_min = ess_min)),
            class = "bamp_diagnostics")
}

#' @method print bamp_diagnostics
#' @export
print.bamp_diagnostics <- function(x, ...) {
  s <- x$summary
  cat(sprintf("bamp convergence diagnostics: %d parameters across %d block(s)\n",
              s$n_params, length(unique(x$by_param$block))))
  cat(sprintf("  max split-Rhat = %.4f (<= %.3f) | min ESS = %.0f (>= %.0f) | flagged = %d\n",
              s$max_rhat, x$thresholds["rhat_max"], s$min_ess, x$thresholds["ess_min"], s$n_flagged))
  cat(if (x$ok) "  OK: no parameter breaches the thresholds.\n" else
      "  WARNING: some parameters breach the thresholds (see $worst).\n")
  if (!x$ok) { cat("  worst:\n"); print(utils::head(x$worst, 5), row.names = FALSE) }
  invisible(x)
}

#' Trace plots for the worst-mixing parameters of a bamp fit
#'
#' @description
#' Diagnostic trace + running-mean plots for the parameters flagged by
#' \code{\link{bamp_diagnostics}} (by default the worst split-\eqn{\hat R}). A quick
#' visual companion to the numeric diagnostics; a well-mixed chain looks like white
#' noise around a flat running mean.
#'
#' @param fit a fitted bamp object (as for \code{\link{bamp_diagnostics}}).
#' @param which \code{"worst"} to plot the most suspect parameters, or a character vector
#'   of \code{"block:parameter"} names to plot specific ones.
#' @param n number of parameters to plot when \code{which = "worst"}.
#' @param diag a precomputed \code{\link{bamp_diagnostics}} result (optional; recomputed if NULL).
#' @return invisibly, the names of the parameters plotted.
#' @seealso \code{\link{bamp_diagnostics}}
#' @export
bamp_traceplot <- function(fit, which = "worst", n = 4, diag = NULL) {
  blocks <- .diag_collect(fit)
  mats <- lapply(blocks, .diag_param_matrix)
  getcol <- function(block, par) {
    M <- mats[[block]]; if (is.null(M) || !(par %in% colnames(M))) return(NULL); M[, par]
  }
  if (identical(which, "worst")) {
    if (is.null(diag)) diag <- bamp_diagnostics(fit)
    sel <- utils::head(diag$worst, n)
    picks <- Map(function(b, p) list(block = b, par = p), sel$block, sel$parameter)
  } else {
    picks <- lapply(which, function(s) { sp <- strsplit(s, ":", fixed = TRUE)[[1]]
      list(block = if (length(sp) > 1) sp[1] else names(blocks)[1], par = utils::tail(sp, 1)) })
  }
  picks <- Filter(function(p) !is.null(getcol(p$block, p$par)), picks)
  if (!length(picks)) stop("none of the requested parameters were found.")
  op <- graphics::par(mfrow = c(length(picks), 1), mar = c(3, 4, 2, 1)); on.exit(graphics::par(op))
  labs <- character(length(picks))
  for (i in seq_along(picks)) {
    p <- picks[[i]]; v <- getcol(p$block, p$par); labs[i] <- paste(p$block, p$par, sep = ":")
    graphics::plot(v, type = "l", col = "grey40", xlab = "", ylab = p$par, main = labs[i])
    graphics::lines(cumsum(v) / seq_along(v), col = "red", lwd = 2)
  }
  invisible(labs)
}
