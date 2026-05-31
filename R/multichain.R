## ===========================================================================
## Multiple independent chains for the coherent / competing-cause / cascade
## samplers (D2). The samplers run one chain; convergence assessment needs
## several from dispersed seeds so bamp_diagnostics() can compute the TRUE
## between-chain split-Rhat (not just the within-chain split). run_chains()
## fits them (optionally forked in parallel); combine_chains() pools the kept
## draws into one fit object for prediction.
##
## Reproducibility does not depend on the parallel backend: each chain seeds the
## sampler's own set.seed(seed), so a given (FUN, args, seeds) is deterministic.
## ===========================================================================

#' Run several independent chains of a bamp sampler
#'
#' @description
#' Fit a coherent / competing-cause / cascade sampler several times from dispersed
#' seeds (optionally forked in parallel) so convergence can be judged with the true
#' between-chain \eqn{\hat R} via \code{\link{bamp_diagnostics}}.
#'
#' @param FUN a fitting function, e.g. \code{\link{bamp_coherent}}, \code{\link{bamp_multicause}}
#'   or \code{\link{bamp_cascade}} (it must take a \code{seed} argument).
#' @param args named list of arguments passed to \code{FUN} (do not include \code{seed}).
#' @param chains number of chains.
#' @param seeds integer seeds, one per chain; default \code{seq_len(chains)}.
#' @param parallel if \code{TRUE}, fork the chains with \code{parallel::mclapply} (no effect on Windows).
#' @param cores number of worker processes; default \code{min(chains, detectCores - 1)}.
#'
#' @return a list of class \code{bamp_multichain} with \code{chains} (the list of fits),
#'   \code{n_chains} and \code{seeds}. Pass it to \code{\link{bamp_diagnostics}} for
#'   between-chain \eqn{\hat R}, or \code{\link{combine_chains}} to pool draws for prediction.
#' @seealso \code{\link{bamp_diagnostics}}, \code{\link{combine_chains}}
#' @export
run_chains <- function(FUN, args = list(), chains = 4, seeds = NULL,
                       parallel = TRUE, cores = NULL) {
  FUN <- match.fun(FUN)
  if (chains < 1) stop("'chains' must be >= 1.")
  if (is.null(seeds)) seeds <- seq_len(chains)
  if (length(seeds) != chains) stop("length(seeds) must equal 'chains'.")
  one <- function(s) do.call(FUN, c(args, list(seed = s)))
  if (parallel && chains > 1 && .Platform$OS.type == "unix") {
    if (is.null(cores)) cores <- max(1L, min(chains, parallel::detectCores() - 1L))
    fits <- parallel::mclapply(seeds, one, mc.cores = cores)
    if (any(vapply(fits, function(f) inherits(f, "try-error") || is.null(f), logical(1))))
      stop("at least one chain failed under mclapply; rerun with parallel = FALSE to see the error.")
  } else {
    fits <- lapply(seeds, one)
  }
  structure(list(chains = fits, n_chains = chains, seeds = seeds), class = "bamp_multichain")
}

#' @method print bamp_multichain
#' @export
print.bamp_multichain <- function(x, ...) {
  cl <- class(x$chains[[1]])[1]
  cat(sprintf("bamp_multichain: %d chains of '%s' (seeds %s)\n",
              x$n_chains, cl, paste(x$seeds, collapse = ", ")))
  cat("  -> bamp_diagnostics() for between-chain Rhat; combine_chains() to pool draws.\n")
  invisible(x)
}

## stack one samples element (vector / matrix / array) along the iteration axis
.mc_stack <- function(parts) {
  a <- parts[[1]]
  if (is.null(dim(a))) do.call(c, parts)
  else if (length(dim(a)) == 2L) do.call(rbind, parts)
  else do.call(function(...) abind::abind(..., along = 1), parts)
}

#' Pool multiple chains into a single fit for prediction
#'
#' @description
#' Concatenate the kept draws of every chain along the iteration axis, returning a single
#' fit object (same class as the individual chains) usable by the \code{predict_*} functions
#' with the combined draws. For \code{\link{bamp_multicause}} the all-cause \code{total}
#' sub-fit is taken from the first chain (the field-share draws are fully pooled).
#'
#' @param mc a \code{bamp_multichain} from \code{\link{run_chains}}.
#' @return a single fit object (e.g. \code{apc_coherent} / \code{apc_multicause}) with pooled draws.
#' @seealso \code{\link{run_chains}}, \code{\link{bamp_diagnostics}}
#' @export
combine_chains <- function(mc) {
  if (!inherits(mc, "bamp_multichain")) stop("'mc' must come from run_chains().")
  fits <- mc$chains; base <- fits[[1]]
  nms <- names(base$samples)
  base$samples <- stats::setNames(lapply(nms, function(nm)
    .mc_stack(lapply(fits, function(f) f$samples[[nm]]))), nms)
  base
}
