## ===========================================================================
## Sparse one-block GMRF draw for the coherent/competing-cause samplers.
##
## The dense reference samplers form Q = X' Omega X + Prec densely and call
## chol(Q) -- O(P^3) per sweep, the wall that stops many causes/strata. Here the
## same draw is done with a SPARSE Cholesky (Matrix package): the design X is
## sparse (a few non-zeros per cell), so crossprod(X, Omega X) is the O(n) STAGE-A
## block assembly for free, and the RW/kronecker prior blocks are banded.
##
## .pg_draw_block_sparse mirrors the dense .pg_draw_block (Rue & Held conditioning
## by Kriging) exactly: x ~ N(Q^{-1} b, Q^{-1}) subject to A x = 0.
##  * REFACTORISES FRESH each call (no cached factor), so a changing Omega cannot
##    invalidate a cached sparsity pattern -- the whole class of update-API bugs
##    disappears, at the cost of redoing the (cheap) symbolic analysis each sweep.
##  * RNG: exactly one rnorm(P). The unconstrained draw is P' L'^{-1} z (system
##    "Lt" then "Pt"), which has covariance Q^{-1} (verified by Monte Carlo) --
##    NOT system "A" on z (that would give Q^{-2}). The fill-reducing permutation
##    reorders which normal hits which coordinate, so draws are not bit-identical
##    to the dense path but share the distribution: parity is by posterior
##    mean/ESS, not by draw.
## ===========================================================================
.pg_draw_block_sparse <- function(Qs, bvec, A, tA = NULL) {
  L <- Matrix::Cholesky(Matrix::forceSymmetric(Qs), perm = TRUE, LDL = FALSE, super = FALSE)
  mu <- as.numeric(Matrix::solve(L, bvec, system = "A"))                  # Q^{-1} b
  z <- stats::rnorm(length(bvec))
  x <- mu + as.numeric(Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "Pt"))
  if (is.null(A)) return(x)
  if (is.null(tA)) tA <- t(A)
  W <- as.matrix(Matrix::solve(L, tA, system = "A"))                      # Q^{-1} A'  (P x k)
  as.numeric(x - W %*% solve(A %*% W, as.numeric(A %*% x)))               # condition on A x = 0
}
