## ===========================================================================
## Polya-Gamma Gibbs engine for the Bayesian age-period-cohort model.
##
## Binomial logit likelihood with RW1/RW2 intrinsic GMRF priors on the age,
## period and cohort effects.  Polya-Gamma data augmentation (Polson, Scott &
## Windle 2013, JASA) renders the whole latent field conditionally Gaussian, so
## the intercept and all three effects are drawn JOINTLY in one exact Gibbs step
## (the "one-block" GMRF sampler, Rue & Held 2005).  Joint sampling removes the
## cross-block autocorrelation that cripples one-at-a-time updates when the data
## are highly informative (large populations) -- the regime typical of
## incidence/mortality data.  There is no Metropolis tuning, no acceptance-rate
## restart and no chain pruning.
##
## Intrinsic priors are scaled to unit generalised variance (Sorbye & Rue 2014)
## so a single hyper-prior is portable across data sets of different size/grid.
##
## Identifiability: each effect is constrained to sum to zero (level).  For RW2,
## the otherwise improper shared linear-trend (drift) direction is pinned by a
## single zero-slope constraint on the period effect.  Constraints are imposed
## jointly by conditioning by Kriging.  All of this is internal; the engine
## returns the same objects as the C++ path.
## ===========================================================================

## --- intrinsic RW structure matrix K (precision = kappa * K) ---------------
.pg_Kmat <- function(L, order) {
  D <- diff(diag(L), differences = order)        # (L-order) x L difference op
  crossprod(D)                                    # D'D
}

## --- Sorbye & Rue (2014) scaling factor: multiply K by this so that the
##     geometric mean of the intrinsic marginal variances equals 1 -----------
.pg_scale <- function(K, order) {
  e <- eigen(K, symmetric = TRUE)
  keep <- e$values > max(e$values) * 1e-9         # drop the (order) null eigenvalues
  V <- e$vectors[, keep, drop = FALSE]
  Sigma <- V %*% diag(1 / e$values[keep], sum(keep)) %*% t(V)   # generalised inverse
  exp(mean(log(diag(Sigma))))
}

## --- vectorised Polya-Gamma draw, PG(b, c).  Exact mean/variance with a
##     Normal approximation (CLT in the count b; near-exact for the large counts
##     of incidence/mortality data).  Analytic limit handles c ~ 0. -----------
.pg_rpg <- function(b, c) {
  ac <- abs(c)
  m <- numeric(length(c)); v <- numeric(length(c))
  small <- ac < 1e-4
  m[small] <- 0.25; v[small] <- 1 / 24
  cs <- ac[!small]; th <- tanh(cs / 2); sech2 <- 1 - th^2
  m[!small] <- th / (2 * cs)
  v[!small] <- ((2 / cs) * th - sech2) / (4 * cs^2)
  out <- rnorm(length(c), b * m, sqrt(pmax(b * v, .Machine$double.eps)))
  pmax(out, 1e-9)                                 # omega must be > 0
}

## --- sample x ~ N(Q^{-1} bvec, Q^{-1}) subject to A x = 0 (conditioning by
##     Kriging, Rue & Held 2005, Alg. 2.6).  A is k x P (or NULL). ------------
.pg_draw_block <- function(Q, bvec, A) {
  R <- chol(Q)                                    # Q = R'R, R upper triangular
  mu <- backsolve(R, forwardsolve(t(R), bvec))
  x  <- mu + backsolve(R, rnorm(length(bvec)))
  if (is.null(A)) return(x)
  W  <- backsolve(R, forwardsolve(t(R), t(A)))    # Q^{-1} A'  (P x k)
  as.numeric(x - W %*% solve(A %*% W, A %*% x))
}

## --- one MCMC chain (joint one-block update) -------------------------------
.bamp_pg_chain <- function(Y, N, ord_a, ord_p, ord_c, ppa,
                           hyper, n_iter, burn_in, thin, seed,
                           prior_scale = TRUE) {
  set.seed(seed)
  I <- nrow(Y); J <- ncol(Y); K <- ppa * (I - 1) + J
  has_a <- ord_a > 0; has_p <- ord_p > 0; has_c <- ord_c > 0
  cohidx <- outer(1:I, 1:J, function(i, j) (I - i) * ppa + j)   # cohort index per cell
  cohv <- as.integer(cohidx)
  ymh <- Y - N / 2
  ## precomputed linear indices for the age<->cohort and period<->cohort
  ## coupling blocks: within a row (or column) each cell hits a distinct cohort,
  ## so these are pure indexed placements (no aggregation needed in the loop)
  ivec <- rep(1:I, times = J); jvec <- rep(1:J, each = I)
  idxT <- (cohv - 1L) * I + ivec        # into an I x K matrix
  idxP <- (cohv - 1L) * J + jvec        # into a J x K matrix

  ## structure matrices + scaling + hyperparameters
  setup <- function(L, ord, hp) {
    Km <- .pg_Kmat(L, ord); s <- if (prior_scale) .pg_scale(Km, ord) else 1
    list(K = Km * s, rank = L - ord, a = hp[1], b = hp[2])
  }
  if (has_a) sa <- setup(I, ord_a, hyper$age)
  if (has_p) sp <- setup(J, ord_p, hyper$period)
  if (has_c) sc <- setup(K, ord_c, hyper$cohort)

  ## parameter index layout in the joint vector beta = (mu, theta, phi, psi)
  ia <- 1L
  ith <- if (has_a) 1L + (1:I)         else integer(0)
  off <- 1L + (if (has_a) I else 0L)
  iph <- if (has_p) off + (1:J)        else integer(0)
  off <- off + (if (has_p) J else 0L)
  ips <- if (has_c) off + (1:K)        else integer(0)
  P   <- 1L + (if (has_a) I else 0L) + (if (has_p) J else 0L) + (if (has_c) K else 0L)

  ## fixed constraint matrix A (sum-to-zero per effect; RW2 period zero-slope)
  Arows <- list()
  mkrow <- function(idx, val) { r <- numeric(P); r[idx] <- val; r }
  if (has_a) Arows <- c(Arows, list(mkrow(ith, 1)))
  if (has_p) Arows <- c(Arows, list(mkrow(iph, 1)))
  if (has_c) Arows <- c(Arows, list(mkrow(ips, 1)))
  ## RW2's second-difference prior does not penalise a linear trend, so in the
  ## full age-period-cohort model the shared drift direction is improper.  Pin
  ## it with a single zero-slope constraint on the period effect (a reporting
  ## convention; fitted rates, curvatures and net drift are unaffected).  RW1
  ## penalises the trend (weakly identified), so it is left free there.
  if (has_a && has_p && has_c && ord_p == 2)
    Arows <- c(Arows, list(mkrow(iph, (1:J) - mean(1:J))))
  A <- if (length(Arows)) do.call(rbind, Arows) else NULL

  ## starting values
  p0 <- sum(Y) / sum(N); mu <- log(p0 / (1 - p0))
  theta <- numeric(I); phi <- numeric(J); psi <- numeric(K)
  kappa <- lambda <- ny <- 1

  nkeep <- length(seq(burn_in + 1, n_iter, by = thin))
  out_theta <- matrix(0, nkeep, I); out_phi <- matrix(0, nkeep, J); out_psi <- matrix(0, nkeep, K)
  out_kap <- out_lam <- out_ny <- out_mu <- out_dev <- numeric(nkeep)
  ksi_sum <- matrix(0, I, J); ksi_n <- 0L; keep <- 0L

  ## constant pieces of b
  b_mu0 <- sum(ymh)
  b_th0 <- if (has_a) rowSums(ymh) else NULL
  b_ph0 <- if (has_p) colSums(ymh) else NULL
  b_ps0 <- if (has_c) { v <- numeric(K); a <- rowsum(as.numeric(ymh), cohv); v[as.integer(rownames(a))] <- a; v } else NULL

  for (it in 1:n_iter) {
    eta <- mu + outer(theta, phi, "+") + matrix(psi[cohidx], I, J)
    omega <- matrix(.pg_rpg(as.numeric(N), as.numeric(eta)), I, J)

    ## ---- assemble joint precision Q (P x P) and rhs b ----
    Q <- matrix(0, P, P); b <- numeric(P)
    omv <- as.numeric(omega)
    rs <- rowSums(omega); cs <- colSums(omega); tot <- sum(omega)
    TPsi <- PPsi <- NULL
    if (has_a && has_c) { TPsi <- matrix(0, I, K); TPsi[idxT] <- omv }
    if (has_p && has_c) { PPsi <- matrix(0, J, K); PPsi[idxP] <- omv }
    if (has_c) csum <- if (!is.null(TPsi)) colSums(TPsi)
                       else if (!is.null(PPsi)) colSums(PPsi)
                       else { v <- numeric(K); ag <- rowsum(omv, cohv); v[as.integer(rownames(ag))] <- ag; v }
    Q[ia, ia] <- tot; b[ia] <- b_mu0
    if (has_a) { Q[ia, ith] <- rs; Q[ith, ia] <- rs
                 Q[ith, ith] <- diag(rs, I) + kappa * sa$K; b[ith] <- b_th0 }
    if (has_p) { Q[ia, iph] <- cs; Q[iph, ia] <- cs
                 Q[iph, iph] <- diag(cs, J) + lambda * sp$K; b[iph] <- b_ph0 }
    if (has_c) { Q[ia, ips] <- csum; Q[ips, ia] <- csum
                 Q[ips, ips] <- diag(csum, K) + ny * sc$K; b[ips] <- b_ps0 }
    if (has_a && has_p) { Q[ith, iph] <- omega; Q[iph, ith] <- t(omega) }
    if (!is.null(TPsi)) { Q[ith, ips] <- TPsi; Q[ips, ith] <- t(TPsi) }
    if (!is.null(PPsi)) { Q[iph, ips] <- PPsi; Q[ips, iph] <- t(PPsi) }
    ## ridge to make the (intrinsically rank-deficient) Q numerically PD; the
    ## constraints below remove the corresponding null directions
    diag(Q) <- diag(Q) + 1e-6 * mean(diag(Q))

    beta <- .pg_draw_block(Q, b, A)
    mu <- beta[ia]
    if (has_a) theta <- beta[ith]
    if (has_p) phi   <- beta[iph]
    if (has_c) psi   <- beta[ips]

    ## ---- precisions: centred (sufficient) Gamma full conditionals ----
    if (has_a) kappa  <- rgamma(1, sa$a + sa$rank / 2, sa$b + 0.5 * as.numeric(theta %*% sa$K %*% theta))
    if (has_p) lambda <- rgamma(1, sp$a + sp$rank / 2, sp$b + 0.5 * as.numeric(phi   %*% sp$K %*% phi))
    if (has_c) ny     <- rgamma(1, sc$a + sc$rank / 2, sc$b + 0.5 * as.numeric(psi   %*% sc$K %*% psi))

    ## ---- ASIS interweaving: non-centred (ancillary) re-draw of each precision
    ## (Yu & Meng 2011).  Rescaling the effect and its precision together breaks
    ## the precision-effect coupling that otherwise slows mixing, especially for
    ## the smoothing of weakly-informed cells in highly informative data. -------
    eta <- mu + outer(theta, phi, "+") + matrix(psi[cohidx], I, J)
    zwork <- ymh / omega
    nc_step <- function(x, prec, a, b, make_deta) {
      xt <- sqrt(prec) * x
      ll <- function(xn) { ek <- eta + make_deta(xn - x); -0.5 * sum(omega * (ek - zwork)^2) }
      lpost <- function(k, xn) ll(xn) + a * log(k) - b * k        # +a*log k: Gamma prior & log-k Jacobian
      propk <- exp(log(prec) + rnorm(1, 0, 0.4))
      xp <- xt / sqrt(propk)
      if (log(runif(1)) < lpost(propk, xp) - lpost(prec, x)) {
        eta <<- eta + make_deta(xp - x); list(x = xp, prec = propk)
      } else list(x = x, prec = prec)
    }
    if (has_a) { r <- nc_step(theta, kappa,  sa$a, sa$b, function(d) matrix(d, I, J));            theta <- r$x; kappa  <- r$prec }
    if (has_p) { r <- nc_step(phi,   lambda, sp$a, sp$b, function(d) matrix(d, I, J, byrow = TRUE)); phi <- r$x; lambda <- r$prec }
    if (has_c) { r <- nc_step(psi,   ny,     sc$a, sc$b, function(d) matrix(d[cohidx], I, J));      psi   <- r$x; ny     <- r$prec }

    if (it > burn_in && ((it - burn_in) %% thin == 0)) {
      keep <- keep + 1L
      out_theta[keep, ] <- theta; out_phi[keep, ] <- phi; out_psi[keep, ] <- psi
      out_kap[keep] <- kappa; out_lam[keep] <- lambda; out_ny[keep] <- ny; out_mu[keep] <- mu
      eta <- mu + outer(theta, phi, "+") + matrix(psi[cohidx], I, J)
      ksi_sum <- ksi_sum + eta; ksi_n <- ksi_n + 1L
      pr <- 1 / (1 + exp(-eta)); yhat <- N * pr
      d1 <- 2 * ((N - Y) * log((N - Y) / (N - yhat)))
      d2 <- 2 * (Y * log(Y / yhat) + (N - Y) * log((N - Y) / (N - yhat)))
      d2[is.nan(d2)] <- d1[is.nan(d2)]
      out_dev[keep] <- sum(d2, na.rm = TRUE)
    }
  }
  list(theta = out_theta, phi = out_phi, psi = out_psi,
       kappa = out_kap, lambda = out_lam, ny = out_ny, my = out_mu,
       deviance = out_dev, ksi = ksi_sum / ksi_n)
}

## --- driver: run chains (optionally in parallel) ---------------------------
.bamp_pg <- function(Y, N, ord_a, ord_p, ord_c, ppa, hyper,
                     n_iter, burn_in, thin, n_chains, parallel = FALSE,
                     prior_scale = TRUE, verbose = FALSE) {
  seeds <- sample.int(.Machine$integer.max, n_chains)
  runner <- function(s) .bamp_pg_chain(Y, N, ord_a, ord_p, ord_c, ppa, hyper,
                                       n_iter, burn_in, thin, s, prior_scale)
  if (isTRUE(parallel) && .Platform$OS.type != "windows") {
    cores <- max(1L, min(n_chains, getOption("mc.cores", 2L)))
    parallel::mclapply(seeds, runner, mc.cores = cores)
  } else {
    lapply(seeds, runner)
  }
}
