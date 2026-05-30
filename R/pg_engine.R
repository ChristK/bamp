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
                           prior_scale = TRUE, init = "empirical",
                           overdisp = FALSE, z_hyper = c(1, 0.05)) {
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
  ## orthonormal basis of the constraint null space {beta : A beta = 0}, used by
  ## the Laplace-MH refinement to sample in the (unconstrained) free coordinates
  Zbasis <- if (is.null(A)) diag(P) else {
    sv <- svd(A, nu = 0, nv = P); sv$v[, (nrow(A) + 1):P, drop = FALSE]
  }

  ## starting values
  p0 <- sum(Y) / sum(N); mu <- log(p0 / (1 - p0))
  theta <- numeric(I); phi <- numeric(J); psi <- numeric(K)
  if (identical(init, "empirical")) {
    ## start age/period at their (centred) marginal empirical log-odds so the
    ## chain begins near the data structure -- short burn-in, and all chains
    ## start in the same region
    sl <- function(num, den) {
      p <- pmin(pmax((num + 0.5) / (den + 1), 1e-8), 1 - 1e-8); log(p / (1 - p))
    }
    if (has_a) { t <- sl(rowSums(Y), rowSums(N)) - mu; theta <- t - mean(t) }
    if (has_p) { t <- sl(colSums(Y), colSums(N)) - mu; phi   <- t - mean(t) }
  }
  kappa <- lambda <- ny <- 1
  ## overdispersion: an additive per-cell effect delta_ij ~ N(0, 1/zeta) in the
  ## linear predictor, with precision zeta ~ Gamma(z_hyper). Under PG both have
  ## closed-form conditionals (Gaussian for delta, Gamma for zeta). Only the
  ## scalar precision is stored (the cell effects are re-sampled at predict time),
  ## matching the iwls output contract.
  has_od <- isTRUE(overdisp)
  zeta <- z_hyper[1] / z_hyper[2]
  delta_ij <- matrix(0, I, J)

  nkeep <- length(seq(burn_in + 1, n_iter, by = thin))
  out_theta <- matrix(0, nkeep, I); out_phi <- matrix(0, nkeep, J); out_psi <- matrix(0, nkeep, K)
  out_kap <- out_lam <- out_ny <- out_mu <- out_dev <- out_zeta <- numeric(nkeep)
  ksi_sum <- matrix(0, I, J); ksi_n <- 0L; keep <- 0L

  ## constant pieces of b
  b_mu0 <- sum(ymh)
  b_th0 <- if (has_a) rowSums(ymh) else NULL
  b_ph0 <- if (has_p) colSums(ymh) else NULL
  b_ps0 <- if (has_c) { v <- numeric(K); a <- rowsum(as.numeric(ymh), cohv); v[as.integer(rownames(a))] <- a; v } else NULL

  ## assemble the joint precision X' diag(w) X + prior for a per-cell weight w
  ## (w = Polya-Gamma omega for the Gibbs draw, or the Fisher weight N p(1-p)
  ## for the Laplace-MH proposal). Returns the P x P precision matrix.
  assemble_prec <- function(w, kap, lam, nyv) {
    Qm <- matrix(0, P, P); wv <- as.numeric(w)
    rs <- rowSums(w); cs <- colSums(w)
    TP <- PP <- NULL
    if (has_a && has_c) { TP <- matrix(0, I, K); TP[idxT] <- wv }
    if (has_p && has_c) { PP <- matrix(0, J, K); PP[idxP] <- wv }
    Qm[ia, ia] <- sum(w)
    if (has_a) { Qm[ia, ith] <- rs; Qm[ith, ia] <- rs; Qm[ith, ith] <- diag(rs, I) + kap * sa$K }
    if (has_p) { Qm[ia, iph] <- cs; Qm[iph, ia] <- cs; Qm[iph, iph] <- diag(cs, J) + lam * sp$K }
    if (has_c) {
      csum <- if (!is.null(TP)) colSums(TP) else if (!is.null(PP)) colSums(PP)
              else { v <- numeric(K); ag <- rowsum(wv, cohv); v[as.integer(rownames(ag))] <- ag; v }
      Qm[ia, ips] <- csum; Qm[ips, ia] <- csum; Qm[ips, ips] <- diag(csum, K) + nyv * sc$K
    }
    if (has_a && has_p) { Qm[ith, iph] <- w; Qm[iph, ith] <- t(w) }
    if (!is.null(TP)) { Qm[ith, ips] <- TP; Qm[ips, ith] <- t(TP) }
    if (!is.null(PP)) { Qm[iph, ips] <- PP; Qm[ips, iph] <- t(PP) }
    Qm
  }
  ## numerically stable log(1 + exp(e))
  softplus <- function(e) ifelse(e > 0, e + log1p(exp(-e)), log1p(exp(e)))

  acc_mh <- 0L; n_mh <- 0L

  for (it in 1:n_iter) {
    eta <- mu + outer(theta, phi, "+") + matrix(psi[cohidx], I, J) +
           (if (has_od) delta_ij else 0)
    omega <- matrix(.pg_rpg(as.numeric(N), as.numeric(eta)), I, J)

    ## ---- Polya-Gamma joint Gibbs draw of (mu, theta, phi, psi) ----
    ## With overdispersion, marginalise the cell effect delta out of the smooth
    ## update (collapsed Gibbs). Given omega and zeta, integrating delta_ij gives
    ## the smooth predictor a per-cell working precision wgt = omega*zeta/(omega+
    ## zeta) and response X'(wgt*z) = X'(zeta/(omega+zeta) * ymh) (z = ymh/omega).
    ## Drawing the smooth block from this marginal, then delta | smooth below, is
    ## an exact joint draw of (smooth, delta) and removes the structured/
    ## unstructured confounding that makes the alternating Gibbs mix very slowly.
    b <- numeric(P)
    if (has_od) {
      wgt <- omega * zeta / (omega + zeta)
      bm  <- ymh * zeta / (omega + zeta)              # = wgt * z
      b[ia] <- sum(bm)
      if (has_a) b[ith] <- rowSums(bm)
      if (has_p) b[iph] <- colSums(bm)
      if (has_c) { v <- numeric(K); ag <- rowsum(as.numeric(bm), cohv); v[as.integer(rownames(ag))] <- ag; b[ips] <- v }
    } else {
      wgt <- omega
      b[ia] <- b_mu0
      if (has_a) b[ith] <- b_th0
      if (has_p) b[iph] <- b_ph0
      if (has_c) b[ips] <- b_ps0
    }
    Q <- assemble_prec(wgt, kappa, lambda, ny)
    ## ridge to make the (intrinsically rank-deficient) Q numerically PD; the
    ## constraints remove the corresponding null directions
    diag(Q) <- diag(Q) + 1e-6 * mean(diag(Q))
    beta <- .pg_draw_block(Q, b, A)
    mu <- beta[ia]
    if (has_a) theta <- beta[ith]
    if (has_p) phi   <- beta[iph]
    if (has_c) psi   <- beta[ips]

    ## ---- overdispersion: cell effect delta_ij | rest ~ N, precision zeta | delta ~ Gamma
    if (has_od) {
      eta0 <- mu + outer(theta, phi, "+") + matrix(psi[cohidx], I, J)   # smooth predictor
      precd <- omega + zeta
      mden <- (ymh - omega * eta0) / precd                             # = omega*(z - eta0)/precd
      delta_ij <- matrix(mden + rnorm(I * J) / sqrt(precd), I, J)
      zeta <- rgamma(1, z_hyper[1] + 0.5 * I * J, z_hyper[2] + 0.5 * sum(delta_ij^2))
    }

    ## ---- precisions: centred (sufficient) Gamma full conditionals ----
    if (has_a) kappa  <- rgamma(1, sa$a + sa$rank / 2, sa$b + 0.5 * as.numeric(theta %*% sa$K %*% theta))
    if (has_p) lambda <- rgamma(1, sp$a + sp$rank / 2, sp$b + 0.5 * as.numeric(phi   %*% sp$K %*% phi))
    if (has_c) ny     <- rgamma(1, sc$a + sc$rank / 2, sc$b + 0.5 * as.numeric(psi   %*% sc$K %*% psi))

    ## ---- ASIS interweaving: non-centred (ancillary) re-draw of each precision
    ## (Yu & Meng 2011).  Rescaling the effect and its precision together breaks
    ## the precision-effect coupling that otherwise slows mixing, especially for
    ## the smoothing of weakly-informed cells in highly informative data. -------
    eta <- mu + outer(theta, phi, "+") + matrix(psi[cohidx], I, J) +
           (if (has_od) delta_ij else 0)
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

    ## ---- Laplace (Newton) Metropolis-Hastings refinement -------------------
    ## Pure Polya-Gamma draws move weakly-informed cells in tiny steps because
    ## the augmented conditional is far tighter than the marginal. This step
    ## proposes a joint Newton move against the TRUE binomial likelihood using
    ## the Fisher weight N p(1-p) (wide where data are sparse), in the free
    ## coordinates gamma (beta = Z gamma), so it mixes those cells properly.
    state <- function(bv) {
      mu_ <- bv[ia]
      th_ <- if (has_a) bv[ith] else 0; ph_ <- if (has_p) bv[iph] else 0
      ps_ <- if (has_c) bv[ips] else NULL
      e <- mu_ + (if (has_a) matrix(bv[ith], I, J) else 0) +
                 (if (has_p) matrix(bv[iph], I, J, byrow = TRUE) else 0) +
                 (if (has_c) matrix(ps_[cohidx], I, J) else 0) +
                 (if (has_od) delta_ij else 0)          # cell effect is a fixed offset here
      pp <- plogis(e); Wt <- N * pp * (1 - pp); ymnp <- Y - N * pp
      g <- numeric(P); g[ia] <- sum(ymnp)
      if (has_a) g[ith] <- rowSums(ymnp) - kappa  * as.numeric(sa$K %*% bv[ith])
      if (has_p) g[iph] <- colSums(ymnp) - lambda * as.numeric(sp$K %*% bv[iph])
      if (has_c) { cg <- numeric(K); ag <- rowsum(as.numeric(ymnp), cohv)
                   cg[as.integer(rownames(ag))] <- ag
                   g[ips] <- cg - ny * as.numeric(sc$K %*% ps_) }
      H <- assemble_prec(Wt, kappa, lambda, ny); diag(H) <- diag(H) + 1e-6 * mean(diag(H))
      Hg <- crossprod(Zbasis, H %*% Zbasis); gg <- as.numeric(crossprod(Zbasis, g))
      R <- chol(Hg)
      mean_g <- as.numeric(crossprod(Zbasis, bv)) + backsolve(R, forwardsolve(t(R), gg))
      lp <- sum(Y * e - N * softplus(e)) -
            0.5 * ((if (has_a) kappa  * as.numeric(bv[ith] %*% sa$K %*% bv[ith]) else 0) +
                   (if (has_p) lambda * as.numeric(bv[iph] %*% sp$K %*% bv[iph]) else 0) +
                   (if (has_c) ny     * as.numeric(ps_     %*% sc$K %*% ps_)     else 0))
      list(R = R, mean_g = mean_g, lp = lp, gamma = as.numeric(crossprod(Zbasis, bv)))
    }
    bcur <- numeric(P); bcur[ia] <- mu
    if (has_a) bcur[ith] <- theta; if (has_p) bcur[iph] <- phi; if (has_c) bcur[ips] <- psi
    cur <- state(bcur)
    gstar <- cur$mean_g + backsolve(cur$R, rnorm(length(cur$mean_g)))
    bstar <- as.numeric(Zbasis %*% gstar)
    prop <- state(bstar)
    logq <- function(R, x, m) sum(log(diag(R))) - 0.5 * sum((R %*% (x - m))^2)
    la <- prop$lp - cur$lp + logq(prop$R, cur$gamma, prop$mean_g) - logq(cur$R, gstar, cur$mean_g)
    n_mh <- n_mh + 1L
    if (is.finite(la) && log(runif(1)) < la) {
      acc_mh <- acc_mh + 1L
      mu <- bstar[ia]
      if (has_a) theta <- bstar[ith]; if (has_p) phi <- bstar[iph]; if (has_c) psi <- bstar[ips]
    }

    if (it > burn_in && ((it - burn_in) %% thin == 0)) {
      keep <- keep + 1L
      out_theta[keep, ] <- theta; out_phi[keep, ] <- phi; out_psi[keep, ] <- psi
      out_kap[keep] <- kappa; out_lam[keep] <- lambda; out_ny[keep] <- ny; out_mu[keep] <- mu
      out_zeta[keep] <- zeta
      eta <- mu + outer(theta, phi, "+") + matrix(psi[cohidx], I, J) +
             (if (has_od) delta_ij else 0)
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
       zeta = if (has_od) out_zeta else NULL,
       deviance = out_dev, ksi = ksi_sum / ksi_n,
       mh_accept = if (n_mh > 0) acc_mh / n_mh else NA_real_)
}

## --- driver: run chains (optionally in parallel) ---------------------------
.bamp_pg <- function(Y, N, ord_a, ord_p, ord_c, ppa, hyper,
                     n_iter, burn_in, thin, n_chains, parallel = FALSE,
                     prior_scale = TRUE, verbose = FALSE,
                     overdisp = FALSE, z_hyper = c(1, 0.05)) {
  seeds <- sample.int(.Machine$integer.max, n_chains)
  runner <- function(s) .bamp_pg_chain(Y, N, ord_a, ord_p, ord_c, ppa, hyper,
                                       n_iter, burn_in, thin, s, prior_scale,
                                       overdisp = overdisp, z_hyper = z_hyper)
  ## Honour a numeric `parallel` as the requested number of cores (matching the
  ## iwls path, where cores <- parallel); a bare TRUE means getOption('mc.cores').
  ## Capped at the number of chains. Previously cores were hard-capped at 2, so
  ## parallel=4 ran 4 chains on only 2 cores and PG wall-clock was ~2x inflated.
  run_par <- (isTRUE(parallel) || (is.numeric(parallel) && parallel > 1)) &&
             .Platform$OS.type != "windows"
  if (run_par) {
    req <- if (is.numeric(parallel)) parallel else getOption("mc.cores", 2L)
    cores <- max(1L, min(n_chains, req))
    parallel::mclapply(seeds, runner, mc.cores = cores)
  } else {
    lapply(seeds, runner)
  }
}
