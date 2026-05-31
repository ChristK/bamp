## ===========================================================================
## PHASE 1 RESEARCH PROTOTYPE -- joint sex-coherent age-period-cohort model.
##
## A single joint posterior for two sexes that (i) borrows strength across sexes
## and (ii) cannot diverge in projection. The model is the design-note Strategy 0
## (Riebler-Held multivariate APC in bamp's Polya-Gamma field):
##
##   logit p_{s,i,j} = mu0 + a_s + theta_{s,i} + phi_j + psi_k + sgn(s) * delta_j
##
##   * phi_j (period), psi_k (cohort)  : SHARED across sexes, intrinsic RW1/RW2
##                                       prior  -> projected as a FREE random walk
##                                       (the total trend is allowed to move).
##   * theta_{s,i} (age)               : sex-specific, intrinsic RW prior.
##   * delta_j (period sex-deviation)  : sgn = +1 (female) / -1 (male), so the
##                                       sexes are a sum-to-zero contrast; given a
##                                       PROPER mean-reverting prior (iid ridge =
##                                       AR1 with rho=0, or AR1 with |rho|<1) so it
##                                       is projected by a STATIONARY process and
##                                       the sex gap cannot drift.
##
## This is a deliberately simple, auditable REFERENCE: an explicit design matrix
## and a dense one-block Polya-Gamma Gibbs draw (reusing the engine helpers
## .pg_rpg / .pg_draw_block / .pg_Kmat / .pg_scale), with NO ASIS / Laplace-MH
## refinement. It is O(P^3) per sweep -- fine for S=2 on moderate grids, not for
## production. Port to the sparse C engine only after this validates (design note
## Phase 3). Cohort sex-deviation and S>2 are documented extensions.
## ===========================================================================

.coh_ord <- function(x) switch(as.character(x),
  "rw1" = 1L, "rw2" = 2L,
  stop("bamp_coherent (prototype) needs age/period/cohort = 'rw1' or 'rw2'."))

## AR1 precision (up to scale): tridiagonal, rho=0 -> identity (iid ridge).
.ar1_prec <- function(L, rho) {
  T <- diag(1 + rho^2, L); T[1, 1] <- 1; T[L, L] <- 1
  if (L >= 2) for (j in 2:L) { T[j, j - 1] <- -rho; T[j - 1, j] <- -rho }
  T
}

## free random-walk extrapolation of a SHARED effect (as in predict_apc)
.cj_predict_rw <- function(vec, lambda, rw, n1, n2) {
  if (n2 > n1) for (i in (n1 + 1):n2)
    vec[i] <- if (rw == 2) 2 * vec[i - 1] - vec[i - 2] + rnorm(1, 0, 1 / sqrt(lambda))
              else                 vec[i - 1]               + rnorm(1, 0, 1 / sqrt(lambda))
  vec
}

## mean-reverting AR1 extrapolation of the sex DEVIATION (rho=0 -> reverts to 0
## immediately; |rho|<1 -> reverts gradually). This is what stops the sex gap
## from diverging, the whole point of the coherent model.
.cj_predict_ar <- function(vec, lambda, rho, n1, n2) {
  if (n2 > n1) for (i in (n1 + 1):n2)
    vec[i] <- rho * vec[i - 1] + rnorm(1, 0, 1 / sqrt(lambda))
  vec
}

#' Joint sex-coherent age-period-cohort model (Phase 1 prototype)
#'
#' @description
#' Fit two sexes in ONE joint posterior so that they borrow strength and -- unlike independent
#' fits -- cannot diverge in projection. Shared period and cohort effects carry the common trend;
#' a sex-specific period \emph{deviation} with a proper mean-reverting prior carries the
#' (non-diverging) sex gap. This is the principled alternative to the Phase 0 \code{\link{bamp_strata}}
#' total-plus-share wrapper; see \code{docs/coherent-forecasting.md}.
#'
#' This is a research-grade REFERENCE implementation (dense one-block Polya-Gamma Gibbs, S = 2,
#' period deviation only); it is correct but not optimised. For production use the design-note
#' roadmap (sparse C engine, cohort deviation, S > 2).
#'
#' @param cases,population named lists of two \code{[periods x agegroups]} matrices (the two sexes),
#'   same dimensions (as for \code{\link{bamp_strata}}).
#' @param age,period,cohort each \code{"rw1"} or \code{"rw2"} (all three are required).
#' @param periods_per_agegroup integer M linking the cohort index, as in \code{\link{bamp}}.
#' @param deviation prior for the sex deviation: \code{"iid"} (ridge; instant reversion) or
#'   \code{"ar1"} (gradual reversion, with \code{rho}). Both are proper and mean-reverting.
#' @param rho AR1 coefficient for \code{deviation = "ar1"}. Now a SAMPLED parameter (the value passed
#'   is the starting point); the posterior draws are returned in \code{samples$rho}. \code{0 <= rho < 1}.
#'   Ignored for \code{"iid"} (fixed at 0).
#' @param deviation_cohort optionally also place a sex deviation on the COHORT axis: \code{"none"}
#'   (default; identical to before), \code{"iid"} or \code{"ar1"} (with \code{rho_c}). Mean-reverting,
#'   sum-to-zero over cohorts.
#' @param rho_c fixed AR1 coefficient for the cohort deviation (\code{0 <= rho_c < 1}); used only when
#'   \code{deviation_cohort = "ar1"}.
#' @param mcmc list with \code{iterations}, \code{burn_in}, \code{thin}.
#' @param hyper list of Gamma hyperparameters \code{c(a, b)} for the precisions \code{age},
#'   \code{period}, \code{cohort}, \code{dev} (period deviation) and \code{dev_cohort}, and the
#'   \code{rho} Beta-prior parameters.
#' @param prior_scale scale intrinsic structure matrices to unit generalised variance
#'   (Sorbye-Rue), as in the main engine.
#' @param seed RNG seed.
#' @param mh_sd_rho proposal standard deviation (logit scale) for the Metropolis update of \code{rho}.
#'
#' @return object of class \code{apc_coherent}: posterior draws of all effects and precisions,
#'   the model/data metadata, needed by \code{\link{predict_coherent}}.
#' @seealso \code{\link{predict_coherent}}, \code{\link{bamp_strata}}, \code{\link{bamp}}
#' @export
bamp_coherent <- function(cases, population, age = "rw1", period = "rw1", cohort = "rw1",
                          periods_per_agegroup,
                          deviation = c("iid", "ar1"), rho = 0,
                          deviation_cohort = c("none", "iid", "ar1"), rho_c = 0,
                          mcmc = list(iterations = 4000, burn_in = 1000, thin = 2),
                          hyper = list(age = c(1, 0.5), period = c(1, 5e-4),
                                       cohort = c(1, 5e-4), dev = c(1, 0.05), rho = c(1, 1),
                                       dev_cohort = c(1, 0.05)),
                          prior_scale = TRUE, seed = 1, mh_sd_rho = 0.3) {
  deviation <- match.arg(deviation)
  if (deviation == "iid") rho <- 0
  if (!(is.numeric(rho) && length(rho) == 1 && rho >= 0 && rho < 1))
    stop("'rho' must satisfy 0 <= rho < 1.")
  deviation_cohort <- match.arg(deviation_cohort)
  if (deviation_cohort != "ar1") rho_c <- 0
  if (!(is.numeric(rho_c) && length(rho_c) == 1 && rho_c >= 0 && rho_c < 1))
    stop("'rho_c' must satisfy 0 <= rho_c < 1.")
  use_dpsi <- deviation_cohort != "none"           # optional cohort-axis sex deviation
  if (!is.list(cases) || length(cases) != 2L || !is.list(population) || length(population) != 2L)
    stop("'cases' and 'population' must each be a list of two [periods x agegroups] matrices (the two sexes).")
  if (is.null(names(cases))) names(cases) <- c("sex1", "sex2")
  names(population) <- names(cases)
  ord_a <- .coh_ord(age); ord_p <- .coh_ord(period); ord_c <- .coh_ord(cohort)

  ## internal orientation: Y[[s]] is [age x period]  (engine convention)
  Y <- lapply(cases, function(x) t(as.matrix(x)))
  Npop <- lapply(population, function(x) t(as.matrix(x)))
  I <- nrow(Y[[1]]); J <- ncol(Y[[1]]); M <- as.integer(periods_per_agegroup)
  K <- M * (I - 1L) + J
  S <- 2L; sgn <- c(1, -1)

  ## structure matrices (+ Sorbye-Rue scaling) and ranks
  scl <- function(Km, ord) if (prior_scale) Km * .pg_scale(Km, ord) else Km
  Ka <- scl(.pg_Kmat(I, ord_a), ord_a)
  Kp <- scl(.pg_Kmat(J, ord_p), ord_p)
  Kc <- scl(.pg_Kmat(K, ord_c), ord_c)

  ## AR1 coherence coefficient rho: sampled (deviation = "ar1") or fixed at 0 ("iid").
  ## The MH ratio for rho needs the CONSTRAINED log-determinant of the deviation
  ## precision over the sum-to-zero null space (Z), the logit-proposal Jacobian and
  ## a Beta(rho) prior -- omitting the determinant biases rho low. (See docs/hardening-plan.md.)
  sample_rho <- deviation == "ar1"
  rho_cur <- if (sample_rho && rho == 0) 0.5 else rho   # rho arg is the starting value for ar1
  Td <- .ar1_prec(J, rho_cur)
  Zd <- svd(matrix(1, 1, J), nu = 0, nv = J)$v[, 2:J, drop = FALSE]   # basis of {x : sum(x)=0}
  ld_constr <- function(rr) 2 * sum(log(diag(chol(crossprod(Zd, .ar1_prec(J, rr) %*% Zd)))))
  Tc <- if (use_dpsi) .ar1_prec(K, rho_c) else NULL                  # cohort-deviation AR1 precision

  ## parameter layout: beta = (mu0, a[S], theta1[I], theta2[I], phi[J], psi[K], delta[J])
  i_mu <- 1L
  i_a  <- 1L + seq_len(S)
  i_t1 <- max(i_a) + seq_len(I)
  i_t2 <- max(i_t1) + seq_len(I)
  i_ph <- max(i_t2) + seq_len(J)
  i_ps <- max(i_ph) + seq_len(K)
  i_de <- max(i_ps) + seq_len(J)
  i_dc <- if (use_dpsi) max(i_de) + seq_len(K) else integer(0)
  P <- if (use_dpsi) max(i_dc) else max(i_de)

  ## explicit design matrix and stacked data (cells ordered j-fast, then i, then s)
  grid <- expand.grid(j = 1:J, i = 1:I, s = 1:S)
  n <- nrow(grid)
  kcell <- (I - grid$i) * M + grid$j
  X <- matrix(0, n, P)
  X[, i_mu] <- 1
  X[cbind(seq_len(n), i_a[grid$s])] <- 1
  X[cbind(seq_len(n), ifelse(grid$s == 1L, i_t1[grid$i], i_t2[grid$i]))] <- 1
  X[cbind(seq_len(n), i_ph[grid$j])] <- 1
  X[cbind(seq_len(n), i_ps[kcell])] <- 1
  X[cbind(seq_len(n), i_de[grid$j])] <- sgn[grid$s]
  if (use_dpsi) X[cbind(seq_len(n), i_dc[kcell])] <- sgn[grid$s]   # cohort-axis sex deviation
  yv <- Nv <- numeric(n)
  for (s in 1:S) {
    sel <- grid$s == s
    yv[sel] <- Y[[s]][cbind(grid$i[sel], grid$j[sel])]
    Nv[sel] <- Npop[[s]][cbind(grid$i[sel], grid$j[sel])]
  }
  Xt <- t(X); bvec <- as.numeric(Xt %*% (yv - Nv / 2))   # b is data-only, constant

  ## sum-to-zero constraints (+ RW2 period zero-slope drift pin, as in the engine)
  mkrow <- function(idx, val) { r <- numeric(P); r[idx] <- val; r }
  Arows <- list(mkrow(i_a, 1), mkrow(i_t1, 1), mkrow(i_t2, 1),
                mkrow(i_ph, 1), mkrow(i_ps, 1), mkrow(i_de, 1))
  if (use_dpsi) Arows <- c(Arows, list(mkrow(i_dc, 1)))   # sum_k dpsi = 0 (Tc proper: no slope row)
  if (ord_p == 2L) Arows <- c(Arows, list(mkrow(i_ph, (1:J) - mean(1:J))))
  A <- do.call(rbind, Arows); tA <- t(A)

  tau_a <- 1e-2                                 # weak normal on the sex main effect
  build_prec <- function(kth, lph, nps, ld, ldc) {
    Pm <- matrix(0, P, P)
    Pm[i_mu, i_mu] <- 1e-6
    Pm[cbind(i_a, i_a)] <- tau_a
    Pm[i_t1, i_t1] <- kth * Ka; Pm[i_t2, i_t2] <- kth * Ka
    Pm[i_ph, i_ph] <- lph * Kp
    Pm[i_ps, i_ps] <- nps * Kc
    Pm[i_de, i_de] <- ld * Td
    if (use_dpsi) Pm[i_dc, i_dc] <- ldc * Tc
    diag(Pm) <- diag(Pm) + 1e-7                 # tiny global ridge for PD chol
    Pm
  }

  hp <- function(x) if (is.null(x)) c(1, 1) else x
  ha <- hp(hyper$age); hpp <- hp(hyper$period); hc <- hp(hyper$cohort); hd <- hp(hyper$dev)
  hr <- hp(hyper$rho); hdc <- hp(hyper$dev_cohort)

  ## init
  set.seed(seed)
  p0 <- sum(yv) / sum(Nv); beta <- numeric(P); beta[i_mu] <- log(p0 / (1 - p0))
  kth <- lph <- nps <- ld <- 1; ldc <- if (use_dpsi) 1 else NA_real_
  Prec <- build_prec(kth, lph, nps, ld, ldc)

  iters <- mcmc$iterations; burn <- mcmc$burn_in; thin <- mcmc$thin
  keep <- seq.int(burn + thin, iters, by = thin); nkeep <- length(keep); store <- 0L
  out <- list(mu0 = numeric(nkeep), a = matrix(0, nkeep, S),
              theta = array(0, c(nkeep, S, I)), phi = matrix(0, nkeep, J),
              psi = matrix(0, nkeep, K), delta = matrix(0, nkeep, J),
              lambda_phi = numeric(nkeep), nu_psi = numeric(nkeep),
              lambda_d = numeric(nkeep), kappa_theta = numeric(nkeep),
              rho = numeric(nkeep))
  if (use_dpsi) { out$dpsi <- matrix(0, nkeep, K); out$lambda_dc <- numeric(nkeep) }

  qform <- function(v, Km) sum(v * (Km %*% v))
  n_acc <- 0L
  for (it in seq_len(iters)) {
    eta <- as.numeric(X %*% beta)
    omega <- .pg_rpg(Nv, eta)
    Q <- Xt %*% (X * omega) + Prec
    beta <- .pg_draw_block(Q, bvec, A, tA)

    th1 <- beta[i_t1]; th2 <- beta[i_t2]; ph <- beta[i_ph]; ps <- beta[i_ps]; de <- beta[i_de]
    kth <- rgamma(1, ha[1]  + (I - ord_a),       ha[2]  + 0.5 * (qform(th1, Ka) + qform(th2, Ka)))
    lph <- rgamma(1, hpp[1] + (J - ord_p) / 2,   hpp[2] + 0.5 * qform(ph, Kp))
    nps <- rgamma(1, hc[1]  + (K - ord_c) / 2,   hc[2]  + 0.5 * qform(ps, Kc))
    ld  <- rgamma(1, hd[1]  + (J - 1) / 2,        hd[2]  + 0.5 * qform(de, Td))

    ## sample rho (ar1 only): logit-RW Metropolis with the constrained determinant.
    ## Td is refreshed at loop scope so the next build_prec() sees the accepted value.
    if (sample_rho) {
      rho_star <- plogis(qlogis(rho_cur) + rnorm(1, 0, mh_sd_rho))
      Td_star <- .ar1_prec(J, rho_star)
      logr <- 0.5 * (ld_constr(rho_star) - ld_constr(rho_cur)) -
              0.5 * ld * (qform(de, Td_star) - qform(de, Td)) +
              (dbeta(rho_star, hr[1], hr[2], log = TRUE) - dbeta(rho_cur, hr[1], hr[2], log = TRUE)) +
              (log(rho_star * (1 - rho_star)) - log(rho_cur * (1 - rho_cur)))
      if (is.finite(logr) && log(runif(1)) < logr) { rho_cur <- rho_star; Td <- Td_star; n_acc <- n_acc + 1L }
    }
    if (use_dpsi) { dpsi <- beta[i_dc]; ldc <- rgamma(1, hdc[1] + (K - 1) / 2, hdc[2] + 0.5 * qform(dpsi, Tc)) }
    Prec <- build_prec(kth, lph, nps, ld, ldc)

    if (it %in% keep) {
      store <- store + 1L
      out$mu0[store] <- beta[i_mu]; out$a[store, ] <- beta[i_a]
      out$theta[store, 1, ] <- th1; out$theta[store, 2, ] <- th2
      out$phi[store, ] <- ph; out$psi[store, ] <- ps; out$delta[store, ] <- de
      out$lambda_phi[store] <- lph; out$nu_psi[store] <- nps
      out$lambda_d[store] <- ld; out$kappa_theta[store] <- kth; out$rho[store] <- rho_cur
      if (use_dpsi) { out$dpsi[store, ] <- dpsi; out$lambda_dc[store] <- ldc }
    }
  }

  structure(list(
    samples = out,
    model = list(age = age, period = period, cohort = cohort,
                 deviation = deviation, rho = rho, rho_sampled = sample_rho,
                 rho_accept = if (sample_rho) n_acc / iters else NA_real_,
                 deviation_cohort = deviation_cohort, rho_c = rho_c, use_dpsi = use_dpsi,
                 ord = c(ord_a, ord_p, ord_c)),
    data = list(cases = cases, population = population, periods_per_agegroup = M,
                I = I, J = J, K = K, strata = names(cases))
  ), class = "apc_coherent")
}


#' Coherent projection from a joint sex-coherent APC model
#'
#' @description
#' Project sex-specific rates (and hazards) from \code{\link{bamp_coherent}}. The shared period and
#' cohort effects are extrapolated as free random walks (the total trend may move); the sex
#' deviation is extrapolated by its mean-reverting process, so the sex gap stays bounded and the
#' sexes aggregate to a coherent total.
#'
#' @param object an \code{apc_coherent} object.
#' @param periods number of future periods to project.
#' @param population optional named list of two \emph{future} \code{[periods x agegroups]} population
#'   matrices (with \code{periods} extra rows) for the population-weighted total; \code{NULL} uses the
#'   in-sample populations (model checking only).
#' @param quantiles quantiles to summarise.
#' @param hazard,period_length as in \code{\link{predict_apc}}: if \code{hazard=TRUE} also return the
#'   per-stratum cumulative hazard \eqn{-\log(1-\mathrm{rate})/\mathrm{period\_length}}.
#'
#' @return list with one entry per sex and a \code{total} entry (quantiles of \code{rate}, and
#'   \code{hazard} if requested, on the \code{[period, agegroup]} grid, plus \code{samples}); plus
#'   \code{deviation} (the projected sex deviation \code{delta}, \code{[period, draw]}) for diagnosing
#'   non-divergence; and the stratum \code{order}.
#' @seealso \code{\link{bamp_coherent}}
#' @export
predict_coherent <- function(object, periods = 0, population = NULL,
                             quantiles = c(0.05, 0.5, 0.95),
                             hazard = FALSE, period_length = 1) {
  if (!inherits(object, "apc_coherent")) stop("'object' must come from bamp_coherent().")
  hazard <- isTRUE(hazard)
  s <- object$samples; md <- object$model; dat <- object$data
  I <- dat$I; J <- dat$J; K <- dat$K; M <- dat$periods_per_agegroup
  ord_p <- md$ord[2]; ord_c <- md$ord[3]
  n1 <- J; n2 <- J + periods; K2 <- (I - 1) * M + n2          # = coh(1, n2)
  D <- length(s$mu0); sgn <- c(1, -1)
  rho_vec <- if (!is.null(s$rho)) s$rho else rep(md$rho, D)   # per-draw rho (old objects: fixed)

  rateF <- rateM <- array(0, c(n2, I, D))
  dev_ext <- matrix(0, n2, D)
  for (d in seq_len(D)) {
    ph <- .cj_predict_rw(c(s$phi[d, ], numeric(n2 - n1)), s$lambda_phi[d], ord_p, n1, n2)
    ps <- .cj_predict_rw(c(s$psi[d, ], numeric(K2 - K)),  s$nu_psi[d],    ord_c, K, K2)
    de <- .cj_predict_ar(c(s$delta[d, ], numeric(n2 - n1)), s$lambda_d[d], rho_vec[d], n1, n2)
    dev_ext[, d] <- de
    dc <- if (isTRUE(md$use_dpsi))
      .cj_predict_ar(c(s$dpsi[d, ], numeric(K2 - K)), s$lambda_dc[d], md$rho_c, K, K2) else NULL
    for (i in 1:I) {
      kk <- (I - i) * M + (1:n2)
      base <- s$mu0[d] + ph + ps[kk]
      cohdev <- if (!is.null(dc)) dc[kk] else 0      # cohort-axis sex deviation (mean-reverting)
      etaF <- base + s$a[d, 1] + s$theta[d, 1, i] + sgn[1] * (de + cohdev)
      etaM <- base + s$a[d, 2] + s$theta[d, 2, i] + sgn[2] * (de + cohdev)
      rateF[, i, d] <- 1 / (1 + exp(-etaF))
      rateM[, i, d] <- 1 / (1 + exp(-etaM))
    }
  }

  if (is.null(population)) {
    popF <- as.matrix(dat$population[[1]]); popM <- as.matrix(dat$population[[2]])  # [period, age]
  } else {
    popF <- as.matrix(population[[dat$strata[1]]]); popM <- as.matrix(population[[dat$strata[2]]])
  }
  np <- min(nrow(popF), n2)
  qf <- function(arr) apply(arr, 1:2, stats::quantile, quantiles)
  pack <- function(rate) {
    e <- list(rate = qf(rate), samples = list(rate = rate))
    if (hazard) { hz <- -log1p(-rate) / period_length; e$hazard <- qf(hz); e$samples$hazard <- hz }
    e
  }
  out <- stats::setNames(list(pack(rateF), pack(rateM)), dat$strata)

  wF <- popF[seq_len(np), , drop = FALSE]; wM <- popM[seq_len(np), , drop = FALSE]
  rt <- array(0, c(np, I, D))
  for (d in seq_len(D)) rt[, , d] <- (wF * rateF[seq_len(np), , d] + wM * rateM[seq_len(np), , d]) / (wF + wM)
  rt[!is.finite(rt)] <- 0
  out$total <- pack(rt)
  out$deviation <- list(samples = dev_ext, quantiles = apply(dev_ext, 1, stats::quantile, quantiles))
  out$order <- dat$strata
  out
}
