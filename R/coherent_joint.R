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

## random-walk extrapolation of a SHARED effect, with optional damped trend (A1)
## and innovation-variance shrinkage (A3). `damp` in [0,1] damps the RW2 drift:
## vec[i] = vec[i-1] + damp*(vec[i-1]-vec[i-2]) + noise, so damp=1 is the usual
## RW2 (linear continuation), damp=0 collapses to RW1 (flat from the last level) --
## the Gardner-McKenzie damped trend, which curbs RW2 over-extrapolation at long
## horizons. `var_damp` in (0,1] geometrically shrinks the per-step innovation sd
## so the predictive bands stop fanning out without bound. Defaults (1,1) reproduce
## the previous free random walk exactly.
.cj_predict_rw <- function(vec, lambda, rw, n1, n2, damp = 1, var_damp = 1) {
  if (n2 > n1) for (i in (n1 + 1):n2) {
    sdk <- (1 / sqrt(lambda)) * var_damp^(i - n1 - 1)
    vec[i] <- if (rw == 2) vec[i - 1] + damp * (vec[i - 1] - vec[i - 2]) + rnorm(1, 0, sdk)
              else                      vec[i - 1]                       + rnorm(1, 0, sdk)
  }
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
#' Fit two or more exhaustive strata (e.g. sexes, or education levels) in ONE joint posterior so that
#' they borrow strength and -- unlike independent fits -- cannot diverge in projection. Shared period
#' and cohort effects carry the common trend; per-stratum \emph{deviations} with proper mean-reverting
#' priors carry the (non-diverging) between-stratum gap. This is the principled alternative to the
#' Phase 0 \code{\link{bamp_strata}} total-plus-share wrapper; see \code{docs/coherent-forecasting.md}.
#'
#' Research-grade REFERENCE implementation (dense one-block Polya-Gamma Gibbs; correct but not
#' optimised -- port to a sparse C engine for production). \code{S = 2} additionally supports a
#' sampled AR1 coefficient (\code{rho}) and a cohort-axis deviation (\code{deviation_cohort});
#' \code{S >= 3} uses a contr.sum period deviation (\eqn{\sum_s d_{s,j}=0}).
#'
#' @param cases,population named lists of \code{S >= 2} \code{[periods x agegroups]} matrices (the
#'   strata), same dimensions (as for \code{\link{bamp_strata}}).
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
#' @param mh_sd_rho proposal standard deviation (logit scale) for the Metropolis update of \code{rho};
#'   the starting value when \code{adapt_rho = TRUE}.
#' @param adapt_rho if \code{TRUE} (default), Robbins-Monro adaptation tunes the \code{rho} proposal sd
#'   toward 0.234 acceptance during burn-in only (the sampling-phase kernel is then fixed, preserving
#'   detailed balance). The tuned value is returned in \code{model$mh_sd_rho}.
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
                          prior_scale = TRUE, seed = 1, mh_sd_rho = 0.3, adapt_rho = TRUE) {
  deviation <- match.arg(deviation)
  if (deviation == "iid") rho <- 0
  if (!(is.numeric(rho) && length(rho) == 1 && rho >= 0 && rho < 1))
    stop("'rho' must satisfy 0 <= rho < 1.")
  deviation_cohort <- match.arg(deviation_cohort)
  if (deviation_cohort != "ar1") rho_c <- 0
  if (!(is.numeric(rho_c) && length(rho_c) == 1 && rho_c >= 0 && rho_c < 1))
    stop("'rho_c' must satisfy 0 <= rho_c < 1.")
  use_dpsi <- deviation_cohort != "none"           # optional cohort-axis sex deviation
  if (!is.list(cases) || !is.list(population) || length(cases) < 2L || length(cases) != length(population))
    stop("'cases' and 'population' must be lists of the same length (>= 2) of [periods x agegroups] matrices.")
  if (is.null(names(cases)))
    names(cases) <- if (length(cases) == 2L) c("sex1", "sex2") else paste0("stratum", seq_along(cases))
  names(population) <- names(cases)
  ## S = 2 keeps the full-featured legacy sampler below (sampled rho, cohort deviation).
  ## S >= 3 dispatches to the general contr.sum period-deviation sampler.
  if (length(cases) >= 3L)
    return(.bamp_coherent_general(cases, population, age, period, cohort, periods_per_agegroup,
                                  deviation, rho, mcmc, hyper, prior_scale, seed))
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
  ## D5: Robbins-Monro adaptation of the logit-RW proposal sd toward 0.234 acceptance,
  ## active only during burn-in (diminishing adaptation -> the sampling-phase kernel is
  ## fixed, so detailed balance holds for the retained draws).
  mh_sd <- mh_sd_rho; lsd <- log(mh_sd_rho); target <- 0.234
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
      rho_star <- plogis(qlogis(rho_cur) + rnorm(1, 0, mh_sd))
      Td_star <- .ar1_prec(J, rho_star)
      logr <- 0.5 * (ld_constr(rho_star) - ld_constr(rho_cur)) -
              0.5 * ld * (qform(de, Td_star) - qform(de, Td)) +
              (dbeta(rho_star, hr[1], hr[2], log = TRUE) - dbeta(rho_cur, hr[1], hr[2], log = TRUE)) +
              (log(rho_star * (1 - rho_star)) - log(rho_cur * (1 - rho_cur)))
      acc <- is.finite(logr) && log(runif(1)) < logr
      if (acc) { rho_cur <- rho_star; Td <- Td_star; n_acc <- n_acc + 1L }
      if (adapt_rho && it <= burn) {                       # diminishing RM step, burn-in only
        lsd <- lsd + min(0.5, 5 / it) * ((if (acc) 1 else 0) - target)
        mh_sd <- exp(max(-8, min(2, lsd)))                 # clamp to a sane sd range
      }
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
                 mh_sd_rho = mh_sd, mh_sd_rho_adapted = isTRUE(adapt_rho && sample_rho),
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
#' @param damping,var_damping trend-damping and innovation-variance shrinkage for the shared period/
#'   cohort extrapolation, as in \code{\link{predict_apc}}. Defaults \code{1, 1} reproduce the free
#'   random walk. The mean-reverting sex deviation is unaffected (it already reverts via \code{rho}).
#'
#' @return list with one entry per sex and a \code{total} entry (quantiles of \code{rate}, and
#'   \code{hazard} if requested, on the \code{[period, agegroup]} grid, plus \code{samples}); plus
#'   \code{deviation} (the projected sex deviation \code{delta}, \code{[period, draw]}) for diagnosing
#'   non-divergence; and the stratum \code{order}.
#' @seealso \code{\link{bamp_coherent}}
#' @export
predict_coherent <- function(object, periods = 0, population = NULL,
                             quantiles = c(0.05, 0.5, 0.95),
                             hazard = FALSE, period_length = 1,
                             damping = 1, var_damping = 1) {
  if (!inherits(object, "apc_coherent")) stop("'object' must come from bamp_coherent().")
  hazard <- isTRUE(hazard)
  if (!(damping >= 0 && damping <= 1)) stop("'damping' must be in [0, 1].")
  if (!(var_damping > 0 && var_damping <= 1)) stop("'var_damping' must be in (0, 1].")
  if (!is.null(object$model$S) && object$model$S > 2L)
    return(.predict_coherent_general(object, periods, population, quantiles, hazard,
                                     period_length, damping, var_damping))
  s <- object$samples; md <- object$model; dat <- object$data
  I <- dat$I; J <- dat$J; K <- dat$K; M <- dat$periods_per_agegroup
  ord_p <- md$ord[2]; ord_c <- md$ord[3]
  n1 <- J; n2 <- J + periods; K2 <- (I - 1) * M + n2          # = coh(1, n2)
  D <- length(s$mu0); sgn <- c(1, -1)
  rho_vec <- if (!is.null(s$rho)) s$rho else rep(md$rho, D)   # per-draw rho (old objects: fixed)

  rateF <- rateM <- array(0, c(n2, I, D))
  dev_ext <- matrix(0, n2, D)
  for (d in seq_len(D)) {
    ph <- .cj_predict_rw(c(s$phi[d, ], numeric(n2 - n1)), s$lambda_phi[d], ord_p, n1, n2, damping, var_damping)
    ps <- .cj_predict_rw(c(s$psi[d, ], numeric(K2 - K)),  s$nu_psi[d],    ord_c, K, K2, damping, var_damping)
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


## ===========================================================================
## General S >= 3 strata: shared phi/psi + per-stratum contr.sum period deviation
## d_{s,j} (sum_s d_{s,j}=0 per period), proper AR1 prior (shared lambda_d, fixed
## rho). This is a DIFFERENT model from the S=2 +/-delta legacy (which pins
## sum_j delta=0); a_s is therefore only identified jointly with the deviation
## mean -- report a_s + mean_j d_{s,j} as the stratum level. Sampled rho and the
## cohort deviation are S=2-only for now. (See docs/hardening-plan.md, coh-Sgt2.)
## ===========================================================================
.bamp_coherent_general <- function(cases, population, age, period, cohort, ppa_arg,
                                   deviation, rho, mcmc, hyper, prior_scale, seed) {
  ord_a <- .coh_ord(age); ord_p <- .coh_ord(period); ord_c <- .coh_ord(cohort)
  S <- length(cases)
  Y <- lapply(cases, function(x) t(as.matrix(x)))
  Npop <- lapply(population, function(x) t(as.matrix(x)))
  I <- nrow(Y[[1]]); J <- ncol(Y[[1]]); M <- as.integer(ppa_arg); K <- M * (I - 1L) + J
  if (deviation == "ar1" && rho > 0.9)
    warning("bamp_coherent: rho > 0.9 barely penalises the common stratum shift; a_s vs deviation-mean will mix poorly. Consider rho <= 0.9.")
  scl <- function(Km, ord) if (prior_scale) Km * .pg_scale(Km, ord) else Km
  Ka <- scl(.pg_Kmat(I, ord_a), ord_a); Kp <- scl(.pg_Kmat(J, ord_p), ord_p); Kc <- scl(.pg_Kmat(K, ord_c), ord_c)
  Td <- .ar1_prec(J, rho)

  ## layout: mu0, a[S], theta[[s]][I], phi[J], psi[K], d[[s]][J]
  i_mu <- 1L; i_a <- 1L + seq_len(S); base <- 1L + S
  i_th <- lapply(seq_len(S), function(s) base + (s - 1L) * I + seq_len(I)); base <- base + S * I
  i_ph <- base + seq_len(J); base <- base + J
  i_ps <- base + seq_len(K); base <- base + K
  i_d  <- lapply(seq_len(S), function(s) base + (s - 1L) * J + seq_len(J)); P <- base + S * J

  grid <- expand.grid(j = 1:J, i = 1:I, s = 1:S); n <- nrow(grid); kcell <- (I - grid$i) * M + grid$j
  X <- matrix(0, n, P); X[, i_mu] <- 1
  X[cbind(seq_len(n), i_a[grid$s])] <- 1
  X[cbind(seq_len(n), vapply(seq_len(n), function(r) i_th[[grid$s[r]]][grid$i[r]], 1L))] <- 1
  X[cbind(seq_len(n), i_ph[grid$j])] <- 1
  X[cbind(seq_len(n), i_ps[kcell])] <- 1
  X[cbind(seq_len(n), vapply(seq_len(n), function(r) i_d[[grid$s[r]]][grid$j[r]], 1L))] <- 1   # +1, constraint makes the contrast
  yv <- Nv <- numeric(n)
  for (s in seq_len(S)) { sel <- grid$s == s
    yv[sel] <- Y[[s]][cbind(grid$i[sel], grid$j[sel])]; Nv[sel] <- Npop[[s]][cbind(grid$i[sel], grid$j[sel])] }
  Xt <- t(X); bvec <- as.numeric(Xt %*% (yv - Nv / 2))

  mkrow <- function(idx, val) { r <- numeric(P); r[idx] <- val; r }
  Arows <- list(mkrow(i_a, 1))
  for (s in seq_len(S)) Arows <- c(Arows, list(mkrow(i_th[[s]], 1)))
  Arows <- c(Arows, list(mkrow(i_ph, 1), mkrow(i_ps, 1)))
  for (j in seq_len(J)) Arows <- c(Arows, list(mkrow(vapply(seq_len(S), function(s) i_d[[s]][j], 1L), 1)))  # sum_s d_{s,j}=0
  if (ord_p == 2L) Arows <- c(Arows, list(mkrow(i_ph, (1:J) - mean(1:J))))
  A <- do.call(rbind, Arows); tA <- t(A)

  hp <- function(x) if (is.null(x)) c(1, 1) else x
  ha <- hp(hyper$age); hpp <- hp(hyper$period); hc <- hp(hyper$cohort); hd <- hp(hyper$dev); tau_a <- 1e-2
  build_prec <- function(kth, lph, nps, ld) {
    Pm <- matrix(0, P, P); Pm[i_mu, i_mu] <- 1e-6; Pm[cbind(i_a, i_a)] <- tau_a
    for (s in seq_len(S)) { Pm[i_th[[s]], i_th[[s]]] <- kth * Ka; Pm[i_d[[s]], i_d[[s]]] <- ld * Td }
    Pm[i_ph, i_ph] <- lph * Kp; Pm[i_ps, i_ps] <- nps * Kc
    diag(Pm) <- diag(Pm) + 1e-7; Pm
  }
  set.seed(seed)
  p0 <- sum(yv) / sum(Nv); beta <- numeric(P); beta[i_mu] <- log(p0 / (1 - p0))
  kth <- lph <- nps <- ld <- 1; Prec <- build_prec(kth, lph, nps, ld)
  iters <- mcmc$iterations; burn <- mcmc$burn_in; thin <- mcmc$thin
  keep <- seq.int(burn + thin, iters, by = thin); nkeep <- length(keep); st <- 0L
  out <- list(mu0 = numeric(nkeep), a = matrix(0, nkeep, S), theta = array(0, c(nkeep, S, I)),
              phi = matrix(0, nkeep, J), psi = matrix(0, nkeep, K), d = array(0, c(nkeep, S, J)),
              lambda_phi = numeric(nkeep), nu_psi = numeric(nkeep),
              lambda_d = numeric(nkeep), kappa_theta = numeric(nkeep))
  qform <- function(v, Km) sum(v * (Km %*% v))
  for (it in seq_len(iters)) {
    eta <- as.numeric(X %*% beta); omega <- .pg_rpg(Nv, eta)
    Q <- Xt %*% (X * omega) + Prec; beta <- .pg_draw_block(Q, bvec, A, tA)
    th_q <- sum(vapply(seq_len(S), function(s) qform(beta[i_th[[s]]], Ka), 0))
    d_q  <- sum(vapply(seq_len(S), function(s) qform(beta[i_d[[s]]], Td), 0))
    kth <- rgamma(1, ha[1]  + S * (I - ord_a) / 2, ha[2]  + 0.5 * th_q)
    lph <- rgamma(1, hpp[1] + (J - ord_p) / 2,     hpp[2] + 0.5 * qform(beta[i_ph], Kp))
    nps <- rgamma(1, hc[1]  + (K - ord_c) / 2,     hc[2]  + 0.5 * qform(beta[i_ps], Kc))
    ld  <- rgamma(1, hd[1]  + (S - 1) * J / 2,      hd[2]  + 0.5 * d_q)
    Prec <- build_prec(kth, lph, nps, ld)
    if (it %in% keep) { st <- st + 1L
      out$mu0[st] <- beta[i_mu]; out$a[st, ] <- beta[i_a]
      for (s in seq_len(S)) { out$theta[st, s, ] <- beta[i_th[[s]]]; out$d[st, s, ] <- beta[i_d[[s]]] }
      out$phi[st, ] <- beta[i_ph]; out$psi[st, ] <- beta[i_ps]
      out$lambda_phi[st] <- lph; out$nu_psi[st] <- nps; out$lambda_d[st] <- ld; out$kappa_theta[st] <- kth
    }
  }
  structure(list(samples = out,
    model = list(age = age, period = period, cohort = cohort, deviation = deviation,
                 rho = rho, S = S, ord = c(ord_a, ord_p, ord_c)),
    data = list(cases = cases, population = population, periods_per_agegroup = M,
                I = I, J = J, K = K, strata = names(cases))), class = "apc_coherent")
}

.predict_coherent_general <- function(object, periods, population, quantiles, hazard, period_length,
                                      damping = 1, var_damping = 1) {
  s <- object$samples; md <- object$model; dat <- object$data
  I <- dat$I; J <- dat$J; K <- dat$K; M <- dat$periods_per_agegroup; S <- md$S
  ord_p <- md$ord[2]; ord_c <- md$ord[3]; rho <- md$rho
  n1 <- J; n2 <- J + periods; K2 <- (I - 1) * M + n2; D <- length(s$mu0); strata <- dat$strata
  rate <- lapply(seq_len(S), function(x) array(0, c(n2, I, D)))
  for (d in seq_len(D)) {
    ph <- .cj_predict_rw(c(s$phi[d, ], numeric(n2 - n1)), s$lambda_phi[d], ord_p, n1, n2, damping, var_damping)
    ps <- .cj_predict_rw(c(s$psi[d, ], numeric(K2 - K)),  s$nu_psi[d],    ord_c, K, K2, damping, var_damping)
    ## per-stratum AR1 projection, then RE-CENTER per period so sum_s d_{s,j}=0 (no common-mode drift)
    dmat <- vapply(seq_len(S), function(st)
      .cj_predict_ar(c(s$d[d, st, ], numeric(n2 - n1)), s$lambda_d[d], rho, n1, n2), numeric(n2))  # [n2 x S]
    dmat <- dmat - rowMeans(dmat)
    for (i in 1:I) {
      kk <- (I - i) * M + (1:n2); base <- s$mu0[d] + ph + ps[kk]
      for (st in seq_len(S))
        rate[[st]][, i, d] <- 1 / (1 + exp(-(base + s$a[d, st] + s$theta[d, st, i] + dmat[, st])))
    }
  }
  if (is.null(population)) pop <- lapply(dat$population, as.matrix)
  else {
    if (!all(strata %in% names(population))) stop("'population' must be a named list with every stratum.")
    pop <- lapply(strata, function(nm) as.matrix(population[[nm]]))
  }
  np <- min(min(vapply(pop, nrow, 1L)), n2)
  qf <- function(arr) apply(arr, 1:2, stats::quantile, quantiles)
  pack <- function(rt) { e <- list(rate = qf(rt), samples = list(rate = rt))
    if (hazard) { hz <- -log1p(-pmin(pmax(rt, 0), 1 - 1e-10)) / period_length; e$hazard <- qf(hz); e$samples$hazard <- hz }; e }
  out <- stats::setNames(lapply(rate, pack), strata)
  wsum <- Reduce(`+`, lapply(seq_len(S), function(st) pop[[st]][seq_len(np), , drop = FALSE]))
  num <- Reduce(`+`, lapply(seq_len(S), function(st) {
    w <- pop[[st]][seq_len(np), , drop = FALSE]
    array(apply(rate[[st]][seq_len(np), , , drop = FALSE], 3, function(r) w * r), c(np, I, D)) }))
  rt_tot <- num / as.vector(wsum); rt_tot[!is.finite(rt_tot)] <- 0
  out$total <- pack(rt_tot); out$order <- strata
  out
}
