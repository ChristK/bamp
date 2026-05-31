## ===========================================================================
## PHASE 1 RESEARCH PROTOTYPE -- joint multinomial age-period-cohort model for
## several COMPETING CAUSES with cross-cause correlated trends.
##
## The C causes partition the deaths in each cell. We stick-break the total
## deaths into C-1 conditional binomial shares (Holmes-Held / Linderman
## multinomial-PG), each an APC field:
##
##   y_{c,i,j} ~ Binomial(R_{c,i,j}, pi_{c,i,j}),  R_1 = total deaths,
##   R_{c} = R_{c-1} - y_{c-1},
##   logit pi_{c,i,j} = mu^c + theta^c_i + phi^c_j + psi^c_k.
##
## Coherence is BY CONSTRUCTION: the implied cause shares sum to one, so the
## cause-specific rates/hazards sum to all-cause. The NEW ingredient over the
## Phase 0 bamp_strata fallback is a CROSS-CAUSE coupling of the period trends:
## the period effects of all causes share a multivariate random walk with
## innovation precision Omega ~ Wishart (prior precision K_p (x) Omega). Omega
## may be NEGATIVELY correlated -- cause replacement (one cause down, another up)
## -- which a shared-factor model cannot represent and independent fits ignore.
## Carried through prediction, it makes the projected cause trends move together.
##
## Auditable dense one-block Polya-Gamma Gibbs (reuses .pg_rpg / .pg_draw_block /
## .pg_Kmat / .pg_scale), no ASIS/Laplace-MH. O(P^3) per sweep -- a reference,
## not production. Age/cohort are cause-specific but not cross-cause coupled
## (period carries the dependence); that and the C++/sparse port are extensions.
## ===========================================================================

## multivariate random-walk extrapolation of the C-cause period field with
## Omega-correlated innovations (this is what carries cross-cause dependence
## into the forecast).
.mvrw_predict <- function(Phi, Omega, rw, n1, n2, damp = 1, var_damp = 1) {
  Cc <- ncol(Phi)
  if (n2 > n1) {
    L <- chol(solve(Omega))                       # cov = Omega^{-1};  z %*% L ~ N(0, cov)
    Phi <- rbind(Phi, matrix(0, n2 - n1, Cc))
    for (j in (n1 + 1):n2) {
      eps <- as.numeric(rnorm(Cc) %*% L) * var_damp^(j - n1 - 1)            # A3 variance shrink
      Phi[j, ] <- if (rw == 2) Phi[j - 1, ] + damp * (Phi[j - 1, ] - Phi[j - 2, ]) + eps  # A1 damped drift
                  else         Phi[j - 1, ] + eps
    }
  }
  Phi
}

## one Gibbs sweep of Bayesian factor analysis on the period increments Y (N x Cm,
## rows ~ N(0, Sigma)), returning the implied LOW-RANK cross-cause precision
## Omega = Sigma^{-1}, Sigma = Lam Lam' + diag(Psi). This replaces the full Wishart
## Omega with a few latent factors -- the shared, cross-cutting (cross-group)
## drivers a disease taxonomy cannot express. Only Sigma is used downstream, and it
## is rotation-invariant, so the loadings need no identification constraint.
.mc_fa_omega <- function(Y, Lam, Psi, cprior = 1, a_psi = 2, b_psi = 1) {
  v <- mean(Y^2) + 1e-12                     # overall increment scale; run FA standardised so the
  Ys <- Y / sqrt(v)                          # O(1) priors are scale-appropriate, then rescale Sigma
  N <- nrow(Ys); C <- ncol(Ys); R <- ncol(Lam); Pinv <- 1 / Psi
  V <- chol2inv(chol(diag(R) + crossprod(Lam, Lam * Pinv)))
  H <- (Ys * rep(Pinv, each = N)) %*% Lam %*% V + matrix(rnorm(N * R), N, R) %*% chol(V)  # scores
  HtH <- crossprod(H)
  for (c in seq_len(C)) {
    Q <- HtH / Psi[c] + diag(R) / cprior
    Lam[c, ] <- solve(Q, crossprod(H, Ys[, c]) / Psi[c]) + backsolve(chol(Q), rnorm(R))
    res <- Ys[, c] - H %*% Lam[c, ]
    Psi[c] <- 1 / rgamma(1, a_psi + N / 2, b_psi + 0.5 * sum(res^2))
  }
  list(Lam = Lam, Psi = Psi,                 # Sigma = v*(LL'+Psi); Omega = Sigma^{-1}
       Omega = chol2inv(chol(tcrossprod(Lam) + diag(Psi))) / v, loadings = Lam * sqrt(v))
}

#' Joint multinomial APC model for competing causes (Phase 1 prototype)
#'
#' @description
#' Fit several competing causes of death in ONE joint posterior, coherent with all-cause by
#' construction and with cross-cause correlated trends. The deaths in each cell are stick-broken into
#' \code{C - 1} conditional binomial-logit APC shares; the cause period trends share a multivariate
#' random walk whose innovation precision \code{Omega} (Wishart prior) captures cross-cause dependence
#' -- including negative correlation (cause replacement). This is the principled alternative to the
#' Phase 0 \code{\link{bamp_strata}} fallback for causes; see \code{docs/coherent-forecasting.md}.
#'
#' Research-grade reference (dense Polya-Gamma Gibbs, period-only cross-cause coupling); correct but
#' not optimised.
#'
#' @param cases named list of \code{C} cause-count matrices (\code{[periods x agegroups]}) that
#'   partition the deaths; their sum is the all-cause total. The list order is the stick-breaking
#'   order (put the most prevalent causes first).
#' @param population a single shared \code{[periods x agegroups]} population matrix (all causes share
#'   the population at risk).
#' @param age,period,cohort each \code{"rw1"} or \code{"rw2"}.
#' @param periods_per_agegroup integer M (cohort index), as in \code{\link{bamp}}.
#' @param order stick-breaking cause order: \code{"prevalence"} (default; most prevalent first, so the
#'   running remainders stay large) or a permutation of the cause names/indices. Coherence holds for
#'   any order; \code{\link{predict_multicause}} returns causes in the original input order.
#' @param mcmc list with \code{iterations}, \code{burn_in}, \code{thin}.
#' @param hyper Gamma hyperparameter \code{age} (per-cause age precision) and the cross-cause Wishart
#'   parameters \code{omega} (period) and \code{omega_c} (cohort), each \code{c(df_extra, v0)}.
#' @param prior_scale Sorbye-Rue scaling of intrinsic structure matrices.
#' @param seed RNG seed.
#' @param factor optional integer \code{R}: use a LOW-RANK factor model for the cross-cause period
#'   covariance (\eqn{\Sigma=\Lambda\Lambda^\top+\Psi}, \code{R} latent factors) instead of the full
#'   Wishart. The factors are shared, \emph{cross-cutting} latent drivers (risk-factor proxies) that a
#'   disease taxonomy/cascade cannot express; they make the coupling identifiable for many causes
#'   (\code{R}\eqn{\times}\code{C} loadings, not \code{C^2/2} correlations). Posterior-mean loadings are
#'   returned by \code{\link{predict_multicause}} as \code{loadings}. \code{NULL} (default) keeps the
#'   full Wishart coupling.
#' @param engine \code{"dense"} (default) or \code{"sparse"}. The sparse engine assembles the
#'   one-block precision sparsely and uses a sparse Cholesky (Matrix package) -- much faster for many
#'   causes (the dense \eqn{O(P^3)} Cholesky is the bottleneck). Same posterior; draws differ from the
#'   dense path (a fill-reducing permutation reorders the RNG), so results match in distribution, not
#'   bit-for-bit.
#'
#' @return object of class \code{apc_multicause}: posterior draws (incl. the cross-cause precision
#'   \code{Omega}), the fitted all-cause total model, and metadata; used by \code{\link{predict_multicause}}.
#' @seealso \code{\link{predict_multicause}}, \code{\link{bamp_strata}}, \code{\link{reconcile_apc}}
#' @export
bamp_multicause <- function(cases, population, age = "rw1", period = "rw1", cohort = "rw1",
                            periods_per_agegroup, order = "prevalence",
                            mcmc = list(iterations = 4000, burn_in = 1000, thin = 2),
                            hyper = list(age = c(1, 0.5), omega = c(2, 1e-4), omega_c = c(2, 1e-4)),
                            prior_scale = TRUE, seed = 1, factor = NULL,
                            engine = c("dense", "sparse")) {
  if (!is.list(cases) || length(cases) < 2L) stop("'cases' must be a list of >= 2 cause matrices.")
  if (is.null(names(cases))) names(cases) <- paste0("cause", seq_along(cases))
  ## stick-breaking is order-dependent; default to most-prevalent-first so the running
  ## remainders stay large (rarest cause = unmodelled reference). Coherence holds for any order.
  causes_orig <- names(cases)
  perm <- if (identical(order, "prevalence"))
            order(vapply(cases, function(x) sum(as.matrix(x)), 0), decreasing = TRUE)
          else if (is.numeric(order)) as.integer(order)
          else if (is.character(order) && all(order %in% causes_orig)) match(order, causes_orig)
          else seq_along(cases)
  if (!setequal(perm, seq_along(cases)))
    stop("'order' must be \"prevalence\" or a permutation of the cause names/indices.")
  cases <- cases[perm]
  ord_a <- .coh_ord(age); ord_p <- .coh_ord(period); ord_c <- .coh_ord(cohort)
  Cn <- length(cases); Cm <- Cn - 1L                  # number of stick-breaking shares
  engine <- match.arg(engine)                         # 'sparse' = sparse Cholesky one-block draw
  if (engine == "sparse" && !requireNamespace("Matrix", quietly = TRUE))
    stop("engine = 'sparse' requires the Matrix package.")
  ## optional LOW-RANK factor coupling on the period covariance (cross-cutting drivers)
  use_factor <- !is.null(factor) && Cm >= 2L
  if (use_factor) {
    Rfac <- as.integer(factor)
    if (Rfac < 1L || Rfac >= Cm)
      stop("'factor' must be an integer in 1..(number of causes - 2) for a low-rank coupling.")
  }
  Y <- lapply(cases, function(x) t(as.matrix(x)))     # [age x period]
  Npop <- t(as.matrix(population))
  I <- nrow(Y[[1]]); J <- ncol(Y[[1]]); M <- as.integer(periods_per_agegroup)
  K <- M * (I - 1L) + J

  ## stick-breaking remainders: R_1 = total deaths; R_c = R_{c-1} - Y_{c-1}
  Ytot <- Reduce(`+`, Y)
  Rl <- vector("list", Cn); Rl[[1]] <- Ytot
  for (c in 2:Cn) Rl[[c]] <- Rl[[c - 1]] - Y[[c - 1]]

  scl <- function(Km, ord) if (prior_scale) Km * .pg_scale(Km, ord) else Km
  Ka <- scl(.pg_Kmat(I, ord_a), ord_a)
  Kp <- scl(.pg_Kmat(J, ord_p), ord_p)
  Kc <- scl(.pg_Kmat(K, ord_c), ord_c)
  if (use_factor) {                                   # increment operator + scale, so Cov(increments)=Sigma
    Dop_p <- diff(diag(J), differences = ord_p)
    sp_p  <- if (prior_scale) .pg_scale(.pg_Kmat(J, ord_p), ord_p) else 1
  }

  ## layout: beta = (mu[Cm], theta^1..^Cm[I], phi (cause-major Cm*J), psi^1..^Cm[K])
  i_mu <- seq_len(Cm)
  i_th <- lapply(seq_len(Cm), function(c) Cm + (c - 1L) * I + seq_len(I))
  base_ph <- Cm + Cm * I
  i_ph <- lapply(seq_len(Cm), function(c) base_ph + (c - 1L) * J + seq_len(J))  # cause-major
  i_ph_all <- unlist(i_ph)
  base_ps <- base_ph + Cm * J
  i_ps <- lapply(seq_len(Cm), function(c) base_ps + (c - 1L) * K + seq_len(K))  # cause-major
  i_ps_all <- unlist(i_ps)
  P <- base_ps + Cm * K

  ## explicit design + stacked data over cells (c, i, j)
  grid <- expand.grid(j = 1:J, i = 1:I, c = 1:Cm)
  n <- nrow(grid); kcell <- (I - grid$i) * M + grid$j
  col_mu <- i_mu[grid$c]                                            # the 4 non-zero design columns/cell
  col_th <- vapply(seq_len(n), function(r) i_th[[grid$c[r]]][grid$i[r]], 1L)
  col_ph <- vapply(seq_len(n), function(r) i_ph[[grid$c[r]]][grid$j[r]], 1L)
  col_ps <- vapply(seq_len(n), function(r) i_ps[[grid$c[r]]][kcell[r]], 1L)
  yv <- Nv <- numeric(n)
  for (c in seq_len(Cm)) {
    sel <- grid$c == c
    yv[sel] <- Y[[c]][cbind(grid$i[sel], grid$j[sel])]
    Nv[sel] <- Rl[[c]][cbind(grid$i[sel], grid$j[sel])]
  }
  if (engine == "sparse") {                                        # sparse design + banded structure matrices
    Xs <- Matrix::sparseMatrix(i = rep(seq_len(n), 4L), j = c(col_mu, col_th, col_ph, col_ps),
                               x = 1, dims = c(n, P))
    bvec <- as.numeric(Matrix::crossprod(Xs, yv - Nv / 2))
    Ka_s <- Matrix::Matrix(Ka, sparse = TRUE); Kp_s <- Matrix::Matrix(Kp, sparse = TRUE)
    Kc_s <- Matrix::Matrix(Kc, sparse = TRUE)
  } else {
    X <- matrix(0, n, P)
    X[cbind(seq_len(n), col_mu)] <- 1; X[cbind(seq_len(n), col_th)] <- 1
    X[cbind(seq_len(n), col_ph)] <- 1; X[cbind(seq_len(n), col_ps)] <- 1
    Xt <- t(X); bvec <- as.numeric(Xt %*% (yv - Nv / 2))
  }

  ## sum-to-zero per effect per cause (+ RW2 period zero-slope per cause)
  mkrow <- function(idx, val) { r <- numeric(P); r[idx] <- val; r }
  Arows <- list()
  for (c in seq_len(Cm)) Arows <- c(Arows, list(mkrow(i_th[[c]], 1), mkrow(i_ph[[c]], 1), mkrow(i_ps[[c]], 1)))
  if (ord_p == 2L) for (c in seq_len(Cm)) Arows <- c(Arows, list(mkrow(i_ph[[c]], (1:J) - mean(1:J))))
  if (ord_c == 2L) for (c in seq_len(Cm)) Arows <- c(Arows, list(mkrow(i_ps[[c]], (1:K) - mean(1:K))))  # pin coupled cohort drift
  A <- do.call(rbind, Arows); tA <- t(A)

  ha <- hyper$age; ho <- if (is.null(hyper$omega)) c(2, 1e-4) else hyper$omega
  hoc <- if (is.null(hyper$omega_c)) c(2, 1e-4) else hyper$omega_c
  nu0   <- Cm + ho[1];  V0inv   <- diag(Cm) * ho[2];  rank_p <- J - ord_p
  nu0_c <- Cm + hoc[1]; V0inv_c <- diag(Cm) * hoc[2]; rank_c <- K - ord_c
  tau_mu <- 1e-2
  build_prec <- function(kth, Omega, Omega_psi) {
    Pm <- matrix(0, P, P)
    Pm[cbind(i_mu, i_mu)] <- tau_mu
    for (c in seq_len(Cm)) Pm[i_th[[c]], i_th[[c]]] <- kth[c] * Ka
    Pm[i_ph_all, i_ph_all] <- kronecker(Omega, Kp)        # cross-cause coupled PERIOD prior
    Pm[i_ps_all, i_ps_all] <- kronecker(Omega_psi, Kc)    # cross-cause coupled COHORT prior
    diag(Pm) <- diag(Pm) + 1e-7
    Pm
  }
  ## sparse Prec: contiguous blocks (mu | theta | phi cause-major | psi cause-major); the
  ## kronecker(Omega, Kp_s) period/cohort blocks are banded-in-K and dense-across-cause.
  build_prec_sparse <- function(kth, Omega, Omega_psi)
    Matrix::bdiag(Matrix::Diagonal(Cm, tau_mu),
                  Matrix::bdiag(lapply(seq_len(Cm), function(c) kth[c] * Ka_s)),
                  kronecker(Matrix::Matrix(Omega, sparse = TRUE), Kp_s),
                  kronecker(Matrix::Matrix(Omega_psi, sparse = TRUE), Kc_s)) +
      Matrix::Diagonal(P, 1e-7)

  set.seed(seed)
  beta <- numeric(P)
  kth <- rep(1, Cm); Omega <- diag(Cm); Omega_psi <- diag(Cm)
  if (use_factor) { Lam_p <- matrix(rnorm(Cm * Rfac, 0, 0.3), Cm, Rfac); Psi_p <- rep(1, Cm) }
  if (engine == "dense") Prec <- build_prec(kth, Omega, Omega_psi)
  iters <- mcmc$iterations; burn <- mcmc$burn_in; thin <- mcmc$thin
  keep <- seq.int(burn + thin, iters, by = thin); nkeep <- length(keep); st <- 0L
  out <- list(mu = matrix(0, nkeep, Cm), theta = array(0, c(nkeep, Cm, I)),
              phi = array(0, c(nkeep, J, Cm)), psi = array(0, c(nkeep, K, Cm)),
              kappa_theta = matrix(0, nkeep, Cm),
              Omega = array(0, c(nkeep, Cm, Cm)), Omega_psi = array(0, c(nkeep, Cm, Cm)))
  if (use_factor) out$Lambda <- array(0, c(nkeep, Cm, Rfac))
  qf <- function(v, Km) sum(v * (Km %*% v))
  for (it in seq_len(iters)) {
    if (engine == "sparse") {
      eta <- as.numeric(Xs %*% beta)
      omega <- .pg_rpg(Nv, eta)
      Qs <- Matrix::crossprod(Xs, Matrix::Diagonal(x = omega) %*% Xs) +
            build_prec_sparse(kth, Omega, Omega_psi)
      beta <- .pg_draw_block_sparse(Qs, bvec, A, tA)
    } else {
      eta <- as.numeric(X %*% beta)
      omega <- .pg_rpg(Nv, eta)
      Q <- Xt %*% (X * omega) + Prec
      beta <- .pg_draw_block(Q, bvec, A, tA)
    }
    Phi <- matrix(beta[i_ph_all], J, Cm)               # [period x cause]
    Psi <- matrix(beta[i_ps_all], K, Cm)               # [cohort x cause]
    for (c in seq_len(Cm))
      kth[c] <- rgamma(1, ha[1] + (I - ord_a) / 2, ha[2] + 0.5 * qf(beta[i_th[[c]]], Ka))
    if (use_factor) {                                  # low-rank factor coupling (period)
      fa <- .mc_fa_omega(sqrt(sp_p) * (Dop_p %*% Phi), Lam_p, Psi_p)
      Lam_p <- fa$Lam; Psi_p <- fa$Psi; Omega <- fa$Omega; Lam_store <- fa$loadings
    } else {                                            # full Wishart coupling (period)
      Omega <- stats::rWishart(1, nu0 + rank_p, solve(V0inv + crossprod(Phi, Kp %*% Phi)))[, , 1]
    }
    Omega_psi <- stats::rWishart(1, nu0_c + rank_c, solve(V0inv_c + crossprod(Psi, Kc %*% Psi)))[, , 1]
    if (engine == "dense") Prec <- build_prec(kth, Omega, Omega_psi)
    if (it %in% keep) {
      st <- st + 1L
      out$mu[st, ] <- beta[i_mu]
      for (c in seq_len(Cm)) out$theta[st, c, ] <- beta[i_th[[c]]]
      out$phi[st, , ] <- Phi; out$psi[st, , ] <- Psi; out$kappa_theta[st, ] <- kth
      out$Omega[st, , ] <- Omega; out$Omega_psi[st, , ] <- Omega_psi
      if (use_factor) out$Lambda[st, , ] <- Lam_store    # loadings on the original increment scale
    }
  }

  message("bamp_multicause: fitting all-cause total model ...")
  total_cases <- Reduce(`+`, lapply(cases, as.matrix))
  fit_total <- bamp(total_cases, population, age = age, period = period, cohort = cohort,
                    periods_per_agegroup = M,
                    mcmc.options = list(number_of_iterations = max(2000, iters %/% 2),
                                        burn_in = burn, step = thin, tuning = max(100, burn %/% 4)),
                    parallel = FALSE, verbose = FALSE)

  structure(list(
    samples = out, total = fit_total,
    model = list(age = age, period = period, cohort = cohort, ord = c(ord_a, ord_p, ord_c),
                 factor = if (use_factor) Rfac else NULL),
    data = list(cases = cases, population = population, periods_per_agegroup = M,
                I = I, J = J, K = K, causes = names(cases),
                causes_orig = causes_orig, perm = perm)
  ), class = "apc_multicause")
}


#' Coherent projection of competing causes with correlated trends
#'
#' @description
#' Project cause-specific rates and additive hazards from \code{\link{bamp_multicause}}. The all-cause
#' total is projected with the ordinary random walk; the cause shares are projected with a
#' multivariate random walk whose innovations are correlated across causes by the posterior \code{Omega}
#' (so the cause trends move together as estimated). Cause rates and hazards are
#' \code{share x total}, hence coherent with all-cause by construction.
#'
#' @param object an \code{apc_multicause} object.
#' @param periods number of future periods to project.
#' @param population optional future \code{[periods x agegroups]} population for the all-cause total;
#'   \code{NULL} uses the in-sample population.
#' @param quantiles quantiles to summarise.
#' @param hazard,period_length if \code{hazard=TRUE} also return cause-specific additive hazards
#'   \code{share x (-log(1-total_rate)/period_length)} (they sum to the all-cause hazard).
#' @param damping,var_damping trend-damping and innovation-variance shrinkage for the forecast
#'   extrapolation (the all-cause total \emph{and} the cause shares), as in \code{\link{predict_apc}}.
#'   Defaults \code{1, 1} reproduce the free random walk; \code{damping < 1} curbs long-horizon
#'   over-extrapolation, \code{var_damping < 1} stops the predictive bands fanning out without bound.
#'
#' @return list with one entry per cause and a \code{total} entry (quantiles of \code{rate}, and
#'   \code{hazard} if requested, on the \code{[period, agegroup]} grid, plus \code{samples}); the
#'   posterior cross-cause correlation \code{cor_omega}; and \code{coherence_maxerr} (largest deviation
#'   of summed cause rates from the total; ~0 by construction).
#' @seealso \code{\link{bamp_multicause}}
#' @export
predict_multicause <- function(object, periods = 0, population = NULL,
                               quantiles = c(0.05, 0.5, 0.95),
                               hazard = FALSE, period_length = 1,
                               damping = 1, var_damping = 1) {
  if (!(damping >= 0 && damping <= 1)) stop("'damping' must be in [0, 1].")
  if (!(var_damping > 0 && var_damping <= 1)) stop("'var_damping' must be in (0, 1].")
  if (!inherits(object, "apc_multicause")) stop("'object' must come from bamp_multicause().")
  hazard <- isTRUE(hazard)
  s <- object$samples; md <- object$model; dat <- object$data
  I <- dat$I; J <- dat$J; K <- dat$K; M <- dat$periods_per_agegroup
  ord_p <- md$ord[2]; ord_c <- md$ord[3]
  Cm <- length(dat$causes) - 1L; Cn <- Cm + 1L; causes <- dat$causes
  n1 <- J; n2 <- J + periods; K2 <- (I - 1) * M + n2

  ## all-cause total rate (separate fit), projected
  pt <- predict_apc(object$total, periods = periods, population = population,
                    damping = damping, var_damping = var_damping)
  totrate <- pt$samples$pr                            # [period, age, Dtot]
  D <- min(dim(s$phi)[1], dim(totrate)[3])
  totrate <- totrate[, , seq_len(D), drop = FALSE]

  rate <- lapply(seq_len(Cn), function(c) array(0, c(n2, I, D)))
  for (d in seq_len(D)) {
    Phi <- .mvrw_predict(matrix(s$phi[d, , ], J, Cm), matrix(s$Omega[d, , ], Cm, Cm), ord_p, n1, n2, damping, var_damping)
    if (!is.null(s$Omega_psi)) {                        # cross-cause coupled cohort projection
      Psi <- .mvrw_predict(matrix(s$psi[d, , ], K, Cm), matrix(s$Omega_psi[d, , ], Cm, Cm), ord_c, K, K2, damping, var_damping)
    } else {                                            # old objects: per-cause free RW (nu_psi)
      Psi <- vapply(seq_len(Cm), function(c)
        .cj_predict_rw(c(s$psi[d, c, ], numeric(K2 - K)), s$nu_psi[d, c], ord_c, K, K2, damping, var_damping), numeric(K2))
    }
    for (i in 1:I) {
      kk <- (I - i) * M + (1:n2)
      pic <- vapply(seq_len(Cm), function(c)
        1 / (1 + exp(-(s$mu[d, c] + s$theta[d, c, i] + Phi[, c] + Psi[kk, c]))), numeric(n2))  # [period x cause]
      rem <- rep(1, n2)                               # stick-breaking -> cause shares
      for (c in seq_len(Cm)) { sh <- pic[, c] * rem; rate[[c]][, i, d] <- sh * totrate[, i, d]; rem <- rem - sh }
      rate[[Cn]][, i, d] <- rem * totrate[, i, d]
    }
  }

  np <- dim(totrate)[1]
  qf <- function(arr) apply(arr, 1:2, stats::quantile, quantiles)
  clip <- function(x) pmin(pmax(x, 0), 1 - 1e-10)
  tot_haz <- if (hazard) -log1p(-clip(totrate)) / period_length else NULL
  pack <- function(rt, sh = NULL) {
    e <- list(rate = qf(rt), samples = list(rate = rt))
    if (hazard) { hz <- if (is.null(sh)) -log1p(-clip(rt)) / period_length else sh * tot_haz
                  e$hazard <- qf(hz); e$samples$hazard <- hz }
    e
  }
  out <- stats::setNames(lapply(seq_len(Cn), function(c) {
    sh <- if (hazard) { z <- rate[[c]] / totrate; z[!is.finite(z)] <- 0; z } else NULL
    pack(rate[[c]], sh)
  }), causes)
  out$total <- pack(totrate)

  summed <- Reduce(`+`, lapply(out[causes], function(e) e$samples$rate))
  out$coherence_maxerr <- max(abs(summed - totrate))
  ## posterior-mean cross-cause trend correlation (drop-robust at Cm = 1)
  cormean <- function(Om) {
    if (Cm == 1L) return(matrix(1, 1, 1))
    M <- vapply(seq_len(dim(Om)[1]),
                function(d) stats::cov2cor(solve(matrix(Om[d, , ], Cm, Cm))), numeric(Cm * Cm))
    matrix(rowMeans(M), Cm, Cm)
  }
  out$cor_omega <- cormean(s$Omega)                   # period trend correlation (stick-breaking order)
  dimnames(out$cor_omega) <- list(causes[seq_len(Cm)], causes[seq_len(Cm)])
  if (!is.null(s$Omega_psi)) {
    out$cor_omega_psi <- cormean(s$Omega_psi)         # cohort trend correlation (stick-breaking order)
    dimnames(out$cor_omega_psi) <- list(causes[seq_len(Cm)], causes[seq_len(Cm)])
  }
  if (!is.null(s$Lambda)) {                           # factor model: posterior-mean loadings on the latent drivers
    out$loadings <- apply(s$Lambda, 2:3, mean)
    dimnames(out$loadings) <- list(causes[seq_len(Cm)], paste0("factor", seq_len(dim(s$Lambda)[3])))
  }
  ## cause rate/hazard entries are accessible by name; report causes in the user's original order
  out$causes <- if (!is.null(dat$causes_orig)) dat$causes_orig else causes
  out
}


#' Stick-breaking order sensitivity for a competing-cause model
#'
#' @description
#' Refit \code{\link{bamp_multicause}} under alternative cause orderings and compare the implied
#' cause-SHARE forecasts (order-invariant in expectation, unlike the raw stick-breaking parameters).
#' A large discrepancy signals the ordering matters for this data and a more order-symmetric model
#' may be preferable. The diagnostic detects, it does not prove, order effects; differences at the
#' Monte-Carlo-noise level are not meaningful (use enough iterations).
#'
#' @param object an \code{apc_multicause} object.
#' @param orders named list of alternative orderings (each a permutation of the original cause names
#'   or indices); default compares the fitted order against its reverse.
#' @param periods horizon for the compared share forecast.
#' @param ... passed to \code{\link{bamp_multicause}} for the refits (use the same \code{mcmc}).
#'
#' @return data.frame with, per alternative ordering, the \code{max_abs_share_diff} and
#'   \code{mean_abs_share_diff} of the posterior-mean cause-share forecast vs the object's ordering.
#' @seealso \code{\link{bamp_multicause}}
#' @export
order_sensitivity <- function(object, orders = NULL, periods = 0, ...) {
  if (!inherits(object, "apc_multicause")) stop("'object' must come from bamp_multicause().")
  dat <- object$data
  cases0 <- dat$cases[order(dat$perm)]; names(cases0) <- dat$causes_orig   # original input order
  n2 <- dat$J + periods
  mean_share <- function(fit) {
    pr <- predict_multicause(fit, periods = periods)
    vapply(dat$causes_orig, function(nm) {
      z <- pr[[nm]]$samples$rate / pr$total$samples$rate; z[!is.finite(z)] <- 0
      apply(z, 1:2, mean)
    }, matrix(0, n2, dat$I))
  }
  base_sh <- mean_share(object)
  if (is.null(orders)) orders <- list(reverse = rev(dat$causes_orig))
  do.call(rbind, lapply(seq_along(orders), function(k) {
    fk <- bamp_multicause(cases0, dat$population, age = object$model$age, period = object$model$period,
                          cohort = object$model$cohort, periods_per_agegroup = dat$periods_per_agegroup,
                          order = orders[[k]], ...)
    d <- abs(mean_share(fk) - base_sh)
    data.frame(order = if (!is.null(names(orders))) names(orders)[k] else as.character(k),
               max_abs_share_diff = max(d), mean_abs_share_diff = mean(d), stringsAsFactors = FALSE)
  }))
}
