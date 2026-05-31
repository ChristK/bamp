#' Prediction for age-period-cohort models
#'
#' @param object apc object
#' @param periods number of periods to predict
#' @param population matrix of (predicted) population, if NULL, population data from original bamp call will be used
#' @param quantiles vector of quantiles to compute
#' @param update boolean. If TRUE, object will be returned with results added to the object
#' @param hazard boolean. If TRUE, additionally return the cumulative cause-specific hazard
#' \eqn{h = -\log(1-p)} (per \code{period_length}) alongside the probabilities. Hazards are additive
#' across competing causes, so these are the quantity to feed to a competing-risk life-table or
#' microsimulation (see Details).
#' @param period_length single positive number: the length of one period in the time units you want
#' the hazard expressed per (e.g. years per period). Default 1 returns the cumulative hazard per period.
#' @param damping trend-damping factor in \[0,1\] for the RW2 drift when extrapolating future period/
#' cohort effects: each future increment is \code{damping} times the previous one (Gardner-McKenzie
#' damped trend). \code{1} (default) is the usual linear RW2 continuation; \code{0} collapses to a flat
#' RW1-style forecast; intermediate values curb long-horizon over-extrapolation. No effect on RW1 point
#' forecasts.
#' @param var_damping innovation-variance shrinkage in (0,1\]: the per-step forecast innovation sd is
#' multiplied by \code{var_damping^(h-1)} at horizon \code{h}, so predictive bands stop fanning out
#' without bound. \code{1} (default) leaves the bands unchanged.
#' @param trend_window for an RW2 period prior, base the forecast drift on the mean increment over the
#' last \code{trend_window} fitted periods (a recent/post-break slope) instead of the single last
#' increment (A5). \code{NULL} (default) keeps the standard RW2 continuation. See
#' \code{\link{changepoint_window}} to choose it from a detected structural break.
#' Only used when \code{hazard=TRUE}.
#'
#' @description Prediction of rates and, if possible, cases from the Bayesian age-period-cohort model
#' using the prior assumptions (random walks) of the model and the estimated variance of the random walk.
#' For example, random walk of first order (rw1) for period effect predicts constant effects for future periods plus noise.
#' 
#' @details This function will return predicted rates for future periods. For this, future period and cohort effects will be predicted.
#' Further age group effects will not be predicted. The rates are random samples from the predictive distribution; number of samples is equal
#' to number of MCMC iterations. Quantiles will be provided for convenience, but all samples are available.
#' If population numbers are given, number of cases will also be predicted. Number of cases
#' will not only be predicted for future periods,
#' but also for the time periods where data are available; this can be used for model assessment.
#'
#' When forecasting several competing causes of death (or disease), fit one model per cause
#' (cause-specific events out of the at-risk population, treating other-cause events as ordinary
#' survivors), call \code{predict_apc(..., hazard=TRUE)} for each, and pass the resulting
#' \emph{hazards} to your life-table or microsimulation. Cause-specific hazards add up
#' (\eqn{H_{all-cause}=\sum_c H_c}); the implied all-cause risk is
#' \eqn{1-\exp(-\sum_c H_c)=1-\prod_c (1-p_c)}, which is \emph{not} \eqn{\sum_c p_c}. Working on the
#' probability scale and summing would double-count the shared population at risk. The mechanical
#' competing-risk coupling (shared survivors) is then handled by the downstream life-table, not here.
#'
#' @return list with quantiles of predicted probabilities (\code{pr}), predicted cases (\code{cases}) and predicted cases per period (\code{cases_period})
#' and a list samples with MCMC samples of pr, cases and cases_period.
#' If \code{hazard=TRUE}, the cumulative cause-specific hazard is added as \code{hazard} (quantiles) and
#' \code{samples$hazard}, on the same \code{[period, age]} grid as \code{pr}.
#' If \code{update=TRUE}, the apc object will be returned with this list (predicted) added.
#' @seealso \code{vignette("prediction", package = "bamp")}
#' @import parallel
#' @importFrom abind abind
#' @export
#'
#' @examples
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' pred <- predict_apc(model, periods=1)
#' plot(pred$pr[2,11,], main="Predicted rate per agegroup", ylab="p")
#' }
predict_apc<-function(object, periods=0, population=NULL, quantiles=c(0.05,0.5,0.95), update=FALSE, hazard=FALSE, period_length=1, damping=1, var_damping=1, trend_window=NULL){
  if (!(is.numeric(damping) && length(damping)==1L && damping>=0 && damping<=1)) stop("'damping' must be in [0, 1].")
  if (!(is.numeric(var_damping) && length(var_damping)==1L && var_damping>0 && var_damping<=1)) stop("'var_damping' must be in (0, 1].")
  if (!is.null(trend_window) && !(is.numeric(trend_window) && length(trend_window)==1L && trend_window>=1)) stop("'trend_window' must be NULL or a positive integer.")
  ksi_prognose <-
    function(prepi, vdb, noa, nop, nop2, noc, zmode){
      my<-prepi[1]
      theta<-prepi[2:(noa+1)]
      phi<-prepi[(noa+2):(noa+nop2+1)]
      psi<-prepi[(noa+nop2+2):(noa+nop2+noc+1)]
      delta<-prepi[noa+nop2+noc+2]
      
      ksi<-array(0, c(nop2,noa))
      for(i in 1:noa){
        for(j in 1:nop2){
          ksi[j,i] <- my + theta[i] + phi[j] + psi[bamp::coh(i,j,noa,vdb)]
          
          if(zmode){
            ksi[j,i] <- ksi[j,i] + (rnorm(1, mean = 0, sd = 1)/sqrt(delta))
          }
        }
      }
      
      return(ksi)
    }
  
  predict_rw <-
    function(prepi, rw, n1, n2, damp=1, vdamp=1, dwin=NULL){
      lambda<-prepi[1]
      phi<-prepi[-1]
      # nothing to extrapolate when there are no future steps (n2==n1). Guard the
      # (n1+1):n2 loops: R's `:` counts DOWN when n2<n1+1, so (n1+1):n1 would be
      # c(n1+1, n1) and wrongly append/overwrite a step. This matters for
      # predict_apc(periods=0) (retrospective model checking), where n2==n1.
      # damp (A1) damps the RW2 drift; vdamp (A3) shrinks the per-step innovation
      # sd geometrically; dwin (A5) bases the RW2 drift on the mean increment over
      # the last dwin periods. damp=1, vdamp=1, dwin=NULL reproduce the free RW.
      if(n2 > n1){
        if(rw == 1){
          for(i in (n1+1):n2){
            phi[i] <- (rnorm(1, mean = 0, sd = 1)/sqrt(lambda))*vdamp^(i-n1-1) + phi[i-1]
          }
        }

        if(rw == 2 && !is.null(dwin)){
          w <- max(1L, min(dwin, n1-1L)); g <- mean(diff(phi[(n1-w):n1]))
          for(i in (n1+1):n2){
            phi[i] <- phi[i-1] + g + (rnorm(1, mean = 0, sd = 1)/sqrt(lambda))*vdamp^(i-n1-1)
            g <- g*damp
          }
        } else if(rw == 2){
          for(i in (n1+1):n2){
            phi[i] <- (rnorm(1, mean = 0, sd = 1)/sqrt(lambda))*vdamp^(i-n1-1) + phi[i-1] + damp*(phi[i-1] - phi[i-2])
          }
        }

        if(rw == 0){
          for(i in (n1+1):n2){
            phi[i] <- phi[i-1]
          }
        }
      }

      return(phi)
    }
  
  hazard <- isTRUE(hazard)
  if (hazard && (!is.numeric(period_length) || length(period_length)!=1L || is.na(period_length) || period_length<=0))
    stop("'period_length' must be a single positive number (length of one period in the time units you want the hazard expressed per).")

  phi<-psi<-NA

  a1<-dim(object$data$cases)[1]
  n1<-dim(object$data$cases)[2]
  n2<-n1+periods
  rwp<-rwc<-0
  
  if (!object$model$period=="")
  {
    rwp<-switch(object$model$period,
                rw1 = 1,
                rw2 = 2
    )
    ch<-length(object$samples$period)
      prep<-parallel::mclapply(1:ch, function(i,samples)cbind(object$samples$period_parameter[[i]],object$samples$period[[i]]), samples)
      phi<-parallel::mclapply(prep, function(prepi, rw, n1, n2, damp, vdamp, dwin){t(apply(prepi, 1, predict_rw, rw, n1, n2, damp, vdamp, dwin))}, rwp, n1, n2, damping, var_damping, trend_window)
  }
  
  if (!object$model$cohort=="")
  {
    rwc<-switch(object$model$cohort,
                rw1 = 1,
                rw2 = 2
    )
    c1<-dim(object$samples$cohort[[1]])[2]
    c2<-bamp::coh(1,n2,a1,object$data$periods_per_agegroup)
    ch<-length(object$samples$cohort)
      prep<-parallel::mclapply(1:ch, function(i,samples)cbind(object$samples$cohort_parameter[[i]],object$samples$cohort[[i]]), samples)
      psi<-parallel::mclapply(prep, function(prepi, rw, n1, n2, damp, vdamp){t(apply(prepi, 1, predict_rw, rw, n1, n2, damp, vdamp))}, rwc, c1, c2, damping, var_damping)
  }
    
  nr.samples<-length(object$samples$intercept[[1]])
    theta<-if(object$model$age!=""){object$samples$age}else{NA}
    
      delta<-if(object$model$overdispersion){object$samples$overdispersion}else{NA}

      prepfx<- function(i, theta, phi, psi, my, delta, a1, n1, c1, nr){
      theta=if(any(is.na(theta))){array(0, c(nr, a1))}else{theta[[i]]}
      phi=if(any(is.na(phi))){array(0, c(nr, a1))}else{phi[[i]]}
      psi=if(any(is.na(psi))){array(0, c(nr, a1))}else{psi[[i]]}
      my=my[[i]]
      delta=if(any(is.na(delta))){rep(0, nr)}else{delta[[i]]}
      return(cbind(my,theta,phi,psi, delta))
    }
    
    prep<- parallel::mclapply(1:ch, prepfx, theta, phi, psi, object$samples$intercept, delta, a1, n1, c1, nr.samples)

  ksi<-parallel::mclapply(prep, function(prepi, vdb, noa, nop, nop2, noc, zmode){
    temp<-apply(prepi, 1, ksi_prognose, vdb, noa, nop, nop2, noc, zmode); return(array(temp,c(nop2,noa,dim(temp)[2])))}, object$data$periods_per_agegroup,
    a1, n1, n2, c2, object$model$overdispersion)

    
  ksi0<-ksi[[1]]
  if (ch>1)
  for (i in 2:ch)
  {
    ksi0<-abind::abind(ksi0,ksi[[i]], along=3)
  }
  
  pr <- exp(ksi0)/(1+exp(ksi0))

  # Optional hazard scale. Convert the per-period event probability p into the
  # cumulative cause-specific hazard h = -log(1-p) over one period, divided by
  # period_length to express it per time unit (e.g. per year). Cause-specific
  # hazards are ADDITIVE across competing causes (H_allcause = sum_c H_c),
  # unlike risks, so these are the correct quantity to feed to a competing-risk
  # life-table / microsimulation. See the Details section of ?predict_apc.
  hz <- if (hazard) -log1p(-pr)/period_length else NULL

  if (is.null(population))population<-object$data$population
  
  n0<-min(dim(population)[1],n2)
  
  predictedcases <- array(apply(pr,3,function(pr1,n,n0){
    return(rbinom(n0*dim(n)[2],n[1:n0,],pr1[1:n0,]))},population,n0),c(n0,dim(pr)[2:3]))
  
  predictperiod<-apply(predictedcases,c(1,3),sum)

  qu_predicted_cases<-apply(predictedcases,1:2,quantile,quantiles)
  qu_predicted_pr<-apply(pr,1:2,quantile,quantiles)
  qu_predictperiod<-apply(predictperiod,1,quantile,quantiles)
  
  
  period<-phi[[1]]
  if (ch>1)
  for (i in 2:ch)
    period<-abind::abind(period,phi[[i]], along=1)
  cohort<-psi[[1]]
  if (ch>1)
  for (i in 2:ch)
    cohort<-abind::abind(cohort,psi[[i]], along=1)

  samples<-list(
    "pr"=pr,
    "cases"=predictedcases,
    "cases_period"=predictperiod,
    "period"=period,
    "cohort"=cohort
  )
  if (hazard) samples[["hazard"]] <- hz

  predicted<-list(
    "pr"=qu_predicted_pr,
    "cases"=qu_predicted_cases,
    "cases_period"=qu_predictperiod,
    "period"=apply(period,2,quantile,quantiles),
    "cohort"=apply(cohort,2,quantile,quantiles),
    "samples"=samples
  )
  if (hazard) predicted[["hazard"]] <- apply(hz, 1:2, quantile, quantiles)

  if (!update){
    return(predicted)}
  else{
    object$predicted=predicted
    return(object)
  }
}
