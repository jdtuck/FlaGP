#' @title FlaGP MCMC with joint proposal
#'
#' @description Addaptive MCMC as defined in Haario et al. 2001 - "An adaptive Metropolis algorithm"
#' @param flagp an \code{flagp} object.
#' @param t.init vector an inital value of the calibration parameters
#' @param prop.sd.init vector of initial proposal distribution sd's
#' @param n.samples integer number of mcmc samples to keep
#' @param n.burn integer number of initial samples to discard
#' @param adapt.par vector of parameters for addaptive metropolis. The first parameter it the first iteration to begin adaptation, the second is the frequency of adaptation, the third is the number of past samples to include for adaptation, and the fourth is the iteration to stop adaptation.
#' @param end.eta integer number of neighbors to use for laGP emulator
#' @param delta.method Model used for discrepancy. 'newGP' is a full GP using laGP default emperical bayes. 'rgasp' is a default RobustGaSP implementation and 'lagp' is an laGP predictor using the ALC criterion
#' @param start.delta initial neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param end.delta final neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param theta.prior prior for calibration parameters, 'beta' or 'unif'
#' @param theta.prior.params parameters for beta prior. Only applicable of theta.prior='beta'
#' @param ssq.prior.params gamma prior parameters for error variance facilitating conjugate updating
#' @param prev.samples Use in accordance with add_samples(). A matrix of samples can be passed in to be used for proposal covariance estimation.
#' @param recalc.llh Should the likelihood be recalculated at every iteration. Generally this should be False and only considered if sample=T. Recalculing the likelihood with a new sample from the emulator (and discrepancy) can be useful if mcmc mixing is a problem. This is only likely to happen if a discrepancy model is being used or if the emulator poorly fits the observations and no discrepancy model is used. Importantly, inference is approximate if recalculation is used, as posterior variance will be artificially increased.
#' @param sample Account for uncertainty in the emulator and discrepancy model by sampling rather than by including the variance in the likelihood. Default is sample=T when bias=F which works well because sampling variability in the emulator is usually small. For a biased calibration we recommend sample=F. This will increase the cost of likelihood calculation, but result in a much more stable mcmc.
#' @param verbose print status updates
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
mcmc = function(flagp,
                t.init=rep(.5,flagp$XT.data$p.t),
                ssq.init=.01,
                prop.cov=diag((.5/3)^2,flagp$XT.data$p.t+1),
                n.samples=10000,n.burn=1000,
                adapt.par = c(100,50,50,1000),
                end.eta=50,
                delta.method='newGP',start.delta=6,end.delta=50,
                theta.prior='beta',theta.prior.params=c(2,2),
                ssq.prior='hcauchy',ssq.prior.params=c(.5),
                prev.samples=NULL,
                recalc.llh=F,
                sample=as.logical(ifelse(flagp$bias,F,T)),
                verbose=T)
{
  ptm = proc.time()
  # constants and initializations
  bias = flagp$bias
  p.t = flagp$XT.data$p.t

  t.curr = t.init
  ssq.curr = ssq.init
  # llt = fit_model(t.curr,flagp,F,sample,end.eta,delta.method,start.delta,end.delta,F,
  #                 theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
  llt = fit_model(t.curr,ssq.curr,flagp,F,sample,end.eta,delta.method,start.delta,end.delta,F,
                  theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
  llt$ll = llt$ll + logit.jacobian(t.curr) + log.jacobian(ssq.curr) # jacobian is added here because it may be hard to move from initial location otherwise (jacobian can be a large decrease in llh)
  # if(sample.type=='logit')
  #   llt$ll = llt$ll
  llt.prev = llt
  ll.prop = list()

  # Storage
  t.store = matrix(nrow=n.samples,ncol=flagp$XT.data$p.t)
  ssq.store = numeric(n.samples)
  eta.store = vector(mode='list',length=n.samples)
  if(bias){
    delta.store = vector(mode='list',length=n.samples)
  } else{
    delta.store = NULL
  }
  ll.store = accept = numeric(n.samples)

  # initialize adaptive parameters
  # s.d = 2.4^2/p.t
  eps = diag(sqrt(.Machine$double.eps),p.t+1)
  if(verbose)
    cat('MCMC Start #--',format(Sys.time(), "%a %b %d %X"),'--#\n')
  for(i in 1:n.samples){
    if(verbose){
      if(i %% 1000 == 0){
        cat('MCMC iteration',i,'#--',format(Sys.time(), "%a %b %d %X"),'--#\n')
      }
    }

    # propose new t
    prop = proposal(t.curr,ssq.curr,p.t,prop.cov+eps)
    t.prop = prop$t.prop
    ssq.prop = prop$ssq.prop

    # if(any(t.prop<0) | any(t.prop>1)){
    #   ll.prop$ll = -Inf
    # } else{
    ll.prop = fit_model(t.prop,ssq.prop,flagp,F,sample,end.eta,delta.method,start.delta,end.delta,F,
                        theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
    ll.prop$ll = ll.prop$ll + logit.jacobian(t.prop) + log.jacobian(ssq.prop)
      # ll.prop = fit_model(t.prop,flagp,F,sample,end.eta,delta.method,start.delta,end.delta,F,
      #                     theta.prior,theta.prior.params,ssq,prior,ssq.prior.params)
    # }
    # if(sample.type=='logit')
    #   ll.prop$ll = ll.prop$ll + logit.jacobian(t.prop) + log.jacobian(ssq.prop)

    # we have to recompute this every time to account for the variability in the likelihood at t.curr,
    # otherwise, the chain does not mix well.
    if(recalc.llh & i>1){
      llt$eta = resample_w(llt$eta$w,flagp$Y.data$obs$trans,flagp$basis$obs$B,ssq.prior.params)
      if(flagp$bias){
        llt$delta = resample_v(llt$delta$v,llt$eta$y.resid,delta.method)
      }
      llt$ll = compute_ll(t.curr,ssq.curr,llt$eta,llt$delta,sample,flagp,theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
      # if(sample.type=='logit')
      llt$ll = llt$ll + logit.jacobian(t.curr) + log.jacobian(ssq.curr)
      if(flagp$bias){
        llt$ssq.hat = llt$delta$ssq.hat
        llt$delta$ssq.hat = NULL
      } else{
        llt$ssq.hat = llt$eta$ssq.hat
        llt$eta$ssq.hat = NULL
      }
    }

    # 4. Accept or reject t
    accept.prob =  exp(ll.prop$ll - llt$ll)
    if(runif(1) < accept.prob){
      t.curr = t.prop
      ssq.curr = ssq.prop
      if(!recalc.llh){
        # accepted, but not recalculating, so update old llt to ll.prop
        llt = ll.prop
      } else{
        # accepted and recalculating, update llt.prev
        llt.prev = ll.prop
      }
      # ssq.store[i] = ll.prop$ssq.hat
      ssq.store[i] = ssq.prop
      eta.store[[i]] = ll.prop$eta
      if(bias){delta.store[[i]] = ll.prop$delta}
      ll.store[i] = ll.prop$ll
      accept[i] = 1
    } else{
      if(recalc.llh){
        llt = llt.prev
      }
      # ssq.store[i] = llt$ssq.hat
      eta.store[[i]] = llt$eta
      if(bias){
        delta.store[[i]] = llt$delta
      }
      ll.store[i] = llt$ll
    }
    t.store[i,] = t.curr
    ssq.store[i] = ssq.curr

    # adapt proposal distribution
    if(i>=adapt.par[1] & i<=adapt.par[4] & i %% adapt.par[2] == 0){
      interval = max(i-adapt.par[2]+1,1):i
      t.interval = t.store[interval,,drop=F]
      ssq.interval = ssq.store[interval]
      lambda = ifelse(i==adapt.par[1],1,update$lambda)
      # if(sample.type=='logit'){
        # update = update_tuning_mv(i, mean(accept[interval]), lambda=lambda, log(t.interval/(1-t.interval)),
        #                           prop.cov/lambda)
      update = update_tuning_mv(i, mean(accept[interval]), lambda=lambda, cbind(log(t.interval/(1-t.interval)),ssq.interval),
                                prop.cov/lambda)
      # } else{
      #   update = update_tuning_mv(i, mean(accept[interval]), lambda=lambda, t.interval,
      #                             prop.cov/lambda)
      # }
      prop.cov = update$lambda * update$Sigma_tune
    }
  }
  mcmc.time = proc.time() - ptm

  return.ids = (n.burn+1):n.samples
  returns = list(t.samp     = t.store[return.ids,,drop=F],
                 ssq.samp   = ssq.store[return.ids],
                 ll.samp    = ll.store[return.ids],
                 eta        = eta.store[return.ids],
                 delta      = delta.store[return.ids],
                 prop.cov   = prop.cov,
                 acpt.ratio = mean(accept[return.ids]),
                 time       = mcmc.time,
                 n.samples  = n.samples,
                 n.burn     = n.burn)
  class(returns) <- c("mcmc", class(returns))
  return(returns)
}
logit.jacobian = function(theta){
  sum(log(theta)+log(1-theta))
}
log.jacobian = function(ssq){
  log(ssq)
}
# proposal = function(t.curr,p.t,sample.type,prop.type,prop.cov){
#   if(prop.type == 'mv'){
#     # multivariate proposal with covariances
#     if(sample.type=='logit'){
#       t.prop = mvtnorm::rmvnorm(1, log(t.curr/(1-t.curr)), prop.cov)
#       t.prop = exp(t.prop)/(1+exp(t.prop))
#     } else{
#       t.prop = mvtnorm::rmvnorm(1, t.curr, (prop.cov))
#     }
#   } else{
#     # diagonal proposal with variance only - sampling may be less efficient than if correlation between the parameters is used for proposal
#     if(sample.type=='logit'){
#       t.prop = rnorm(p.t,log(t.curr/(1-t.curr)), sqrt(diag(prop.cov)))
#     } else{
#       t.prop = rnorm(p.t,t.curr, sqrt(diag((prop.cov))))
#     }
#   }
#   return(t.prop)
# }
# propose both theta and sigma^2
proposal = function(t.curr,ssq.curr,p.t,prop.cov){
  # propose theta on logit scale and error variance on log scale
  par.prop = mvnfast::rmvn(1,c(log(t.curr/(1-t.curr)),log(ssq.curr)), prop.cov)
  return(list(t.prop=exp(par.prop[1:p.t])/(1+exp(par.prop[1:p.t])),
              ssq.prop=exp(par.prop[p.t+1])))
}
# Function to add samples to existing mcmc object, starting from the last interation.
# mcmc() has been updated, this function needs to be as well to reflect that
add_samples = function(flagp,mcmc,n.samples=100,
                       end.eta=50,sample=F,recalc.llh=F,
                       delta.method='newGP',start.delta=6,end.delta=50,
                       adapt.par = c(100,50,50,1000),
                       theta.prior='beta',theta.prior.params=c(2,2),
                       ssq.prior='hcauchy',ssq.prior.params=c(.5),
                       verbose=T){

  # flagp,
  # t.init=rep(.5,flagp$XT.data$p.t),
  # ssq.init=.01,
  # prop.cov=diag((.5/3)^2,flagp$XT.data$p.t+1),
  # n.samples=10000,n.burn=1000,
  # adapt.par = c(100,50,50,1000),
  # end.eta=50,
  # delta.method='newGP',start.delta=6,end.delta=50,
  # theta.prior='beta',theta.prior.params=c(2,2),
  # ssq.prior='hcauchy',ssq.prior.params=c(.5),
  # prev.samples=NULL,
  # recalc.llh=F,#sample.type='logit',prop.type='mv',
  # sample=as.logical(ifelse(flagp$bias,F,T)),
  # verbose=T)

  add.samples = mcmc(flagp,
                     t.init=mcmc$t.samp[nrow(mcmc$t.samp),],
                     ssq.init=mcmc$ssq.samp[length(mcmc$ssq.samp)],
                     prop.cov=mcmc$prop.cov,
                     n.samples=n.samples,n.burn=0,adapt.par=adapt.par,
                     end.eta=end.eta,delta.method=delta.method,start.delta=start.delta,end.delta=end.delta,
                     theta.prior=theta.prior,
                     ssq.prior=ssq.prior,ssq.prior.params=ssq.prior.params,
                     prev.samples=mcmc$t.samp,recalc.llh=recalc.llh,
                     sample=sample,verbose=verbose)

  mcmc$t.samp = rbind(mcmc$t.samp,add.samples$t.samp)
  mcmc$ssq.samp = c(mcmc$ssq.samp,add.samples$ssq.samp)
  mcmc$ll.samp = c(mcmc$ll.samp,add.samples$ll.samp)
  mcmc$eta = append(mcmc$eta,add.samples$eta)
  mcmc$delta = append(mcmc$delta,add.samples$delta)

  # acceptance ratio is weighted avg of number of samples
  prev.samples = mcmc$n.samples-mcmc$n.burn
  mcmc$acpt.ratio = (prev.samples*mcmc$acpt.ratio + n.samples*add.samples$acpt.ratio)/(prev.samples+n.samples)
  mcmc$time = mcmc$time + add.samples$time
  mcmc$n.samples = mcmc$n.samples + n.samples
  return(mcmc)
}

#' @title Parallel wrapper for mcmc function
#'
#' @description Addaptive MCMC as defined in Haario et al. 2001 - "An adaptive Metropolis algorithm"
#' @param flagp an \code{flagp} object.
#' @param t.init vector an inital value of the calibration parameters
#' @param prop.sd.init matrix initial covariance matrix for the proposal distribution
#' @param n.samples integer number of mcmc samples to keep
#' @param n.burn integer number of initial samples to discard
#' @param update.every integer sample frequency that proposal covariance should be updated, only working for update.every=1 currently
#' @param stop.update integer sample index where proposal covariance updates should stop
#' @param end.eta integer number of neighbors to use for laGP emulator
#' @param delta.method Model used for discrepancy. 'newGP' is a full GP using laGP default emperical bayes. 'rgasp' is a default RobustGaSP implementation and 'lagp' is an laGP predictor using the ALC criterion
#' @param start.delta initial neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param end.delta final neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param theta.prior prior for calibration parameters, 'beta' or 'unif'
#' @param ssq.prior.params gamma prior parameters for error variance facilitating conjugate updating
#' @param prev.samples Use in accordance with add_samples(). A matrix of samples can be passed in to be used for proposal covariance estimation.
#' @param verbose print status updates
#' @param n.chains number of parallel chains to run
#' @param n.cores number of cores for parallel processing
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
parallel_mcmc = function(flagp,t.init=rep(.5,flagp$XT.data$p.t),
                         prop.cov=diag((.5/3)^2,flagp$XT.data$p.t),
                          n.samples=10000,n.burn=1000,
                          update.every=50,stop.update=1000,
                          end.eta=50,
                          delta.method='newGP',start.delta=6,end.delta=50,
                          theta.prior='beta',theta.prior.params=c(2,2),
                          ssq.prior.params=c(1,1e-3),prev.samples=NULL,
                          recalc.llh=F,
                          verbose=T,n.chains=2,n.cores=2)
{
  ptm = proc.time()

  cl = parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)

  # constants and initializations
  bias = flagp$bias
  p.t = flagp$XT.data$p.t
  if(is.null(t.init) | nrow(t.init)!=n_chains){
    t.init = tgp::lhs(n.chains,matrix(rep(c(0,1),n.chains),nrow=n.chains,byrow=T))
  }

  # do parallel mcmc chains
  results = foreach::foreach(i=1:n.chains) %dopar% mcmc(flagp, t.init[i,], prop.cov, n.samples, n.burn, update.every, stop.update,
                                                        end.eta, delta.method, start.delta, end.delta,
                                                        theta.prior, theta.prior.params,
                                                        ssq.prior.params, prev.samples, recalc.llh, verbose)

  parallel::stopCluster(cl)
  # combine results
  results$t.samp = c()
  results$ssq.samp = c()
  for(i in 1:length(results)){
    results$t.samp = rbind(results$t.samp,results[[i]]$t.samp)
    results$ssq.samp = c(results$ssq.samp,results[[i]]$ssq.samp)
  }
  class(results) <- c("mcmc", class(results))
  return(results)
}

#' this function updates a block Gaussian random walk proposal Metropolis-Hastings tuning parameters
#' @param k is a positive integer that is the current MCMC iteration.
#' @param accept is a positive number between 0 and 1 that represents the batch accpetance rate. The target of this adaptive tuning algorithm is an acceptance rate of between 0.44 (a univariate update) and 0.234 (5+ dimensional upate).
#' @param lambda is a positive scalar that scales the covariance matrix Sigma_tune and diminishes towards 0 as k increases.
#' @param batch_samples is a \eqn{batch.size \times d}{batch.size x d} dimensional matrix that consists of the batch samples of the d-dimensional parameter being sampled.
#' @param Sigma_tune is a \eqn{d \times d}{d x d} positive definite covariance matrix of the batch samples used to generate the multivariate Gaussian random walk proposal.
#'
#' @keywords internal
update_tuning_mv <- function(k, accept, lambda, batch_samples,
                             Sigma_tune) {
  # target acceptance rate by dimension of parameter vector, min of .234
  arr <- c(0.44, 0.35, 0.32, 0.25, 0.234)
  d = ncol(batch_samples)
  dimension <- min(d,5)
  batch_size <- nrow(batch_samples)
  optimal_accept <- arr[dimension]
  times_adapted <- floor(k / batch_size)
  gamma1 <- 1.0 / ((times_adapted + 3.0)^0.8)
  gamma2 <- 10.0 * gamma1
  adapt_factor <- exp(gamma2 * (accept - optimal_accept))
  ## update the MV scaling parameter
  lambda_out <- lambda * adapt_factor
  ## center the batch of MCMC samples
  batch_samples_tmp <- batch_samples
  for (j in 1:d) {
    mean_batch = mean(batch_samples[, j])
    for (i in 1:batch_size) {
      batch_samples_tmp[i, j] <- batch_samples[i, j] - mean_batch
    }
  }
  Sigma_tune_out <- Sigma_tune + gamma1 *
    (t(batch_samples_tmp) %*% batch_samples_tmp / (batch_size-1.0) - Sigma_tune)
  accept_out <- 0.0
  batch_samples_out <- matrix(0, batch_size, d)
  return(
    list(
      accept          = accept_out,
      lambda          = lambda_out,
      batch_samples   = batch_samples_out,
      Sigma_tune      = Sigma_tune_out
    )
  )
}
