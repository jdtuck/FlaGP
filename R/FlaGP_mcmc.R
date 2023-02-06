# Single site MCMC with uniform proposal steps. Option to adapt proposal steps to target acceptance rate
mcmc_ss = function(flagp,t.init,n.samples=100,n.burn=0,
                   prop.step=rep(.2,flagp$XT.data$p.t),
                   end.eta=50,
                   lagp.delta=F,start.delta=6,end.delta=50,
                   adapt=T,target.accept=.3,
                   theta.prior='beta',ssq.prior='gamma',verbose=T)
{
  ptm = proc.time()[3]
  bias = flagp$bias
  p.t = flagp$XT.data$p.t
  n.pc = flagp$basis$sim$n.pc

  t.curr = t.init
  #llt = fit_model(t.curr,flagp,F,T,end.eta,lagp.delta,start.delta,end.delta,F,theta.prior,ssq.prior)
  #ssq = llt$ssq.hat
  #w = llt$eta$w$sample
  #v = llt$delta$v$sample
  ll.prop = list()
  accept = matrix(0,nrow=n.samples,ncol=p.t)#; accept[1,] = 1

  # Things to store
  t.store = matrix(nrow=n.samples,ncol=flagp$XT.data$p.t)#;t.store[1,] = t.curr
  ssq.store = numeric(n.samples)#;ssq.store[1] = ssq
  #w.store = array(dim=c(n.samples,dim(w)));w.store[1,,] = w
  #v.store = array(dim=c(n.samples,dim(v)))
  y.resid.store = array(dim=c(n.samples,length(flagp$Y.data$obs$ind),flagp$Y.data$n))

  GPs.store = vector(mode='list',length=n.samples); GPs.store[[1]] = llt
  ll.store = numeric(n.samples); ll.store[1] = llt$ll

  #if(bias){v.store[1,,] = v}
  cat('MCMC Start #--',format(Sys.time(), "%a %b %d %X"),'--#\n')
  for(i in 1:n.samples){
    if(verbose){
      if(i %% 1000 == 0){
        cat('MCMC iteration',i,'#--',format(Sys.time(), "%a %b %d %X"),'--#\n')
      }
    }
    for(j in 1:p.t){
      t.prop = t.curr # initialize t.prop to t so that the parameter we are not updating is fixed at its previous value
      t.prop[j] = t.curr[j] + (prop.step[j] * runif(1,-0.5, 0.5))
      if(t.prop[j]<0 | t.prop[j]>1){
        ll.prop$ll = -Inf
      } else{
        ll.prop = fit_model(t.prop,flagp,F,T,end.eta,lagp.delta,start.delta,end.delta,F,theta.prior,ssq.prior)
      }

      # we have to recompute this to account for the variability in the likelihood at t.curr,
      # otherwise, the chain does not mix well.
      llt = fit_model(t.curr,flagp,F,T,end.eta,lagp.delta,start.delta,end.delta,F,theta.prior,ssq.prior)

      # 4. Accept or reject t
      if(runif(1) < exp(ll.prop$ll - llt$ll)){
        t.curr[j] = t.prop[j]
        ssq = ll.prop$ssq.hat
        w = ll.prop$eta$w$sample
        v = ll.prop$delta$v$sample
        accept[i,j] = accept[i,j] + 1
      }
      t.store[i,j] = t.curr[j]
    }
    # store w,v,ssq after both updates have been made
    ssq.store[i] = ssq
    #w.store[i,,] = w
    GPs.store[[i]] = llt
    ll.store[i] = llt$ll

    if(adapt & i<=n.burn){
      # adaptive proposal step
      if( i %% 100 == 0 ){
        acpt.ratio = colSums(accept[(i-99):i,,drop=F])/100
        for(jj in 1:p.t){
          if(acpt.ratio[jj]<(target.accept-.02) & prop.step[jj]>.05){
            prop.step[jj]=.9*prop.step[jj]
          } else if(acpt.ratio[jj]>(target.accept+.02)){
            prop.step[jj]=1.1*prop.step[jj]
          }
        }
      }
    }
  }
  mcmc.time = proc.time()[3] - ptm

  returns = list(t.samp=t.store[(n.burn+1):n.samples,,drop=F],
                 ssq.samp=ssq.store[(n.burn+1):n.samples],
                 ll.samp=ll.store[(n.burn+1):n.samples],
                 GPs = GPs.store[(n.burn+1):n.samples],
                 acpt.ratio=colSums(accept[(n.burn+1):n.samples,,drop=F])/(n.samples-n.burn),
                 time=mcmc.time,
                 prop.step=prop.step,
                 n.samples = n.samples,
                 n.burn = n.burn,
                 lagp.delta = lagp.delta)
  class(returns) <- c("mcmc", class(returns))
  return(returns)
}

# Estimate optimal step size for mcmc proposals to target a specific acceptance rate
tune_step_sizes = function(flagp,n.samples,n.levels,target.accept.rate=NULL,min.step=.01,max.step=1,
                           end.eta=50,lagp.delta=F,start.delta=6,end.delta=50)
{
  logit = function(x){
    log(x/(1-x))
  }

  p.t = flagp$XT.data$p.t

  step.sizes = matrix(logseq(min.step,max.step,n=n.levels),nrow=n.levels,ncol=p.t,byrow = F)

  acc = matrix(nrow=n.levels,ncol=p.t)
  for(i in 1:n.levels){
    acc[i,] = do_mcmc(flagp,rep(.5,p.t),n.samples,0,step.sizes[i,],F,end.eta,lagp.delta,start.delta,end.delta)$acpt.ratio
  }

  # Compute GLM for each parameter
  step.tuned = numeric(p.t)
  if(is.null(target.accept.rate)){
    target.logit = rep(log(1 / (exp(1) - 1)),p.t)
  } else{
    target.logit = rep(logit(target.accept.rate),p.t)
  }
  for(j in 1:p.t){
    y = cbind(acc[,j]*n.samples,n.samples-(acc[,j]*n.samples))
    #x = cbind(log(step.sizes[,j]))
    x = log(step.sizes[,j])
    glm.model = glm(y~x,family=binomial(link='logit'))
    coefs = as.numeric(coef(glm.model))
    if(target.logit[j]<min(x) | target.logit[j]>max(x)){
      cat('extrapolating')
    }
    #step.tuned[j] = exp((target.logit[j]-coefs[1])/coefs[2])
    step.tuned[j] = predict(glm.model,data.frame(x=target.logit[j]),type='response')
  }
  return(step.tuned)
}

#' @title FlaGP MCMC with joint proposal
#'
#' @description Addaptive MCMC as defined in Haario et al. 2001 - "An adaptive Metropolis algorithm"
#' @param flagp an \code{flagp} object.
#' @param t.init vector an inital value of the calibration parameters
#' @param cov.init matrix initial covariance matrix for the proposal distribution
#' @param n.samples integer number of mcmc samples to keep
#' @param n.burn integer number of initial samples to discard
#' @param update.every integer sample frequency that proposal covariance should be updated
#' @param stop.update integer sample index where proposal covariance updates should stop
#' @param end.eta integer number of neighbors to use for laGP emulator
#' @param lagp.delta T/F laGP bias model. If F, a full GP is used
#' @param start.delta initial neighborhood size for bias model. Only applicable if lagp.delta=T
#' @param end.delta final neighborhood size for bias model. Only applicable if lagp.delta=T
#' @param theta.prior prior for calibration parameters, 'beta' or 'unif'
#' @param ssq.prior prior for error variance, 'gamma' is a non-informative gamma on the precision
#' @param prev.samples Use in accordance with add_samples(). A matrix of samples can be passed in to be used for proposal covariance estimation.
#' @param verbose print status updates
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
mcmc = function(flagp,t.init=NULL,cov.init=NULL,
                      n.samples=10000,n.burn=1000,
                      update.every=1,stop.update=1000,
                      end.eta=50,
                      lagp.delta=F,start.delta=6,end.delta=50,
                      theta.prior='beta',ssq.prior='gamma',prev.samples=NULL,
                      verbose=T)
{
  ptm = proc.time()[3]

  # constants and initializations
  bias = flagp$bias
  p.t = flagp$XT.data$p.t
  if(is.null(t.init)){
    t.init = rep(.5,p.t)
  }
  if(is.null(cov.init)){
    # diagonal with sd s.t. 99% mass between [0,1]
    cov.init = diag((.5/3)^2,p.t)
  }
  t.curr = t.init
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
  ll.store = numeric(n.samples)
  accept = numeric(n.samples)

  # initialize adaptive parameters
  s.d = 2.4^2/p.t
  eps = diag(sqrt(.Machine$double.eps),p.t)
  C = cov.init
  # If we aren't updating C, let C.store = cov.init which will be returned
  # if(update.every>0){
  #   C.store = array(dim = c(n.samples/update.every, p.t, p.t))
  # } else{
  #   C.store = NULL
  # }

  cat('MCMC Start #--',format(Sys.time(), "%a %b %d %X"),'--#\n')
  for(i in 1:n.samples){
    if(verbose){
      if(i %% 1000 == 0){
        cat('MCMC iteration',i,'#--',format(Sys.time(), "%a %b %d %X"),'--#\n')
      }
    }

    # propose new t
    t.prop = mvtnorm::rmvnorm(1, t.curr, s.d * C + eps)
    if(any(t.prop<0) | any(t.prop>1)){
      ll.prop$ll = -Inf
    } else{
      ll.prop = fit_model(t.prop,flagp,F,T,end.eta,lagp.delta,start.delta,end.delta,F,theta.prior,ssq.prior)
    }

    # we have to recompute this every time to account for the variability in the likelihood at t.curr,
    # otherwise, the chain does not mix well.
    llt = fit_model(t.curr,flagp,F,T,end.eta,lagp.delta,start.delta,end.delta,F,theta.prior,ssq.prior)

    # 4. Accept or reject t
    accept.prob =  exp(ll.prop$ll - llt$ll)
    if(runif(1) < accept.prob){
      t.curr = t.prop
      ssq.store[i] = ll.prop$ssq.hat
      eta.store[[i]] = ll.prop$eta
      if(bias){delta.store[[i]] = ll.prop$delta}
      ll.store[i] = ll.prop$ll
      accept[i] = accept[i] + 1
    } else{
      ssq.store[i] = llt$ssq.hat
      eta.store[[i]] = llt$eta
      if(bias){delta.store[[i]] = llt$delta}
      ll.store[i] = llt$ll
    }
    t.store[i,] = t.curr
    if(i==2)
      t.mean = apply(t.store[1:2,,drop=F],2,mean)
    # adapt proposal distribution
    if(i>2 & update.every>0 & i<= stop.update){# & i %% update.every == 0){
      # if previous samples were passed to the function, use them for estimation
      # C = cov(rbind(prev.samples,t.store[1:i,])) + diag(eps,p.t)
      # C.store[as.integer(i/update.every),,] = C

      # single sample updates

      # center new sample with old mean
      tmp = t.store[i,,drop=F] - t.mean
      # use covariance recursion formula
      C = ((i-2)/(i-1))*C + (1/i)*t(tmp)%*%tmp
      # C.store[i,,] = C
      # update mean to include new sample
      t.mean = (1/i)*((i-1)*t.mean + t.curr)
    }
  }
  mcmc.time = proc.time()[3] - ptm

  return.ids = (n.burn+1):n.samples
  returns = list(t.samp     = t.store[return.ids,,drop=F],
                 ssq.samp   = ssq.store[return.ids],
                 ll.samp    = ll.store[return.ids],
                 eta        = eta.store[return.ids],
                 delta      = delta.store[return.ids],
                 prop.cov   = C,
                 acpt.ratio = mean(accept[return.ids]),
                 time       = mcmc.time,
                 n.samples  = n.samples,
                 n.burn     = n.burn)
  class(returns) <- c("mcmc", class(returns))
  return(returns)
}

# Function to add samples to existing mcmc object, starting from the last interation.
add_samples = function(flagp,mcmc,type='joint',n.samples=100,
                       end.eta=50,
                       lagp.delta=F,start.delta=6,end.delta=50,
                       update.every=1,stop.update=10000,
                       adapt=T,target.accept=.3,
                       theta.prior='beta',ssq.prior='gamma',verbose=T){

  if(type=='joint'){
    add.samples = mcmc_joint(flagp,mcmc$t.samp[nrow(mcmc$t.samp),],mcmc$prop.cov[dim(mcmc$prop.cov)[1],,],
                             n.samples,0,update.every,stop.update,
                             end.eta,lagp.delta,start.delta,end.delta,
                             theta.prior,ssq.prior,mcmc$t.samp,verbose)
    mcmc$prop.cov = abind::abind(mcmc$prop.cov,add.samples$prop.cov,along=1)
  } else{
    add.samples = mcmc_ss(flagp,mcmc$t.samp[nrow(mcmc$t.samp),],n.samples,0,mcmc$prop.step,
                          end.eta,lagp.delta,start.delta,end.delta,adapt=adapt,
                          target.accept=target.accept,verbose)
    mcmc$prop.step = add.samples$prop.step
  }


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
