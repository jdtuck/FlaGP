
# Returns predictions at X.pred.orig using the calibration parameters theta. To do this
# the emulator must be 'fit' at theta.
predict_w = function(flagp,X.pred.orig=NULL,theta=NULL,end=50,sample=F,n.samples=1,ann=T)
{
  n.pc = flagp$basis$sim$n.pc
  y=NULL
  if(!is.null(X.pred.orig)){
    n.x.pred = nrow(X.pred.orig)
  } else{
    n.x.pred = 1
  }

  if(!is.null(theta)){
    X = transform_xt(X.sim = flagp$XT.data$sim$X$orig,
                     T.sim = flagp$XT.data$sim$T$orig,
                     X.obs = X.pred.orig,
                     T.obs = matrix(rep(theta,n.x.pred),nrow=n.x.pred,byrow=T))

    X = get_SC_inputs(flagp$lengthscales,X,n.pc)

    XX = lapply(1:n.pc,function(i) cbind(X$X.obs[[i]],X$T.obs[[i]]))
    w = aGPsep_SC_mv(X=X$XT.sim,
                     Z=flagp$basis$sim$V.t,
                     XX=XX,
                     start=start,
                     end=end)
  } else{
    X = transform_xt(X.sim = flagp$XT.data$sim$X$orig,
                     X.obs = X.pred.orig)

    X = get_SC_inputs(flagp$lengthscales,X,n.pc)
    XX = lapply(1:n.pc,function(i) X$X.obs[[i]])
    X = X$XT.sim

    w = aGPsep_SC_mv(X=X,
                     Z=flagp$basis$sim$V.t,
                     XX=XX,
                     start=start,
                     end=end,bias=flagp$bias,
                     ann = ann)
  }

  if(sample){
    w$sample = mvtnorm::rmvnorm(n.samples,as.numeric(w$mean),diag(as.numeric(w$var)))
    dim(w$sample) = c(n.samples,dim(w$mean))
  } else{
    w$sample=w$mean
  }
  return(w)
}

# This function is for predicting an laGP delta model
predict_v = function(flagp,delta,X.pred.orig=NULL,start=6,end=50,sample=F,n.samples=1)
{
  n.pc = ncol(flagp$basis$obs$D)
  y=NULL
  if(!is.null(X.pred.orig)){
    n.x.pred = nrow(X.pred.orig)
  } else{
    n.x.pred = 1
  }

  X = transform_xt(X.sim = flagp$XT.data$sim$X$orig,
                   X.obs = X.pred.orig)
  v = aGPsep_SC_mv(X=X$obs$X$trans,
                   Z=delta$V.t,
                   XX=X$obs$X$trans,
                   start=start,
                   end=end,bias=T)

  if(sample){
    v$sample = mvtnorm::rmvnorm(n.samples,as.numeric(v$mean),diag(as.numeric(v$var)))
    dim(v$sample) = c(n.samples,dim(v$mean))
  } else{
    v$sample=v$mean
  }
  return(v)
}

mv_delta_predict = function(X.pred.orig,delta,flagp,sample=F,n.samples=1, start=6, end=50)
{
  n.pc = ncol(flagp$basis$obs$D)
  v = NULL
  # standardize X.pred.orig
  X.pred.std = unit_xform(X.pred.orig,X.min=flagp$XT.data$sim$X$min,X.range=flagp$XT.data$sim$X$range)$trans

  if(!delta$lagp){
    # build full GPs
    delta.GPs = lapply(1:n.pc, function(k) laGP::newGPsep(X=flagp$XT.data$obs$X$trans,
                                                          Z=delta$V.t[k,],
                                                          d=delta$mle[[k]]$d, g=delta$mle[[k]]$g))
    # predict from GPs
    pred = lapply(1:n.pc, function(k) laGP::predGPsep(delta.GPs[[k]],X.pred.std,lite=T))
    for(i in 1:n.pc){
      v$mean = rbind(v$mean,pred[[i]]$mean)
      v$var = rbind(v$var,pred[[i]]$s2)
    }
  } else{
    # laGP model:
    # X - training inputs are observed X's
    # Z - training outputs are from the delta object which used residuals y - eta and basis D
    # XX - are the prediction locations scaled appropriately
    pred = aGPsep_SC_mv(X=flagp$XT.data$obs$X$trans,
                          Z=delta$V.t,
                          XX=X.pred.std,
                          start=start,
                          end=end,bias=T)
    v$mean = pred$mean
    v$var = pred$var
  }

  if(sample){
    v$sample = mvtnorm::rmvnorm(n.samples,as.numeric(v$mean),diag(as.numeric(v$var)))
    dim(v$sample) = c(n.samples,dim(v$mean)); v$sample = drop(v$sample)
  } else{
    v$sample=v$mean
  }
  return(v)
}

#' @title FlaGP Prediction
#'
#' @description Prediction with \code{mcmc} or \code{map} object
#' @param model an \code{mcmc} or \code{map} object.
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
predict.flagp = function(flagp,model=NULL,X.pred.orig=NULL,n.samples=1,samp.ids=NULL,return.samples=F,support='obs',
                         end.eta=50,lagp.delta=F,start.delta=6,end.delta=50,return.eta=F,return.delta=F, native = T, conf.int=F, ann=T)
{
  if(!is.null(model)){
    if(class(model)[1] == 'mcmc'){
      pred = mcmc_predict(flagp,model,X.pred.orig,samp.ids,n.samples,return.samples,support,end.eta,start.delta,end.delta,return.eta,return.delta)
    } else if(class(model)[1] == 'map'){
      pred = map_predict(flagp,model,X.pred.orig,n.samples,return.samples,support,end.eta,start.delta,end.delta)
    } else{
      stop('model must be of class mcmc or map')
    }
  } else{
    cat('No calibration model, emulation only prediction.')
    if(is.null(X.pred.orig))
      stop('must give X.pred.orig')
    pred = em_only_predict(flagp,X.pred.orig,n.samples,return.samples,support,end.eta,native,conf.int,ann=ann)
  }

  return(pred)
}

map_predict = function(flagp,map,X.pred.orig=NULL,n.samples=1,return.samples=F,support='obs',
                       end.eta=50,start.delta=6,end.delta=50)
{
  start.time = proc.time()[3]
  returns = list()
  if(support=='obs'){
    B = flagp$basis$obs$B
    D = flagp$basis$obs$D
    ym = flagp$Y.data$obs$mean
    ysd = flagp$Y.data$obs$sd
    n.y = nrow(flagp$Y.data$obs$orig)
  } else{
    B = flagp$basis$sim$B
    D = flagp$basis$sim$D
    ym = flagp$Y.data$sim$mean
    ysd = flagp$Y.data$sim$sd
    n.y = nrow(flagp$Y.data$sim$orig)
  }
  # n = 1 if no X model
  n.pred = ifelse(!is.null(X.pred.orig),nrow(X.pred.orig),1)
  # transform_theta
  theta = map$theta.hat * flagp$XT.data$sim$T$range + flagp$XT.data$sim$T$min

  # emulator predictions
  w = predict_w(flagp,X.pred.orig,theta,end=end.eta,sample=T,n.samples)

  if(!flagp$bias){
    # unbiased prediction
    returns$y.samp = array(dim=c(n.samples,n.y,n.pred))
    for(i in 1:n.pred){
      if(n.samples>1){
        returns$y.samp[,,i] = t(t(mvtnorm::rmvnorm(n.samples,mean=B%*%w$mean[,i,drop=F],
                                        sigma = map$ssq.hat*diag(1,n.y) + B%*%diag(w$var[,i])%*%t(B))) *
          ysd + ym)
      } else{
        returns$y.samp = B%*%w$mean[,i,drop=F] * ysd + ym
      }
    }
    returns$time = proc.time()[3] - start.time
    returns$y.mean = apply(returns$y.samp,c(2,3),mean)
    returns$y.conf.int = apply(returns$y.samp,2:3,quantile,c(.025,.975))
    if(!return.samples)
      returns$y.samp = NULL
  } else{
    # biased prediction
    eta.samp = array(dim=c(n.y,n.samples,n.pred))
    delta.samp = array(dim=c(n.y,n.samples,n.pred))
    v = mv_delta_predict(X.pred.orig,map$delta,flagp,F,start=start.delta,end=end.delta)
    for(j in 1:n.pred){
      eta.samp[,,j] = t(mvtnorm::rmvnorm(n.samples,mean=B%*%w$mean[,j,drop=F],
                                       sigma = B%*%diag(w$var[,j])%*%t(B))) *
        ysd + ym
      delta.samp[,,j] = t(mvtnorm::rmvnorm(n.samples,mean=D%*%v$mean[,j,drop=F],
                                         sigma = D%*%diag(v$var[,j],nrow=length(v$var[,j]))%*%t(D) +
                                           map$ssq.hat*diag(n.y) )) * ysd
    }
    returns$y.samp = eta.samp + delta.samp
    returns$time = proc.time()[3] - start.time
    returns$y.mean = apply(returns$y.samp,c(1,3),mean) # do i need to compute this or just add eta.mean to delta.mean?
    returns$eta.mean = apply(eta.samp,c(1,3),mean)
    returns$delta.mean = apply(delta.samp,c(1,3),mean)
    returns$y.conf.int = apply(returns$y.samp,c(1,3),quantile,c(.025,.975))
    if(!return.samples)
      returns$y.samp = NULL
  }
  return(returns)
}

mcmc_predict = function(flagp ,mcmc, X.pred.orig=NULL, samp.ids=NULL, n.samples = 1, return.samples=F, support='obs',
                        end.eta = 50, start.delta = 6, end.delta = 50, return.eta = F, return.delta = F)
{
  returns = list()
  start.time = proc.time()
  if(support=='obs'){
    B = flagp$basis$obs$B
    D = flagp$basis$obs$D
    ym = flagp$Y.data$obs$mean
    ysd = flagp$Y.data$obs$sd
    n.y = flagp$Y.data$obs$n.y
  } else{
    B = flagp$basis$sim$B
    D = flagp$basis$sim$D
    ym = flagp$Y.data$sim$mean
    ysd = flagp$Y.data$sim$sd
    n.y = flagp$Y.data$sim$n.y
  }
  if(is.null(samp.ids)){
    samp.ids = seq(1,(mcmc$n.samples - mcmc$n.burn), length.out = n.samples)
  } else{
    n.samples = length(samp.ids)
  }
  n.pred = ifelse(!is.null(X.pred.orig),nrow(X.pred.orig),1)
  t.pred = t(t(mcmc$t.samp[samp.ids,,drop=F]) * flagp$XT.data$sim$T$range + flagp$XT.data$sim$T$min)
  ssq.samp = mcmc$ssq.samp[samp.ids]
  returns$y.samp = array(0,dim=c(n.samples,n.y,n.pred))
  returns$eta.samp = array(0,dim=c(n.samples,n.y,n.pred))
  if(flagp$bias)
    returns$delta.samp = array(0,dim=c(n.samples,n.y,n.pred))

  for(i in 1:n.samples){
    w = predict_w(flagp,X.pred.orig,t.pred[i,],sample=T,end=end.eta)
    returns$eta.samp[i,,] = B%*%drop(w$sample)
    sigma = ssq.samp[i]*diag(n.y)
    if(flagp$bias){
      # Biased prediction add delta model
      v = FlaGP:::mv_delta_predict(X.pred.orig,mcmc$delta[[i]],flagp,sample=T,n.samples=1,start=start.eta,end=end.eta)
      returns$delta.samp[i,,] = D%*%drop(v$sample)
      for(j in 1:n.pred){
        # sigma noise is on standardized scale, so add noise before scaling back to native with ysd and ym
        returns$y.samp[i,,j] = mvtnorm::rmvnorm(1,mean=returns$eta.samp[i,,j]+returns$delta.samp[i,,j],sigma=sigma) * ysd + ym
      }
    } else{
      # Unbiased prediction eta only
      for(j in 1:n.pred){
        returns$y.samp[i,,j] = mvtnorm::rmvnorm(1,mean=returns$eta.samp[i,,j],sigma=sigma) * ysd + ym
      }
    }
  }

  returns$pred.time = proc.time() - start.time
  returns$y.mean = apply(returns$y.samp,2:3,mean)
  returns$y.conf.int = apply(returns$y.samp,2:3,quantile,c(.025,.975))

  if(return.eta){
    # convert to standard scale first
    for(i in 1:n.samples){
      returns$eta.samp[i,,] = returns$eta.samp[i,,] * ysd + ym
    }
    returns$eta.mean = apply(returns$eta.samp,2:3,mean)
    returns$eta.conf.int = apply(returns$eta.samp,2:3,quantile,c(.025,.975))
  }
  if(flagp$bias & return.delta){
    # convert to standard scale first
    for(i in 1:n.samples){
      returns$delta.samp[i,,] = returns$delta.samp[i,,] * ysd
    }
    returns$delta.mean = apply(returns$delta.samp,2:3,mean)
    returns$delta.conf.int = apply(returns$delta.samp,2:3,quantile,c(.025,.975))
  }
  if(!return.samples)
    returns$y.samp = NULL
    returns$eta.samp = NULL
    if(flagp$bias)
      returns$delta.samp = NULL

  return(returns)
}

em_only_predict = function(flagp, X.pred.orig, n.samples = 1, return.samples=F, support='obs', end.eta = 50, native = T, conf.int = F, ann = NULL)
{
  start.time = proc.time()
  returns = list()
  n.pred = nrow(X.pred.orig)
  # get predictive samples of w at X.pred.orig
  w = predict_w(flagp,X.pred.orig,end = end.eta,sample = T,n.samples = n.samples, ann = ann)
  returns$w = w
  # convert samples of w to samples of y on native scale
  returns$y.samp = array(0,dim=c(n.samples,flagp$Y.data$sim$n.y,n.pred))
  for(i in 1:n.samples){
    returns$y.samp[i,,] = flagp$basis$sim$B %*% drop(w$sample[i,,])
    if(native){
      for(j in 1:n.pred){
        returns$y.samp[i,,j] = returns$y.samp[i,,j] * flagp$Y.data$sim$sd + flagp$Y.data$sim$mean
      }
    }
  }

  returns$pred.time = proc.time() - start.time
  returns$y.mean = apply(returns$y.samp,2:3,mean)
  if(conf.int)
    returns$y.conf.int = apply(returns$y.samp,2:3,quantile,c(.025,.975))
  if(!return.samples)
    returns$y.samp = NULL
  return(returns)
}
