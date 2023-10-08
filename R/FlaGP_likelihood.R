
# w,v are needed if sample=F
compute_ll = function(theta,eta,delta,sample,flagp,
                      theta.prior='beta',theta.prior.params=c(2,2),
                      ssq.prior.params=c(1,1e-3)){
  if(!is.null(delta)){
    ssq.hat = delta$ssq.hat
    y.resid = delta$y.resid
  } else{
    ssq.hat = eta$ssq.hat
    y.resid = eta$y.resid
  }
  if(sample){
    # compute simple llh - Sigma_w and Sigma_v are accounted for in calc of y.resid
    ldetS = flagp$Y.data$obs$n.y*log(ssq.hat)
    # Don't compute S inverse directly - bad for huge d_y
    ySinvty = apply(y.resid^2,2,sum)/ssq.hat
    ll = -.5*(ldetS + ySinvty)
  } else{
    lambda = 1/ssq.hat
    B.hat = flagp$precomp$BDtBDinvtBD%*%y.resid
    SigB.hat = lapply(1:flagp$Y.data$n, function(i) as.matrix(diag(c(eta$w$var[,i],delta$v$var[,i]),nrow=length(c(eta$w$var[,i],delta$v$var[,i]))) + ssq.hat*flagp$precomp$BDtBDinv))
    l_beta_hat = sapply(1:flagp$Y.data$n, function(i) mvtnorm::dmvnorm(x=as.numeric(B.hat[,i]),sigma = SigB.hat[[i]], log = T))
    ll = sapply(1:flagp$Y.data$n, function(i) -.5*((flagp$Y.data$obs$n.y-flagp$precomp$rankBD)*log(2*pi*ssq.hat) + flagp$precomp$ldetBDtBD +
                                        lambda*t(y.resid[,i])%*%flagp$precomp$LLHmat%*%y.resid[,i]) + l_beta_hat[i])
  }
  if(theta.prior=='beta'){
    # beta prior - default
    p.theta = sum(dbeta(theta,theta.prior.params[1],theta.prior.params[2],log=T))
  } else if(theta.prior=='unif'){
    p.theta = sum(dunif(theta,0,1,log=T))
  }
  return(sum(ll) + p.theta)
}
