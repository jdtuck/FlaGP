#' @title FlaGP MAP calibration
#'
#' @description MAP estimate of calibration parameters
#' @param flagp x an \code{flagp} object.
#' @param n.restarts size of random sample used to initialize the optimizer
#' @param opts list of options for crs::snomadr
#' @details Returns FlaGP data object
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
map = function(flagp,n.restarts=1,init=NULL,seed=1,
                opts=list("MAX_BB_EVAL" = 10000, "INITIAL_MESH_SIZE" = .1,
                          "MIN_POLL_SIZE" = "r0.001", "DISPLAY_DEGREE" = 0),
                end.eta=50,lagp.delta=F,start.delta=6,end.delta=50,
                theta.prior='beta',theta.prior.params=c(2,2),
                ssq.prior.params=c(1,1e-3)){

  p.t = flagp$XT.data$p.t
  if(!is.null(init)){
    # if init passed, check for correct dimensions
    if(ncol(init) != p.t){
      stop("init should be a matrix with p.t columns")
    }
    n.restarts = nrow(init)
  } else{
    # if init not passed generate with LHS
    if(n.restarts>1){
      set.seed(seed)
      init = rbind(init,lhs::maximinLHS(max(2,n.restarts),p.t))
    } else{
      init = matrix(rep(.5,p.t),nrow=1)
    }
  }

  out = NULL
  time = proc.time()[3]
  # try doParallel::foreach or parLapply
  if(n.restarts>1){
    cores = min(n.restarts,parallel::detectCores()-4)
    cl <- parallel::makeCluster(cores)
    outi = parallel::mclapply(1:n.restarts, function(i) crs::snomadr(fit_model, n = p.t, bbin = rep(0, p.t), bbout = 0,
                                                                     x0 = init[i,],
                                                                     lb = rep(0, p.t),
                                                                     ub = rep(1, p.t),
                                                                     random.seed = seed,
                                                                     opts = opts,
                                                                     print.output=F,
                                                                     flagp = flagp,
                                                                     lite=T,negll=T,
                                                                     sample=F,
                                                                     end.eta=end.eta,
                                                                     lagp.delta=lagp.delta,
                                                                     start.delta=start.delta,
                                                                     end.delta=end.delta,
                                                                     theta.prior=theta.prior,
                                                                     theta.prior.params=theta.prior.params,
                                                                     ssq.prior.params=ssq.prior.params,
                                                                     map=T),mc.cores=min(parallel::detectCores()-1,n.restarts))
    obj = Inf
    parallel::stopCluster(cl)
    for(i in 1:n.restarts){
      if(outi[[i]]$objective<obj){
        out = outi[[i]]
        obj = out$objective
      }
    }
  } else{
    out = crs::snomadr(fit_model, n = p.t, bbin = rep(0, p.t), bbout = 0,
                       x0 = init,
                       lb = rep(0, p.t),
                       ub = rep(1, p.t),
                       random.seed = seed,
                       opts = opts,
                       print.output=F,
                       flagp = flagp,
                       lite=T,negll=T,
                       sample=F,
                       end.eta=end.eta,
                       lagp.delta=lagp.delta,
                       start.delta=start.delta,
                       end.delta=end.delta,
                       theta.prior=theta.prior,
                       theta.prior.params=theta.prior.params,
                       ssq.prior.params=ssq.prior.params,
                       map=T)
  }

  opt.time = proc.time()[3] - time
  return = fit_model(out$solution,flagp,lite=F,sample=F,end.eta=end.eta,
                     lagp.delta=lagp.delta,start.delta=start.delta,end.delta=end.delta,
                     theta.prior=theta.prior,theta.prior.params=theta.prior.params,ssq.prior.params=ssq.prior.params,map=T)
  return$delta$lagp = lagp.delta
  return$snomadr = out
  return$theta.hat = out$solution; return$theta = NULL
  return$time = opt.time
  return$init = init[1:n.restarts,]
  class(return) = c('map',class(return))
  return(return)
}
