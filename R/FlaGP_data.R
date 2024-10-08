w_to_y = function(w,B,mean=NULL,sd=NULL){
  Y = B%*%w
  if(!is.null(sd)){
    Y = Y * sd
  }
  if(!is.null(mean)){
    Y = Y + mean
  }
  return(Y)
}

# FUNCTION: get_basis
{### Accepts a matrix of response values Y, and returns the basis decomposition from SVD
## Parameters:
# Y (numeric)       : matrix of response values to decompose into basis representation. Y should be
#                     formatted such that each row represents the multivariate response for a single experiment
# n.pc (int)        : number of principal components desired from basis decomposition
# pct.var (numeric) : desired percentage of variability explained by PC's. Overrides n.pc if specified
## Returns: List containings the following
# B (numeric)       : matrix of n.pc basis vectors for Ytrans
# V.t (numeric)      : matrix of n.pc response vectors for Ytrans
# n.pc (integer)    : number of PC's used for the basis decomposition
## Calls:
## Called by:
# user when generating sim basis
# mv_calib.bias for generating discrepancy basis
}
get_basis = function(Y, n.pc = NULL, pct.var = .95, full.basis=F, B=NULL, V.t = NULL, bias=FALSE, D=NULL,
                     rsvd=FALSE, k = NULL, nu = NULL, nv = NULL, p = 10, q = 2, sdist = "normal"){

  return.list = list()

  # Compute SVD
  nExp = ncol(Y) # number of experiments or simulations
  if(is.null(B)){
    ptm = proc.time()[3]
    if(!rsvd){
      svdY = svd(Y)
    } else{
      if(!is.null(n.pc)){
        svdY = rsvd::rsvd(Y,k=n.pc)
      } else{
        stop('n.pc must be specified for rsvd.')
      }
    }
    return.list$svd.time = proc.time()[3] - ptm
    return.list$B = svdY$u %*% diag(svdY$d) / sqrt(nExp)
    return.list$V.t = t(svdY$v) * sqrt(nExp)

    # Return only the required number of basis vectors
    percent.variance = zapsmall((svdY$d^2) / sum(svdY$d^2))
    if(!is.null(pct.var) & is.null(n.pc)){
      # percentage of variance explained specified
      n.pc = which.max(cumsum(percent.variance)>=pct.var)
    }

    return.list$n.pc = n.pc
    if(full.basis){
      return.list$pct.var = percent.variance
      return.list$B = as.matrix(return.list$B)
      return.list$V.t = as.matrix(return.list$V.t)
    } else{
      return.list$pct.var = percent.variance[1:return.list$n.pc]
      return.list$B = as.matrix(return.list$B[,1:return.list$n.pc,drop=F])
      return.list$V.t = as.matrix(return.list$V.t[1:return.list$n.pc,,drop=F])
    }
  } else{
    # manual basis matrix was passed
    return.list$B = B
    return.list$n.pc = ncol(B)
    return.list$V.t = V.t
    if(is.null(V.t))
      return.list$V.t = solve(t(B)%*%B)%*%t(B)%*%Y
  }

  if(bias & is.null(D)){
    stop('must manually specify bias basis matrix')
  }
  return.list$D = D
  return(return.list)
}

# FUNCTION: get_obs_basis
{### Generated obs basis vectors by interpolating sim basis vectors
## Parameters:
# sim.basis       : output from get_basis
# Y.obs        : matrix (n.y x n) field data
# y.ind.sim : matrix of field data indices of observation
## Returns: List containings the following
# B (numeric)       : matrix of n.pc basis vectors
# V.t (numeric)      : matrix of n.pc response vectors
# n.pc (integer)    : number of PC's used for the basis decomposition
## Calls:
## Called by:
# user when generating obs basis
}
get_obs_basis = function(sim.basis,Y.obs,y.ind.sim,y.ind.obs,sigma.y=NULL){
  obs.basis = list()
  n.pc = sim.basis$n.pc
  obs.basis$n.pc = n.pc
  bias = ifelse(is.null(sim.basis$D),F,T)
  if(!isTRUE(all.equal(y.ind.obs,y.ind.sim))){
    B = matrix(nrow=nrow(y.ind.obs),ncol=obs.basis$n.pc)
    if(bias){
      D = matrix(nrow=nrow(y.ind.obs),ncol=ncol(sim.basis$D))
    }
    if(ncol(y.ind.sim)==1){
      # 1d interpolation
      for(j in 1:n.pc){
        B[,j] = pracma::interp1(as.numeric(y.ind.sim),as.numeric(sim.basis$B[,j]),as.numeric(y.ind.obs))
      }
      if(bias){
        for(j in 1:ncol(sim.basis$D)){
          D[,j] = pracma::interp1(as.numeric(y.ind.sim),as.numeric(sim.basis$D[,j]),as.numeric(y.ind.obs))
        }
      }
    } else if(ncol(y.ind.sim)==2){
      # 2d interpolation
      y = unique(y.ind.sim[,1]); x = unique(y.ind.sim[,2])
      yp = y.ind.obs[,1]; xp = y.ind.obs[,2]
      for(i in 1:n.pc){
        B[,i] = pracma::interp2(x,y,matrix(sim.basis$B[,i],length(y),length(x),byrow = F),xp,yp)
      }
      if(bias){
        for(j in 1:ncol(sim.basis$D)){
          D[,j] = pracma::interp2(x,y,matrix(sim.basis$D[,j],length(y),length(x),byrow = F),xp,yp)
        }
      }
    }
    obs.basis$B = B
    if(bias){obs.basis$D=D}
  } else{
    cat('setting Bobs==Bsim \n')
    obs.basis$B = sim.basis$B
    if(bias){
      cat('setting Dobs==Dsim \n')
      obs.basis$D = sim.basis$D
    }
  }
  return(obs.basis)
}

# FUNCTION: mv_lengthscales
# TO DO: this function needs to be updated for only using a subset of the data
### Computes estimated length scale parameters for use in stretching and compressing the inputs space
# Parameters:
# X (numeric) : untransformed matrix of inputs
# V.t (numeric) : matrix
# ls_prior_ab (numeric) : 2 vector of parameters for gamma prior on length-scales
mv_lengthscales = function(XT,p.x,p.t,V.t,g,subsample,m,K,seed,ls.prior,ls.parallel,make.cluster,ls.subsample.size,verbose){
  ptm = proc.time()[3]
  n.pc = nrow(V.t)
  dConfig = laGP::darg(NULL,XT)
  if(!ls.prior){dConfig$ab = c(0,0)}
  if(ls.parallel & make.cluster){
    cl = parallel::makeCluster(min(n.pc,parallel::detectCores()-1))
    doParallel::registerDoParallel(cl)
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(subsample=='blhs'){
    # check size of m
    m.max = ncol(V.t)^(1/(ncol(XT)-1))
    if(m>m.max)
      stop(paste0('Must have ls.m < ',m.max,'. See documentation for laGP::blhs()'))

    if(ls.parallel & n.pc>1){
      tryCatch({
        estLenscales = foreach::foreach(i=1:n.pc) %dopar% blhs.loop.lite(y=V.t[i,],X=XT,m=m,K=K,da=dConfig,g=g)
      },error=function(e){
        if(make.cluster)
          parallel::stopCluster(cl)
        stop("\nError in parallel processing. This could be due to a memory issue.\nTry ls.parallel=F. This will be slower because the estimation will be done in serial. If that does not work, try ls.subsample='strat' or 'rand' for stratified or simple random sampling with ls.subsample.size<m.\n")
      })
    } else{
      estLenscales = lapply(1:n.pc, function(i) blhs.loop.lite(y=V.t[i,],X=XT,m=m,K=K,da=dConfig,g=g))
    }
  } else{
    if(subsample=='strat' & ls.subsample.size<nrow(XT)){
      sample.time = proc.time()[3]
      samp.id = strat.sample(data.frame(XT),ls.subsample.size,m,K)
      sample.time = proc.time()[3] - sample.time
    } else if(subsample=='rand' & ls.subsample.size<nrow(XT)){
      sample.time = proc.time()[3]
      samp.id = lapply(1:K,function(i) sample(1:ncol(V.t),ls.subsample.size))
      sample.time = proc.time()[3] - sample.time
    } else{
      if(verbose)
        message("Estimating with all data.\n")
      K = 1
      samp.id = list(1:nrow(XT))
      sample.time = 0
    }
    get_ls = function(XT,V.t,d,g,K,i){
      gp = laGP::newGPsep(X = XT,Z = V.t,d = d$start,g = g,dK = TRUE)
      mle = list(conv=1)
      i = 0
      maxit=100
      while(mle$conv==1 & i<=1000){
        mle = laGP::mleGPsep(gp,
                             param='d',
                             tmin = d$min,
                             tmax = 10*d$max,
                             ab=d$ab,
                             maxit=maxit)
        i = i+maxit
        if(mle$conv==1){
          message(paste('Warning: mleGPsep for bootstrap replicate',K,'and basis component',i,'has not converged after',i*100,'iterations.'))
        }
      }
      return(mle$d)
    }
    tmp.ls = array(dim=c(n.pc,ncol(XT),K))
    for(k in 1:K){
      if(ls.parallel){
        tmp.ls[,,k] = foreach::foreach(i=1:n.pc,.combine='rbind') %dopar% get_ls(XT[samp.id[[k]],],V.t[i,samp.id[[k]]],dConfig,g,k,i)
      } else{
        tmp.ls[,,k] = t(sapply(1:n.pc, function(i) get_ls(XT[samp.id[[k]],],V.t[i,samp.id[[k]]],dConfig,g,k,i)))
        ##### DEV #####
        # try using mlegp - didn't seem to change much
        # tmp.ls[,,k] = t(sapply(1:n.pc, function(i) 1/nugget=g,nugget.known=1,verbose=0,simplex.ntries = 1)$beta))
      }
    }
    invisible(utils::capture.output(laGP::deleteGPseps()))
    estLenscales = apply(tmp.ls,1:2,median) # median lengthscale over subsamples
    estLenscales = lapply(seq_len(n.pc), function(i) estLenscales[i,]) # convert to list of row vectors
    if(ls.parallel & make.cluster)
      parallel::stopCluster(cl)
  }

  return = list()
  return$XT = estLenscales
  if(p.x>0){
    return$X = lapply(1:n.pc, function(i) estLenscales[[i]][1:p.x])
  } else{
    return$X = list()
  }
  if(p.t>0){
    return$T = lapply(1:n.pc, function(i) estLenscales[[i]][(p.x+1):(p.x+p.t)])
  } else{
    return$T = list()
  }
  return$time = ptm-proc.time()[3]
  return$subsample = subsample
  if(subsample %in% c('strat','rand')){
    return$sample.time = sample.time
    return$samp.id = samp.id
  }
  return$m = m
  return$k = K
  return$g = g
  return$prior = ls.prior
  return$subsample.size = ls.subsample.size
  return$parallel = ls.parallel
  return(return)
}

# FUNCTION: unit_xform
### Transform columns to lie in unit interval
# Parameters
# X (numeric) : matrix to be transformed
# Returns
# Xtrans (numeric) : transformed matrix
unit_xform = function(X,X.min=NULL,X.range=NULL){
  X = as.matrix(X)
  if(is.null(X.min)){
    X.min = apply(X,2,min)
  }
  if(is.null(X.range)){
    Xmax = apply(X,2,max)
    X.range = Xmax - X.min
  }
  if(!isTRUE(all.equal(X.range,rep(0,length(X.range))))){
    Xtrans = t( (t(X) - X.min) / X.range )
  } else{
    # just return X if dummyx
    Xtrans = X
  }
  return(list(orig=X,
              trans=Xtrans,
              min=X.min,
              range=X.range))
}

transform_y = function(Y.sim,y.ind.sim=NULL,Y.obs=NULL,y.ind.obs=NULL,center=T,scale=T,scaletype='scalar'){
  sim.list = list()
  if(!is.null(Y.obs)){
    obs.list = list()
  } else{
    obs.list = NULL
  }


  n.y.sim = nrow(Y.sim)
  if(is.null(y.ind.sim))
    y.ind.sim = matrix(seq(1:n.y.sim))
  dim = ncol(y.ind.sim)

  # scale Y.sim
  sim.list$orig = Y.sim
  if(center){
    sim.list$mean = rowMeans(Y.sim)
    sim.list$trans = sim.list$orig - sim.list$mean
  } else{
    sim.list$mean = rep(0,n.y.sim)
    sim.list$trans = sim.list$orig
  }

  # check
  if(center){
    if(!all(abs(rowMeans(sim.list$trans)) < 1e-8)){
      stop('Y.sim not properly centered \n')
    }
  }

  if(scale){
    if(scaletype=='rowwise'){
      sim.list$sd = apply(sim.list$trans,1,sd)
      # dont scale places with var 0
      sim.list$sd[sim.list$sd==0]=1
    } else if(scaletype=='scalar') {
      if(center){
        # we already centered data so we can compute the sd faster
        sd = sum(sim.list$trans^2)/(prod(dim(sim.list$trans))-1)
      } else{
        # this is breaking for huge matrices
        sd = sd(sim.list$trans)
      }
      sim.list$sd = rep(sd,n.y.sim)
    }
    sim.list$trans = sim.list$trans / sim.list$sd
  } else{
    sim.list$sd = rep(1,n.y.sim)
  }

  # scale Y.obs
  if(!is.null(Y.obs)){
    obs.list$orig = Y.obs
    obs.list$n.y = nrow(Y.obs)
    if(is.null(y.ind.obs))
      y.ind.obs = matrix(seq(1:obs.list$n.y))
    if(!isTRUE(all.equal(y.ind.obs,y.ind.sim))){
      # interpolation is needed
      if(dim==1){
        # 1d interpolation of sim mean onto yobs support
        obs.list$mean = pracma::interp1(as.numeric(y.ind.sim),sim.list$mean,as.numeric(y.ind.obs))
        obs.list$sd = pracma::interp1(as.numeric(y.ind.sim),sim.list$sd,as.numeric(y.ind.obs))
        obs.list$trans = t(scale(t(Y.obs),obs.list$mean,obs.list$sd))
      } else if(dim==2){
        # 2d interpolation of sim mean onto yobs support
        x = unique(y.ind.sim[,1]); y = unique(y.ind.sim[,2])
        xp = y.ind.obs[,1]; yp = y.ind.obs[,2]
        obs.list$mean = pracma::interp2(x,y,matrix(sim.list$mean,length(y),length(x),byrow = T),xp,yp)
        if(scale){
          if(scaletype=='rowwise'){
            obs.list$sd = pracma::interp2(x,y,matrix(sim.list$sd,length(y),length(x),byrow = T),xp,yp)
          } else if(scaletype=='scalar'){
            obs.list$sd = rep(sd,nrow(Y.obs))
          }
        } else{
          obs.list$sd = rep(1,nrow(Y.obs))
        }
        obs.list$trans = t(scale(t(Y.obs),obs.list$mean,obs.list$sd))
      }
    } else{
      obs.list$trans = t(scale(t(Y.obs),sim.list$mean,sim.list$sd))
      obs.list$mean = sim.list$mean
      obs.list$sd = sim.list$sd
    }
  }
  sim.list$ind = y.ind.sim
  sim.list$n.y = nrow(Y.sim)
  obs.list$ind = y.ind.obs
  return(list(sim=sim.list,obs=obs.list,
              n=ncol(Y.obs),m=ncol(Y.sim)))
}

# FUNCTION: transform_xt
### Transforms inputs to lie on the unit hypercube
# T.obs are optional known calibration parameters for the observed data
# transforming these is useful for checking predictions at obs locations
transform_xt = function(X.sim=NULL,T.sim=NULL,X.obs=NULL,T.obs=NULL){
  sim.list = list()
  obs.list = list()
  if(is.null(X.sim) & is.null(T.sim)){
    stop('One of X.sim or T.sim must be specified')
  }
  if(!is.null(X.sim)){
    p.x = ncol(X.sim)
    sim.list$X = unit_xform(X.sim)
  } else{
    p.x = 0
  }
  if(!is.null(T.sim)){
    p.t = ncol(T.sim)
    sim.list$T = unit_xform(T.sim)
  } else{
    p.t = 0
  }

  if(!is.null(X.obs)){
    obs.list$X = unit_xform(X.obs,X.min=sim.list$X$min,X.range=sim.list$X$range)
  }
  if(!is.null(T.obs)){
    obs.list$T = unit_xform(T.obs,X.min=sim.list$T$min,X.range=sim.list$T$range)
  }
  return(list(sim=sim.list,obs=obs.list,p.x=p.x,p.t=p.t))
}

# FUNCTOIN: mean_sd_xform
### Centers and/or scales a matrix
mean_sd_xform = function(X,center=T,scale=T){
  tmp = t(scale(t(X),center,scale))
  if(scale){
    sd = attr(tmp,'scaled:scale')
  }else{
    sd = rep(1,nrow(X))
  }
  return(list(orig=X,
              trans=as.matrix(tmp),
              mean=attr(tmp,'scaled:center'),
              sd=sd))
}

get_SC_inputs = function(ls,XT.data,n.pc){
  return = list()
  if(length(ls$X)>0){
    X.sim = lapply(1:n.pc, function(i) sc_inputs(XT.data$sim$X$trans,ls$X[[i]]))
  } else{
    X.sim = NULL
  }
  if(length(ls$T)>0){
    T.sim =  lapply(1:n.pc, function(i) sc_inputs(XT.data$sim$T$trans,ls$T[[i]]))
  } else{
    T.sim = NULL
  }
  return$XT.sim = lapply(1:n.pc, function(i) cbind(X.sim[[i]],T.sim[[i]]))

  if(!is.null(XT.data$obs$X)){
    return$X.obs = lapply(1:n.pc, function(i) sc_inputs(XT.data$obs$X$trans,ls$X[[i]]))
  }
  if(!is.null(XT.data$obs$T)){
    return$T.obs = lapply(1:n.pc, function(i) sc_inputs(XT.data$obs$T$trans,ls$T[[i]]))
  }
  return$p.x = XT.data$p.x
  return$p.t = XT.data$p.t
  return(return)
}
# FUNCTION: sc_inputs
### Function to apply stretching and compressing to inputs
# Parameters
# X (numeric) : matrix of inputs to be stretched and compressed by lengthscale parameters
# ls: vector of lengthscale parameters length(ls)==ncol(X)
# Returns
# Xtrans (numeric) : transformed matrix
# sc_inputs = function(X,ls){
#   return(t(t(X) / sqrt(ls))) # t(t(X)) is for the case where X is a vector (1d input space)
# }
sc_inputs = function(X,ls){
  return(sweep(X, 2, sqrt(ls), FUN = '/'))
}

#' @title FlaGP data object constructor
#'
#' @description Builds FlaGP data object and does necessary precomputing
#' @param X.sim experimental inputs for simulator output data
#' @param T.sim calibration inputs for simulator output data
#' @param X.obs observed inputs for experimental data
#' @param T.obs optional "known" calibration parameters for observations. This parameter does nothing in the model, but the data object will store scaled copies of it so it can be compared to calibration results
#' @param Y.sim simulator output matrix
#' @param y.ind.sim data indices for functional response simulations
#' @param Y.obs observed data matrix
#' @param y.ind.obs data indices for functional response observations
#' @param center center simulations to mean zero
#' @param scale scale simulations to unit variance
#' @param scaletype "scalar" or "rowwise" scaling to unit variance. Functional outputs should likely use "scalar" and multivariate outputs with different units should use "rowwise"
#' @param n.pc number of basis components to use for emulation. Defaults to 95% variance explained.
#' @param pct.var choose number of basis components s.t. this proportion of variation is accounted for in simulations. Defaults to 95% variance explained.
#' @param B optional matrix of basis vectors, if null B is computed via SVD(Y.sim)
#' @param V.t optional precomputed matrix of simulation weights for B
#' @param sigma.y assumed standard error of observations (not currently implemented)
#' @param ls.subsample type of subsampling procedure for lengthscale estimation. default 'blhs' is bootstrapped block latin hypercube sampling, 'strat' is stratified random sampling, and 'rand' is random sampling. Otherwise no subsampling is used.
#' @param ls.nugget gp nugget for lengthscale estimation
#' @param ls.m m parameter in laGP::blhs controlling data blocking for LHS subsample. Larger m will result in faster lengthscale estimation but is more memory intensive.
#' @param ls.K Number of bootstrap replicate estimations to do when ls.subsample is 'blhs','strat', or 'rand'. Estimates can be highly variable for small K.
#' @param ls.prior use detault prior in laGP::newGP for MAP lengthscale estimation, if FALSE, MLE estimation
#' @param ls.parallel do lengthscale estimation in parallel
#' @param ls.subsample.size size of subsamples when ls.subsample is 'strat' or 'rand'
#' @param bias calibration with a discrepancy model
#' @param D matrix of basis vectors for discrepancy model
#' @param small return small (memory) data object
#' @param verbose print status updates while building data object and doing precomputing
#' @param rsvd use random svd for fast svd calculations on massive Y.sim matrices. Must specify n.pc.
#' @details Returns FlaGP data object
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
flagp = function(X.sim=NULL,T.sim=NULL,X.obs=NULL,T.obs=NULL,                                           # X and T data
                    Y.sim,y.ind.sim=NULL,Y.obs=NULL,y.ind.obs=NULL,center=T,scale=T,scaletype='scalar', # Y data
                    n.pc = NULL, pct.var = .95, B = NULL, V.t = NULL, sigma.y=NULL,                     # sim basis
                    ls.subsample = 'strat', ls.nugget=1e-7, ls.m = 1, ls.K = 1, ls.prior=T, ls.parallel=T, make.cluster=T, ls.subsample.size = 250, # length scale estimation
                    bias=F,D=NULL,                                                                      # discrepancy
                    small=F,seed=NULL,verbose=T,
                    rsvd = F){                                                       # additional flags
  data = list(bias=bias,num=list(m=ncol(Y.sim),
                                 n=max(0,ncol(Y.obs)),
                                 p.x = max(0,ncol(X.sim)),
                                 p.t = max(0,ncol(T.sim)),
                                 g = ls.nugget))
  # print information about data
  if(verbose){
    cat('Building FlaGP data object.\n')
    cat('m:', data$num$m,'\n')
    cat('n:', data$num$n,'\n')
    if(is.null(n.pc)){
      cat('pct.var:', pct.var,'\n')
    } else{
      cat('n.pc:', n.pc,'\n')
    }
    cat('p.x:',data$num$p.x,'\n')
    cat('p.t:',data$num$p.t,'\n')
    cat('dim Y.sim:', dim(Y.sim),'\n')
    cat('dim Y.obs:', dim(Y.obs),'\n')
  }
  if(bias & data$num$n<5)
      warning(paste0('Discrepancy model will be a GP fit to only n=',n,' data points.\n'))
  # precomputing and data building including lengthscale estimation
  start.time = proc.time()[3]
  if(verbose){cat('transforming X,T... ')}
  data$XT.data = transform_xt(X.sim,T.sim,X.obs,T.obs)
  if(verbose){cat('done.\n')}
  if(verbose){cat('transforming Y... ')}
  data$Y.data = transform_y(Y.sim,y.ind.sim,Y.obs,y.ind.obs,center,scale,scaletype)
  if(verbose){cat('done.\n')}
  if(verbose){cat('computing sim basis... ')}
  data$basis = list(); class(data$basis) = c('basis',class(data$basis))
  data$basis$sim = get_basis(data$Y.data$sim$trans,n.pc,pct.var,F,B,V.t,bias=bias,D=D,
                        rsvd=rsvd)
  if(verbose){cat('done.\n')}
  if(!is.null(Y.obs)){
    precomp = TRUE
    if(verbose){cat('computing obs basis... ')}
    data$basis$obs = get_obs_basis(data$basis$sim,data$Y.data$obs$trans,y.ind.sim,y.ind.obs,sigma.y)
    if(verbose){cat('done.\n')}
  } else{
    # em only
    data$basis$obs = NULL
    precomp = FALSE
  }
  if(verbose){cat('estimating lengthscale parameters... ')}
  data$lengthscales = mv_lengthscales(cbind(data$XT.data$sim$X$trans,data$XT.data$sim$T$trans),data$num$p.x,data$num$p.t,
                                      data$basis$sim$V.t,data$num$g,ls.subsample,ls.m,ls.K,seed,ls.prior,ls.parallel,make.cluster,ls.subsample.size,verbose=verbose)
  if(verbose){cat('done.\n')}
  if(verbose){cat('stretching and compressing inputs... ')}
  data$SC.inputs = get_SC_inputs(data$lengthscales, data$XT.data, data$basis$sim$n.pc)
  if(verbose){cat('done.\n')}

  if(precomp){
    precomp.start.time = proc.time()[3]
    if(verbose){cat('precomputing for fast calibration... ')}
    # do precomputing necessary for calibration
    precomp = list()
    BD = cbind(data$basis$obs$B,
               data$basis$obs$D)
    tBD = t(BD)
    precomp$rankBD = pracma::Rank(BD)
    BDtBD = tBD%*%BD
    precomp$BDtBDinv = chol2inv(chol(BDtBD + 1e-8*diag(1,nrow(BDtBD))))
    precomp$ldetBDtBD = determinant(BDtBD)$modulus
    precomp$BDtBDinvtBD = precomp$BDtBDinv%*%tBD
    precomp$LLHmat = diag(1,data$Y.data$obs$n.y)-BD%*%precomp$BDtBDinv%*%tBD
    data$precomp = precomp
    data$precomp$time = proc.time()[3] - precomp.start.time
    if(verbose){cat('done.\n')}
  } else{
    data$precomp = NULL
  }
  if(small){
    # remove potentially huge sim data matrices
    data$Y.data$sim$orig = NULL
    data$Y.data$sim$trans = NULL
  }
  data$time = proc.time()[3] - start.time
  # how can I return the objects in data without making a copy? I need to give it a class.
  class(data) = c('flagp',class(data))
  return(data)
}

# A more lightweight version of the laGP::blhs.loop function that does not return any additional information
# other than the lengthscale estimates
blhs.loop.lite = function (y, X, m, K, da, g = 0.001, maxit = 100)
{
  boot_theta_loop <- matrix(NA, nrow = K, ncol = ncol(X))
  failed=0
  for (kk in 1:K) {
    #sub <- laGP::blhs(y = y, X = X, m = m)
    #######
    # blhs code to avoid another copy
    if (length(y) != nrow(X))
      stop("dimension mismatch")
    if (m <= 0)
      stop("m must be positive")
    if (m > nrow(X)^(1/(ncol(X) - 1)))
      stop("must have <= N^[1/(d-1]")
    D <- as.data.frame(cbind(y, X))
    d <- ncol(D)
    k <- d - 1
    names(D)[1:d] <- c("y", paste("x", 1:k, sep = ""))
    xblock <- paste("block_x", 1:k, sep = "")
    D[xblock] <- 1
    for (i in 1:m) {
      for (j in 1:k) {
        D[, d + j] <- ifelse(D[, j + 1] > (i - 1)/m & D[,
                                                        j + 1] < i/m, i, D[, d + j])
      }
    }
    global_block <- D[, d + 1]
    for (u in 2:k) {
      new <- (D[, d + u] - 1) * m^(u - 1)
      global_block <- cbind(global_block, new)
    }
    D$Global_block <- rowSums(global_block[, 1:k])
    Label <- matrix(NA, nrow = m, ncol = k)
    for (h in 1:k) {
      Label[, h] <- sample(1:m, m, replace = FALSE)
    }
    su <- Label[, 1]
    for (f in 2:k) {
      nn <- (Label[, f] - 1) * m^(f - 1)
      su <- cbind(su, nn)
    }
    sub <- rowSums(su[, 1:k])
    ysub <- D[D$Global_block %in% sub, 1]
    Xsub <- D[D$Global_block %in% sub, 2:d]
    #######
    if(length(ysub)==0){
      if(failed==0){
        cat("\nblhs failed for bootstrap sample",kk,',')
        failed = 1
      } else if(failed==1 & kk<K){
        cat(kk,',')
      } else{
        cat(kk,'\n')
      }
      boot_theta_loop[kk, ] <- NA
    } else{
      gpsepi <- laGP::newGPsep(as.matrix(Xsub), ysub, d = da$start,
                               g = g, dK = TRUE)
      # gpsepi <- laGP::newGPsep(as.matrix(sub$xs), sub$ys, d = da$start,
      #                    g = g, dK = TRUE)
      mle <- laGP::mleGPsep(gpsepi, tmin = da$min, tmax = 10 * da$max,
                            ab = da$ab, maxit = maxit)
      i = 0
      while(mle$conv==1){
        mle <- laGP::mleGPsep(gpsepi, tmin = da$min, tmax = 10 * da$max,
                              ab = da$ab, maxit = maxit)
        i = i+maxit
        if(mle$conv==1){
          message(paste('Warning: mleGPsep for bootstrap replicate',K,'and basis component',i,'has not converged after',i*100,'iterations.'))
        }
      }
      laGP::deleteGPsep(gpsepi)
      boot_theta_loop[kk, ] <- mle$d
    }

  }
  theta.hat_median <- apply(boot_theta_loop, 2, median, na.rm=T)
  if(is.na(theta.hat_median[1]))
    stop('blhs failed.\n')
  return(theta.hat_median)
}
blhs.lite = function (D, m)
{
  d <- ncol(D)
  k <- d - 1
  if (m > nrow(D)^(1/k))
    stop("must have <= N^[1/(d-1]")
  names(D)[1:d] <- c("y", paste("x", 1:k, sep = ""))
  xblock <- paste("block_x", 1:k, sep = "")
  D[xblock] <- 1
  for (i in 1:m) {
    for (j in 1:k) {
      D[, d + j] <- ifelse(D[, j + 1] > (i - 1)/m & D[,
                                                      j + 1] < i/m, i, D[, d + j])
    }
  }
  global_block <- D[, d + 1]
  for (u in 2:k) {
    new <- (D[, d + u] - 1) * m^(u - 1)
    global_block <- cbind(global_block, new)
  }
  D$Global_block <- rowSums(global_block[, 1:k])
  Label <- matrix(NA, nrow = m, ncol = k)
  for (h in 1:k) {
    Label[, h] <- sample(1:m, m, replace = FALSE)
  }
  su <- Label[, 1]
  for (f in 2:k) {
    nn <- (Label[, f] - 1) * m^(f - 1)
    su <- cbind(su, nn)
  }
  sub <- rowSums(su[, 1:k])
  ysub <- D[D$Global_block %in% sub, 1]
  Xsub <- D[D$Global_block %in% sub, 2:d]
  return(list(xs = Xsub, ys = ysub))
}

# wrapper function to call Strata
# generates strata by cutting then calls Strata
strat.sample = function(X.df,n,n.div,K){
  # make a grouping variables
  dx = ncol(X.df)
  nstrata = n.div^dx
  intervals = seq(0,1,length.out=n.div+1)
  intervals[1] = -Inf
  for(i in 1:ncol(X.df)){
    X.df[,ncol(X.df)+1] = cut(X.df[,i],breaks=intervals)
  }
  samples = Strata(X.df,names(X.df)[(dx+1):ncol(X.df)],size=n/nstrata,n=n,K=K)
  return(samples)
}
# Stratified random sampling function. If a strata is empty, random strata are chosen for extra samples
# so that the total number of samples is n.
# DescTools::Strata function returns the wrong ids due to a merge call which adds rows to X and seemingly changes the order
Strata = function (x, stratanames, size, n, K)
{
  factor_fg <- unlist(lapply(x[, stratanames, drop = FALSE],
                             is.factor))
  lvl <- c(lapply(lapply(x[, names(which(!factor_fg)), drop = FALSE],
                         factor), levels), lapply(x[, names(which(factor_fg)),
                                                    drop = FALSE], levels))
  strat <- expand.grid(lvl[stratanames])
  strat$stratum <- factor(1:nrow(strat))
  x$id <- 1:nrow(x)
  x <- merge(x, strat)
  n.strat <- table(x$stratum)
  size = rep(ceiling(size),length(n.strat))
  size[n.strat==0] = 0
  x = x[,c('id','stratum')]
  x = split(x, x$stratum)
  idx = array(dim=c(length(x),max(size),K))

  size.orig = size
  for(k in 1:K){
    # for each K, fix size to sum to n
    size = size.orig
    if(sum(size)<n){
      # randomly add some
      diff = sum(size) - n
      while(diff<0){
        cand.add = which(size>0)
        add = sample(cand.add,1)
        size[add] = size[add] + 1
        diff = sum(size) - n
      }
    } else if(sum(size)>n){
      # randomly remove some
      diff = sum(size) - n
      while(diff>0){
        cand.rm = which(size>0)
        rm = sample(cand.rm,1)
        size[rm] = size[rm] - 1
        diff = sum(size) - n
      }
    }
    for(i in 1:length(size)){
      # sample the i'th strata
      if(size[i]>0){
          tmp <- sample(nrow(x[[i]]),size[i])
          idx[i,1:size[i],k] = x[[i]]$id[tmp]
      }
    }
  }
  # collapse idx
  dim(idx) = c(length(x)*max(size),K)
  # convert to list
  idx = lapply(seq_len(K), function(i) idx[,i]) # convert to list of column vectors
  # remove NA's that indicated an empty strata
  idx = lapply(idx, function(z) z[!is.na(z)]) # convert to list of column vectors
  return(idx)
}
