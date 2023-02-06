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
get_basis = function(Y, n.pc = 1, pct.var = NULL, full.basis=F, B=NULL, bias=F, D=NULL,
                     rsvd=FALSE, k = NULL, nu = NULL, nv = NULL, p = 10, q = 2, sdist = "normal"){

  return.list = list()

  # Compute SVD
  nExp = ncol(Y)
  if(is.null(B)){
    ptm = proc.time()[3]
    if(!rsvd){
      svdY = svd(Y)
    } else{
      svdY = rsvd(Y,k,nu,nv,p,q,sdist)
    }
    return.list$svd.time = proc.time()[3] - ptm
    return.list$B = svdY$u %*% diag(svdY$d) / sqrt(nExp)
    return.list$V.t = t(svdY$v) * sqrt(nExp)

    # Return only the required number of basis vectors
    percent.variance = zapsmall((svdY$d) / sum(svdY$d))
    if(!is.null(pct.var)){
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
  obs.basis$V.t = solve(t(obs.basis$B)%*%obs.basis$B)%*%t(obs.basis$B)%*%Y.obs
  return(obs.basis)
}

# FUNCTION: mv_lengthscales
# TO DO: this function needs to be updated for only using a subset of the data
### Computes estimated length scale parameters for use in stretching and compressing the inputs space
# Parameters:
# X (numeric) : untransformed matrix of inputs
# V.t (numeric) : matrix
# ls_prior_ab (numeric) : 2 vector of parameters for gamma prior on length-scales
mv_lengthscales = function(XT.data,V.t,g=1e-7,subsample=F,subsample.size=100,seed=NULL,ls.prior=F){
  ptm = proc.time()
  XT = cbind(XT.data$sim$X$trans,XT.data$sim$T$trans)
  if(subsample){
    if(!is.null(seed)){
      set.seed(seed)
    }
    samp.id = sample(1:nrow(XT),subsample.size)
    XT = XT[samp.id,]
    V.t = V.t[,samp.id]
  }
  p.x = XT.data$p.x
  p.t = XT.data$p.t
  n.pc = nrow(V.t)

  dConfig = laGP::darg(d=list(mle = TRUE, max = 100),X = XT)
  if(!ls.prior){dConfig$ab = c(0,0)}
  GPlist = vector(mode='list',length=n.pc)
  for(i in 1:n.pc){
    GPlist[[i]] = laGP::newGPsep(X = XT,
                           Z = V.t[i, ],
                           d = rep(dConfig$start, ncol(XT)),
                           g = g,
                           dK = TRUE)
  }
  # try doParallel::foreach (not sure how to use it with pointers to C objects)
  cores = min(n.pc,parallel::detectCores()-4)
  cl <- parallel::makeCluster(cores)
  estLenscales = parallel::mclapply(1:n.pc, function(i) laGP::mleGPsep(GPlist[[i]],param='d',
                                                      tmin = dConfig$min,
                                                      tmax = dConfig$max,
                                                      ab=dConfig$ab,
                                                      maxit=100)$d,
                         mc.cores = parallel::detectCores())
  parallel::stopCluster(cl)
  for(i in n.pc){
    laGP::deleteGPsep(GPlist[[i]])
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
  return$time = ptm-proc.time()
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

transform_y = function(Y.sim,y.ind.sim=NULL,Y.obs=NULL,y.ind.obs=NULL,center=T,scale=T, scaletype='scalar'){
  sim.list = list()
  obs.list = list()

  n.y.sim = nrow(Y.sim)
  if(is.null(y.ind.sim))
    y.ind.sim = matrix(seq(1:n.y.sim))
  dim = ncol(y.ind.sim)

  # scale Y.sim
  sim.list$orig = Y.sim
  if(center){
    sim.list$mean = rowMeans(Y.sim)
  } else{
    sim.list$mean = rep(0,n.y.sim)
  }

  # center
  sim.list$trans = sim.list$orig - sim.list$mean
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
      sd = sd(sim.list$trans)
      sim.list$sd = rep(sd,n.y.sim)
    }
  } else{
    sim.list$sd = rep(1,n.y.sim)
  }

  # scale
  sim.list$trans = sim.list$trans / sim.list$sd

  # check
  if(scale){
    if(scaletype=='rowwise'){
      if(!all(abs(apply(sim.list$trans,1,sd)-1) < 1e8)){
        stop('Y.sim not properly scaled to unit variance in rows \n')
      }
    } else if(scaletype=='scalar'){
      if(!abs(sd(sim.list$trans)-1) < 1e-8){
        stop('Y.sim not properly scaled to unit variance \n')
      }
    }
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
#' @details Returns FlaGP data object
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
flagp = function(X.sim=NULL,T.sim=NULL,X.obs=NULL,T.obs=NULL,                              # X and T data
                    Y.sim,y.ind.sim=NULL,Y.obs=NULL,y.ind.obs=NULL,center=T,scaletype='scalar', # Y data
                    n.pc = NULL, pct.var = .95, B = NULL, sigma.y=NULL,                           # sim basis
                    sc.nugget=1e-7, sc.subsample = F, sc.subsample.size = 250, ls.prior=T,      # length scale estimation
                    bias=F,D=NULL,small=F,precomp=T,verbose=F){
  start.time = proc.time()
  if(verbose){cat('transforming X, T \n')}
  XT.data = transform_xt(X.sim,T.sim,X.obs,T.obs)
  if(verbose){cat('transforming Y \n')}
  Y.data = transform_y(Y.sim,y.ind.sim,Y.obs,y.ind.obs,center,T,scaletype)
  if(verbose){cat('computing sim basis \n')}
  basis = list(); class(basis) = c('basis',class(basis))
  basis$sim = get_basis(Y.data$sim$trans,n.pc,pct.var,F,B,bias=bias,D=D)
  if(!is.null(Y.obs)){
    if(verbose){cat('computing obs basis \n')}
    basis$obs = get_obs_basis(basis$sim,Y.data$obs$trans,y.ind.sim,y.ind.obs,sigma.y)
  } else{
    # em only
    basis$obs = NULL
    precomp = FALSE
  }
  if(verbose){cat('estimating lengthscale parameters \n')}
  ls = mv_lengthscales(XT.data,basis$sim$V.t,sc.nugget,sc.subsample,sc.subsample.size,ls.prior=ls.prior)
  if(verbose){cat('stretching and compressing inputs \n')}
  SC.inputs = get_SC_inputs(ls,XT.data,basis$sim$n.pc)

  data = list(XT.data=XT.data,
              Y.data=Y.data,
              basis = basis,
              lengthscales=ls,
              SC.inputs=SC.inputs,
              bias=bias)

  if(precomp){
    if(verbose){cat('precomputing for fast calibration \n')}
    # do precomputing necessary for calibration
    precomp = list()
    precomp$I.n.y = diag(1,Y.data$obs$n.y)
    BD = cbind(basis$obs$B,basis$obs$D)
    tBD = t(BD)
    precomp$rankBD = pracma::Rank(BD)
    BDtBD = tBD%*%BD
    precomp$BDtBDinv = solve(BDtBD + 1e-8*diag(1,nrow(BDtBD)))
    precomp$ldetBDtBD = determinant(BDtBD)$modulus
    precomp$BDtBDinvtBD = precomp$BDtBDinv%*%tBD
    precomp$LLHmat = precomp$I.n.y-BD%*%precomp$BDtBDinv%*%tBD
    data$precomp = precomp
  } else{
    data$precomp = NULL
  }
  if(small){
    # remove potentially huge sim data matrices
    data$Y.data$sim$orig = NULL
    data$Y.data$sim$trans = NULL
  }
  data$time = proc.time() - start.time
  class(data) = c('flagp',class(data))
  return(data)
}
get_calib_data = function(flagp){

  # store the objects needed for calibration in an object smaller than flagp, which can be inefficient to pass around through functions
  flagp = flagp
  # remove the things from flagp that are expensive to pass around during MCMC and aren't needed
  flagp$Y.data$sim = NULL
  flagp$Y.data$obs$orig = NULL
  flagp$XT.data$sim = NULL; flagp$XT.data$obs = NULL
  flagp$basis$sim$V.t = NULL
  flagp$basis$obs$V.t = NULL
  flagp$SC.inputs$X.sim = NULL; flagp$SC.inputs$T.sim = NULL; flagp$SC.inputs$T.obs = NULL

  flagp$n = flagp$Y.data$n
  flagp$n.y = flagp$Y.data$obs$n.y
  flagp$I.n.y = diag(1,flagp$n.y)
  flagp$BD = cbind(flagp$basis$obs$B,flagp$basis$obs$D)
  flagp$tBD = t(flagp$BD)
  flagp$rankBD = rankMatrix(flagp$BD)
  flagp$BDtBD = flagp$tBD%*%flagp$BD
  if(flagp$rankBD != ncol(flagp$BD)){
    flagp$BDtBDinv = solve(flagp$BDtBD + 1e-8*diag(1,nrow(flagp$BDtBD)))
  } else{
    flagp$BDtBDinv = solve(flagp$BDtBD) # need to get rid of solves to make code faster
  }
  flagp$ldetBDtBD = determinant(flagp$BDtBD)$modulus
  flagp$BDtBDinvtBD = flagp$BDtBDinv%*%flagp$tBD
  flagp$LLHmat = flagp$I.n.y-flagp$BD%*%flagp$BDtBDinv%*%flagp$tBD
  flagp$Z = lapply(1:flagp$basis$sim$n.pc,function(jj) flagp$basis$sim$V.t[jj,])
  ####
  return(flagp)
}
