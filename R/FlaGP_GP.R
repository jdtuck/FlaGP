GP_fit_isotropic = function(X,d,g,y=NULL,lite=F)
{
  n = nrow(X)
  dx = ncol(X)

  out = .C("covar_symm_R",
           col = as.integer(dx),
           X = as.double(t(X)),
           n = as.integer(n),
           d = as.double(d),
           g = as.double(g),
           K = as.double(matrix(0,n,n)),
           PACKAGE = "FlaGP")
  dim(out$K) = c(n,n)
  returns = list(K = out$K,
                 Kchol = chol(out$K))
  if(!lite){
    if(is.null(y))
      stop('Must give y if lite=F')
    returns$Ki = chol2inv(returns$Kchol)
    returns$phi = t(y)%*%returns$Ki%*%y
  }
  returns$d = d
  returns$g = g
  returns$n = n
  returns$dx = dx
  return(returns)
}
# Given the cholesky of the covariance matrix draw a random sample
# from the mean zero with stderr 'sd'
GP_sample = function(gp_fit,sd)
{
  sd^2*drop(rnorm(nrow(gp_fit$Kchol)) %*% gp_fit$Kchol)
}

# With covariance K defined on points xy predict from the GP
# at XX using training points y
GP_predict = function(gp_fit,X,XX,y,predvar=F)
{
  # cross-covariance of training locations and prediction locations
  n1 = nrow(XX)
  n2 = nrow(X)
  dx = ncol(X)
  out = .C("covar_R",
           col = as.integer(dx),
           X1 = as.double(t(XX)),
           n1 = as.integer(n1),
           X2 = as.double(t(X)),
           n2 = as.integer(n2),
           d = as.double(gp_fit$d),
           K = as.double(matrix(0,n1,n2)),
           PACKAGE = "FlaGP")
  k = out$K
  dim(k) = c(n1,n2)
  # mean prediction
  ktKi = k%*%gp_fit$Ki
  returns = list(mean=ktKi%*%y)

  if(predvar){
    phidf = gp_fit$phi/n2
    ktKik = ktKi%*%t(k)
    returns$s2 = phidf*(1+gp_fit$g-ktKik)
  }
  return(returns)
}

# compute power exponential covariance matrix at points x,y
gauss_cov = function(x1,x2=NULL,d=1,g=1e-8)
{
  if(is.null(x2)){
    # symmetric case
    dist = as.matrix(distances::distances(x1))
    K = exp(-dist^2/d) + diag(g,nrow(dist))
  } else{
    # nonsymmetric case
    dist = sqrt(plgp::distance(x1,x2))
    K = exp(-dist^2/d)
  }
  return(K)
}
# compute matern covariance matrix at points x,y
matern_cov = function(rho,x1,x2=NULL,nu=7/2)
{
  if(is.null(x2)){
    d = distances::distances(x1)
  } else{
    d = sqrt(plgp::distance(x1,x2))
  }
  fields::Matern(d,range=rho,smoothness=nu)
}
