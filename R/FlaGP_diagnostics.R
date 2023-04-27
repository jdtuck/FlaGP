# ypred_samp n.samples x n.y x n
energy_score = function(y.samp,y,terms=F){
  n.samples = dim(y.samp)[1]
  n.y = dim(y.samp)[2]
  n = dim(y.samp)[3]

  # make sure the field data shape matches the samples shape
  if(!all(dim(y)==c(n.y,n)))
    dim(y) = c(n.y,n)

  get_es = function(i,terms){
    diff = t(t(y.samp[,,i])-y[,i])
    error_term = 2*mean(apply(diff,1,norm,type='2')) # mean over samples of norm of difference
    uq_term = (1/(n.samples^2 - n.samples))*sum(sapply(1:n.samples, function(k)
      sum(sapply(1:n.samples, function(j)
        norm(y.samp[k,,i]-y.samp[j,,i],type='2')))))
    if(terms){
      return(list(error_term=error_term,uq_term=uq_term))
    } else{
      return(error_term - uq_term)
    }

  }
  if(terms){
    es = matrix(unlist(parallel::mclapply(1:n, function(i) get_es(i,terms))),ncol=2,byrow=T)
  } else{
    es = unlist(
      parallel::mclapply(1:n, function(i) get_es(i,terms))
    )
  }

  return(es)
}
# y is a matrix of dim (n.y x n.pred)
# y.samp is a matrix of dim (n.samp x n.y x n.pred)
# y.conf is an array of dim (2 x n.y x n.pred)
interval_score = function(y,alpha=.05,y.samp=NULL,y.conf.int=NULL,terms=F){
  if(!is.null(y.samp)){
    n = dim(y.samp)[3]
    n.y = dim(y.samp)[1]
    y.conf.int = apply(y.samp,c(1,3),quantile,c(.025,.975))
  } else{
    n = dim(y.conf.int)[3]
    n.y = dim(y.conf.int)[2]
  }
  # make sure the field data shape matches the samples shape
  if(!all(dim(y)==c(n.y,n)))
    dim(y) = c(n.y,n)

  width = array(dim=c(n.y,n))
  penalty = array(dim=c(n.y,n))
  i_score = array(dim=c(n.y,n))

  for(k in 1:n){ # loop over pred locations
    # i_score = (u-l)+2/alpha(l-x)I(x<l)+2/alpha(x-u)I(x>u)
    width[,k] = (y.conf.int[2,,k]-y.conf.int[1,,k])
    penalty[,k] = (2/alpha)*(y.conf.int[1,,k] - y[,k])*as.integer(y[,k]<y.conf.int[1,,k]) +
      (2/alpha)*(y[,k] - y.conf.int[2,,k])*as.integer(y[,k]>y.conf.int[2,,k])
    i_score[,k] = width[,k] + penalty[,k]
  }
  if(terms){
    return(list(i_score=i_score,width=width,penalty=penalty))
  } else{
    return(i_score)
  }
}

plot_basis = function(b,y.ind.sim,legend=T,xlab='x')
{

  w = b$sim$V.t
  z = b$obs$V.t
  n.pc = nrow(w)
  m = ncol(w)

  # plot simulation basis vectors
  matplot(y.ind.sim,b$sim$B,type='l',xlab=xlab,ylab='Sim Basis Vectors')
  if(legend){
    legend("bottomleft", inset=c(.05,.05), legend = paste0('PC ',1:n.pc),col=seq_len(n.pc),cex=1,lty=seq_len(n.pc))
  }

  # plot simulation and observed basis weights
  if(!is.null(z)){
    n = ncol(z)
    if(n.pc>1){
      pairs(rbind(t(w),t(z)),col=c(rep('blue',m),rep('orange',n)),pch=c(rep(1,m),rep(16,n)),labels=paste0('PC ',1:n.pc),
            oma=c(3,3,3,15))
      if(legend){
        par(xpd = TRUE)
        legend("bottomright", pch=1, col=c('blue','orange'), legend=c('w','z'))
      }
    } else{
      hist(t(w),col='blue');abline(v=t(z),col='orange')
      if(legend){
        legend("topright", pch=1, col=c('blue','orange'), legend=c('w','z'))
      }
    }
  } else{
    if(n.pc>1){
      labs = paste0('PC',1:n.pc)
      pairs(t(w),col=rep('blue',m),pch=1,labels=labs,
            oma=c(3,3,3,15),asp=1)
      if(legend){
        par(xpd = TRUE)
        legend("bottomright", pch=1, col=c('blue'), legend=c('w'))
      }
    } else{
      hist(t(w),col='blue');
      if(legend){
        legend("topright", pch=1, col=c('blue'), legend=c('w'))
      }
    }
  }
}

#' @title FlaGP data plots
#'
#' @description Plots X,Y data as well as computed bases
#' @param flagp an \code{flagp} object.
#' @param basis T/F indicating if the bases should be plotted
#' @param legend T/F indicating if legends should be added to the plots
#' @param xlab an optional descriptive label for the data x axes
#' @param ylab an optional descriptive label for the data y axes
#' @param ... additional graphical parameters.
#' @details Plots the X,Y data, both simulated and observed and optionally plots the bases as well
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
plot.flagp = function(flagp,basis=T,legend=T,xlab='x',ylab='y',...){
  if(ncol(flagp$Y.data$sim$ind == 1)){
    par(mfrow=c(1,2),mar=c(4,4,4,.3))
    matplot(flagp$Y.data$sim$ind,flagp$Y.data$sim$orig,type='l',col='orange',xlab=xlab,ylab=ylab,main='Original Scale',...)
    if(!is.null(flagp$Y.data$obs)){
      if(ncol(flagp$Y.data$obs$ind) == 1){
        matplot(flagp$Y.data$obs$ind,flagp$Y.data$obs$orig,pch=1,add=T,col='black',...)
        matplot(flagp$Y.data$obs$ind,flagp$Y.data$obs$orig,type='l',lty=1,add=T,col='black',...)
      }
    }
    legend('topleft',inset=c(.05,.05),legend=c('simulations','field obs.'),lty=1,col=c('orange','black'),...)

    matplot(flagp$Y.data$sim$ind,flagp$Y.data$sim$trans,type='l',col='orange',xlab=xlab,ylab=ylab,main='Standardized Scale',...)
    if(!is.null(flagp$Y.data$obs)){
      if(ncol(flagp$Y.data$obs$ind) == 1){
        matplot(flagp$Y.data$obs$ind,flagp$Y.data$obs$trans,pch=1,add=T,col='black',...)
        matplot(flagp$Y.data$obs$ind,flagp$Y.data$obs$trans,type='l',lty=1,add=T,col='black',...)
      }
    }

    if(!is.null(flagp$Y.data$obs)){
      legend('topleft',inset=c(.05,.05),legend=c('simulations','field obs.'),lty=1,col=c('orange','black'),...)
    } else{
      legend('topleft',inset=c(.05,.05),legend=c('simulations'),lty=1,col=c('orange'),...)
    }
  }else{
    stop('Default data visualization not implemented for 2-D functions.')
  }

  if(basis){
    plot_basis(flagp$basis,flagp$Y.data$sim$ind,legend,xlab,...)
  }
}

#' @title MCMC Plot
#'
#' @description Generate plots to assess MCMC
#' @param x an \code{mcmc} object.
#' @param labels an optional vector of labels for each calibration parameter
#' @param nrow number of rows in display grid of plots
#' @param ncol number of columns in display grid of plots
#' @param ... additional graphical parameters.
#' @details A pairs plot and trace plots are displayed in an nrow x ncol grid
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
plot.mcmc = function(x, labels=NULL, nrow=2, ncol=2, ...){
  if(class(x)[1]!='mcmc')
    stop('x must be an object of class mcmc')
  p.t = ncol(x$t.samp)
  if(is.null(labels)){labels = paste0('t ',1:p.t)}
  if(p.t>1){
    panel.hist <- function(x,...) {
      usr <- par("usr")
      on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5))
      his <- hist(x, plot = FALSE)
      breaks <- his$breaks
      nB <- length(breaks)
      y <- his$counts
      y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1],y, ...)
      # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
    }
    pairs(x$t.samp,diag.panel = panel.hist, labels = labels, ...)
    par(mfrow=c(nrow,ncol))
    for(i in 1:p.t){
      plotlims = c(max(0,min(x$t.samp[,i])-.05),min(1,max(x$t.samp[,i])+.05))
      plot(x$t.samp[,i],type='l',ylab=labels[i],xlab='MCMC iteration (post-burn)',ylim=plotlims, ...)
    }
    plot(x$ssq.samp,type='l',ylab='error variance',xlab='MCMC iteration (post-burn)', ...)
  } else{
    par(mfrow = c(2,1))
    plot(x$ssq.samp,type='l',ylab='error variance',xlab='MCMC iteration (post-burn)', ...)
    plotlims = c(max(0,min(x$t.samp[,1])-.05),min(1,max(x$t.samp[,1])+.05))
    plot(x$t.samp[,1],type='l',ylab=labels[1],xlab='MCMC iteration (post-burn)',ylim=plotlims, ...)
  }
}

