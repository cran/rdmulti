###################################################################
# rdmcplot: RD plots with multiple cutoffs
# !version 0.4 07-Jan-2020
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' RD plots with multiple cutoffs.
#'
#' \code{rdmc()} RD plots with multiple cutoffs.
#'
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}
#'
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}
#'
#' Gonzalo Vazquez-Bare, UC Santa Barbara. \email{gvazquez@econ.ucsb.edu}
#'
#' @references
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2019). \href{https://sites.google.com/site/rdpackages/rdmulti/Cattaneo-Titiunik-VazquezBare_2019_rdmulti.pdf}{Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}. \emph{Working paper}.
#'
#'
#' @param Y outcome variable.
#' @param X running variable.
#' @param C cutoff variable.
#' @param pvec vector of cutoff-specific polynomial orders. See \code{rdplot()} for details.
#' @param hmat matrix of cutoff-specific bandwidths. See \code{rdplot()} for details.
#' @param nbinsmat matrix of cutoff-specific number of bins. See \code{rdplot()} for details.
#' @param binselectvec vector of cutoff-specific bins selection method. See \code{rdplot()} for details.
#' @param scalevec vector of cutoff-specific scale factors. See \code{rdplot()} for details.
#' @param kernelvec vector of cutoff-specific kernels. See \code{rdplot()} for details.
#' @param weightsvec vector of cutoff-specific weights. See \code{rdplot()} for details.
#' @param supportmat matrix of cutoff-specific support conditions. See \code{rdplot()} for details..
#' @param covsvec vector of cutoff-specific covariates. See \code{rdplot()} for details.
#' @param nobins omits bins plot.
#' @param nopoly omits polynomial curve plot.
#' @param nodraw omits plot.
#'
#'
#' @return
#' \item{clist}{list of cutoffs}
#' \item{cnum}{number of cutoffs}
#' \item{X0}{matrix of X values for control units}
#' \item{X1}{matrix of X values for treated units}
#' \item{Yhat0}{estimated polynomial for control units}
#' \item{Yhat1}{estimated polynomial for treated units}
#' \item{Xmean}{bin average of X values}
#' \item{Ymean}{bin average for Y values}
#'
#'
#' @examples
#' # Toy dataset
#' X <- runif(1000,0,100)
#' C <- c(rep(33,500),rep(66,500))
#' Y <- (1 + X + (X>=C))*(C==33)+(.5 + .5*X + .8*(X>=C))*(C==66) + rnorm(1000)
#' # rdmcplot with standard syntax
#' tmp <- rdmcplot(Y,X,C)
#'
#'
#' @export

rdmcplot = function(Y,X,C,
                    pvec=NULL,
                    hmat=NULL,
                    nbinsmat=NULL,
                    binselectvec=NULL,
                    scalevec=NULL,
                    kernelvec=NULL,
                    weightsvec=NULL,
                    supportmat=NULL,
                    covsvec=NULL,
                    nobins=FALSE,
                    nopoly=FALSE,
                    nodraw=FALSE){

  #################################################################
  # Setup and error checking
  #################################################################

  if (!is.numeric(C)){stop('C has to be numeric')}
  if (max(C,na.rm=TRUE)>=max(X,na.rm=TRUE) | min(C,na.rm=TRUE)<=min(X,na.rm=TRUE)){stop('cutoff variable outside range of running variable')}

  clist = sort(unique(C))
  cnum = length(clist)

  D = as.numeric(X>=C)

  if (is.null(pvec)){pvec = rep(4,cnum)}
  if (!is.null(hmat)){
    if (is.null(dim(hmat))){
      hmat = matrix(hmat,nrow=cnum,ncol=2)
    }
    haux = hmat
  }
  else{
    haux = matrix(Inf,ncol=2,nrow=cnum)
  }
  if (!is.null(nbinsmat)){
    if (is.null(dim(nbinsmat))){
      nbinsmat = matrix(nbinsmat,nrow=cnum,ncol=2)
    }
  }
  if (is.null(binselectvec)){binselectvec = rep('esmv',cnum)}
  if (is.null(scalevec)){scalevec = rep(1,cnum)}
  if (is.null(kernelvec)){kernelvec = rep('uni',cnum)}

  X0 = NULL
  X1 = NULL
  YHAT0 = NULL
  YHAT1 = NULL
  XMEAN = NULL
  YMEAN = NULL

  colorlist = c('darkblue','darkred','darkgreen','darkorange','gray50','khaki4','brown3','blue','darkgoldenrod4','cyan4')


  #################################################################
  # Construct variables for plots
  #################################################################

  count = 1
  for (c in clist){

    yc = Y[C==c & X<=c+haux[count,2] & X>=c-haux[count,1]]
    xc = X[C==c & X<=c+haux[count,2] & X>=c-haux[count,1]]
    dc = D[C==c & X<=c+haux[count,2] & X>=c-haux[count,1]]
    yc0 = yc[dc==0]
    yc1 = yc[dc==1]
    xc0 = xc[dc==0]
    xc1 = xc[dc==1]

    aux = rdrobust::rdplot(yc,xc,c=c,
                           p=pvec[count],
                           h=hmat[count,],
                           nbins=nbinsmat[count,],
                           binselect=binselectvec[count],
                           scale=scalevec[count],
                           kernel=kernelvec[count],
                           weights=weightsvec[count],
                           covs=covsvec[count],
                           support=supportmat[count,],
                           hide=TRUE)

    xmean = aux$vars_bins[,2]
    ymean = aux$vars_bins[,3]
    x0 = aux$vars_poly[aux$vars_poly[,1]<c,1]
    yhat0 = aux$vars_poly[aux$vars_poly[,1]<c,2]
    x1 = aux$vars_poly[aux$vars_poly[,1]>c,1]
    yhat1 = aux$vars_poly[aux$vars_poly[,1]>c,2]

    length(xmean) = length(Y)
    length(ymean) = length(Y)
    length(x0) = length(Y)
    length(yhat0) = length(Y)
    length(x1) = length(Y)
    length(yhat1) = length(Y)

    XMEAN = cbind(XMEAN,xmean)
    YMEAN = cbind(YMEAN,ymean)
    X0 = cbind(X0,x0)
    X1 = cbind(X1,x1)
    YHAT0 = cbind(YHAT0,yhat0)
    YHAT1 = cbind(YHAT1,yhat1)

    count = count + 1

  }

  #################################################################
  # Plots
  #################################################################

  if (nodraw==FALSE){

    xlim = c(min(X),max(X))
    ylim = c(min(YHAT0,YHAT1,na.rm=TRUE),max(YHAT0,YHAT1,na.rm=TRUE))

    plot(NA,xlim=xlim,ylim=ylim,xlab='Running variable',ylab='Outcome')

    if (nobins==FALSE & nopoly==FALSE){
      for (c in 1:cnum){
        lines(X0[,c],YHAT0[,c],col=colorlist[c])
        lines(X1[,c],YHAT1[,c],col=colorlist[c])
        points(XMEAN[,c],YMEAN[,c],col=colorlist[c])
        abline(v=clist[c],col=colorlist[c],lty='dotted')
      }
    }
    else if (nobins==TRUE & nopoly==FALSE){
      for (c in 1:cnum){
        lines(X0[,c],YHAT0[,c],col=colorlist[c])
        lines(X1[,c],YHAT1[,c],col=colorlist[c])
        abline(v=clist[c],col=colorlist[c],lty='dotted')
      }
    }
    else if (nobins==FALSE & nopoly==TRUE){
      for (c in 1:cnum){
        points(XMEAN[,c],YMEAN[,c],col=colorlist[c])
        abline(v=clist[c],col=colorlist[c],lty='dotted')
      }
    }

  }

  #################################################################
  # Return values
  #################################################################

  output = list(clist = clist,
                cnum = cnum,
                X0 = X0,
                X1 = X1,
                Yhat0 = YHAT0,
                Yhat1 = YHAT1,
                Xmean = XMEAN,
                Ymean = YMEAN)

  return(output)
}
