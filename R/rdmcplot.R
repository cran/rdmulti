###################################################################
# rdmcplot: RD plots with multiple cutoffs
# !version 0.2 12-Jul-2018
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' RD plots with multiple cutoffs.
#'
#' \code{rdmc()} RD plots with multiple cutoffs.
#'
#'
#' @author
#' Matias Cattaneo, University of Michigan. \email{cattaneo@umich.edu}
#'
#' Rocio Titiunik, University of Michigan. \email{titiunik@umich.edu}
#'
#' Gonzalo Vazquez-Bare, University of Michigan. \email{gvazquez@umich.edu}
#'
#' @references
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2018). \href{https://sites.google.com/site/rdpackages/rdmulti/Cattaneo-Titiunik-VazquezBare_2018_rdmulti.pdf}{Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}. \emph{Working paper, University of Michigan}.
#'
#'
#' @param Y outcome variable.
#' @param X running variable.
#' @param C cutoff variable.
#' @param hvec bandwidths to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param pvec order of the polynomials to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param noscatter omits scatter plot.
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
                    hvec=NULL,
                    pvec=NULL,
                    noscatter=FALSE,
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
  if (is.null(hvec)){hvec = rep(Inf,cnum)}

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

  c = clist[1]
  h = hvec[1]
  p = pvec[1]

  yc = Y[C==c & X<=c+h & X>=c-h]
  xc = X[C==c & X<=c+h & X>=c-h]
  dc = D[C==c & X<=c+h & X>=c-h]
  yc0 = yc[dc==0]
  yc1 = yc[dc==1]
  xc0 = xc[dc==0]
  xc1 = xc[dc==1]

  xseq0 = seq(min(xc0,na.rm=TRUE),max(xc0,na.rm=TRUE),length.out=length(xc0))
  XP0 = poly(xseq0-c,degree=p,raw=TRUE)
  XP0 = cbind(rep(1,nrow(XP0)),XP0)

  xseq1 = seq(min(xc1,na.rm=TRUE),max(xc1,na.rm=TRUE),length.out=length(xc1))
  XP1 = poly(xseq1-c,degree=p,raw=TRUE)
  XP1 = cbind(rep(1,nrow(XP1)),XP1)

  aux = rdrobust::rdplot(yc,xc,c=c,p=p,hide=TRUE)
  aux1 = aux$genvars
  xmean = aux1$rdplot_mean_x
  xmean[xmean<c & !is.na(xmean)] = sort(xmean[xmean<c & !is.na(xmean)]) # remove if rdplot is fixed
  ymean = aux1$rdplot_mean_y

  yhat0 = XP0%*%aux$coef[,1]
  yhat1 = XP1%*%aux$coef[,2]

  length(xseq0) = length(Y)
  length(xseq1) = length(Y)
  length(yhat0) = length(Y)
  length(yhat1) = length(Y)
  length(xmean) = length(Y)
  length(ymean) = length(Y)

  X0 = cbind(X0,xseq0)
  YHAT0 = cbind(YHAT0,yhat0)
  X1 = cbind(X1,xseq1)
  YHAT1 = cbind(YHAT1,yhat1)
  XMEAN = cbind(XMEAN,xmean)
  YMEAN = cbind(YMEAN,ymean)


  count = 2
  for (c in clist[-1]){

    h = hvec[count]
    p = pvec[count]

    yc = Y[C==c & X<=c+h & X>=c-h]
    xc = X[C==c & X<=c+h & X>=c-h]
    dc = D[C==c & X<=c+h & X>=c-h]
    yc0 = yc[dc==0]
    yc1 = yc[dc==1]
    xc0 = xc[dc==0]
    xc1 = xc[dc==1]

    xseq0 = seq(min(xc0,na.rm=TRUE),max(xc0,na.rm=TRUE),length.out=length(xc0))
    XP0 = poly(xseq0-c,degree=p,raw=TRUE)
    XP0 = cbind(rep(1,nrow(XP0)),XP0)

    xseq1 = seq(min(xc1,na.rm=TRUE),max(xc1,na.rm=TRUE),length.out=length(xc1))
    XP1 = poly(xseq1-c,degree=p,raw=TRUE)
    XP1 = cbind(rep(1,nrow(XP1)),XP1)

    aux = rdrobust::rdplot(yc,xc,c=c,p=p,hide=TRUE)
    aux1 = aux$genvars
    xmean = aux1$rdplot_mean_x
    xmean[xmean<c & !is.na(xmean)] = sort(xmean[xmean<c & !is.na(xmean)]) # remove if rdplot is fixed
    ymean = aux1$rdplot_mean_y

    yhat0 = XP0%*%aux$coef[,1]
    yhat1 = XP1%*%aux$coef[,2]

    length(xseq0) = length(Y)
    length(xseq1) = length(Y)
    length(yhat0) = length(Y)
    length(yhat1) = length(Y)
    length(xmean) = length(Y)
    length(ymean) = length(Y)

    X0 = cbind(X0,xseq0)
    YHAT0 = cbind(YHAT0,yhat0)
    X1 = cbind(X1,xseq1)
    YHAT1 = cbind(YHAT1,yhat1)
    XMEAN = cbind(XMEAN,xmean)
    YMEAN = cbind(YMEAN,ymean)

    count = count + 1

  }


  #################################################################
  # Plots
  #################################################################

  if (nodraw==FALSE){

    xlim = c(min(X),max(X))
    ylim = c(min(YHAT0,YHAT1,na.rm=TRUE),max(YHAT0,YHAT1,na.rm=TRUE))

    plot(NA,xlim=xlim,ylim=ylim,xlab='Running variable',ylab='Outcome')

    for (c in 1:cnum){
      lines(X0[,c],YHAT0[,c],col=colorlist[c])
      lines(X1[,c],YHAT1[,c],col=colorlist[c])
      abline(v=clist[c],col=colorlist[c],lty='dotted')
    }

    if (noscatter==FALSE){
      for (c in 1:cnum){
        points(XMEAN[,c],YMEAN[,c],col=colorlist[c])
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
