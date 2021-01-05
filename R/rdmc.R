###################################################################
# rdmc: analysis of RD designs with multiple cutoffs
# !version 0.7 30-Dec-2020
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Analysis of RD designs with multiple cutoffs
#'
#' \code{rdmc()} analyzes RD designs with multiple cutoffs.
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
#' Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2020). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf}{Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}. \emph{Stata Journal}, forthcoming.
#'
#'
#' @param Y outcome variable.
#' @param X running variable.
#' @param C cutoff variable.
#' @param fuzzy specifies a fuzzy design. See \code{rdrobust()} for details.
#' @param derivvec vector of cutoff-specific order of derivatives. See
#'   \code{rdrobust()} for details.
#' @param pooled_opt options to be passed to \code{rdrobust()} to calculate
#'   pooled estimand.
#' @param verbose displays the output from \code{rdrobust} for estimating the
#'   pooled estimand.
#' @param pvec vector of cutoff-specific polynomial orders. See
#'   \code{rdrobust()} for details.
#' @param qvec vector of cutoff-specific polynomial orders for bias estimation.
#'   See \code{rdrobust()} for details.
#' @param hmat matrix of cutoff-specific bandwidths. See \code{rdrobust()} for
#'   details.
#' @param bmat matrix of cutoff-specific bandwidths for bias estimation. See
#'   \code{rdrobust()} for details.
#' @param rhovec vector of cutoff-specific values of rho. See \code{rdrobust()}
#'   for details.
#' @param covsvec vector of cutoff-specific covariates. See \code{rdrobust()}
#'   for details.
#' @param covs_dropvec vector indicating whether collinear covariates should be
#'   dropped at each cutoff. See \code{rdrobust()} for details.
#' @param kernelvec vector of cutoff-speficif kernels. See \code{rdrobust()} for
#'   details.
#' @param weightsvec vector of cutoff-speficif weights. See \code{rdrobust()}
#'   for details.
#' @param bwselectvec vector of cutoff-speficif bandwidth selection methods. See
#'   \code{rdrobust()} for details.
#' @param scaleparvec vector of cutoff-speficif scale parameters. See
#'   \code{rdrobust()} for details.
#' @param scaleregulvec vector of cutoff-speficif scale regularization
#'   parameters. See \code{rdrobust()} for details.
#' @param masspointsvec vector indicating how to handle repeated values at each
#'   cutoff. See \code{rdrobust()} for details.
#' @param bwcheckvec vector indicating the value of bwcheck at each cutoff. See
#'   \code{rdrobust()} for details.
#' @param bwrestrictvec vector indicating whether computed bandwidths are
#'   restricted to the range or runvar at each cutoff. See \code{rdrobust()} for
#'   details.
#' @param stdvarsvec vector indicating whether variables are standardized at
#'   each cutoff. See \code{rdrobust()} for details.
#' @param vcevec vector of cutoff-speficif variance-covariance estimation
#'   methods. See \code{rdrobust()} for details.
#' @param nnmatchvec vector of cutoff-speficif nearestneighbors for variance
#'   estimation. See \code{rdrobust()} for details.
#' @param cluster cluster ID variable. See \code{rdrobust()} for details.
#' @param level confidence level for confidence intervals. See \code{rdrobust()}
#'   for details.
#' @param plot plots cutoff-specific estimates and weights.
#'
#'
#' @return
#' \item{tau}{pooled estimate}
#' \item{se.rb}{robust bias corrected standard error for pooled estimate}
#' \item{pv.rb}{robust bias corrected p-value for pooled estimate}
#' \item{ci.rb.l}{left limit of robust bias corrected CI for pooled estimate}
#' \item{ci.rb.r}{right limit of robust bias corrected CI for pooled estimate}
#' \item{hl}{bandwidth to the left of the cutoff for pooled estimate}
#' \item{hr}{bandwidth to the right of the cutofffor pooled estimate}
#' \item{Nhl}{sample size within bandwidth to the left of the cutoff for pooled estimate}
#' \item{Nhr}{sample size within bandwidth to the right of the cutoff for pooled estimate}
#' \item{B}{vector of bias-corrected estimates}
#' \item{V}{vector of robust variances of the estimates}
#' \item{Coefs}{vector of conventional estimates}
#' \item{W}{vector of weights for each cutoff-specific estimate}
#' \item{Nh}{vector of sample sizes within bandwidth}
#' \item{CI}{robust bias-corrected confidence intervals}
#' \item{H}{matrix of bandwidths}
#' \item{Pv}{vector of robust p-values}
#' \item{rdrobust.results}{results from rdrobust for pooled estimate}
#' \item{cfail}{Cutoffs where rdrobust() encountered problems}
#'
#'
#' @examples
#' # Toy dataset
#' X <- runif(1000,0,100)
#' C <- c(rep(33,500),rep(66,500))
#' Y <- (1 + X + (X>=C))*(C==33)+(.5 + .5*X + .8*(X>=C))*(C==66) + rnorm(1000)
#' # rdmc with standard syntax
#' tmp <- rdmc(Y,X,C)
#'
#'
#' @export

rdmc <- function(Y,X,C,fuzzy=NULL,derivvec=NULL,pooled_opt=NULL,verbose=FALSE,
                pvec=NULL,qvec=NULL,hmat=NULL,bmat=NULL,rhovec=NULL,
                covsvec=NULL,covs_dropvec=NULL,kernelvec=NULL,weightsvec=NULL,
                bwselectvec=NULL,scaleparvec=NULL,scaleregulvec=NULL,
                masspointsvec=NULL,bwcheckvec=NULL,bwrestrictvec=NULL,
                stdvarsvec=NULL,vcevec=NULL,nnmatchvec=NULL,cluster=NULL,
                level=95,plot=FALSE){

  #################################################################
  # Setup and error checking
  #################################################################

  if (!is.numeric(C)){stop('C has to be numeric')}
  if (max(C,na.rm=TRUE)>=max(X,na.rm=TRUE) | min(C,na.rm=TRUE)<=min(X,na.rm=TRUE)){stop('cutoff variable outside range of running variable')}

  clist <- sort(unique(C))
  cnum <- length(clist)

  Xc <- X - C

  if (!is.null(hmat)){
    if (is.null(dim(hmat))){
      hmat <- matrix(hmat,nrow=cnum,ncol=2)
    }
  }
  if (!is.null(bmat)){
    if (is.null(dim(bmat))){
      bmat <- matrix(bmat,nrow=cnum,ncol=2)
    }
  }
  if (is.null(covs_dropvec)) covs_dropvec <- rep(TRUE,cnum)
  if (is.null(kernelvec)) kernelvec <- rep('tri',cnum)
  if (is.null(bwselectvec)) bwselectvec <- rep('mserd',cnum)
  if (is.null(vcevec)) vcevec <- rep('nn',cnum)
  if (is.null(nnmatchvec)) nnmatchvec <- rep(3,cnum)
  if (is.null(scaleparvec)) scaleparvec <- rep(1,cnum)
  if (is.null(scaleregulvec)) scaleregulvec <- rep(1,cnum)
  if (is.null(masspointsvec)) masspointsvec <- rep('adjust',cnum)
  if (is.null(bwrestrictvec)) bwrestrictvec <- rep(TRUE,cnum)
  if (is.null(stdvarsvec)) stdvarsvec <- rep(FALSE,cnum)

  B <- matrix(NA,nrow=1,ncol=cnum+2)
  V <- matrix(NA,nrow=1,ncol=cnum+2)
  Coefs <- matrix(NA,nrow=1,ncol=cnum+2)
  Nh <- matrix(NA,nrow=2,ncol=cnum+2)
  CI <- matrix(NA,nrow=2,ncol=cnum+2)
  Pv <- matrix(NA,nrow=1,ncol=cnum+2)
  H <- matrix(NA,nrow=2,ncol=cnum+2)
  W <- matrix(NA,nrow=1,ncol=cnum)
  Cfail <- numeric()

  #################################################################
  # Calculate pooled estimate
  #################################################################

  aux1 <- paste0('rdrobust::rdrobust(Y,Xc,,fuzzy=fuzzy,',pooled_opt,')')

  rdr <- eval(parse(text=aux1))

  hl <- rdr$bws[1,1]
  hr <- rdr$bws[1,2]
  Nhl <- rdr$N_h[1]
  Nhr <- rdr$N_h[2]

  B[1,cnum+2] <- rdr$Estimate[2]
  V[1,cnum+2] <- rdr$se[3]^2
  Coefs[1,cnum+2] <- rdr$Estimate[1]
  CI[,cnum+2] <- rdr$ci[3,]
  H[,cnum+2] <- rdr$bws[1,]
  Nh[,cnum+2] <- rdr$N_h
  Pv[1,cnum+2] <- rdr$pv[3]

  #################################################################
  # Calculate cutoff-specific estimates and weights
  #################################################################

  count <- 1
  count_fail <- 0
  for (c in clist){

    yc <- Y[abs(C-c)<=.Machine$double.eps]
    xc <- Xc[abs(C-c)<=.Machine$double.eps]

    rdr.tmp <- try(rdrobust::rdrobust(yc,xc,
                                      fuzzy=fuzzy,
                                      deriv=derivvec[count],
                                      p=pvec[count],
                                      q=qvec[count],
                                      h=hmat[count,],
                                      b=bmat[count,],
                                      rho=rhovec[count],
                                      covs=covsvec[count],
                                      covs_drop=covs_dropvec[count],
                                      kernel=kernelvec[count],
                                      weights=weightsvec[count],
                                      bwselect=bwselectvec[count],
                                      scalepar=scaleparvec[count],
                                      scaleregul=scaleregulvec[count],
                                      masspoints=masspointsvec[count],
                                      bwcheck=bwcheckvec[count],
                                      bwrestrict=bwrestrictvec[count],
                                      stdvars=stdvarsvec[count],
                                      vce=vcevec[count],
                                      nnmatch=nnmatchvec[count],
                                      cluster=cluster,
                                      level=level),
                   silent=TRUE)

    if (class(rdr.tmp)!="try-error"){
      B[1,count] <- rdr.tmp$Estimate[2]
      V[1,count] <- rdr.tmp$se[3]^2
      Coefs[1,count] <- rdr.tmp$Estimate[1]
      CI[,count] <- rdr.tmp$ci[3,]
      H[,count] <- rdr.tmp$bws[1,]
      Nh[,count] <- rdr.tmp$N_h
      Pv[1,count] <- rdr.tmp$pv[3]
    } else{
      Cfail <- c(Cfail,c)
      count_fail <- count_fail + 1
    }

    count <- count + 1

  }

  # Weights

  W[1,] <- colSums(Nh[,1:cnum],na.rm=TRUE)/sum(Nh[,1:cnum],na.rm=TRUE)
  W[is.na(W)] <- 0

  # Weighted estimate

  Baux <- B
  Vaux <- V
  Coefsaux <- Coefs
  Baux[is.na(Baux)] <- 0
  Vaux[is.na(Vaux)] <- 0
  Coefsaux[is.na(Coefsaux)] <- 0

#  B[1,cnum+1] <- B[1,1:cnum]%*%t(W)
#  V[1,cnum+1] <- V[1,1:cnum]%*%t(W^2)
#  Coefs[1,cnum+1] <- Coefs[1,1:cnum]%*%t(W)

  B[1,cnum+1] <- Baux[1,1:cnum]%*%t(W)
  V[1,cnum+1] <- Vaux[1,1:cnum]%*%t(W^2)
  Coefs[1,cnum+1] <- Coefsaux[1,1:cnum]%*%t(W)

#  Nh[,cnum+1] <- rowSums(Nh[,1:cnum])
  Nh[,cnum+1] <- rowSums(Nh[,1:cnum],na.rm=TRUE)

  CI[1,cnum+1] <- B[1,cnum+1]-sqrt(V[1,cnum+1])*qnorm(1-(1-level/100)/2)
  CI[2,cnum+1] <- B[1,cnum+1]+sqrt(V[1,cnum+1])*qnorm(1-(1-level/100)/2)
  Pv[1,cnum+1] <- 2*(1-pnorm(abs(B[1,cnum+1]/sqrt(V[1,cnum+1]))))

  #################################################################
  # Display results
  #################################################################

  if (verbose==TRUE){
    cat(summary(rdr))
    cat('\n')
  }

  cat('\n')
  cat(paste0(rep('=',80),collapse='')); cat('\n')
  cat(format('Cutoff',  width=11))
  cat(format('Coef.',   width=8))
  cat(format('P-value', width=16))
  cat(format('95% CI',  width=16))
  cat(format('hl',      width=9))
  cat(format('hr',      width=9))
  cat(format('Nh',      width=5))
  cat(format('Weight',  width=5))
  cat('\n')
  cat(paste0(rep('=',80),collapse=''))
  cat('\n')

  for (k in 1:cnum){
    cat(format(sprintf('%4.3f',clist[k]),        width=11))
    cat(format(sprintf('%7.3f',Coefs[k]),        width=9))
    cat(format(sprintf('%1.3f',Pv[k]),           width=9))
    cat(format(sprintf('%4.3f',CI[1,k]),         width=10))
    cat(format(sprintf('%4.3f',CI[2,k]),         width=10))
    cat(format(sprintf('%4.3f',H[1,k]),          width=9))
    cat(format(sprintf('%4.3f',H[2,k]),          width=9))
    cat(format(sprintf('%4.0f',Nh[1,k]+Nh[2,k]), width=8))
    cat(format(sprintf('%1.3f',W[k]),            width=5))
    cat('\n')
  }
  cat(paste0(rep('-',80),collapse='')); cat('\n')

  cat(format('Weighted',                                 width=11))
  cat(format(sprintf('%7.3f',Coefs[1,cnum+1]),           width=9))
  cat(format(sprintf('%1.3f',Pv[1,cnum+1]),              width=9))
  cat(format(sprintf('%4.3f',CI[1,cnum+1]),              width=10))
  cat(format(sprintf('%4.3f',CI[2,cnum+1]),              width=13))
  cat(format('  .',                                      width=9))
  cat(format('  .',                                      width=6))
  cat(format(sprintf('%4.0f',Nh[1,cnum+1]+Nh[2,cnum+1]), width=11))
  cat(format(' .',                                       width=5))
  cat('\n')

  cat(format('Pooled',                                   width=11))
  cat(format(sprintf('%7.3f',Coefs[cnum+2]),             width=9))
  cat(format(sprintf('%1.3f',Pv[cnum+2]),                width=9))
  cat(format(sprintf('%4.3f',CI[1,cnum+2]),              width=10))
  cat(format(sprintf('%4.3f',CI[2,cnum+2]),              width=10))
  cat(format(sprintf('%4.3f',H[1,cnum+2]),               width=9))
  cat(format(sprintf('%4.3f',H[2,cnum+2]),               width=9))
  cat(format(sprintf('%4.0f',Nh[1,cnum+2]+Nh[2,cnum+2]), width=11))
  cat(format(' .',                                       width=5))
  cat('\n')
  cat(paste0(rep('=',80),collapse='')); cat('\n')

  if (count_fail>0){
    warning("rdrobust() could not run in one or more cutoffs.")
  }

  #################################################################
  # Plots
  #################################################################

  if (plot==TRUE){

    ylim <- c(min(CI,na.rm=TRUE)*1.3,max(CI,na.rm=TRUE)*1.3)
    xlim <- c(min(clist),max(clist))

    par(mfrow=c(1,2))

    plot(NA,ylim=ylim,xlim=xlim,ylab='Treatment effect',xlab='Cutoff')

    polygon(x=c(min(clist),min(clist),max(clist),max(clist)),y=c(CI[1,1],CI[2,1],CI[2,1],CI[1,1]),
            border = NA,col='gray87')
    points(clist,Coefs[1:cnum],col='darkblue',pch=16)
    arrows(clist,CI[1,-1],clist,CI[2,-1],length=0.05,angle=90,code=3)
    abline(h=rdr$Estimate[1],col='gray34')
    abline(h=0,lty='dotted')

    legend('bottomright',legend=c('Estimates','Pooled estimate'),pch=c(16,NA),lty=c(NA,1),
           col=c('darkblue','gray34'),bty='n',cex=0.75)

    mtext('Bars are 95% CIs for estimates. \nShaded area is the 95% CI for pooled.',cex=0.8)

    barplot(W,xlab='Cutoff',ylab='Weight',names.arg=clist,space=1.5)

  }

  #################################################################
  # Return values
  #################################################################

  colnames(B) <- c(1:cnum,"weighted","pooled")
  colnames(V) <- c(1:cnum,"weighted","pooled")
  colnames(Coefs) <- c(1:cnum,"weighted","pooled")
  colnames(CI) <- c(1:cnum,"weighted","pooled")
  colnames(Nh) <- c(1:cnum,"weighted","pooled")
  colnames(H) <- c(1:cnum,"weighted","pooled")
  colnames(Pv) <- c(1:cnum,"weighted","pooled")

  rownames(Nh) <- c("left","right")
  rownames(H) <- c("left","right")

  output <- list(tau = rdr$Estimate[1],
                se.rb = rdr$se[3],
                pv.rb = rdr$pv[3],
                ci.rb.l = rdr$ci[3,1],
                ci.rb.r = rdr$ci[3,2],
                hl = rdr$bws[1,1],
                hr = rdr$bws[1,2],
                Nhl = rdr$N_h[1],
                Nhr = rdr$N_h[2],
                B = B,
                V = V,
                Coefs = Coefs,
                W = W,
                Nh = Nh,
                CI = CI,
                H = H,
                Pv = Pv,
                rdrobust.results = rdr,
                cfail = Cfail)

  return(output)

}
