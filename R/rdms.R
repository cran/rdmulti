###################################################################
# rdms: analysis of RD designs with multiple scores
# !version 0.3 21-Oct-2019
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Analysis of RD designs with cumulative cutoffs or two running variables
#'
#' \code{rdms()} analyzes RD designs with cumulative cutoffs or two running variables.
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
#' @param C vector of cutoffs.
#' @param X2 if specified, second running variable.
#' @param zvar if X2 is specified, treatment indicator.
#' @param C2 if specified, second vector of cutoffs.
#' @param rangemat matrix of cutoff-specific ranges for the running variable.
#' @param xnorm normalized running variable to estimate pooled effect.
#' @param pooled.opt options to be passed to \code{rdrobust()} to calculate pooled estimand.
#' @param derivvec vector of cutoff-specific order of derivatives. See \code{rdrobust()} for details.
#' @param pvec vector of cutoff-specific polynomial orders. See \code{rdrobust()} for details.
#' @param qvec vector of cutoff-specific polynomial orders for bias estimation. See \code{rdrobust()} for details.
#' @param hmat matrix of cutoff-specific bandwidths. See \code{rdrobust()} for details.
#' @param bmat matrix of cutoff-specific bandwidths for bias estimation. See \code{rdrobust()} for details.
#' @param rhovec vector of cutoff-specific values of rho. See \code{rdrobust()} for details.
#' @param covsvec vector of cutoff-specific covariates. See \code{rdrobust()} for details.
#' @param kernelvec vector of cutoff-speficif kernels. See \code{rdrobust()} for details.
#' @param weightsvec vector of cutoff-speficif weights. See \code{rdrobust()} for details.
#' @param bwselectvec vector of cutoff-speficif bandwidth selection methods. See \code{rdrobust()} for details.
#' @param vcevec vector of cutoff-speficif variance-covariance estimation methods. See \code{rdrobust()} for details.
#' @param cluster cluster ID variable. See \code{rdrobust()} for details.
#' @param nnmatchvec vector of cutoff-speficif nearestneighbors for variance estimation. See \code{rdrobust()} for details.
#' @param scaleparvec vector of cutoff-speficif scale parameters. See \code{rdrobust()} for details.
#' @param scaleregulvec vector of cutoff-speficif scale regularization parameters. See \code{rdrobust()} for details.
#' @param fuzzy specifies a fuzzy design. See \code{rdrobust()} for details.
#' @param level confidence level for confidence intervals. See \code{rdrobust()} for details.
#' @param plot plots cutoff-specific and pooled estimates.
#'
#'
#' @return
#' \item{B}{vector of bias-corrected coefficients}
#' \item{V}{variance-covariance matrix of the estimators}
#' \item{Coefs}{vector of conventional coefficients}
#' \item{Nh}{vector of sample sizes within bandwidth at each cutoff}
#' \item{CI}{bias corrected confidence intervals}
#' \item{H}{bandwidth used at each cutoff}
#' \item{Pv}{vector of robust p-values}
#'
#'
#' @examples
#' # Toy dataset: cumulative cutoffs
#' X <- runif(1000,0,100)
#' C <- c(33,66)
#' Y <- (1+X)*(X<C[1])+(0.8+0.8*X)*(X>=C[1]&X<C[2])+(1.2+1.2*X)*(X>=C[2]) + rnorm(1000)
#' # rmds: basic syntax
#' tmp <- rdms(Y,X,C)
#'
#'
#' @export

rdms = function(Y,X,C,X2=NULL,zvar=NULL,C2=NULL,
                #range.l=NULL,
                #range.r=NULL,
                rangemat=NULL,
                xnorm=NULL,
                pooled.opt=NULL,
                derivvec=NULL,
                pvec=NULL,
                qvec=NULL,
                hmat=NULL,
                bmat=NULL,
                rhovec=NULL,
                covsvec=NULL,
                kernelvec=NULL,
                weightsvec=NULL,
                bwselectvec=NULL,
                vcevec=NULL,
                cluster=NULL,
                nnmatchvec=NULL,
                scaleparvec=NULL,
                scaleregulvec=NULL,
                level=95,
                fuzzy=NULL,
                plot=FALSE){

  #################################################################
  # Setup and error checking
  #################################################################

  if (!is.null(X2) & is.null(zvar)){stop('Need to specify zvar when X2 is specified')}
  if (!is.null(X2) & is.null(C2)){stop('Need to specify C2 if X2 is specified')}

  cnum = length(C)

  if (!is.null(C2)){
    if (cnum!=length(C2)){stop('cutoff coordinates incorrectly specified')}
  }

  if (!is.null(rangemat)){
    if (is.null(dim(rangemat))){
      rangemat = matrix(rangemat,nrow=cnum,ncol=2)
    }
  }
  else {
    rangemat = cbind(rep(-Inf,cnum),rep(Inf,cnum))
  }

  if (!is.null(hmat)){
    if (is.null(dim(hmat))){
      hmat = matrix(hmat,nrow=cnum,ncol=2)
    }
  }

  if (is.null(kernelvec)){kernelvec = rep('tri',cnum)}
  if (is.null(bwselectvec)){bwselectvec = rep('mserd',cnum)}
  if (is.null(vcevec)){vcevec = rep('nn',cnum)}
  if (is.null(nnmatchvec)){nnmatchvec = rep(3,cnum)}
  if (is.null(scaleparvec)){scaleparvec = rep(1,cnum)}
  if (is.null(scaleregulvec)){scaleregulvec = rep(1,cnum)}

  B = matrix(NA,nrow=1,ncol=cnum+1)
  V = matrix(NA,nrow=1,ncol=cnum+1)
  Coefs = matrix(NA,nrow=1,ncol=cnum+1)
  Nh = matrix(NA,nrow=2,ncol=cnum+1)
  CI = matrix(NA,nrow=2,ncol=cnum+1)
  Pv = matrix(NA,nrow=1,ncol=cnum+1)
  H = matrix(NA,nrow=2,ncol=cnum+1)

  c.disp = NULL

  #################################################################
  # Calculate cutoff-specific estimates
  #################################################################

  if (is.null(X2)){

    for (c in 1:cnum){

      xc = X - C[c]
      Rc = rangemat - C

      yc = Y[xc>=Rc[c,1] & xc<=Rc[c,2]]
      xc = xc[xc>=Rc[c,1] & xc<=Rc[c,2]]

      rdr.tmp = rdrobust::rdrobust(yc,xc,deriv=derivvec[c],
                                   p=pvec[c],
                                   q=qvec[c],
                                   h=hmat[c,],
                                   b=bmat[c,],
                                   rho=rhovec[c],
                                   covs=covsvec[c],
                                   kernel=kernelvec[c],
                                   weights=weightsvec[c],
                                   bwselect=bwselectvec[c],
                                   vce=vcevec[c],
                                   cluster=cluster,
                                   nnmatch=nnmatchvec[c],
                                   scalepar=scaleparvec[c],
                                   scaleregul=scaleregulvec[c],
                                   level=level,
                                   fuzzy=fuzzy)

      B[1,c] = rdr.tmp$Estimate[2]
      V[1,c] = rdr.tmp$se[3]^2
      Coefs[1,c] = rdr.tmp$Estimate[1]
      CI[,c] = rdr.tmp$ci[3,]
      H[,c] = rdr.tmp$bws[1,]
      Nh[,c] = rdr.tmp$Nh
      Pv[1,c] = rdr.tmp$pv[3]

      c.disp = c(c.disp,round(C[c],2))

    }

  } else {

    for (c in 1:cnum){

      xc = sqrt((X-C[c])^2+(X2-C2[c])^2)*(2*zvar-1)

      yc = Y[xc>=rangemat[c,1] & xc<=rangemat[c,2]]
      xc = xc[xc>=rangemat[c,1] & xc<=rangemat[c,2]]

      rdr.tmp = rdrobust::rdrobust(yc,xc,
                                   p=pvec[c],
                                   q=qvec[c],
                                   h=hmat[c,],
                                   b=bmat[c,],
                                   rho=rhovec[c],
                                   covs=covsvec[c],
                                   kernel=kernelvec[c],
                                   weights=weightsvec[c],
                                   bwselect=bwselectvec[c],
                                   vce=vcevec[c],
                                   cluster=cluster,
                                   nnmatch=nnmatchvec[c],
                                   scalepar=scaleparvec[c],
                                   scaleregul=scaleregulvec[c],
                                   level=level,
                                   fuzzy=fuzzy)

      B[1,c] = rdr.tmp$Estimate[2]
      V[1,c] = rdr.tmp$se[3]^2
      Coefs[1,c] = rdr.tmp$Estimate[1]
      CI[,c] = rdr.tmp$ci[3,]
      H[,c] = rdr.tmp$bws[1,]
      Nh[,c] = rdr.tmp$Nh
      Pv[1,c] = rdr.tmp$pv[3]

      c.disp = c(c.disp,paste0('(',round(C[c],2),',',round(C2[c],2),')'))

    }

  }


  #################################################################
  # Calculate pooled estimates
  #################################################################

  if (!is.null(xnorm)){

    aux1 = paste0('rdrobust::rdrobust(Y,xnorm,fuzzy=fuzzy,',pooled.opt,')')

    rdr = eval(parse(text=aux1))

    B[1,cnum+1] = rdr$Estimate[2]
    V[1,cnum+1] = rdr$se[3]^2
    Coefs[1,cnum+1] = rdr$Estimate[1]
    CI[,cnum+1] = rdr$ci[3,]
    H[,cnum+1] = rdr$bws[1,]
    Nh[,cnum+1] = rdr$Nh
    Pv[1,cnum+1] = rdr$pv[3]

  }


  #################################################################
  # Display results
  #################################################################

  cat('\n')
  cat(paste0(format('Cutoff',  width=14),
             format('Coef.',   width=8),
             format('P-value', width=15),
             format('95% CI',  width=16),
             format('hl',      width=9),
             format('hr',      width=9),
             format('Nh',      width=10))); cat('\n')


  for (k in 1:cnum){
    cat(paste0(format(toString(c.disp[k]),         width=12),
               format(toString(round(Coefs[k],3)), width=12),
               format(toString(round(Pv[k],3)),    width=8),
               format(toString(round(CI[1,k],3)),  width=10),
               format(toString(round(CI[2,k],3)),  width=10),
               format(toString(round(H[1,k],3)),   width=9),
               format(toString(round(H[2,k],3)),   width=9),
               format(toString(Nh[1,k]+Nh[2,k]),   width=10)))
    cat('\n')
  }

  if (!is.null(xnorm)){
    cat(paste0(format('Pooled',                            width=12),
               format(toString(round(Coefs[cnum+1],3)),    width=12),
               format(toString(round(Pv[cnum+1],3)),       width=8),
               format(toString(round(CI[1,cnum+1],3)),     width=10),
               format(toString(round(CI[2,cnum+1],3)),     width=10),
               format(toString(round(H[1,cnum+1],3)),      width=9),
               format(toString(round(H[2,cnum+1],3)),      width=9),
               format(toString(Nh[1,cnum+1]+Nh[2,cnum+1]), width=10)))
  }
  cat('\n')

  #################################################################
  # Plots
  #################################################################


  #################################################################
  # Return values
  #################################################################

  colnames(B) = c(1:cnum,"pooled")
  colnames(V) = c(1:cnum,"pooled")
  colnames(Coefs) = c(1:cnum,"pooled")
  colnames(CI) = c(1:cnum,"pooled")
  colnames(Nh) = c(1:cnum,"pooled")
  colnames(H) = c(1:cnum,"pooled")
  colnames(Pv) = c(1:cnum,"pooled")

  rownames(Nh) = c("left","right")
  rownames(H) = c("left","right")

  output = list(B = B,
                V = V,
                Coefs = Coefs,
                Nh = Nh,
                CI = CI,
                H = H,
                Pv = Pv)

  return(output)

}
