###################################################################
# rdms: analysis of RD designs with multiple scores
# !version 0.2 12-Jul-2018
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Analysis of RD designs with cumulative cutoffs or two running variables
#'
#' \code{rdms()} analyzes RD designs with cumulative cutoffs or two running variables.
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
#' @param C vector of cutoffs.
#' @param X2 if specified, second running variable.
#' @param zvar if X2 is specified, treatment indicator.
#' @param C2 if specified, second vector of cutoffs.
#' @param range.l range of the running variable to be used for estimation around the cutoff from the left.
#' @param range.r range of the running variable to be used for estimation around the cutoff from the right.
#' @param xnorm normalized running variable to estimate pooled effect.
#' @param pooled.opt options to be passed to \code{rdrobust()} to calculate pooled estimand.
#' @param hvec bandwidths to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param bvec bandwidths for the bias to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param pvec order of the polynomials to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param kernelvec kernels to be passed to \code{rdrobust()} to calculate cutoff-specific estimates.Should be a vector of length equal to the number of different cutoffs.
#' @param fuzzy specifies a fuzzy design.
#' @param plot plots cutoff-specific estimates and weights.
#'
#'
#' @return
#' \item{B}{vector of bias-corrected coefficients}
#' \item{V}{variance-covariance matrix of the estimators}
#' \item{Coefs}{vector of conventional coefficients}
#' \item{Nh}{vector of sample sizes within bandwidth at each cutoff}
#' \item{CI}{bias corrected confidence intervals}
#' \item{H}{bandwidth used at each cutoff}
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
                range.l=NULL,
                range.r=NULL,
                xnorm=NULL,
                pooled.opt=NULL,
                hvec=NULL,
                bvec=NULL,
                pvec=NULL,
                kernelvec=NULL,
                fuzzy=NULL,
                plot=FALSE){

  #################################################################
  # Setup and error checking
  #################################################################

  if (!is.null(X2) & is.null(zvar)){stop('Need to specify zvar when X2 is specified')}
  if (!is.null(X2) & is.null(C2)){stop('Need to specify C2 if X2 is specified')}

  ncutoffs = length(C)

  if (!is.null(C2)){
    if (ncutoffs!=length(C2)){stop('cutoff coordinates incorrectly specified')}
  }

  if (!is.null(range.l)){
    if (is.null(range.r)){range.r = range.l}
  }

  B = NULL
  V = matrix(0,nrow=ncutoffs+1,ncol=ncutoffs + !is.null(xnorm))
  Nh = NULL
  Coefs = NULL
  CI = matrix(NA,nrow=2,ncol=ncutoffs + !is.null(xnorm))
  Pv = NULL
  H = NULL

  c.disp = NULL


  #################################################################
  # Calculate cutoff-specific estimates
  #################################################################

  if (is.null(X2)){

    for (c in 1:ncutoffs){

      xc = X - C[c]

      if (is.null(range.l)){
        rc = Inf
        rt = Inf
      } else {
        rc = range.l[c]
        rt = range.r[c]
      }

      h = hvec[c]
      b = bvec[c]
      p = pvec[c]
      if (!is.null(kernelvec)){
        kernel = kernelvec[c]
      } else{
        kernel = 'tri'
      }

      yc = Y[xc>=-rc & xc<=rt]
      xc = xc[xc>=-rc & xc<=rt]

      rdr.tmp = rdrobust::rdrobust(yc,xc,h=h,b=b,p=p,kernel=kernel,fuzzy=fuzzy)

      Nh = c(Nh,rdr.tmp$Nh[1]+rdr.tmp$Nh[2])
      B = c(B,rdr.tmp$Estimate[2])
      V[c,c] = (rdr.tmp$se[3])^2
      Coefs = c(Coefs,rdr.tmp$Estimate[1])
      CI[,c] = c(rdr.tmp$ci[3,])
      Pv = c(Pv,rdr.tmp$pv[3])
      H = c(H,rdr.tmp$bws[1,1])

      c.disp = c(c.disp,round(C[c],2))

    }

  } else {

    for (c in 1:ncutoffs){

      xc = sqrt((X-C[c])^2+(X2-C2[c])^2)*(2*zvar-1)

      if (is.null(range.l)){
        rc = Inf
        rt = Inf
      } else {
        rc = range.l[c]
        rt = range.r[c]
      }

      h = hvec[c]
      b = bvec[c]
      p = pvec[c]
      if (!is.null(kernelvec)){
        kernel = kernelvec[c]
      } else{
        kernel = 'tri'
      }

      yc = Y[xc>=-rc & xc<=rt]
      xc = xc[xc>=-rc & xc<=rt]

      rdr.tmp = rdrobust::rdrobust(yc,xc,h=h,b=b,p=p,kernel=kernel,fuzzy=fuzzy)

      Nh = c(Nh,rdr.tmp$Nh[1]+rdr.tmp$Nh[2])
      B = c(B,rdr.tmp$Estimate[2])
      V[c,c] = (rdr.tmp$se[3])^2
      Coefs = c(Coefs,rdr.tmp$Estimate[1])
      CI[,c] = c(rdr.tmp$ci[3,])
      Pv = c(Pv,rdr.tmp$pv[3])
      H = c(H,rdr.tmp$bws[1,1])

      c.disp = c(c.disp,paste0('(',round(C[c],2),',',round(C2[c],2),')'))

    }

  }


  #################################################################
  # Calculate pooled estimates
  #################################################################

  if (!is.null(xnorm)){

    aux1 = paste0('rdrobust::rdrobust(Y,xnorm,',pooled.opt,'fuzzy=fuzzy)')

    rdr = eval(parse(text=aux1))

    Nh = c(Nh,rdr$Nh[1]+rdr$Nh[2])
    B = c(B,rdr$Estimate[2])
    V[ncutoffs+1,ncutoffs+1] = (rdr$se[3])^2
    Coefs = c(Coefs,rdr$Estimate[1])
    CI[,ncutoffs+1] = c(rdr$ci[3,])
    Pv = c(Pv,rdr$pv[3])
    H = c(H,rdr$bws[1,1])

  }


  #################################################################
  # Display results
  #################################################################

  cat('\n')
  cat(paste0(format('Cutoff',  width=15),
             format('Coef.',   width=11),
             format('P-value', width=19),
             format('95% CI',  width=20),
             format('h',       width=10),
             format('Nh',      width=10))); cat('\n')


  for (k in 1:ncutoffs){
    cat(paste0(format(toString(c.disp[k]),         width=15),
               format(toString(round(Coefs[k],3)), width=12),
               format(toString(round(Pv[k],3)),    width=12),
               format(toString(round(CI[1,k],3)),  width=12),
               format(toString(round(CI[2,k],3)),  width=12),
               format(toString(round(H[k],3)),     width=12),
               format(toString(Nh[k]),             width=11)))
    cat('\n')
  }

  if (!is.null(xnorm)){
    cat(paste0(format('Pooled', width=15),
               format(toString(round(Coefs[ncutoffs+1],3)),width=12),
               format(toString(round(Pv[ncutoffs+1],3)),   width=12),
               format(toString(round(CI[1,ncutoffs+1],3)), width=12),
               format(toString(round(CI[2,ncutoffs+1],3)), width=12),
               format(toString(round(H[ncutoffs+1],3)),    width=12),
               format(toString(Nh[ncutoffs+1]),            width=11)))
  }
  cat('\n')

  #################################################################
  # Plots
  #################################################################


  #################################################################
  # Return values
  #################################################################

  output = list(B = B,
                V = V,
                Coefs = Coefs,
                Nh = Nh,
                CI = CI,
                H = H)

  return(output)

}
