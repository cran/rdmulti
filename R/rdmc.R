###################################################################
# rdmc: analysis of RD designs with multiple cutoffs
# !version 0.1 11-Apr-2018
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Analysis of RD designs with multiple cutoffs
#'
#' \code{rdmc()} analyzes RD designs with multiple cutoffs.
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
#' @param pooled.opt options to be passed to \code{rdrobust()} to calculate pooled estimand.
#' @param hvec bandwidths to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param bvec bandwidths for the bias to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param pvec order of the polynomials to be passed to \code{rdrobust()} to calculate cutoff-specific estimates. Should be a vector of length equal to the number of different cutoffs.
#' @param kernelvec kernels to be passed to \code{rdrobust()} to calculate cutoff-specific estimates.Should be a vector of length equal to the number of different cutoffs.
#' @param fuzzy specifies a fuzzy design.
#' @param plot plots cutoff-specific estimates and weights.
#' @param verbose displays the output from \code{rdrobust} for estimating the pooled estimand.
#'
#'
#' @return
#' \item{tau}{pooled estimate}
#' \item{se.rb}{robust bias corrected standard error for pooled estimate}
#' \item{pv.rb}{robust bias corrected p-value for pooled estimate}
#' \item{ci.rb.l}{left limit of robust bias corrected CI for pooled estimate}
#' \item{ci.rb.r}{right limit of robust bias corrected CI for pooled estimate}
#' \item{hl}{bandwidth to the left of the cutoff for pooled estimate}
#' \item{hr}{bandwidth to the right of the cutoff for pooled estimate}
#' \item{Nhl}{sample size within bandwidth to the left of the cutoff for pooled estimate}
#' \item{Nhr}{sample size within bandwidth to the right of the cutoff for pooled estimate}
#' \item{B}{vector of bias-corrected coefficients}
#' \item{V}{variance-covariance matrix of the estimators}
#' \item{Coefs}{vector of conventional coefficients}
#' \item{W}{vector of weights for each cutoff-specific estimate}
#' \item{Nh}{vector of sample sizes within bandwidth at each cutoff}
#' \item{CI}{bias corrected confidence intervals}
#' \item{H}{bandwidth used at each cutoff}
#' \item{rdrobust.results}{results from rdrobust for pooled estimate}
#'
#'
#' @examples
#' # Toy dataset
#' X <- runif(1000,0,100)
#' C <- c(rep(33,500),rep(66,500))
#' Y <- (1 + X + (X>=C))*(C==33)+(.5 + .5*X + .8*(X>=C))*(C==66) + rnorm(1000)
#' # rdmc with standard syntax
#' tmp <- rdmc(Y,X,C)
#' # rdmc with cutoff-specific bandwidths
#' tmp <- rdmc(Y,X,C,hvec=c(9,13))
#'
#'
#' @export

rdmc = function(Y,X,C,pooled.opt=NULL,
                hvec=NULL,
                bvec=NULL,
                pvec=NULL,
                kernelvec=NULL,
                fuzzy=NULL,
                plot=FALSE,
                verbose=FALSE){

  #################################################################
  # Setup and error checking
  #################################################################

  if (!is.numeric(C)){stop('C has to be numeric')}
  if (max(C,na.rm=TRUE)>=max(X,na.rm=TRUE) | min(C,na.rm=TRUE)<=min(X,na.rm=TRUE)){stop('cutoff variable outside range of running variable')}

  clist = sort(unique(C))
  cnum = length(clist)

  Xc = X - C

  B = NULL
  V = matrix(0,nrow=cnum+1,ncol=cnum+1)
  Nh = NULL
  W = NULL
  Coefs = NULL
  CI = matrix(NA,nrow=2,ncol=cnum+1)
  Pv = NULL
  H = NULL

  #################################################################
  # Calculate pooled estimate
  #################################################################

  aux1 = paste0('rdrobust::rdrobust(Y,Xc,',pooled.opt,',fuzzy=fuzzy)')

  rdr = eval(parse(text=aux1))

  hl = rdr$h_l
  hr = rdr$h_r
  Nhl = rdr$N_h_l
  Nhr = rdr$N_h_r

  Nh = c(Nh,Nhl+Nhr)
  B = c(B,rdr$Estimate[2])
  V[1,1] = (rdr$se[3])^2
  Coefs = c(Coefs,rdr$Estimate[1])
  CI[,1] = c(rdr$ci[3,])
  Pv = c(Pv,rdr$pv[3])
  H = c(H,rdr$h_l)


  #################################################################
  # Calculate cutoff-specific estimates and weights
  #################################################################

  count = 1
  for (c in clist){

    n.aux = length(Y[C==c & Xc<=hr & Xc>=-hl])
    weight = n.aux / (Nhl + Nhr)

    W = c(W,weight)

    yc = Y[C==c]
    xc = Xc[C==c]

    h = hvec[count]
    b = bvec[count]
    p = pvec[count]
    if (!is.null(kernelvec)){
      kernel = kernelvec[count]
    } else{
      kernel = 'tri'
    }

    rdr.tmp = rdrobust::rdrobust(yc,xc,h=h,b=b,p=p,kernel=kernel,fuzzy=fuzzy)

    Nh = c(Nh,rdr.tmp$N_h_l+rdr.tmp$N_h_r)
    B = c(B,rdr.tmp$Estimate[2])
    V[count+1,count+1] = (rdr.tmp$se[3])^2
    Coefs = c(Coefs,rdr.tmp$Estimate[1])
    CI[,count+1] = c(rdr.tmp$ci[3,])
    Pv = c(Pv,rdr.tmp$pv[3])
    H = c(H,rdr.tmp$h_l)

    count = count + 1
  }


  #################################################################
  # Display results
  #################################################################

  if (verbose==TRUE){
    cat(summary(rdr))
    cat('\n')
  }

  cat('\n')
  cat(paste0(format('Cutoff',  width=11),
             format('Coef.',   width=11),
             format('P-value', width=19),
             format('95% CI',  width=20),
             format('h',       width=10),
             format('Nh',      width=10),
             format('Weight',  width=10))); cat('\n')


  for (k in 1:cnum){
    cat(paste0(format(toString(round(clist[k],3)),   width=11),
               format(toString(round(Coefs[k+1],3)), width=12),
               format(toString(round(Pv[k+1],3)),    width=12),
               format(toString(round(CI[1,k+1],3)),  width=12),
               format(toString(round(CI[2,k+1],3)),  width=12),
               format(toString(round(H[k+1],3)),     width=12),
               format(toString(Nh[k+1]),             width=11),
               format(toString(round(W[k],3)),       width=12)))
    cat('\n')
  }

  cat(paste0(format('Pooled', width=11),
             format(toString(round(Coefs[1],3)), width=12),
             format(toString(round(Pv[1],3)),    width=12),
             format(toString(round(CI[1,1],3)),  width=12),
             format(toString(round(CI[2,1],3)),  width=12),
             format(toString(round(H[1],3)),     width=12),
             format(toString(Nh[1]),             width=11),
             format(' .',                        width=12)))
  cat('\n')


  #################################################################
  # Plots
  #################################################################

  if (plot==TRUE){

    ylim = c(min(CI*1.3),max(CI)*1.3)
    xlim = c(min(clist),max(clist))

    par(mfrow=c(1,2))


    plot(NA,ylim=ylim,xlim=xlim,ylab='Treatment effect',xlab='Cutoff')

    polygon(x=c(min(clist),min(clist),max(clist),max(clist)),y=c(CI[1,1],CI[2,1],CI[2,1],CI[1,1]),
            border = NA,col='gray87')

    points(clist,Coefs[-1],col='darkblue',pch=16)
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

  output = list(tau = rdr$Estimate[1],
                se.rb = rdr$se[3],
                pv.rb = rdr$pv[3],
                ci.rb.l = rdr$ci[3,1],
                ci.rb.r = rdr$ci[3,2],
                hl = rdr$h_l,
                hr = rdr$h_r,
                Nhl = rdr$N_h_l,
                Nhr = rdr$N_h_r,
                B = B,
                V = V,
                Coefs = Coefs,
                W = W,
                Nh = Nh,
                CI = CI,
                H = H,
                rdrobust.results = rdr)

  return(output)

}
