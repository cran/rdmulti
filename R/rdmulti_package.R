###################################################################
# rdmulti: analysis of RD designs with multiple cutoffs or scores
# !version 0.4 07-Jan-2020
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################


#' rdmulti: analysis of RD Designs with multiple cutoffs or scores
#'
#' The regression discontinuity (RD) design is a popular quasi-experimental design
#' for causal inference and policy evaluation. The \code{'rdmulti'} package provides tools
#' to analyze RD designs with multiple cutoffs or scores: \code{\link{rdmc}()} estimates
#' pooled and cutoff-speficif effects in multi-cutoff designs, \code{\link{rdmcplot}()}
#' draws RD plots for multi-cutoff RD designs and \code{\link{rdms}()} estimates effects in
#' cumulative cutoffs or multi-score designs. For more details, and related \code{Stata} and
#' \code{R} packages useful for analysis of RD designs, visit \url{https://sites.google.com/site/rdpackages}.
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
#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @importFrom graphics legend
#' @importFrom graphics barplot
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom stats poly
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#'
#'
#' @aliases rdmulti_package
"_PACKAGE"
