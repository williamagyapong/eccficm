#' @title ECCFIC Sample Estimates
#' @name  eccfic
#' @description Computes Expected Conditional Characteristic Function-based
#' Independence Criterion (ECCFIC) statistics, which are new multivariate measures of
#' independence.
#' @param x Data of first sample serving as the response
#' @param y Data of second sample conditioned on
#' @param est The estimation method for the conditional mean. Set
#'   to 'slice' for slicing estimation or to 'kernel' for kernel.
#'   regression estimation. See 'Details' for additional information
#' @param kernel A kernel embedding to use for x, 'gaussian' or 'distance
#' @param index Exponent on distance in (0,2], if a 'distance' kernel is used
#'   for x. Defaults to 1
#' @param ns Number of slices. Defaults to NULL and must be specified when
#'   method is set to 'slice' for a continuous y
#' @param bw Bandwidth of the gaussian smoothing kernel applied on y, bandwidths
#'   suggested by Silverman (1986) are used unless otherwise specified.
#' @param sigma Bandwidth, if a 'gaussian' kernel is used for x, default is
#'   heuristic median pairwise distances of x.
#' @param wt Use weight function or not, if TRUE, f^2 is used by default
#' @return \code{eccfic} returns a list with components \describe{
#'   \item{eccfic}{sample ECCFIC estimate of \eqn{H_k^2(X|Y)}}
#'    \item{eccficVar}{sample ECCFIC variance estimate}
#'   \item{eccficCor}{sample correlation type statistic } }
#'
#' @details
#' When estimation method (est) is set to 'slice', the ECCFIC measure is computed by slicing
#' on Y in which case Y can originally be categorical or continuous. Estimating ECCFIC
#' is straightforward if Y is categorical.  For a continuous Y, we change Y into finite
#' categories with 'ns' levels. Slicing in multivariate and high-dimensional situations
#' is very challenging. So, for a continuous multivariate Y, we recommend using the kernel approach
#' or slicing Y before hand to be able to use the slicing approach. For help on slicing, please see
#' the works by \insertCite{zhu2010sufficient}{eccficm},
#' \insertCite{li2008projective}{eccficm}, and
#' \insertCite{cook2014fused}{eccficm}.
#'
#' Kernel regression is an alternative approach to use for a continuous Y. To use
#' the kernel regression method simply set \code{est} to 'kernel'.
#'
#' Refer to section 4 of \insertCite{yin2019new}{eccficm} for more information on
#' the estimation approaches.
#'
#' Though ECCFIC is not very sensitive to the number of slices, it is suggested that
#'  each slice should have at least 5 and at most n/2 data points.
#'
#' ECD \insertCite{yin2019new}{eccficm}, \code{\link{ecd}}, belongs to the family of generalized
#' ECCFIC measures. In fact, \code{ecdCov^2} = 2*\code{eccfic}, when both measures
#' are estimated by the slicing method and a distance kernel is used for the ECCFIC.
#'
#' @note See \code{\link{eccfic.test}} for a test of multivariate independence based on the
#' ECCFIC statistic as well as the correlation type statistic.
#'
#' @author Chenlu Ke (\email{chenlu.ke@@uky.edu}) and Xiangrong Yin
#'
#' @references \insertRef{ke2019expected}{eccficm}
#' \url{https://www.tandfonline.com/doi/abs/10.1080/01621459.2019.1604364}
#'
#'\insertRef{zhu2010sufficient}{eccficm}
#'
#'\insertRef{li2008projective}{eccficm}
#'
#'\insertRef{cook2014fused}{eccficm}
#'
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100)
#'
#' # Using kernel regression estimation method
#' eccfic(x,y, est="kernel")
#'
#' # Using slicing estimation method for a continuous Y
#' eccfic(x,y, est="slice", ns=5)
#'
#' # Using slicing method for a categorical Y
#' y <- rep(1:4, 25)
#' x <- y + rnorm(100)
#' # don't need to provide 'ns' because Y is already categorical
#' eccfic(x,as.factor(y),est='slice',kernel='distance')
#'
#' @seealso \code{\link{ecd}} \code{\link{ecd.test}}
#'
#' @importFrom stats dist sd cov dnorm median
#' @export


eccfic <- function(x,y, est="kernel",ns=NULL, kernel="gaussian", bw="default", sigma="default",index=1.0, wt=F) {

  est <- tolower(est) # make method case insensitive
  kernel <- tolower(kernel) # make case insensitive

  if(est=="slice") {
    # use slicing estimation method
    if((is.null(ns)) & !(is.factor(y) | is.character(y)))
      stop("Using slicing method, \"ns\" cannot be NULL for a continuous Y")

      # compute the eccfic measures
      HXY <- eccfic.slice(x, y, ns, kernel, sigma, index) # estimated H^2(X|Y)
      HXX <- eccfic.slice(x, x, ns, kernel, sigma, index) # estimated H^2(X|X)


  } else if (est=="kernel")  {
    # use the kernel regression estimation method
    # compute the eccfic measures
    HXY <- eccfic.kernreg(x, y, kernel, sigma, bw, index, wt) # estimated H^2(X|Y)
    HXX <- eccfic.kernreg(x, x, kernel, sigma, bw, index, wt) # estimated H^2(X|X)
  } else {
    stop("Unknown estimation method. Available options, 'slice' or 'kernel'")
  }

  rhok <- HXY/HXX # Correlation type statistic

  return(list(eccfic=HXY,eccficCor=rhok))
}

