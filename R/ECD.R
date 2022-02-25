#'@title ECD Sample Estimates
#'@name  ecd
#'@description Computes Expectation of the Conditional Difference (ECD) statistics.
#' These covariance and correlation type statistics are multivariate measures of
#' independence.
#'@param x Data of first sample serving as the response
#'@param y Data of second sample conditioned on
#'@param est The estimation method (options: 'slice', 'kernel.epa', 'kernel.gau').
#'See 'Details' for additional information
#'@param ns Number of slices. Defaults to NULL and must be specified when method
#'  is set to 'slice' for a continuous y
#'@param bw Bandwidth of the kernel function, default bandwidths  are used
#'  unless otherwise specified
#'@param index Exponent on Euclidean distance in (0,2], with default value of 1
#'
#'@details
#' When estimation method (est) is set to 'slice', the ECD measure is computed by slicing
#' on Y in which case Y can originally be categorical or continuous. Estimating ECD
#' is straightforward if Y is categorical. For a continuous Y, we first change Y into
#' finite categories with 'ns' levels. Slicing in multivariate and high-dimensional situations
#' is very challenging. So, for a continuous multivariate Y, we recommend using the kernel approach
#' or slicing Y before hand to be able to use the slicing approach. For help on slicing,
#'  please see the works by
#' \insertCite{zhu2010sufficient}{eccficm}, \insertCite{li2008projective}{eccficm},
#' and \insertCite{cook2014fused}{eccficm}.
#'
#' Kernel estimation is an alternative approach to use for a continuous Y. The
#' available method options, 'kernel.epa' and 'kernel.gau', make use of the
#' Epanechnikov kernel function and Gaussian kernel function, respectively.
#'
#' Refer to section 4 of \insertCite{yin2019new}{eccficm} for more information on
#' the estimation approaches.
#'
#' ECD \insertCite{yin2019new}{eccficm}, \code{\link{ecd}}, belongs to the family of generalized
#' ECCFIC measures. In fact, \code{ecdCov^2} = 2*\code{eccfic}, when both measures
#' are estimated by the slicing method and a distance kernel is used for the ECCFIC.
#'
#'@author Qinqcong Yuan (\email{yuanq3@@miamioh.edu}) and Xiangrong Yin
#'@references
#'\insertRef{yin2019new}{eccficm}
#'  \url{http://www3.stat.sinica.edu.tw/ss_newpaper/SS-2017-0538_na.pdf}
#'
#'\insertRef{zhu2010sufficient}{eccficm}
#'
#'\insertRef{li2008projective}{eccficm}
#'
#'\insertRef{cook2014fused}{eccficm}
#'
#'@return \code{ecd} returns a list with components
#'\describe{
#'  \item{ecdCov}{sample ECD covariance type statistic denoted \eqn{C_n(X|Y)}}
#'  \item{ecdCor}{sample ECD correlation type statistic}
#'  \item{ecdVarX}{sample ECD variance type statistic of X}
#'  }
#' @note See \code{\link{ecd.test}} for a test of multivariate independence based on the
#' covariance and correlation type statistics.
#' @examples
#' library(MASS)
#' n <- 30; p <- 5; q <- 1
#' mu <- rep(0,(p+q))
#' sigma <- diag((p+q))
#' DATA <- mvrnorm(n, mu, sigma)
#' X <- DATA[,1:p]
#' Y <- DATA[,(p+1):(p+q)]
#'
#' # using slicing estimation method on a continuous Y
#' ecd(x = X, y = Y, est = "slice", ns = 5)
#'
#' # using Gaussian kernel estimation method
#' ecd(x = X, y = Y, est = "kernel.gau")
#'
#' # using Epanechnikov kernel estimation method
#' ecd(x = X, y = Y, est = "kernel.epa")
#'
#' #This example shows how to directly compute the ECD measures if Y is a
#' #categorical variable
#' library(boot)
#' #gravity data is in the boot package,there are 81 measurements in a series of
#' #eight experiments
#' #gravity$series shows the categories 1-8 of the observations.
#'str(gravity)
#' ecd(gravity$g, gravity$series, est="slice")
#'
#'@seealso \code{\link{eccfic}}  \code{\link{eccfic.test}}
#'@importFrom stats dist sd
#'@export

ecd <- function(x,y,est="slice",ns=5,bw="default",index=1.0) {

  est <- tolower(est) # make method case insensitive

  # choosing an estimation method
  if(est=="slice") {

    if((is.null(ns)) & !(is.factor(y) | is.character(y)))
      stop("Using slicing method, \"ns\" cannot be NULL for a continuous Y")
    # use slicing estimation method to compute the ecd measure
    estimates <- ecd.slice(x, y, ns, index)

  } else if ((est=="kernel.epa") | (est=="kernel.gau"))  {
    # use the kernel estimation method to compute the ecd measure
    estimates <- ecd.kernel(x, y, est, bw, index)
  } else
    stop("Unknown method. Available options are 'slice', 'kernel.epa', and 'kernel.gau'")

  ecdcov <- estimates$ecdCov # covariance type statistic
  ecdcor <- estimates$ecdCor # correlation type statistic
  ecdvarX <- estimates$ecdVarX # variance type statistic
  return(list(ecdCov=ecdcov,ecdCor=ecdcor, ecdVarX=ecdvarX))
}

