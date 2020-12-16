#' @title ECD Sample Estimates
#' @name  ecd
#' @description Computes distance covariance and distance correlation type statistics, which are multivariate measures of dependence.
#' @param x data of first sample serving as the response
#' @param y data of second sample conditioned on
#' @param method ECD estimation method (options: slice, kernel.epa, kernel.gau)
#' @param ns number of slices. Defaults to NULL and must be specified when method is set to 'slice' for a continuous y
#' @param bw bandwidth of the kernel function, default
#'           bandwidths  are used unless otherwise specified
#' @param index exponent on Euclidean distance, in (0,2]. Defaults to 1
#' @return \code{ecd} returns a list with components
#' \describe{
#' \item{ecdCov}{sample covariance type statistic denoted \eqn{C_n(X|Y)}}
#' \item{ecdCor}{sample correlation type statistic}
#' \item{ecdVarX}{ECD variance of x sample}
#' }
#' @examples
#' library(MASS)
#' n=30
#' p=5
#' q=1
#' mu=rep(0,(p+q))
#' sigma=diag((p+q))
#' DATA=mvrnorm(n, mu, sigma)
#' X=DATA[,1:p]
#' Y=DATA[,(p+1):(p+q)]
#' # using slicing estimation method on a continuous Y
#' ecd(x = X, y = Y, method = "slice", ns = 5)
#' # using Gaussian kernel estimation method
#' ecd(x = X, y = Y, method = "kernel.gau")
#' # using Epanechnikov kernel estimation method
#' ecd(x = X, y = Y, method = "kernel.epa")
#'
#' #This example shows how to directly compute the ECD measures if Y is a
#' #categorical variable
#' library(boot)
#' #gravity data is in the boot package,there are 81 measurements in a series of
#' #eight experiments
#' #gravity$series shows the categories 1-8 of the observations. 999 permutations
#' #are used to compute p-value.
#'str(gravity)
#' ecd(gravity$g, gravity$series, method="slice")
#' @seealso \code{\link{ecdcor.test}} \code{\link{ecdcov.test}}
#' @importFrom stats dist sd
#' @export

ecd <- function(x,y,method,ns=NULL,bw="default",index=1.0) {

  # get the appropriate function to compute the ecd measure based on the value in 'method'

  third.arg <- NULL # Third argument of ECD methods. It is the number of slices for ecd.slice() and it is the kernel method for ecd.kernel()

  # choosing an estimation method
  method <- tolower(method) # make method case insensitive

  if(method=="slice") {
    # use slicing method
    if((is.null(ns)) & !(is.factor(y) | is.character(y)))
      stop("Using slicing method, \"ns\" cannot be NULL for a continuous Y")

    ecd.method <- ecd.slice
    third.arg <- ns # this time third.arg becomes the number of slices

  } else if ((method=="kernel.epa") | (method=="kernel.gau"))  {
    # use the kernel method
    ecd.method <- ecd.kernel
    third.arg <- method # this time third.arg becomes the kernel method
  } else {
    stop("Unknown method. Please check the help documentation.")
  }

  # ecd.method now becomes a generic function for both slicing and kernel estimation

  # compute the ecd measure
  estimates <- ecd.method(x, y, third.arg,bw, index)
  ecdcov <- estimates$ecdCov
  ecdcor <- estimates$ecdCor
  ecdvarX <- estimates$ecdVarX
  return(list(ecdCov=ecdcov,ecdCor=ecdcor, ecdVarX=ecdvarX))
}

