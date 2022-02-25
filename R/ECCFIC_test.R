#' @title ECCFIC Independence Test
#' @name eccfic.test
#' @description Performs a permutation test of independence between
#' two random vectors based on the ECCFIC statistic or the correlation type statistic.
#' @param B  number of permutation replicates
#' @inheritParams eccfic
#' @param ts the test statistic, 'eccfic' or 'eccficCor'
#' @details
#' For the permutation test, if no value is given for 'B', (200+500/n) replicates
#' are used by default, where n is the sample size. At least 99 and at most 999 random permutations
#' should be sufficient practically.
#'
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
#' @author Chenlu Ke (\email{chenlu.ke@@uky.edu}) and Xiangrong Yin
#' @references \insertRef{ke2019expected}{eccficm}
#' \url{https://www.tandfonline.com/doi/abs/10.1080/01621459.2019.1604364}
#'
#' @return \code{eccfic.test} returns a list containing the following components:
#' \describe{
#'  \item{estimation.method}{a character string indicating the method used to compute the ECCFIC statistic}
#'  \item{method}{a character string describing the type of ECCFIC test performed}
#'  \item{statistic}{the value of the ECCFIC test statistic}
#'  \item{p.value}{approximate p-value of the test based on a permutation test}
#'  \item{alternative }{a character string describing the alternative hypothesis}
#'  \item{estimates}{vector containing estimates of H2(X|Y), H2(X|X), and Rho(X|Y)}
#'  \item{data.name}{a character string containing the names of the data used}
#' }
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100)
#'
#' # Using Gaussian Kernel Regression estimation method
#' eccfic.test(x,y) # based on default permutations
#' eccfic.test(x,y,B=199) # based on 199 permutations
#'
#' # Using slicing estimation method on a continuous Y
#' eccfic.test(x,y,est='slice',ns=5,kernel='distance') # based on default permutations
#' eccfic.test(x,y,est='slice',ns=5,kernel='distance',B=199) # based on 199 permutations
#'
#' # Using slicing method for a categorical Y
#' y <- rep(1:4, 25)
#' x <- y + rnorm(100)
#' # don't need to provide 'ns' because Y is already categorical
#' eccfic.test(x,as.factor(y),est='slice',kernel='distance')
#'
#'@seealso \code{\link{eccfic}} \code{\link{ecd.test}}
#' @export

eccfic.test <-
function(x, y, est="kernel", ns=NULL,kernel="gaussian", B="default", ts="eccficCor",bw="default",sigma="default", index=1.0, wt=F){
  # x <- X
  # y <- Y
  est <- tolower(est); kernel <- tolower(kernel);  ts <- tolower(ts)
  if(ts=='eccficcor') ts <- 'eccficCor'
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(y)

  if(B=="default") {
    # Set default number of permutations
    B <- floor((200+5000/n)) # based on a method suggested by Szekely, et. al. (2007)
    if(B < 99) warning("Insufficient number of permutations! Consider changing \"B\"
                       from default to a value greater than 99 but not exceeding 999")
  }

  # Set index of the test statistic to be used (1:eccficCov, 2:eccficCor)
  tstats.index <- ifelse(ts=='eccficCor',2, ifelse(ts=='eccfic', 1, 2))

  # compute sample statistics
  estimates <- eccfic(x,y,est,ns=ns,kernel=kernel, bw=bw, sigma=sigma,index=index, wt=wt)
  eccfic <- estimates$eccfic
  eccficCor <- estimates$eccficCor
  test.statistic <- estimates[[tstats.index]]

  # Compute p-value via permutation test
  reps <- replicate(B,eccfic(x,y[sample(1:n),],est,ns,kernel,bw,sigma,index, wt)[[tstats.index]])
  pvalue <- (sum(reps>=test.statistic)+1)/(B+1)

  estimates <- c(eccfic, eccficCor)
  names(estimates) <- c("eccfic", "eccficCor")
  names(test.statistic) <- ts

  est.method <- paste("Estimation Method: ", ifelse(est=='slice', 'Slicing on Y', 'Kernel Regression'))
  test.method <- paste(ts," permutation test of independence - ", est.method)
  alternative <- "X is not independent of Y"
  dataname <- ifelse((xname=="X")&(yname=="Y"), paste("X, Y, Replicates = ", B),
                     paste("X -",xname, ", Y -", yname, ", Replicates = ", B))

  # Put together all output elements
  output <- list(
    method = test.method,
    estimation.method = est.method,
    estimates = estimates,
    statistic = test.statistic,
    alternative = alternative,
    p.value = pvalue,
    data.name = dataname
  )

  class(output) <- 'htest'
  return(output)
}



