#' @title ECD Independence Test
#' @name ecd.test
#' @description Performs a permutation test of independence between
#' two random vectors based on the ECD covariance and correlation type statistics.
#' @param B  number of permutation replicates
#' @inheritParams ecd
#' @param ts the test statistic, 'ecdcor' or 'ecdcov'
#' @details
#' For the permutation test, if no value is given for 'B', (200+500/n) replicates
#' are used by default, where n is the sample size. At least 99 and at most 999 random permutations
#' should be sufficient practically.
#'
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
#' @author Qinqcong Yuan (\email{yuanq3@@miamioh.edu}) and Xiangrong Yin
#' @references \insertRef{yin2019new}{eccficm}
#'
#' @return \code{ecd.test} returns a list containing the following components:
#' \describe{
#'  \item{estimation.method}{a character string indicating the method used to compute the ECD measure}
#'  \item{method}{a character string describing the type of ECD test performed}
#'  \item{statistic}{the value of the ecd test statistic}
#'  \item{p.value}{approximate p-value of the test based on a permutation test}
#'  \item{alternative }{a character string describing the alternative hypothesis}
#'  \item{estimates}{vector containing ECD statistics, ecdCov, ecdCor, and ecdVarX}
#'  \item{data.name}{a character string containing the names of the data used}
#' }
#' @examples
#' library(MASS)
#' n <- 30; p <- 5; q <- 1
#' mu <- rep(0,(p+q))
#' sigma <- diag((p+q))
#' DATA <- mvrnorm(n, mu, sigma)
#' X <- DATA[,1:p]
#' Y <- DATA[,(p+1):(p+q)]
#'
#' # Using slicing estimation method on a continuous Y
#' ecd.test(X, Y, est = "slice", ns = 5,ts="ecdcor", B = "default")
#'
#' # Using Gaussian kernel estimation method
#' ecd.test(X, Y, est = "kernel.gau", ts="ecdcor", B = "default")
#'
#' # Using Epanechnikov kernel estimation method
#' ecd.test(X, Y, est = "kernel.epa", ts="ecdcor", B = "default")
#'
#' # This example shows how to directly perform the test if Y is a categorical variable
#' library(boot)
#' #gravity data is in the boot package,there are 81 measurements in a series
#' # of eight experiments
#' #gravity$series shows the categories 1-8 of the observations. 999 permutations
#' # are used to compute p-value.
#'str(gravity)
#' ecd.test(gravity$g, gravity$series, est="slice", ts="ecdcor", B=999)
#'
#'@seealso \code{\link{ecd}} \code{\link{eccfic.test}}
#' @export

ecd.test <-
function(x, y, est="slice", ns=5, ts="ecdcor", B="default",bw="default", index=1.0){

  est <- tolower(est); ts <- tolower(ts) # make case insensitive
  ts <- ifelse(ts=='ecdcov','ecdCov', 'ecdCor')
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

  # Set index of the test statistic to be used (1:ecdCov, 2:ecdCor)
  tstats.index <- ifelse(ts=='ecdCor',2, ifelse(ts=='ecdCov', 1, 2))

  # Compute sample statistics
  estimates <- ecd(x,y,est,ns,bw,index)
  ecdcov <- estimates$ecdCov
  ecdcor  <- estimates$ecdCor
  ecdvarx <- estimates$ecdVarX
  test.statistic <- estimates[[tstats.index]]

  # compute p-value using permutation test
  reps <- replicate(B,ecd(x,y[sample(1:n),],est,ns,bw,index)[[tstats.index]])
  pvalue <- (sum(reps>=test.statistic)+1)/(B+1)

  estimates <- c(ecdcor, ecdcov, ecdvarx)
  names(estimates) <- c("ecdCor", "ecdCov", "ecdVar(X)")
  names(test.statistic) <- ts
  kernel.used <- ""
  if(est=="kernel.epa") {kernel.used<-"Kernel estimation (Epanechnikov)"}
  else if(est=="kernel.gau"){kernel.used <-"Kernel estimation (Gaussian)"}

  est.method <- paste("Estimation Method: ", ifelse(est=='slice', 'Slicing', kernel.used))
  test.method <- paste(ts," permutation test of independence - ", est.method)
  alternative <- "X is not independent of Y"
  dataname <- ifelse((xname=="X")&(yname=="Y"), paste("X, Y, Replicates", B),
                     paste("X -",xname, ", Y -", yname, ", Replicates", B))

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



