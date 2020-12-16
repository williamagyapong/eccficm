#' @title ECD Independence Test
#' @name ecdcov.test
#' @description Performs a permutation test of independence based on
#' the distance covariance type of measure
#' @param B number of permutation replicates
#' @inheritParams ecd
#' @return \code{ecdcov.test} returns a list containing the following components:
#' \describe{
#'  \item{estimation.method}{a character string indicating the method used to
#'  compute the ECD measure}
#'  \item{method}{a character string describing the type of ECD test performed}
#'  \item{statistic}{the value of the ecd test statistic}
#'  \item{p.value}{approximate p-value of the test based on a permutation test}
#'  \item{alternative }{a character string describing the alternative hypothesis}
#'  \item{estimates}{vector containing ecdCov, ecdCor, and ecdVarX}
#'  \item{data.name}{a character string containing the names of the data used}
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
#' ecdcov.test(x = X, y = Y, method = "slice", B = "default", ns = 5)
#' # using Gaussian kernel estimation method
#' ecdcov.test(x = X, y = Y, method = "kernel.gau", B = "default")
#' # using Epanechnikov kernel estimation method
#' ecdcov.test(x = X, y = Y, method = "kernel.epa", B = "default")
#'
#' #This example shows how to directly perform the test if Y is a categorical
#' #variable
#' library(boot)
#' #gravity data is in the boot package,there are 81 measurements in a series of
#' #eight experiments
#' #gravity$series shows the categories 1-8 of the observations. 999 permutations
#' # are used to compute p-value.
#'str(gravity)
#' ecdcor.test(gravity$g, gravity$series, method="slice", B=999)
#'@seealso @seealso \code{\link{ecd}} \code{\link{ecdcor.test}}
#' @export

ecdcov.test <-
function(x, y, method, ns=NULL, B="default",bw="default", index=1.0){
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))

  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(y)

  if(B=="default") {
    # first check the value of B
    B <- floor((200+5000/n)) # based on a method suggested by Szekely, et. al. (2007)
    if(B < 99) warning("Insufficient number of permutations! Consider changing \"B\" from default to a number greater than 99.")
  }
  permute.vec <- rep(0,B)
  dataname <- ifelse((xname=="X")&(yname=="Y"), paste("X, Y, replicates = ", B), paste("X ->",xname, ", Y ->", yname, ", replicates = ", B))


  #set.seed(456)
  Yind=matrix(rep(seq(1,n),B),ncol=B)
  Ysample <- apply(Yind,2,sample)
  for(s in 1:B){
    X <- x[Ysample[,s],]
    permute.vec[s] <- ecd(X,y,method,ns,bw,index)$ecdCov
  }

  # compute sample statistics
  estimates <- ecd(x,y,method,ns,bw,index)
  ecdcov <- estimates$ecdCov
  ecdcor <- estimates$ecdCor
  ecdvarx <- estimates$ecdVarX

  pvalue <- (1+sum(permute.vec>=ecdcov))/(B+1)
  #pval <- (sum(replicate(B,ecd(x,sample(y),method,ns,bw,index)$ecdCov)>=rhoc.observed)+1)/(B+1)

  estimates <- c(ecdcov, ecdcor, ecdvarx)
  names(estimates) <- c("ecdCov", "ecdCor", "ecdVarX")
  teststat <- ecdcov # not very sure (need to confirm)
  names(teststat) = 'Cn(X|Y)'
  kernel.used <- ""
  if(method=="kernel.epa") {kernel.used<-"Epanechnikov kernel"}
  else if(method=="kernel.gau"){kernel.used <-"Gaussian kernel"}
  est.method <- paste("Estimation approach: ", ifelse(method=='slice', 'Slicing', kernel.used))
  test.method <- paste("ecdCov test of independence (permutation test) - ", est.method)
  alternative <- "the random vector X is not independent of the random vector Y"

  output <- list(
    method = test.method,
    estimation.method = est.method,
    estimates = estimates,
    statistic = teststat,
    alternative = alternative,
    p.value = pvalue,
    data.name = dataname
  )

  class(output) <- 'htest'
  #return(list(Rhoc=rhoc.observed,p_value=pvalue, p_value2=pval))
  return(output)
}



