# @title ECCFIC Slicing Estimation Method
# @name eccfic.slice
# @description computes the ECCFIC measure between x and y by slicing on y,
# where y is categorical or continuous. It is straightforward when y is
# originally categorical, otherwise y is first sliced.
# @inheritParams eccfic
# @author Chenlu Ke, Qinqcong Yuan and William Agyapong
# @references \insertRef{ke2019expected}{eccficm}
# @usage eccfic.slice(x, y,ns=NULL, kernel='gaussian', sigma='default', index=1)
# @examples
# y <- rep(1:4, 25)
# x <- y + rnorm(100)
# eccfic.slice(x,y,ns=5)
#@seealso @seealso \code{\link{eccfic.kernreg}} \code{\link{eccfic.test}}

eccfic.slice <-
function(x, y, ns=NULL, kernel='gaussian', sigma='default', index=1){

  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  q <- ncol(y)

  # validate inputs
  if (n != m) stop("Sample sizes for x and y must agree")

  if (index <= 0 || index > 2) {
    warning("Index must be in (0,2], using default index=1")
    index=1.0 # reset index to default
  }

  #n=1, trivial
  if(n==1) return(0)


  # compute heuristic sigma
  if(kernel=='gaussian'){
    if(sigma=='default') {
      sigma <- sqrt(0.5*median(dist(x)^2))
    }
  }

  ## compute H2(X|X), see the natural estimator under Thm 7 (page 15)
  if(identical(x,y)){
    if(kernel=='gaussian'){
      K <- dnorm(as.matrix(dist(x,diag=T,upper=T)), mean=0, sd=sigma)
      return((dnorm(0,0,sigma)-mean(K))*sqrt(2*pi)*sigma)
    }else{
      if(kernel=='distance'){
        K <- 0.5*(as.matrix(dist(x,diag=T,upper=T))^index)
        return(mean(K))
      }
    }
  }

  ## compute ECCFIC: H2(X|Y)

  # slice the continuous Y variable
  if(!(is.factor(y) | is.character(y)) & !is.null(ns)) {
    if(q>1)  stop("No support for slicing on a multivariate Y. See the 'Details' section of the help documentation")

    y <- sliceY(y, ns) # see function in the 'sliceY.R' file
  }

  if(kernel=='gaussian'){
    K <- dnorm(as.matrix(dist(x,diag=T,upper=T)), mean=0, sd=sigma)
    L <- ifelse(as.matrix(dist(y,diag = T,upper = T))==0,1,0)
    L <- n*L/rowSums(L)-1
    return(sum(K*L)/(n^2)*sqrt(2*pi)*sigma)
  }else{
    if(kernel=='distance'){
      K <- -0.5*(as.matrix(dist(x,diag=T,upper=T))^index)
      L <- ifelse(as.matrix(dist(y,diag = T,upper = T))==0,1,0)
      L <- n*L/rowSums(L)-1
      return(sum(K*L)/(n^2))
    }
  }
}
