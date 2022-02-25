# @title ECCFIC Kernel Regression Estimation Method
# @name eccfic.kernreg
# @description computes the ECCFIC between x and a continuous y only,
# @inheritParams eccfic
# @author Chenlu Ke, Qinqcong Yuan and William Agyapong
# @references \insertRef{ke2019expected}{eccficm}
# @usage eccfic.kernreg(x, y, kernel='gaussian', sigma='default', bw='default', index=1, wt=F)
# @examples
# x <- rnorm(100)
# y <- x + rnorm(100)
# eccfic.kernreg(x,y)
# @seealso @seealso \code{\link{eccfic}} \code{\link{eccfic.slice}}  \code{\link{eccfic.test}}
#
eccfic.kernreg <-
function(x, y, kernel='gaussian', sigma='default', bw='default', index=1, wt=F){

  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  q <- ncol(y)

  # validate inputs
  if (n != m) stop("Sample sizes for x and y must agree")

  if (index <= 0 | index > 2) {
    warning("Index must be in (0,2], using default index=1")
    index=1.0 # reset index to default
  }

  #n=1, trivial
  if (n==1){
    return(0)
  }

  # Compute H2(X|X) or H2(Y|Y)
  # compute heuristic sigma
  if(kernel=='gaussian'){
    if(sigma=='default') {
      sigma <- sqrt(0.5*median(dist(x)^2))
    }
  }

  # compute H2(X|X) or H2(Y|Y)
  ## see the natural estimator under Thm 7 (page 15)
  if(identical(x, y)){
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

  #compute ECCFIC: H2(X|Y)
  if(bw=='default'){
    if(q==1){
      bw <- 1.06*sd(y)*n^(-1/5)
    }else{
      bw <- (4/n/(q+2))^(1/(q+4))*sum(diag(cov(y)))/q
    }
  }

  H <- diag(n) - matrix(1, n, n)/n
  if(kernel=='gaussian'){
    if(sigma =='default'){
      sigma <- sqrt(0.5*median(dist(x)^2))
      if(sigma==0){sigma <- 0.001} # ensure sigma isn't 0
    }
    K <- dnorm(as.matrix(dist(x,diag=T,upper=T)), mean=0, sd=sigma)*sqrt(2*pi)*sigma
  }else{
    if(kernel=='distance'){
      normx <- matrix(apply(x*x,1,sum)^(index/2),n,n)
      K <- 0.5*(normx+t(normx)-as.matrix(dist(x,diag=T,upper=T))^index)
    }
  }

  if(wt==T){
    G <- dnorm(as.matrix(dist(y,diag=T,upper=T)), mean=0, sd=bw)
    return(sum(diag(K%*%H%*%G%*%G%*%H))/n^3)
  }else{
    G <- dnorm(as.matrix(dist(y,diag=T,upper=T)), mean=0, sd=bw)
    Gstar <- t(G/rowSums(G))
    return(sum(diag(K%*%H%*%Gstar%*%t(Gstar)%*%H))/n)
  }
}
