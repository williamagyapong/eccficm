#----------------------------------------------------------------------------
# This function depends on the kern.fun() and mhat() functions defined in the
# helper section that follows to compute ECD measures for both the kernel method
# with Epanechnikov kernel, and kernel method with Gaussian kernel
#____________________________________________________________________________

# @title ECD kernel Estimation Methods
# @description Implements ECD kernel estimation methods
# @inheritparams ecd

ecd.kernel <-
function(x,y, est="kernel.gau", bw="default", index=1.0){
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  q <- ncol(y)

  if (n != m) stop("Sample sizes for x and y must agree")

  if (index <= 0 || index > 2) {
    warning("Index must be in (0,2], using default index=1")
    index=1.0 # reset index to default
  }

  xdist <- (as.matrix(dist(x)))^index
  term1 <- sum(xdist)/(n^2)

  # computing Cn(X|X)
  CnXX <- CnXY <- sqrt(term1) # a measure analogous to variance of X

  if(!identical(x,y)) {
    term2 <- as.numeric(mhat(x,y, est, bw, index)) # see mhat() at the end of this file
    CnXY <- sqrt(abs(term1 - term2)) # Cn(X|Y), a measure analogous to covariance between X and Y
    # print(term1-term2)
  }

  rhoc <- CnXY/CnXX # this is the correlation type statistic
  return(list(ecdCov=CnXY, ecdCor=rhoc, ecdVarX=CnXX))
}

#-----------------------------------------------------------
# Helper functions
#__________________________________________________________

# kernel function
kern.fun <-
  function(t, type) {
    if(type=="kernel.epa") {
      ifelse(abs(t)<=1, 3/4*(1-t^2), 0)
    }else if(type=="kernel.gau") {
      1/sqrt(2*pi)*exp(-1/2*t^2) # Gaussian kernel
    } else {
      stop("kernel estimation method needs to be specified (kernel.epa for
        Epanechnikov kernel estimation or kernel.gau for Gaussian kernel estimation")
    }
  }

# Estimating the second term
mhat <-
  function(x,y, est, bw,index){

    n <- nrow(x)
    sigma <- sd(y)

    if(est=="kernel.epa") {
      bw <- ifelse(bw=="default", 1.5*sigma*n^(-1/(4+1/4)), bw)
    } else if(est=="kernel.gau") {
      bw <- ifelse(bw=="default", 1.06*sigma*n^(-1/5), bw)
    }

    # mat1 <- matrix(rep(y,n),ncol=n)
    # mat2 <- (mat1-t(mat1))/bw
    # Kamat <- matrix(mapply(kern.fun,as.vector(mat2), est),nrow=n,byrow=FALSE)
    Kamat <- kern.fun(as.matrix(dist(y))/bw, est) # this replaces lines 72 - 74 to accommodate multivariate Y
    Wa <- Kamat/(apply(Kamat,1,sum))
    Gmat <- as.matrix(dist(x)^index)
    tmp <- Wa%*%Gmat%*%t(Wa)
    result <- sum(diag(tmp))/n # find the trace and divide by n
    return(result)
  }

# profvis({
#   n <- 30; p <- 5; q <- 1
#   mu <- rep(0,(p+q))
#   sigma <- diag((p+q))
#   DATA <- mvrnorm(n, mu, sigma)
#   X <- DATA[,1:p]
#   Y <- DATA[,(p+1):(p+q)]
#
#    ecd.kernel(X,Y)
# })




