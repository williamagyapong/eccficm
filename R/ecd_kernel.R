#----------------------------------------------------------------------------
# The function depends on the kernel.fun() and mhatY() functions defined in the
# helper section that follows to compute ECD measures for both the kernel method with Epanechnikov
# kernel, and kernel method with Gaussian kernel
#____________________________________________________________________________

# @title ECD kernel Estimation Methods
# @description Implements ECD kernel estimation methods
# @param x data of first sample serving as the response
# @param y data of second sample conditioned on
# @param method ECD estimation method (options:kernel.epa, kernel.gau)
# @param bw bandwidth
# @param index exponent on Euclidean distance, in (0,2]

ecd.kernel <-
function(x,y, method, bw="default", index=1.0){
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  q <- ncol(y)

  # validate inputs
  if(q>1)  stop("Sorry, this current implementation does not support a multivariate Y. Check out for updates.")

  if (n != m) stop("Sample sizes for x and y must agree")

  if (index <= 0 || index > 2) {
    warning("Index must be in (0,2], using default index=1")
    index=1.0 # reset index to default
  }

  xdist <- (as.matrix(dist(x)))^index
  term1 <- sum(xdist)/(n^2)

  # computing Cn(X|X)
  CnXX <- CnXY <- sqrt(term1) # a measure analogous to variance of X (not very sure)

  if(!identical(x,y)) {
    term2 <- as.numeric(mhatY(x,y, method, bw))
    CnXY <- sqrt(term1 - term2) # Cn(X|Y), a measure analogous to covariance between X and Y (not very sure)
  }

  rhoc <- CnXY/CnXX
  return(list(ecdCov=CnXY, ecdCor=rhoc, ecdVarX=CnXX))
}

#-----------------------------------------------------------
# Helper functions
#__________________________________________________________

kernel.fun <-
  function(t, type) {
    if(type=="kernel.epa") {

      kernel <- ifelse(abs(t)<=sqrt(1), 3/4*(1-t^2), 0) # want to replace sqrt() by 1
    }else if(type=="kernel.gau") {
      kernel <- 1/sqrt(2*pi)*exp(-1/2*t^2) # Gaussian kernel
    } else {
      stop("Type of kernel needs to be specified (kernel.epa for Epanechnikov kernel estimation or kernel.gau for Gaussian kernel estimation")
    }
    return(kernel)
  }


mhatY <-
  function(x,y, method, bw){
    n <- nrow(x)
    sigma <- sd(y)

    if(method=="kernel.epa") {
      bw <- ifelse(bw=="default", 1.5*sigma*n^(-1/(4+1/4)), bw)
    } else if(method=="kernel.gau") {
      bw <- ifelse(bw=="default", 1.06*sigma*n^(-1/5), bw)
    }

    mat1 <- matrix(rep(y,n),ncol=n)
    mat2 <- (mat1-t(mat1))/bw
    Kamat <- matrix(mapply(kernel.fun,as.vector(mat2), method),nrow=n,byrow=FALSE)
    Wa <- Kamat/(apply(Kamat,1,sum))
    Gmat <- as.matrix(dist(x))
    tmp <- Wa%*%Gmat%*%t(Wa)
    result <- sum(diag(tmp))/n # find the trace and divide by n
    return(result)
  }
