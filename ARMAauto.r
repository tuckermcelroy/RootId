ARMAauto <- function(phi,theta,maxlag)
{
	#### Function computes autocovariances of ARMA (p,q) from lag zero
	#	to maxlag, with inputs phi and theta.
	#	(1 - phi[1]z ... - phi[p]z^p) X_t = (1 + theta[1]z ...+ theta[q]z^q) WN
	#  output: autocovariance string of length maxlag
	#  for absent AR or MA portions, pass in NULL

polymul <- function(a,b)
{
	u <- convolve(Re(a),rev(Re(b)),type="open") - convolve(Im(a),rev(Im(b)),type="open")
	v <- convolve(Im(a),rev(Re(b)),type="open") + convolve(Re(a),rev(Im(b)),type="open")
	return(Re(u + 1i*v))
}

p <- length(phi)
q <- length(theta)
gamMA <- polymul(c(1,theta),rev(c(1,theta)))
gamMA <- gamMA[(q+1):(2*q+1)]

if (p > 0) 
{
Amat <- matrix(0,nrow=(p+1),ncol=(2*p+1))
for(i in 1:(p+1))
{
	Amat[i,i:(i+p)] <- c(-1*rev(phi),1)
}
Amat <- cbind(Amat[,(p+1)],as.matrix(Amat[,(p+2):(2*p+1)]) + t(matrix(apply(t(matrix(Amat[,1:p],p+1,p)),2,rev),p,p+1)))

Bmat <- matrix(0,nrow=(q+1),ncol=(p+q+1))
for(i in 1:(q+1))
{
	Bmat[i,i:(i+p)] <- c(-1*rev(phi),1)
}
Bmat <- t(matrix(apply(t(Bmat),2,rev),p+q+1,q+1))
Bmat <- matrix(apply(Bmat,2,rev),q+1,p+q+1)
Bmat <- Bmat[,1:(q+1)]
Binv <- solve(Bmat)

gamMix <- Binv %*% gamMA
if (p <= q) gamMix <- matrix(gamMix[1:(p+1),],p+1,1) else 
{ gamMix <- matrix(c(gamMix,rep(0,(p-q))),p+1,1) }
gamARMA <- solve(Amat) %*% gamMix 
} else gamARMA <- gamMA[1]

gamMA <- as.vector(gamMA)
if (maxlag <= q) gamMA <- gamMA[1:(maxlag+1)] else gamMA <- c(gamMA,rep(0,(maxlag-q)))
gamARMA <- as.vector(gamARMA)
if (maxlag <= p) gamARMA <- gamARMA[1:(maxlag+1)] else {
for(k in 1:(maxlag-p))
{
	len <- length(gamARMA)
	acf <- gamMA[p+1+k]
	if (p > 0) acf <- acf + sum(phi*rev(gamARMA[(len-p+1):len]))
	gamARMA <- c(gamARMA,acf)
}
}
return(gamARMA)
}
