seas.proxim <- function(psi,phi,theta,rho,omega)
{
	p.order <- length(phi)
	q.order <- length(theta)
	zeta <- rho + (1-rho)*exp(psi)/(1+exp(psi))
	zeta <- zeta*exp(1i*omega)
	val <- 1
	if(p.order > 0) 
	{
		xi.p <- zeta^(-seq(1,p.order))
		val <- val*(Mod(1 - sum(phi*xi.p)))^2
	}
	if(q.order > 0)
	{
		xi.q <- zeta^(-seq(1,q.order))
		val <- val/(Mod(1 + sum(theta*xi.q)))^2
	}
	return(val)
}
