flipIt <- function(psi,q)
{
		# invertibility adjustment
		#  psi is a vector with q ma parameters, followed by ar parameters, with last
		#	entry the logged innovation variance
		#  uses ARMA convention, with plus for MA and minus for AR
		#  returns new psi vector of same form, but where the new parameters correspond
		#	to ar and ma invertible polynomials, i.e., with roots outside the unit circle
		#   the (log) innovation variance is changed appropriately

		polyroot2 <- function(poly)
		{
			ar.rep <- -1*poly[-1]/poly[1]
			p.order <- length(ar.rep)
			A.mat <- NULL
			if(p.order > 1) { A.mat <- diag(p.order)[1:(p.order-1),] }
			A.mat <- rbind(ar.rep,A.mat)
			eig.A <- eigen(A.mat)
			ar.roots <- 1/eig.A$values
			return(ar.roots)
		}

		ceps2wold <- function(ceps,q)
		{
			m <- length(ceps)
			if(q > m) { ceps <- c(ceps,rep(0,q-m)) }
			wold <- 1
			wolds <- wold
			for(j in 1:q)
			{
				wold <- sum(seq(1,j)*ceps[1:j]*wolds[j:1])/j
				wolds <- c(wolds,wold)
			}
			return(wolds)
		}

		roots2ceps <- function(roots,m)
		{
			p <- length(roots)
			ceps <- rep(0,m)
			for(k in 1:m)
			{	
				ceps[k] <- -1*sum(roots^(-k))/k
			}
			return(ceps)
		}

		# Slight modification to ensure no real zeros are passed to polyroot2
		r <- length(psi) - 1
		if (r > 0) psi[1:r] <- psi[1:r] + .000001
			
		p <- r - q
		if (p==0) { ar <- NULL } else { ar <- psi[(q+1):(q+p)] }
		if (q==0) { ma <- NULL } else { ma <- psi[1:q] }

		factor <- 1
		if (p > 0)
		{
			newAR <- 1
			arRoots <- polyroot2(c(1,-1*ar))
			for(k in 1:length(arRoots))
			{
				if (Mod(arRoots[k]) < 1)
				{
					factor <- factor/arRoots[k]^2
					arRoots[k] <- 1/arRoots[k]
				}
			}
			my.ceps <- roots2ceps(arRoots,length(arRoots))
			newAR <- ceps2wold(Re(my.ceps),length(arRoots))
			ar <- -1*newAR[-1]
			ar <- Re(ar)
			factor <- Re(factor)
		}

		if (q > 0)
		{
			newMA <- 1
			maRoots <- polyroot2(c(1,ma))
			for(k in 1:length(maRoots))
			{
				if (Mod(maRoots[k]) < 1)
				{
					factor <- factor*maRoots[k]^2
					maRoots[k] <- 1/maRoots[k]
				}
			}
			my.ceps <- roots2ceps(maRoots,length(maRoots))
			newMA <- ceps2wold(Re(my.ceps),length(maRoots))
			ma <- newMA[-1]
			ma <- Re(ma)
			factor <- Re(factor)
		}

		newPsi <- c(ma,ar,psi[r+1] - log(factor))
		return(newPsi)
}
