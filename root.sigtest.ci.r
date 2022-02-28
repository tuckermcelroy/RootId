root.sigtest.ci <- function(signal,p.order,diffs,omega,monte,rhos)
{

      ##########################################################################
      #
      #       root.sigtest.ci
      #           Copyright (C) 2019  Tucker McElroy
      #
      #    This program is free software: you can redistribute it and/or modify
      #    it under the terms of the GNU General Public License as published by
      #    the Free Software Foundation, either version 3 of the License, or
      #    (at your option) any later version.
      #
      #    This program is distributed in the hope that it will be useful,
      #    but WITHOUT ANY WARRANTY; without even the implied warranty of
      #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      #    GNU General Public License for more details.
      #
      #    You should have received a copy of the GNU General Public License
      #    along with this program.  If not, see <https://www.gnu.org/licenses/>.
      #
      ############################################################################

      ################# Documentation #####################################
      #
      #       Purpose: compute Root seasonality diagnostic for input signal
      #       Inputs:
      #                 signal: time series signal to be tested  
	#			p.order: order of AR model to be fitted; leave as NULL
	#				if AIC will be used to determine order
	#			diffs: number of trend differences (1-B) to be applied to signal
	#			omega: argument of a vector of complex numbers (inside unit circle) 
	#				corresponding to the null hypothesis	
	#			monte: number of monte carlo sims used to compute distributions
 	#			rhos: various values of rho
      #       Outputs:
	#			pvals: p-values for each rho
      #	  Requires: ARMAauto.r, seas.proxim.r, flipIt.r
	#
      ############################################

	arp.fityw <- function(data,p)
	{
		gamma.hat <- acf(data,lag=p,type="covariance",plot=FALSE)$acf[,,1]
		phi <- solve(toeplitz(gamma.hat[1:p]),gamma.hat[2:(p+1)])
		sigma2 <- gamma.hat[1] - sum(phi*gamma.hat[2:(p+1)])
		hess <- sigma2*diag(solve(toeplitz(gamma.hat[1:p])))
		return(list(phi,sigma2,hess))
	}

	ar.data <- signal 
	if(diffs > 0) { ar.data <- diff(ar.data,differences = diffs) }
	T <- length(ar.data)

	# fit AR model with OLS and AIC for order
	if(length(p.order)>0)
	{
		my.arfit <- ar.ols(ar.data,order.max=p.order,aic=FALSE)
	} else
	{
		max.ord <- floor(T/3)
		my.vars <- NULL
		for(p in 1:max.ord) {
			my.vars <- c(my.vars,arp.fityw(ar.data,p)[[2]]) }
#		bics <- T*(1+log(my.vars)) + (1+seq(1,max.ord))*log(T)
		aics <- T*(1+log(my.vars)) + (1+seq(1,max.ord))*2
		p.order <- which.min(aics)
		my.arfit <- ar.ols(ar.data,order.max=p.order,aic=FALSE)
#		p.order <- my.arfit$order
	}
 
	mesh <- length(rhos)
	pvals <- rep(0,mesh)
	if(p.order > 0) {

	# uncertainty of AR coefficient estimates
	new.ar <- flipIt(c(my.arfit$ar,log(my.arfit$var)),0)
	gamma.ar <- ARMAauto(new.ar[1:p.order],NULL,p.order)*exp(new.ar[(p.order+1)])
	gamma.mat <- toeplitz(gamma.ar[1:p.order])
	sigma2 <- (gamma.ar[1] - t(gamma.ar[2:(p.order+1)]) %*% 
		solve(gamma.mat) %*% gamma.ar[2:(p.order+1)])[1,1]
	v.half <- solve(chol(gamma.mat))*sqrt(sigma2)

	# get test statistics and confidence interval
	gauss.draws <- matrix(rnorm(p.order*monte),nrow=p.order)
	m <- length(omega)
	for(k in 1:mesh)
	{
		test.stats <- NULL
		dist.draws <- NULL
		rho <- rhos[k]
		for(j in 1:m)
		{
			omega.j <- omega[j]
			test.stat.j <- T*seas.proxim(-Inf,phi=new.ar[1:p.order],theta=NULL,
				rho=rho,omega=omega.j)
			test.stats <- c(test.stats,test.stat.j)
			xi.j <- rho*exp(1i*omega.j)^(-seq(1,p.order))
 
			# uncertainty of signal test statistic
			z.draw.j <- t(Conj(xi.j)) %*% v.half %*% gauss.draws	
			dist.draw.j <- (Mod(z.draw.j))^2
			dist.draws <- rbind(dist.draws,dist.draw.j)
		}
		test.stat <- min(test.stats)
		dist.draw <- apply(dist.draws,2,min)
		pvals[k] <- length(dist.draw[dist.draw > test.stat])/monte
#		if(k %% 100 == 0) print(k)
	}
	
	}  

	return(pvals)
}

