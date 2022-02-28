#Rfunctions for the project 'Identification of the  Differencing Operator 
#of a Non-stationary Time Series  via Testing for Zeroes in the Spectral Density'

#load the following packages; install if needed first
library('data.table')#data manipulation
library('ggplot2')#plotting
library('itsmr')#Rfunction 'hannan' (Gomez (2013))

######################################################
#auxilary Rfunctions
my.ccf<-function(a, b) 
#function to calculate sample cross-covariance from Rpackage "waveslim" (without de-meaaning)
#a=series1, b=series2  
{
  n <- length(a)#=T
  #a<-a-mean(a)
  #b<-b-mean(b) 
  a <- c(a, rep(0, n))
  b <- c(b, rep(0, n))
  x <- Re(fft(fft(a) * Conj(fft(b)), inverse = TRUE))/2/n^2
  ccf<-x[c((n + 2):(2 * n), 1:n)]
  return(ccf)#denom=T; for lag h from -(T-1) to T-1
}

polymult <- function(a,b)
#function to perform polynomial multiplication
{
  n <- length(a)
  m <- length(b)
  a <- c(a,rep(0,m))
  b <- c(b,rep(0,n))
  out <- rep(0,n+m-1)
  for(i in 1:(n+m-1))
  {
    out[i] <- a[1:i] %*% b[i:1]
  }
  return(out)
}

Bartlett<-function(x,fixedb=0.3) 
#Bartlett taper at x, using "fixed-b" fraction
{
  y<-x/fixedb
  value<-ifelse(y<=1 & y>=-1,1-abs(y),0)
  return(value)
}  

Parzen<-function(x,fixedb=0.3)
#Parzen taper at x, using "fixed-b" fraction
{
  y<-x/fixedb
  value<-ifelse(y<=0.5 & y >=-0.5, 1-6*y^2+6*abs(y)^3, 2*(1-abs(y))^3)  
  value<-ifelse(y<=1 & y>=-1, value, 0)
  return(value)
}  
######################################################
#functions for the estimator and subsampling test
mySchurest<-function(x,lambda=0,fixedb=0.3,taper=Bartlett)
#function that implements estimator hat{f}_b(\lambda), b=fixedb fraction
#input:
#x: sample data (vector of length T)
#lambda: frequency of interest  
#fixedb: fixed-b fraction
#output: the estimate
{
  T<-length(x)
  GammaHat<-my.ccf(a=x,b=x)

  hh<- -(T-1):(T-1) #sequence of lags
  xx<- hh/T #rescaled lags, the argument of the taper function
  
  Taperb<-taper(xx,fixedb)
  fhat<-sum(Taperb*GammaHat*exp(-1i*lambda*hh)) 
  est<-Re(fhat)#\hat{d}_1 from McElroy, T. and Jach, A. (2019), The Econometrics Journal,  22, 97-116
  return(est)
}  

mytest<-function(lambda,x,b,nominal=0.01,param=0,pow=1,myest,...)
#function that implements subsampling test of the current article, 
#subsampling test for estimator 'myest' with rate 'T^pow'
#lambda: frequency of interest
#x: sample data (vector of length T)
#b: subsample size 
#nominal: nominal significance level (can be a vector)
#param: value of the parameter under H0
#pow: power of T in the rate T^pow
#myest: name of a function that computes the estimator
#...: other values passed to 'myest'  
#output: list with p-value of the test, test stat, subsampling quantile, subsampling statistics  
{ 
  #create T-b+1 subsamples of size b based on x and tag them with l=1,2,...,T-b+1
  T<-length(x)
  mylist<-shift(x=x, n=0:(T-b), fill=NA, type="lead")
  mylist<-lapply(mylist,function(x)head(x,b))
  df<-reshape2::melt(mylist)
  names(df)<-c("x","l")#l variable indexes subsamples, l=1,2,...,T-b+1
  dt<-data.table(df)#convert to a data.table format 
  
  #large-sample estimate
  est<-myest(x=x,lambda=lambda,...)

  temp<-dt[,.(stat=(b^pow)*(myest(x=x,lambda=lambda,...)-param)),by=.(l)]#use param to center or 
  #temp<-dt[,.(stat=(b^pow)*(myest(x=x,...)-est)),by=.(l)]#use a large-sample estimate to center
  
  #subsample values of the test stat
  substat<-as.vector(unlist(temp[,"stat"]))
  
  #pvalue of the test
  stat<-(T^pow)*(est-param)
  pval<-sum(substat>=stat)/length(substat)
  
  list(pval=pval,stat=stat,q=quantile(substat,1-nominal,names=FALSE),substat=substat)
}

myjointtest<-function(lambdas,x,b,nominal=0.01,param=0,pow=1,myest,...)
#function that implements joint subsampling test of the current article  
#lambdas: frequencies to be tested (vector)
#x: sample data (vector) 
#b: subsample size (vector of length of at least 3, for adaptive selection), otherwise the first value
#nominal: nominal significance level
#fixedb: fixed-b fraction
#param=0,pow=1,fixedb,taper: other arguments used in 'myest'  
#output: list with the p-value of the joint test, subsamplng stats, optimal_b
{ 
  d<-length(lambdas)

  if(length(b)>=3)#if at least 3 candidate block sizes
  {
    mydist<-mystat<-vector('list',length=length(b))
    for(i in 1:length(b))
    {
      stat<-substat<-NULL#variable to store results for a given b and all lambdas
      for(l in 1:d)
      {
        out<-mytest(lambda=lambdas[l],x=x,b=b[i],nominal=nominal,param=param,pow=pow,myest,...)
        #myest will use lambda and x for its first arguments

        stat<-c(stat,out$stat)#stat is a vector of length d
        substat<-rbind(substat,out$substat)#substat is a matrix d x (T-b+1)
        
      }#end of l/lambdas loop

      mydist[[i]]<-apply(substat,2,min)# min wrt j (ie, columns) yields a vector of length T-b+1
      mystat[[i]]<-min(stat)
       
    }#end of i loop over block sizes  
    
    ds<-NULL
    #find KS distances between distribution functions for b_i and b_{i+1}
    for(i in 1:(length(b)-1))
    {
      temp<-ks.test(mydist[[i]],mydist[[i+1]])
      ds<-c(ds,temp$statistic)
    }  
    
    iopt<-which.min(ds)#which slot/index 'i' of b_1,b_2,... corresponds to the smallest D
    bopt<-b[iopt]#get the corresponding bopt
    minsubstat<-mydist[[iopt]]#extract the corresponding subsampling dist
    minstat<-mystat[[iopt]]#extract the corresponding large sample stat

  } else {  
    b<-bopt<-b[1]#if 1 or 2 b values, take the first one and call it bopt 
    stat<-substat<-NULL
    for(l in 1:d)
    {
      out<-mytest(lambda=lambdas[l],x=x,b=b,nominal=nominal,param=param,pow=pow,myest,...)
                  #myest will use lambda and x for its first arguments

      stat<-c(stat,out$stat)#stat is a vector of length d
      substat<-rbind(substat,out$substat)#substat is a matrix d x (T-b+1)
      
    }#end of lambdas loop
    
    minsubstat<-apply(substat,2,min)# min wrt j (ie, columns) yields a vector of length T-b+1
    minstat<-min(stat)
  }
  #find p-value
  pvaljoint<-sum(minsubstat>=minstat)/length(minsubstat)
  
  list(pval=pvaljoint,minstat=minstat,minsubstat=minsubstat,bopt=bopt)
}

myprocedureBackward<-function(filters,lambdas,x,b,nominal=0.05,param=0,pow=1,myest,...)
#function that implements Stage2 of the current article
#filters: filters associated with lambdas (list)
#lambdas: frequencies to be tested (vector)
#x: data to be differenced (vector)
#b: fixed-b fraction
#nominal: nominal significance level
#output: list with p-value, final set J, test stat, subsampling stats, b_opt (if adaptive b selection used)
{ 
  d<-length(lambdas)
  
  poly<-1
  for(i in 1:d)
    poly<-polymult(poly,filters[[i]])
  
  q.order<-2*sum(lambdas>0 & lambdas!=pi)+sum(lambdas==pi)+sum(lambdas==0)
  T<-length(x)-q.order
  
  w <- stats::filter(x=x,filter=poly,method="convolution",sides=1)[q.order+(1:T)]#filter x to obtain w
  
  out<-myjointtest(lambdas=lambdas,x=w,b=b,nominal=nominal,param=param,pow=pow,myest,...)
  
  pval<-out$pval
  J<-1:d#this will be returned if no over-differencing is present, ie, null is rejected, pval<=nominal
  stat<-out$minstat
  substat<-out$minsubstat
  bopt<-out$bopt
  
  if(pval>nominal)#if (full) joint hypothesis is not rejected, do BD to identify over-differenced freqs
  {
    J<-integer(length=0)#this will be returned if over-differencing is present at all frequencies
    #that is if the series was stationary (hence there should be no differencing done at all)
    j<-d-1
    while(j>=1)
    {
      combinations<-combn(d,j) 
      
      mylist<-NULL#for collecting combinations for which H0 is rejected
      for(k in 1:ncol(combinations))
      {
        Jj<-combinations[,k]
        poly<-1
        for(i in Jj)
          poly<-polymult(poly,filters[[i]])
        
        q.order<-2*sum(lambdas[Jj]>0 & lambdas[Jj]!=pi)+sum(lambdas[Jj]==pi)+sum(lambdas[Jj]==0)
        
        T<-length(x)-q.order
        
        w <- stats::filter(x=x,filter=poly,method="convolution",sides=1)[q.order+(1:T)]#filter x to obtain w
        
        out<-myjointtest(lambdas=lambdas[Jj],x=w,b=b,nominal=nominal,param=param,pow=pow,myest,...)
        
        #update all info
        pval<-out$pval
        stat<-out$minstat
        substat<-out$minsubstat
        bopt<-out$bopt
        
        if(pval<=nominal) #collect only pvals <= nominal
        {
          temp<-list(pval=pval,J=Jj,stat=stat,substat=substat,bopt=bopt)
          mylist<-c(mylist,temp)
        }  
      }# end of for(k in 1:ncol(combinations))  
      if(!is.null(mylist))#if after examining all Jj sets (doing k loop), there is at least one that rejects the null, return the info
      {
        pvals<-mylist[seq(1,length(mylist),5)]
        
        jmin<-seq(1,length(mylist),5)[which.min(pvals)]#take all p-values and find the index of the smallest
                                                       #(first smallest in case of ties)   
        
        pval<-mylist[[jmin]]
        
        #########################in case pvalue = 0, there might be lack of uniqueness      
        
        if(pval==0)#if the smallest = 0, check if there is more than one equal to 0
        {
          j0<-seq(1,length(mylist),5)[pvals==0]#find all equal to 0; there is at least 1 due to if(pval==0)
          
          pval<-mylist[[j0[1]]]
          J<-mylist[[j0[1]+1]]
          stat<-mylist[[j0[1]+2]]
          substat<-mylist[[j0[1]+3]]
          bopt<-mylist[[j0[1]+4]]
          
          j0winner<-j0[1]
          j0dist<-abs(stat-max(substat))
          
          if(length(j0)>1)
          #if there are more than one p-values equal to 0, 
          #choose the one where the stat is farthest away from the subs distribution
          {
            for(ii in 2:length(j0))
            {
              pval<-mylist[[j0[ii]]]
              J<-mylist[[j0[ii]+1]]
              stat<-mylist[[j0[ii]+2]]
              substat<-mylist[[j0[ii]+3]]
              bopt<-mylist[[j0[ii]+4]]
              if(abs(stat-max(substat))>j0dist) 
              {j0winner<-j0[ii];j0dist<-abs(stat-max(substat))}#update the current winner
            }#end of ii loop  
            
          }#end of if length(j0)>1  
          jmin<-j0winner                
        }
        ####################################
        
        pval<-mylist[[jmin]]
        Jj<-mylist[[jmin+1]]
        stat<-mylist[[jmin+2]]
        substat<-mylist[[jmin+3]]
        bopt<-mylist[[jmin+4]]
        
        #overwrite mylist and return it
        
        mylist<-list(pval=pval,J=Jj,stat=stat,substat=substat,bopt=bopt)

        return(mylist)
      } else {j<-j-1}#if p-value is large, reduce j and continue
    }#end  while(j>=1)  
  }#end if(pval>nominal)
  
  return(list(pval=pval,J=J,stat=stat,substat=substat,bopt=bopt))
}
################################################################################
Gomez2013SeasonalPlotting<-function(x,k=6,c=0.11,b=0.55,s=12) 
#function that implements Stage1 of the current article  
#and returns a ggplot2 object and some other useful info for plotting
#input:
#x: data (vector)
#k,c,b: tuning constants recommended by Gomez 2013
#s: the number of periods to consider
#if s<2, no seasonal frequencies are considered, only trend frequency;
#for monthly data s=12, the frequencies are 2*pi*j/12, j=1,2,...,5, 
#frequency pi (6th frequency) and frequency 0 (7th frequency)
#output: list with several slots; freq gives the frequencies considered, 
#poly the corresponding differencing polynomials (NULL indicates that such polynomial should be omitted)
#plus some additional info
{
  n<-length(x)
  
  if(s>=2)
  {  
    p<-ifelse(s%%2==0,s/2-1,(s-1)/2)
  } else p<-0
  
  #polynomials corresponding to frequencies: sesonal, pi, zero
  poly<-vector('list',length=2+p)#for collecting filters; NULL value corresponds to no uroot
  
  order<-k+(2*p)#order of AR
  ar.fit <- ar.ols(x=x,order=order,aic=FALSE)
  ar.roots <- polyroot(c(1,-1*ar.fit$ar))#AR poly as in (original) Gomez 2013
  
  inv.ar.roots<-1/ar.roots
  re.inv.ar.roots<-Re(inv.ar.roots)
  im.inv.ar.roots<-Im(inv.ar.roots)
  
  alphan<-0.5-1/n#for Step1 of (original) Gomez 2013
  hn<-1/n^alphan
  
  #zero frequency 
  
  #Step 1 of (original) Gomez 2013
  #freq=0
  temp<-(re.inv.ar.roots>1-hn)&(im.inv.ar.roots<hn)#temp = TRUE or FALSE
  uindex<-(1:order)[temp]#index of the uroot at freq 0 among the AR roots
  
  unitroot<-NULL
  #Step2 of (original) Gomez 2013
  if(length(uindex)==0)#if no uroot at freq=0, do another check in this step
  {
    arma.fit<-hannan(x=x,p=1,q=1)
    ph<- -arma.fit$phi
    th<- arma.fit$theta
    lambda<- -(ph)
    
    betan<-0.5-1/n^b
    hn2<-1/n^betan
    
    temp<-((lambda>1-hn2)&(abs(ph-th)>c)) #TRUE or FALSE
    if(temp) poly[[p+2]]<-c(1,-1)#put zero-freq polynomial after seasonal ones and pi
  }  else {
    poly[[p+2]]<-c(1,-1)
    uindex<-uindex[1]#in case there are two 
    unitroot<-c(unitroot,re.inv.ar.roots[uindex]+1i*im.inv.ar.roots[uindex])
  } 
  
  
  freq<-NULL#for collecting all frequencies
  unitroottemp<-NULL#for collecting only unitroot values
  
  #seasonal frequencies, adaptation of Gomez 2013 
  if(s>=2)
  {  
    omega<-c(2*pi*(1:p)/s)
    
    arg.inv.ar.roots<-Arg(inv.ar.roots)
    mod.inv.ar.roots<-Mod(inv.ar.roots)
    
    #drop the freq=0 and do pi separately below
    zeroindex<-uindex#uindex comes from Step1&2 above
    
    if(length(zeroindex)>0)#only if zero freq is relevant
    {
      mod.inv.ar.roots<-mod.inv.ar.roots[-c(zeroindex)]    
      arg.inv.ar.roots<-arg.inv.ar.roots[-c(zeroindex)]    
    }
    
    order<-length(mod.inv.ar.roots)#update order
    
    for(j in 1:p)
    {
      #drop the freq=0 frequency
      tempp<-(mod.inv.ar.roots>1-hn)&(abs(arg.inv.ar.roots-omega[j])<hn)#large Mod and close to omega_j
      tempn<-(mod.inv.ar.roots>1-hn)&(abs(omega[j]+arg.inv.ar.roots)<hn)#large Mod and close to -omega_j
      
      upindex<-(1:order)[tempp]
      unindex<-(1:order)[tempn]
      
      if((length(upindex)==1)&(length(unindex==1)))#identify indices i and i'
      {
        #use slot j as we put pi and zero frequencies at the end
        poly[[j]]<-c(1,-2*cos(omega[j]),1)
        
        z<-mod.inv.ar.roots[upindex]*exp(1i*arg.inv.ar.roots[upindex])
        zstar<-Conj(z)
        unitroottemp<-c(unitroottemp,c(z,zstar))
      }
      #(1-rho*exp(i*w))*(1-rho*exp(-i*w)), so rho*exp(i*w)=inverse root (with +w),
      #take mod of inverse root to obtain rho, so the filter is 1-2*rho*cos(omega)*B+rho^2*B^2
    }#end of j loop  
    
    
    piindex<-(1:order)[(mod.inv.ar.roots>1-hn)&abs(abs(arg.inv.ar.roots)-pi)<hn]
    
    if(length(piindex)>0)#only if pi frequency is relevant
    {
      piindex<-piindex[1]#in case there are two 
      
      poly[[p+1]]<-c(1,1)#1+z
      
      zpi<-mod.inv.ar.roots[piindex]*exp(1i*arg.inv.ar.roots[piindex])
      unitroottemp<-c(unitroottemp,zpi)#update uniroottemp
    }  
    
    freq<-c(freq,omega,pi)#update freq (initially, NULL) with seasonal freq, incl.pi
    
  }#end of if for seasonal
  
  #combine seasonal uroots and trend uroot and update
  #combine seasonal frequencies and trend frequency and update
  unitroottemp<-c(unitroottemp,unitroot)
  freq<-c(freq,0)

  #auxiliary variables for plotting  
  xvals <- seq(-1000,1000)/1000
  yvals <- sqrt(1-xvals^2)
  df0<-data.frame(x=rep(xvals,2),y=c(yvals,-yvals))
  
  df<-data.frame(x=re.inv.ar.roots,y=im.inv.ar.roots)
  
  dfunitroot<-data.frame(x=Re(unitroottemp),y=Im(unitroottemp))
  
  pp<-ggplot(df,aes(x=x,y=y))+geom_point()
  pp<-pp+geom_point(data=dfunitroot,aes(x=x,y=y),col="red",shape=1,size=5)
  pp<-pp+geom_point(data=df0,aes(x=x,y=y),col="black",size=0.1)+theme_bw()+coord_cartesian(xlim=c(-2,2),ylim=c(-2,2))

  list(freq=freq,poly=poly,pp=pp,unitroot=unitroottemp,df=df,df0=df0)
}  

Stage1and2<-function(x,fixedb)
#functions that carries out Stage1&2 of the current article 
#input:
#x: [num] sample data (vector of length n)
#fixedb: [num] fixed-b fraction (can be a vector of fractions)  
#output: matrix whose columns contain:
#   [1] Current procedure number of frequency zero roots
#   [2] Current procedure number of frequency pi/6 roots
#   [3] Current procedure number of frequency 2*pi/6 roots
#   [4] Current procedure number of frequency 3*pi/6 roots
#   [5] Current procedure number of frequency 4*pi/6 roots
#   [6] Current procedure number of frequency 5*pi/6 roots
#   [7] Current procedure number of frequency 6*pi/6 roots
#   [8] Current procedure number of frequency 6*pi/6 roots
#   [9] fixed-b fraction used
#   [10] subsample size used 
{
  output<-Gomez2013SeasonalPlotting(x=x,k=6,c=0.11,b=0.55,s=12)#checks zero and all seasonal frequencies

  temp<-!unlist(lapply(output$poly,is.null))
  lambda0<-output$freq[temp]#take only non-NULL frequencies (extracted using 'temp')
  delta0<-output$poly[temp]#here only non-NULL ones are extracted (using 'temp')
  candidates<-(1:7)[temp]#six seasonal ones and a zero-freq, in this order
  
  n<-length(x)
  
  #a sequence of subsample sizes for the adaptive rule
  BASE<-0.75
  LOB<-0.03
  UPB<-0.20
  j.min<-ceiling( log(c(UPB),base=BASE) )
  j.max<-floor( log(c(LOB),base=BASE) )
  
  b.values<-rev(round(BASE^(j.min:j.max)*n))
  b.values<-unique(b.values)
  b.values<-b.values[b.values>=5]# at least 5 observations in the subsample
  
  b.values<-tail(b.values,4)[1]#take the 4th subsample from the end (or disable this line)
  
  fixedb.values<-fixedb
  
  output.matrix<-NULL

  for(b in 1:length(b.values))
  {  
    for(f in 1:length(fixedb.values))
    {  
      
      out<-myprocedureBackward(filters=delta0,
                               lambdas=lambda0,x=x,b=b.values[b],
                               nominal=0.05,param=0,pow=1,
                               myest=mySchurest,fixedb=fixedb.values[f],
                               taper=Parzen)
      
        k<-1
        pval<-round(out[[k]],4)#k because the p-value is first
        J<-out[[k+1]]
        stat<-out[[k+2]]
        substat<-out[[k+3]]
        bopt<-out[[k+4]]
        
        selected<-integer(0)
        if(length(J)>0) selected<-candidates[J]
        outputtemp<-1*(c(7,1:6)%in%selected)
        outputtemp<-c(outputtemp,b.values[b],fixedb.values[f])#add fixedb and subsize to the info
        output.matrix<-rbind(output.matrix,outputtemp)
        
    }#end of f loop over fixedb values loop
  }#end of the loop over subsamaple values
  return(output.matrix)
}  

