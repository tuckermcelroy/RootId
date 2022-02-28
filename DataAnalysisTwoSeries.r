#Rscript for peforming data analysis of the article/project
#"Identification of the Differencing Operator of a Non-stationary
#Time Series via Testing for Zeroes in the Spectral Density"

#Please, save this file and all the other files 
#of the project in a folder called RcodeAndData and 
#set working directory to the folder '.../RcodeAndData/', ie, the location of this file

#setwd('.../RcodeAndData/')#modify, uncomment, run

rm(list=ls())

Sys.time()

#load packages
library('data.table')#data manipulation
library('seasonal')#TRAMO

checkX13()

source("Functions.R")#Rfunctions for the current article

#Functions for comparison with:
#McElroy, T. (2021). A diagnostic for seasonality based upon polynomial roots of
#ARMA models. Journal of Official Statistics 37 (2), 1-28.
source("ARMAauto.r")
source("flipIt.r")
source("seas.proxim.r")
source("root.sigtest.ci.r")

#auxiliary function to remove NAs in data processing
strip.na <- function(x)
{
  n <- length(x)
  period <- 12
  inds <- seq(1,n)[is.na(x)==FALSE]
  if(max(diff(inds)) == 1) 
  { 
    start.year <- floor(time(x)[min(inds)])
    start.season <- 1 + period*(time(x)[min(inds)] - start.year)
    start.date <- c(start.year,start.season)
    x <- ts(x[min(inds):max(inds)],start=start.date,frequency=period)
  }
  return(x)
}

#open a graphical window
par(mfrow=c(2,3),mar=c(2,2,3,2))

#load data
x.ts<-read.table(file='Series1Series2.txt',header=TRUE)#'Series1Series2.txt' from the current folder
x.ts <- ts(x.ts,start=1959,frequency=12)

#################applying the method of the current article (ZT=Stage1and2 and MG=Stage1)
x.delta<-NULL
for(k in 1:ncol(x.ts)) 
{
  x <- x.ts[,k]
  x <- strip.na(x)
  T <- length(x)

  m <- seas(x = x,transform.function = "auto",x11 = "")
  z <- series(m,"b1")
  if(udg(m)$aictrans == "Log(y)") { z <- log(z) }
  #x is the original time series
  #z is the 'linearized' time series
  
  idnon<- udg(m)$idnonseasonaldiff
  ids <- udg(m)$idseasonaldiff

  MAIN<-ifelse(k==1,"SERIES 1","SERIES 2")

  plot.ts(as.vector(x),main=paste(MAIN));grid()
  acf(as.vector(z),lag.max=100,main="")

  #MG method
  output<-Gomez2013SeasonalPlotting(x=z,k=6,c=0.11,b=0.55,s=12)#checks zero and all seasonal frequencies
  print(output$poly)#MG method info

  plot(output$df0$x,output$df0$y,type='l',xlab='Re',ylab='Im')
  points(output$df$x,output$df$y,cex=1,pch=20)
  
  z<-as.vector(z)
  #ZT method
  outputdelta<-Stage1and2(x=z,fixedb=0.5)#returns a vector with 7 entries containing 0s and 1s 
  
  temp<-cbind(idnon,ids,outputdelta,k)
  #print(temp)
  
  x.delta<-rbind(x.delta,temp)
}

#x.delta contains (in columns):
#   [1] TRAMO number of nonseasonal differences
#   [2] TRAMO number of seasonal differences
#   [3] Current procedure number of frequency zero roots
#   [4] Current procedure number of frequency pi/6 roots
#   [5] Current procedure number of frequency 2*pi/6 roots
#   [6] Current procedure number of frequency 3*pi/6 roots
#   [7] Current procedure number of frequency 4*pi/6 roots
#   [8] Current procedure number of frequency 5*pi/6 roots
#   [9] Current procedure number of frequency 6*pi/6 roots
#   [9] Current procedure number of frequency 6*pi/6 roots
#   [10] subsample size used
#   [11] fixed-b fraction used
#   [12] index of a series (1,2,...)

print(x.delta)

#########################################comparison with 
#McElroy, T. (2021). A diagnostic for seasonality based upon polynomial roots of
#ARMA models. Journal of Offical Statistics 37 (2), 1-28.

RNGkind(sample.kind = 'Rejection')

dfrho<-NULL

for(k in 1:ncol(x.ts)) 
{
  print(Sys.time())
  print(k)
  x <- x.ts[,k]
  x <- strip.na(x)
  T <- length(x)
  
  m <- seas(x = x,transform.function = "auto",x11 = "")
  z <- series(m,"b1")
  if(udg(m)$aictrans == "Log(y)") { z <- log(z) }
  z<-as.numeric(z)
  #x is the original time series
  #z is the 'linearized' time series
  
 
  if(k==1)
  {  
    freq.values<-vector('list',length=8)
    for(f in 1:7)
      freq.values[[f]]<-f-1
    
    freq.values[[8]]<-0:6#several frequencies based on MG
  }
  
  if(k==2)
  {  
    freq.values<-vector('list',length=8)
    for(f in 1:7)
      freq.values[[f]]<-f-1
    
    freq.values[[8]]<-c(0,1,2,4,5,6)#several frequencies based on MG
  }
  
  
  rho.values<-seq(0.85,0.999,0.001)
  
  
  for(f in 1:length(freq.values))
  {
    set.seed(123)#to control the randomness in Monte Carlo
    myfreq<-freq.values[[f]]*pi/6
    rd<-root.sigtest.ci(z,NULL,0,myfreq,monte=10000,rho.values)
    dftemp<-data.frame(series=k,freq=toString(freq.values[[f]]),rho=rho.values,pvalue=rd)
    dfrho<-rbind(dfrho,dftemp)
  }
}#end of loop that goes over the ts

dt<-data.table(dfrho)

dt01<-dt[ ,.(cutoffrho=rho[which(pvalue>0.01)[1]]),by=.(freq,series)]
dt05<-dt[ ,.(cutoffrho=rho[which(pvalue>0.05)[1]]),by=.(freq,series)]

print(dt01)
print(dt05)

#################forecasting and computing (squared) forecast errors
H<-24#forecast horizon
#MYFILTERS is a list with 9 slots containing 
#trend&seasonal polynomials corresponding to 
#frequencies described in x.delta above, respectively

MYFILTERS<-vector('list',length=9)
MYFILTERS[[1]]<-c(1,-1)#(1-z)
MYFILTERS[[2]]<-c(1,rep(0,11),-1)#(1-z^12)

MYFILTERS[[3]]<-c(1,-1)#(1-z)
for(j in 1:5)
{
  MYFILTERS[[3+j]]<-c(1,-2*cos(2*pi*j/12),1)#1-2cos(2*pi*j/12)B+1*B^2, j=1,2,...,5
}  
MYFILTERS[[9]]<-c(1,1)#(1+z), frequency pi

dfSSE<-NULL#for collecting sum of squared forecast errors 

for(k in 1:ncol(x.ts))
{
  x <- x.ts[,k]
  x <- strip.na(x)
  T <- length(x)
  
  m <- seas(x = x,transform.function = "auto",x11 = "")
  z <- series(m,"b1")
  if(udg(m)$aictrans == "Log(y)") { z <- log(z) }
  z<-as.numeric(z)
    
  ####################################################filters and filtered data y1,y2,y3

  #Stage1 (MG method); each series is treated separately
  tempfilter1<-1
  if(k==1)#Series1
    for(i in 3:9)#all frequencies are identified (columns 3-9) in Stage1 (MG method)
        tempfilter1<-polymult(tempfilter1,MYFILTERS[[i]])
  
  if(k==2)#Series2
    for(i in c(3:5,7:9) )#all except for j=3 frequency are identified(columns 3-5,7-9) in Stage1 (MG method)
      tempfilter1<-polymult(tempfilter1,MYFILTERS[[i]])
  
  y <- stats::filter(x=z,tempfilter1,method="convolution",sides=1)
  y1<-strip.na(y)
  y1<-as.numeric(y1)

  acf(y1,lag.max=100,main="STAGE 1")

  #Stage1and2 (ZT method); the output from Stage1and2 is stored in Rvariable 'x.delta' (obtained above);  
  #both series (1 and 2) are analyzed here, ie, for k=1 and 2, tempfilter2 is computed accordingly
  tempfilter2<-1
  for(i in 3:9)#columns 3-9
  {  
    if(x.delta[k,i]>0)#x.delta[k,i] will be 0 or 1; we use only those with 1's
    {  
      tempfilter2<-polymult(tempfilter2,MYFILTERS[[i]])
    }
  }

  y <- stats::filter(x=z,tempfilter2,method="convolution",sides=1)
  y2<-strip.na(y)
  y2<-as.numeric(y2)

  acf(y2,lag.max=100,main="STAGE 1 and STAGE 2")
  
  #ROOT method (McElroy, T. (2021));
  #both series (1 and 2) are analyzed here, ie, for k=1 and 2, tempfilter3 is the same (includes all frequencies)
  tempfilter3<-1
  for(i in 3:9)
  {
      tempfilter3<-polymult(tempfilter3,MYFILTERS[[i]])
  }
  
  y <- stats::filter(x=z,tempfilter3,method="convolution",sides=1)
  y3<-strip.na(y)
  y3<-as.numeric(y3)
  
  #acf(y3,lag.max=100,main="ROOT")#no plotting

  ####################################################forecasting
  
  method<-c('Stage1','Stage1and2','Root')
  
  for(i in 1:3)
  { 

    text<-paste0('y<-y',i)#to y assign y1 (MG=Stage1), y2 (ZT=Stage1and2) or y3 (Root=McElroy (2021))
    eval(parse(text=text))
    
    #to tempfilter assign tempfilter1 (MG=Stage1), tempfilter2 (ZT=Stage1and2) or tempfilter3 (Root=McElroy (2021))   
    text<-paste0('tempfilter<-tempfilter',i)
    eval(parse(text=text))
  
    ytrain<-y[1:(length(y)-H)]#training data (part of differenced data); leave the last H obs. for forecasting
  
    method.fit <- ar.ols(ytrain)
  
    # get pseudo-AR polynomials for forecasting
    phi.method <- c(1,-1*method.fit$ar[1:method.fit$order])
    phi.method <- polymult(phi.method,tempfilter)
  
    yhat<-z[1:(length(z)-H)]#z is undifferenced ('linearized') data
    for(h in 1:H)
    {
      forecast <- -sum(phi.method[-1]*yhat[length(yhat):(length(yhat)-length(phi.method)+2)])
      yhat <- c(yhat,forecast)
    }
      
    #forecast errors 
    forecasterrors<-tail(z-yhat,H)#z is undifferenced ('linearized') data
  
    #plot squared forecast errors
    if(i==1)
    {plot(0,0,type='n',ylim=c(0,0.1),xlim=c(1,H),xlab='Forecast horizon',ylab='Squared forecast error');grid();
      abline(h=0);
    lines(1:H,forecasterrors^2,col='green')} else {
           if(i==2) lines(1:H,forecasterrors^2,col='red') else lines(1:H,forecasterrors^2,col='blue')}
  
    #sum of squared forecast errors
    SSE<-round(sum(forecasterrors^2),6)
    
    dftemp<-data.frame(series=k,
                       method=paste(method[i]),
                       SumSqErrors=SSE)
    
    dfSSE<-rbind(dfSSE,dftemp)

  }#end of i loop over the methods
}#end of k loop over the series  

print(dfSSE)




