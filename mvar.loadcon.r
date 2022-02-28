mvar.loadcon <- function(data.mat,type,geo,sa)
{
	# reads in data.mat in format from Census Bureau
	#  for construction, and returns extracted time series
	#  based on user choices.
	# 
	# type: integer 1 through 3
	# sa: 0 for non-seasonally adjusted, 1 for seasonally adjusted

n.all <- dim(data.mat)[1]
n.time <- max(data.mat[,1])
inds.series <- unique(data.mat[,2])
n.cross <- length(inds.series)

inds.units <- seq(1,n.all)[data.mat[,3]==type]
#length(inds.units)
inds.values <- seq(1,n.all)[data.mat[,4]==0]
#length(inds.values)
inds.geo <- seq(1,n.all)[data.mat[,5]==geo]
#length(inds.geo)
inds.nsa <- seq(1,n.all)[data.mat[,6]==sa]
#length(inds.nsa)
inds.database <- intersect(intersect(intersect(inds.units,inds.values),inds.nsa),inds.geo)
#length(inds.database)

data.ts <- NULL
if(length(inds.database) > 0) {
data.ts <- matrix(NA,nrow=n.time,ncol=n.cross)
series.na <- NULL
for(k in 1:n.cross)
{
	inds.cross <- seq(1,n.all)[data.mat[,2]==inds.series[k]]
	inds.cross <- intersect(inds.cross,inds.database)
	if(length(inds.cross)==0) { series.na <- c(series.na,k) }
	# get time points where data is non-missing
	inds.times <- data.mat[inds.cross,1]
	data.ts[inds.times,k] <- data.mat[inds.cross,7]
} 
if(length(series.na)>0) data.ts <- data.ts[,-series.na]
}

return(data.ts)
}


