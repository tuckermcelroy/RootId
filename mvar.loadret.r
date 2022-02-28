mvar.loadret <- function(data.mat,type,sa)
{
	# reads in data.mat in format from Census Bureau
	#  for retail trade, and returns extracted time series
	#  based on user choices.
	# 
	# type: integer 1 through 5
	# sa: 0 for non-seasonally adjusted, 1 for seasonally adjusted

n.all <- dim(data.mat)[1]
n.time <- max(data.mat[,1])
inds.series <- unique(data.mat[,2])
n.cross <- length(inds.series)

inds.sales <- seq(1,n.all)[data.mat[,3]==type]
#length(inds.sales)
inds.values <- seq(1,n.all)[data.mat[,4]==0]
#length(inds.values)
inds.nsa <- seq(1,n.all)[data.mat[,6]==sa]
#length(inds.nsa)
inds.database <- intersect(intersect(inds.sales,inds.values),inds.nsa)
#length(inds.database)

data.ts <- matrix(NA,nrow=n.time,ncol=n.cross)
for(k in 1:n.cross)
{
	inds.cross <- seq(1,n.all)[data.mat[,2]==inds.series[k]]
	inds.cross <- intersect(inds.cross,inds.database)
	# get time points where data is non-missing
	inds.times <- data.mat[inds.cross,1]
	data.ts[inds.times,k] <- data.mat[inds.cross,7]
}

return(data.ts)
}


