nvals <- c(5000,10000,100000,500000,1000000)
for (link in c("logit","probit"))
{
	for (d in c(2,5,10))
	{
		par(mfrow=c(2,2), lwd=2, cex.axis=2, cex=2)
		res <- list()
		for (n in nvals)
		{
			load(paste("multirun_",n,"_",d,"_",link,".RData",sep=""))
			res <- c(res, mr)
		}
		for (i in 1:2) #our, fm
		{
			for (j in 1:2) #beta, p,b
			{
				if (j == 1) #beta
				{
					for (dim in 1:d)
					{
						ypts <- c()
						for (k in 1:5)
							ypts <- c(ypts, mean(sapply(1:length(res[[k]][[i]]), function(x) mean((res[[k]][[i]][[x]][2:(d+1),dim] - mr_params$beta[,dim])^2))))
						plot(nvals, ypts)
						if (dim < d)
							par(new=TRUE)
					}
				}
				else #p + b
				{
					for (rowidx in c(1,d+2))
					{
						ypts <- c()
						ref <- if (rowidx==1) { mr_params$p } else { mr_params$b }
						for (k in 1:5)
							ypts <- c(ypts, mean(sapply(1:length(res[[k]][[i]]), function(x) mean((res[[k]][[i]][[x]][rowidx,] - ref)^2))))
						plot(nvals, ypts)
						if (rowidx==1)
							par(new=TRUE)
					}
				}
			}
		}
	}
}
