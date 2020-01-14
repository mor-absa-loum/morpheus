# flexmix optimization to get beta
fmOptim <- function(X, Y, K, link)
{
	dat <- as.data.frame( cbind(Y,X) )
	fm <- flexmix( cbind(Y, 1-Y) ~ .-Y, data=dat, k=K,
		model = FLXMRglm(family = binomial(link = link)) )
	p <- mean(fm@posterior[["scaled"]][,1])
	out <- refit(fm)
	beta_b <- sapply( seq_len(K), function(i) as.double( out@components[[1]][[i]][,1] ) )
	list("p"=p, "beta"=beta_b[2:nrow(beta_b),], "b"=beta_b[1,])
	NULL
}

# Our package optimization for beta (using mu as a starting point)
ourOptim <- function(X, Y, K, link)
{
	M <- computeMoments(X, Y)
	mu <- computeMu(X, Y, list(K=K,M=M))
	x_init = list(p=rep(1/K,K-1), beta=mu, b=rep(0,K))
	optimParams(X, Y, K, link, M, 1)$run(x_init)
	NULL
}

# Get timings for both methods with the same beta matrix
getTimings <- function(link)
{
	timings <- list('fm'=matrix(0,nrow=10,ncol=7),'our'=matrix(0,nrow=10,ncol=7))
	K <- 2
	for (d in c(2,5,10))
	{
		beta <- matrix(runif(d*K,min=-5,max=5),ncol=K)
		for (logn in 4:6)
		{
			n <- 10^logn
			io <- generateSampleIO(n, rep(1/K,K-1), beta, runif(K), link)
			timings[['fm']][d,logn] <- system.time(fmOptim(io$X,io$Y,K,link))[3]
			timings[['our']][d,logn] <- system.time(ourOptim(io$X,io$Y,K,link))[3]
		}
	}
	timings
}

#model = binomial
link <- "logit"
ncores <- 1
N <- 100

cmd_args <- commandArgs()
for (arg in cmd_args)
{
	if (substr(arg,1,1)!='-')
	{
		spl <- strsplit(arg,'=')[[1]]
		if (spl[1] == "link") {
			link <- spl[2]
		} else if (spl[1] == "nc") {
			ncores <- as.integer(spl[2])
		} else if (spl[1] == "N") {
			N <- as.integer(spl[2])
		}
	}
}

library(morpheus)
library(flexmix)
source("../patch_Bettina/FLXMRglm.R")

tm <-
	if (ncores == 1) {
		lapply(1:N, function(i) {
			print(paste("Run",i))
			getTimings(link)
		})
	} else {
		library(parallel)
		mclapply(1:N, function(i) {
			print(paste("Run",i))
			getTimings(link)
		},
		mc.preschedule=FALSE, mc.cores=ncores)
	}
tm_params <- list("link"=link, "N"=N, "nc"=ncores)

save("tm", "tm_params", file="timings.RData")
