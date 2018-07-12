library(morpheus)

#model = binomial
K <- 2
p <- .5
b <- c(-.2, .5)
# Default values:
link = "logit"
N <- 100
d <- 2
n <- 1e4
ncores <- 1
nstart <- 3 #nstart-1 random starting points for each MC run

cmd_args <- commandArgs()
for (arg in cmd_args)
{
	if (substr(arg,1,1)!='-')
	{
		spl <- strsplit(arg,'=')[[1]]
		if (spl[1] == "nc") {
			ncores <- as.integer(spl[2])
		} else if (spl[1] == "N") {
			N <- as.integer(spl[2])
		} else if (spl[1] == "n") {
			n <- as.integer(spl[2])
		} else if (spl[1] == "d") {
			d <- as.integer(spl[2])
		} else if (spl[1] == "link") {
			link <- spl[2]
		} else if (spl[1] == "nstart") {
			nstart <- spl[2]
		}
	}
}
betas <- list(
	matrix( c(1,-2, 3,1), ncol=K ), #d=2
	matrix( c(1,2,-1,0,3, 2,-3,0,1,0), ncol=K ), #d=5
	matrix( c(1,2,-1,0,3,4,-1,-3,0,2, 2,-3,0,1,0,-1,-4,3,2,0), ncol=K ) ) #d=10
beta <- betas[[ ifelse( d==2, 1, ifelse(d==5,2,3) ) ]]

ms <- multiRun(
	list(n=n,p=p,beta=beta,b=b,optargs=list(K=K,link=link,nstart=nstart)), list(
		function(fargs) {
			# 1 start
			library(morpheus)
			K <- fargs$optargs$K
			op <- optimParams(K, fargs$optargs$link, fargs$optargs)
			x_init <- c(rep(1/K,K-1), as.double(fargs$mu), rep(0,K))
			do.call(rbind,op$run(x_init))
		},
		function(fargs) {
			# B starts
			library(morpheus)
			K <- fargs$optargs$K
			op <- optimParams(K, fargs$optargs$link, fargs$optargs)
			best_val <- Inf
			best_par <- list()
			for (i in 1:fargs$optargs$nstart)
			{
				x_init <- c(rep(1/K,K-1), as.double(i*fargs$mu), rep(0,K))
				par <- op$run(x_init)
				val <- op$f( op$linArgs(par) )
				if (val < best_val)
				{
					best_par <- par
					best_val <- val
				}
			}
			do.call(rbind,best_par)
		}),
		prepareArgs = function(fargs) {
			library(morpheus)
			io = generateSampleIO(fargs$n, fargs$p, fargs$beta, fargs$b, fargs$optargs$link)
			fargs$optargs$M <- computeMoments(io$X, io$Y)
			mu <- computeMu(io$X, io$Y, fargs$optargs)
			fargs$mu <- mu
		}, N=N, ncores=ncores, verbose=TRUE)
for (i in 1:2)
	ms[[i]] <- alignMatrices(ms[[i]], ref=rbind(p,beta,b), ls_mode="exact")

ms_params <- list("N"=N, "nc"=ncores, "n"=n, "K"=K, "link"=link,
	"p"=p, "beta"=beta, "b"=b, "nstart"=nstart)

save(ms, ms_params, file="multistart.RData")
