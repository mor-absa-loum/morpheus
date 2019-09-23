optimBeta <- function(N, n, K, p, beta, b, link, ncores)
{
	library(morpheus)
	res <- multiRun(
		list(n=n, p=p, beta=beta, b=b, optargs=list(K=K, link=link)),
		list(
			# morpheus
			function(fargs) {
				library(morpheus)
				K <- fargs$optargs$K
				M <- computeMoments(fargs$X, fargs$Y)
				fargs$optargs$M <- M
				mu <- computeMu(fargs$X, fargs$Y, fargs$optargs)
				res2 <- NULL
				tryCatch({
					op <- optimParams(K,fargs$optargs$link,fargs$optargs)
					x_init <- list(p=rep(1/K,K-1), beta=mu, b=rep(0,K))
					res2 <- do.call(rbind, op$run(x_init))
				}, error = function(e) {
					res2 <- NA
				})
				res2
			}
#			,
#			# flexmix
#			function(fargs) {
#				library(flexmix)
#				source("../patch_Bettina/FLXMRglm.R")
#				K <- fargs$optargs$K
#				dat <- as.data.frame( cbind(fargs$Y,fargs$X) )
#				res2 <- NULL
#				tryCatch({
#					fm <- flexmix( cbind(V1, 1-V1) ~ .-V1, data=dat, k=K,
#						model = FLXMRglm(family = binomial(link = link)) )
#					p <- mean(fm@posterior[["scaled"]][,1])
#					out <- refit(fm)
#					beta_b <- sapply( seq_len(K), function(i) {
#						as.double( out@components[[1]][[i]][,1] )
#					} )
#					res2 <- rbind(p, beta_b[2:nrow(beta_b),], beta_b[1,])
#				}, error = function(e) {
#					res2 <- NA
#				})
#				res2
#			}
		),
		prepareArgs = function(fargs, index) {
			library(morpheus)
			io = generateSampleIO(fargs$n, fargs$p, fargs$beta, fargs$b, fargs$optargs$link)
			fargs$X = io$X
			fargs$Y = io$Y
			fargs$optargs$K = ncol(fargs$beta)
			fargs$optargs$M = computeMoments(io$X,io$Y)
			fargs
		}, N=N, ncores=ncores, verbose=TRUE)
	p <- c(p, 1-sum(p))
	for (i in 1:length(res)) {
		for (j in N:1) {
			if (is.null(res[[i]][[j]]) || is.na(res[[i]][[j]]))
				res[[i]][[j]] <- NULL
		}
		print(paste("Count valid runs for ",i," = ",length(res[[i]]),sep=""))
		res[[i]] <- alignMatrices(res[[i]], ref=rbind(p,beta,b), ls_mode="exact")
	}
	res
}

#model = binomial; default values:
link = "logit"
N <- 10
d <- 2
n <- 1e4
ncores <- 1

cmd_args <- commandArgs()
for (arg in cmd_args)
{
	if (substr(arg,1,1)!='-') {
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
		}
	}
}

if (d == 2) {
	K <- 2
	p <- .5
	b <- c(-.2, .5)
	beta <- matrix( c(1,-2, 3,1), ncol=K )
} else if (d == 5) {
	K <- 2
	p <- .5
	b <- c(-.2, .5)
	beta <- matrix( c(1,2,-1,0,3, 2,-3,0,1,0), ncol=K )
} else if (d == 10) {
	K <- 3
	p <- c(.3, .3)
	b <- c(-.2, 0, .5)
	beta <- matrix( c(1,2,-1,0,3,4,-1,-3,0,2, 2,-3,0,1,0,-1,-4,3,2,0, -1,1,3,-1,0,0,2,0,1,-2), ncol=K )
} else if (d == 20) {
	K <- 3
	p <- c(.3, .3)
	b <- c(-.2, 0, .5)
	beta <- matrix( c(1,2,-1,0,3,4,-1,-3,0,2,2,-3,0,1,0,-1,-4,3,2,0, -1,1,3,-1,0,0,2,0,1,-2,1,2,-1,0,3,4,-1,-3,0,2, 2,-3,0,1,0,-1,-4,3,2,0,1,1,2,2,-2,-2,3,1,0,0), ncol=K )
}

mr <- optimBeta(N, n, K, p, beta, b, link, ncores)
mr_params <- list("N"=N, "n"=n, "K"=K, "d"=d, "link"=link,
	"p"=c(p,1-sum(p)), "beta"=beta, "b"=b)

save("mr", "mr_params", file=paste("multirun_",n,"_",d,"_",link,".RData",sep=""))
