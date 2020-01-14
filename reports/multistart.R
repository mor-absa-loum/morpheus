library(morpheus)

testMultistart <- function(N, n, p, beta, b, link, nstart, ncores)
{
  res <- multiRun(
    list(n=n, p=p, beta=beta, b=b, link=link, nstart=nstart),
    list(
      function(fargs) {
        # 1 start
        library(morpheus)
        K <- ncol(fargs$beta)
        mu <- computeMu(fargs$X, fargs$Y, list(K=K, M=fargs$M))
        op <- optimParams(fargs$X, fargs$Y, K, fargs$link, fargs$M, 1)
        x_init <- list(p=rep(1/K,K-1), beta=mu, b=rep(0,K))
				res2 <- NULL
				tryCatch({
          res2 <- do.call(rbind, op$run(x_init))
				}, error = function(e) {})
				res2
      },
      function(fargs) {
        # B starts
        library(morpheus)
        K <- ncol(fargs$beta)
				d <- nrow(fargs$beta)
        op <- optimParams(fargs$X, fargs$Y, K, fargs$link, fargs$M)
        best_val <- Inf
        best_par <- list()
        for (i in 1:fargs$nstart)
        {
          #x_init <- list(p=rep(1/K,K-1), beta=i*fargs$mu, b=rep(0,K))
          M <- matrix(rnorm(d*K), nrow=d, ncol=K)
          M <- t(t(M) / sqrt(colSums(M^2)))
          x_init <- list(p=rep(1/K,K-1), beta=M, b=rep(0,K))
          par <- NULL
					tryCatch({
            par <- op$run(x_init)
          }, error = function(e) {})
          if (!is.null(par))
          {
            val <- op$f( op$linArgs(par) )
            if (val < best_val)
            {
              best_par <- par
              best_val <- val
            }
          }
        }
        do.call(rbind,best_par) #return NULL on empty list
      }
    ),
    prepareArgs = function(fargs, index) {
      library(morpheus)
      io = generateSampleIO(fargs$n, fargs$p, fargs$beta, fargs$b, fargs$link)
      fargs$M <- computeMoments(io$X, io$Y)
      fargs$X <- io$X
      fargs$Y <- io$Y
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

# Default values:
link = "logit"
N <- 10
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

if (d == 2) {
  p <- .5
  b <- c(-.2, .5)
  beta <- matrix( c(1,-2, 3,1), ncol=2 )
} else if (d == 5) {
  p <- .5
  b <- c(-.2, .5)
  beta <- matrix( c(1,2,-1,0,3, 2,-3,0,1,0), ncol=2 )
} else if (d == 10) {
  p <- c(.3, .3)
  b <- c(-.2, 0, .5)
  beta <- matrix( c(1,2,-1,0,3,4,-1,-3,0,2, 2,-3,0,1,0,-1,-4,3,2,0, -1,1,3,-1,0,0,2,0,1,-2), ncol=3 )
}

mr <- testMultistart(N, n, p, beta, b, link, nstart, ncores)
mr_params <- list("N"=N, "nc"=ncores, "n"=n, "link"=link,
	"p"=c(p,1-sum(p)), "beta"=beta, "b"=b, "nstart"=nstart)

save("mr", "mr_params", file=paste("res_",n,"_",d,"_",link,"_",nstart,".RData",sep=""))
