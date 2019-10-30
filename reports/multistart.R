library(morpheus)

testMultistart <- function(N, n, d, K, p, beta, b, link, nstart, ncores)
{
  res <- multiRun(
    list(n=n,p=p,beta=beta,b=b,optargs=list(K=K,d=d,link=link,nstart=nstart)),
    list(
      function(fargs) {
        # 1 start
        library(morpheus)
        K <- fargs$optargs$K
        op <- optimParams(K, fargs$optargs$link, fargs$optargs)
        x_init <- list(p=rep(1/K,K-1), beta=fargs$mu, b=rep(0,K))
				res2 <- NULL
				tryCatch({
          res2 <- do.call(rbind, op$run(x_init))
				}, error = function(e) {})
				res2
      },
      function(fargs) {
        # B starts
        library(morpheus)
        K <- fargs$optargs$K
				d <- fargs$optargs$d
        op <- optimParams(K, fargs$optargs$link, fargs$optargs)
        best_val <- Inf
        best_par <- list()
        for (i in 1:fargs$optargs$nstart)
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
        # Bet that at least one run succeded:
        do.call(rbind,best_par)
      }
    ),
    prepareArgs = function(fargs, index) {
      library(morpheus)
      io = generateSampleIO(fargs$n, fargs$p, fargs$beta, fargs$b, fargs$optargs$link)
      fargs$optargs$M <- computeMoments(io$X, io$Y)
      mu <- computeMu(io$X, io$Y, fargs$optargs)
      fargs$mu <- mu
			fargs
    }, N=N, ncores=ncores, verbose=TRUE)
  for (i in 1:2)
    res[[i]] <- alignMatrices(ms[[i]], ref=rbind(p,beta,b), ls_mode="exact")
  res
}

#model = binomial
K <- 2
p <- .5
b <- c(-.2, .5)
# Default values:
link = "logit"
N <- 10
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

mr <- testMultistart(N, n, d, K, p, beta, b, link, nstart, ncores)
mr_params <- list("N"=N, "nc"=ncores, "n"=n, "K"=K, "d"=d, "link"=link,
	"p"=c(p,1-sum(p)), "beta"=beta, "b"=b, "nstart"=nstart)

save("mr", "mr_params", file=paste("res_",n,"_",d,"_",link,"_",nstart,".RData",sep=""))
