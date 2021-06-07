optimBeta <- function(N, n, p, beta, b, link, ncores)
{
  library(morpheus)
  res <- multiRun(
    list(n=n, p=p, beta=beta, b=b, link=link),
    list(
      # morpheus
      function(fargs) {
        library(morpheus)
        K <- ncol(fargs$beta)
        M <- computeMoments(fargs$X, fargs$Y)
        mu <- computeMu(fargs$X, fargs$Y, list(K=K, M=M))
        op <- optimParams(fargs$X, fargs$Y, K, fargs$link, M, 1) #only 1 OpenMP core
        x_init <- list(p=rep(1/K,K-1), beta=mu, b=rep(0,K))
        res2 <- NULL
        tryCatch({
          res2 <- do.call(rbind, op$run(x_init))
        }, error = function(e) {})
        res2
      }
			,
      # flexmix
      function(fargs) {
        library(flexmix)
        source("../patch_Bettina/FLXMRglm.R")
        K <- ncol(fargs$beta)
        dat <- as.data.frame( cbind(fargs$Y,fargs$X) )
        res2 <- NULL
        tryCatch({
          fm <- flexmix( cbind(V1, 1-V1) ~ ., data=dat, k=K,
            model = FLXMRglm(family = binomial(link = link)) )
          pf <- colMeans(fm@posterior[["scaled"]])
          out <- refit(fm)
          beta_b <- sapply( seq_len(K), function(i) {
            as.double( out@components[[1]][[i]][,1] )
          } )
          res2 <- rbind(pf, beta_b[2:nrow(beta_b),], beta_b[1,])
        }, error = function(e) {
          res2 <- NA
        })
        res2
      }
    ),
    prepareArgs = function(fargs, index) {
      library(morpheus)
      io = generateSampleIO(fargs$n, fargs$p, fargs$beta, fargs$b, fargs$link)
      fargs$X = io$X
      fargs$Y = io$Y
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

mr <- optimBeta(N, n, p, beta, b, link, ncores)
mr_params <- list("N"=N, "nc"=ncores, "n"=n, "link"=link,
  "p"=c(p,1-sum(p)), "beta"=beta, "b"=b)

save("mr", "mr_params", file=paste("res_",n,"_",d,"_",link,".RData",sep=""))
