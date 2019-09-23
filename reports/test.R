library(morpheus)
morph <- function(fargs) {
	K <- fargs$optargs$K
	M <- computeMoments(fargs$X, fargs$Y)
	fargs$optargs$M <- M
	mu <- computeMu(fargs$X, fargs$Y, fargs$optargs)
	res2 <- NULL
	tryCatch({
		op <- optimParams(K,link,fargs$optargs)
		x_init <- list(p=rep(1/K,K-1), beta=mu, b=rep(0,K))
		res2 <- do.call(rbind, op$run(x_init))
	}, error = function(e) {
		res2 <- NA
	})
	res2
}

#model = binomial; default values:
link = "probit"
N <- 10
d <- 2
n <- 1e4
ncores <- 1

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

fargs = list(n=n, p=p, beta=beta, b=b)
fargs$optargs = list(link=link)

io = generateSampleIO(fargs$n, fargs$p, fargs$beta, fargs$b, fargs$optargs$link)
fargs$X = io$X
fargs$Y = io$Y
fargs$optargs$K = ncol(fargs$beta)
fargs$optargs$M = computeMoments(io$X,io$Y)

res2 <- morph(fargs)

save("res2", file="test.RData")
