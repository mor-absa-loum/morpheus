# NOTE: discard top 2% of highest values
prms <- function(name, idx)
{
  load(name)
  d <- nrow(mr[[1]][[1]])-2
	if (idx > length(mr))
		mr[[idx]] = mr[[1]]
  p <- colMeans(do.call(rbind, lapply(mr[[idx]], function(m) m[1,])))
	bVects <- lapply(mr[[idx]], function(m) m[2+d,])
	q98 <- Inf #quantile(sapply(bVects, function(bv) sum(abs(bv))), 0.98)
	bFiltered <- Filter(function(bv) sum(abs(bv)) < q98, bVects)
  b <- colMeans(do.call(rbind, bFiltered))
	betaMatrices <- lapply(mr[[idx]], function(m) m[2:(d+1),])
	q98 <- Inf #quantile(sapply(betaMatrices, function(bm) sum(abs(bm))), 0.98)
	bmFiltered <- Filter(function(bm) sum(abs(bm)) < q98, betaMatrices)
	beta <- (1/length(bmFiltered)) * Reduce("+", bmFiltered)
  list(p, beta, b, mr_params)
}

pprms <- function(link, prefix="./")
{
  toprint <- matrix(nrow=0, ncol=13) #13=1+2+1 + 1+2+1 + 1+3+1
	for (n in c("5000", "10000", "100000", "500000", "1000000"))
  {
    for (method in 1:2)
    {
			row <- c()
      for (d in c(2,5,10))
      {
        name <- paste0(prefix, "res_", n, "_", d, "_", link, ".RData")
        params <- prms(name, method)
        row <- c( row,
          sum(abs(params[[1]] - params[[4]]$p)),
          colSums(abs(params[[2]] - params[[4]]$beta)),
          sum(abs(params[[3]] - params[[4]]$b)) )
      }
      toprint <- rbind(toprint, row)
    }
  }
	print(formatC(toprint, format="e", digits=1)) #for reporting
	return (toprint)
}
