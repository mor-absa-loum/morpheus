prms <- function(name, idx)
{
  load(name)
  d <- nrow(mr[[1]][[1]])-2
  p <- colMeans(do.call(rbind, lapply(mr[[idx]], function(m) m[1,])))
  b <- colMeans(do.call(rbind, lapply(mr[[idx]], function(m) m[2+d,])))
  L <- length(mr[[1]])
  beta <- (1/L) * Reduce("+", lapply(mr[[idx]], function(m) m[2:(d+1),]))
  list(p, beta, b, mr_params)
}

pprms <- function(link)
{
  for (n in c("5000", "10000", "100000", "500000", "1000000"))
  {
    method  =1
    #for (method in 1:2)
    #{
      toprint <- c()
      for (d in c(2,5,10))
      {
        name <- paste0("res_", n, "_", d, "_", link, "_6,3,1.RData")
        params <- prms(name, method)
        toprint <- c(toprint, c(
          sum(abs(params[[1]] - params[[4]]$p)),
          colSums(abs(params[[2]] - params[[4]]$beta)),
          sum(abs(params[[3]] - params[[4]]$b))
        ))
      }
      print(toprint, digits=2)
    #}
  }
}
