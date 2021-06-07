K <- list("2"=2, "5"=2, "10"=3)
p <- list("2"=c(.5,.5), "5"=c(.5,.5), "10"=c(.3, .3, .4))
b <- list("2"=c(-.2, .5), "5"=c(-.2, .5), "10"=c(-.2, 0, .5))
beta <- list(
  "2"=matrix( c(1,-2, 3,1), ncol=2 ),
  "5"=matrix( c(1,2,-1,0,3, 2,-3,0,1,0), ncol=2 ),
  "10"=matrix( c(1,2,-1,0,3,4,-1,-3,0,2, 2,-3,0,1,0,-1,-4,3,2,0, -1,1,3,-1,0,0,2,0,1,-2), ncol=3 ))

for (n in c("5000", "10000", "100000")) {
  for (d in c("2", "5", "10")) {
    load(paste0("res_",n,"_",d,"_logit.RData"))

    # p
    for (m in 1:2) {
      err <- c()
      for (i in 1:K[[d]])
        err <- c(err, mean(abs(morpheus:::.extractParam(mr, 1, i)[[1]] - p[[d]][i])))
      print(paste0("p", m, " ", mean(err)))
    }

    # b
    for (m in 1:2) {
      err <- c()
      for (i in 1:K[[d]])
        err <- c(err, mean(abs(morpheus:::.extractParam(mr, 2+d, i)[[1]] - b[[d]][i])))
      print(paste0("b", m, " ", mean(err)))
    }

    # beta
    for (m in 1:2) {
      for (i in 1:K[[d]]) {
        err <- c()
        for (j in 1:d)
          err <- c(err, mean(abs(morpheus:::.extractParam(mr, 1+j, i)[[1]] - beta[[d]][1+j,i])))
        print(paste0("beta", m, "_", i, " ", mean(err)))
      }
    }
  }
}
