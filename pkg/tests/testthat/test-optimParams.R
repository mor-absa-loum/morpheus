context("OptimParams")

naive_f = function(link, M1,M2,M3, p,β,b)
{
  d = length(M1)
  K = length(p)
  λ <- sqrt(colSums(β^2))

  # Compute β x2,3 (self) tensorial products
  β2 = array(0, dim=c(d,d,K))
  β3 = array(0, dim=c(d,d,d,K))
  for (k in 1:K)
  {
    for (i in 1:d)
    {
      for (j in 1:d)
      {
        β2[i,j,k] = β[i,k]*β[j,k]
        for (l in 1:d)
          β3[i,j,l,k] = β[i,k]*β[j,k]*β[l,k]
      }
    }
  }

  res = 0
  for (i in 1:d)
  {
    term = 0
    for (k in 1:K)
      term = term + p[k]*.G(link,1,λ[k],b[k])*β[i,k]
    res = res + (term - M1[i])^2
    for (j in 1:d)
    {
      term = 0
      for (k in 1:K)
        term = term + p[k]*.G(link,2,λ[k],b[k])*β2[i,j,k]
      res = res + (term - M2[i,j])^2
      for (l in 1:d)
      {
        term = 0
        for (k in 1:K)
          term = term + p[k]*.G(link,3,λ[k],b[k])*β3[i,j,l,k]
        res = res + (term - M3[i,j,l])^2
      }
    }
  }
  res
}

test_that("naive computation provides the same result as vectorized computations",
{
  h <- 1e-7 #for finite-difference tests
  tol <- 1e-3 #large tolerance, necessary in some cases... (generally 1e-6 is OK)
  for (dK in list( c(2,2), c(5,3)))
  {
    d = dK[1]
    K = dK[2]

    M1 = runif(d, -1, 1)
    M2 = matrix(runif(d*d,-1,1), ncol=d)
    M3 = array(runif(d*d*d,-1,1), dim=c(d,d,d))

    for (link in c("logit","probit"))
    {
      op = new("OptimParams", "li"=link, "M1"=as.double(M1),
        "M2"=as.double(M2), "M3"=as.double(M3), "K"=as.integer(K))

      for (var in seq_len((2+d)*K-1))
      {
        p = runif(K, 0, 1)
        p = p / sum(p)
        β <- matrix(runif(d*K,-5,5),ncol=K)
        b = runif(K, -5, 5)
        x <- c(p[1:(K-1)],as.double(β),b)

        # Test functions values
        expect_equal( op$f(x), naive_f(link,M1,M2,M3, p,β,b) )

        # Test finite differences ~= gradient values
        dir_h <- rep(0, (2+d)*K-1)
        dir_h[var] = h

        expect_equal( op$grad_f(x)[var], ( op$f(x+dir_h) - op$f(x) ) / h, tol )
      }
    }
  }
})

# TODO: test computeW
#    computeW = function(θ)
#    {
#      require(MASS)
#      dd <- d + d^2 + d^3
#      M <- Moments(θ)
#      Id <- as.double(diag(d))
#      E <- diag(d)
#      v1 <- Y * X
#      v2 <- Y * t( apply(X, 1, function(Xi) Xi %o% Xi - Id) )
#      v3 <- Y * t( apply(X, 1, function(Xi) { return (Xi %o% Xi %o% Xi
#        - Reduce('+', lapply(1:d, function(j) as.double(Xi %o% E[j,] %o% E[j,])), rep(0, d*d*d))
#        - Reduce('+', lapply(1:d, function(j) as.double(E[j,] %o% Xi %o% E[j,])), rep(0, d*d*d))
#        - Reduce('+', lapply(1:d, function(j) as.double(E[j,] %o% E[j,] %o% Xi)), rep(0, d*d*d))) } ) )
#      Wtmp <- matrix(0, nrow=dd, ncol=dd)
#      
#
#g <- matrix(nrow=n, ncol=dd); for (i in 1:n) g[i,] = c(v1[i,], v2[i,], v3[i,]) - M
#
#
#
#
#
#
#  p <- θ$p
#  β <- θ$β
#  b <- θ$b
#
#
#
#
##  # Random generation of the size of each population in X~Y (unordered)
##  classes <- rmultinom(1, n, p)
##
##  #d <- nrow(β)
##  zero_mean <- rep(0,d)
##  id_sigma <- diag(rep(1,d))
##  X <- matrix(nrow=0, ncol=d)
##  Y <- c()
##  for (i in 1:ncol(β)) #K = ncol(β)
##  {
##    newXblock <- MASS::mvrnorm(classes[i], zero_mean, id_sigma)
##    arg_link <- newXblock %*% β[,i] + b[i]
##    probas <-
##      if (li == "logit")
##      {
##        e_arg_link = exp(arg_link)
##        e_arg_link / (1 + e_arg_link)
##      }
##      else #"probit"
##        pnorm(arg_link)
##    probas[is.nan(probas)] <- 1 #overflow of exp(x)
##    X <- rbind(X, newXblock)
##    Y <- c( Y, vapply(probas, function(p) (rbinom(1,1,p)), 1) )
##  }
#
#
#
#
#
#
#
#
#  Mhatt <- c(
#    colMeans(Y * X),
#    colMeans(Y * t( apply(X, 1, function(Xi) Xi %o% Xi - Id) )),
#    colMeans(Y * t( apply(X, 1, function(Xi) { return (Xi %o% Xi %o% Xi
#      - Reduce('+', lapply(1:d, function(j) as.double(Xi %o% E[j,] %o% E[j,])), rep(0, d*d*d))
#      - Reduce('+', lapply(1:d, function(j) as.double(E[j,] %o% Xi %o% E[j,])), rep(0, d*d*d))
#      - Reduce('+', lapply(1:d, function(j) as.double(E[j,] %o% E[j,] %o% Xi)), rep(0, d*d*d))) } ) ) ))
#  λ <- sqrt(colSums(β^2))
#  β2 <- apply(β, 2, function(col) col %o% col)
#  β3 <- apply(β, 2, function(col) col %o% col %o% col)
#  M <- c(
#    β  %*% (p * .G(li,1,λ,b)),
#    β2 %*% (p * .G(li,2,λ,b)),
#    β3 %*% (p * .G(li,3,λ,b)) )
#  print(sum(abs(Mhatt - M)))
#
#save(list=c("X", "Y"), file="v2.RData")
#
#
#
#
#browser()
#      for (i in 1:n)
#      {
#        gi <- t(as.matrix(c(v1[i,], v2[i,], v3[i,]) - M))
#        Wtmp <- Wtmp + t(gi) %*% gi / n
#      }
#      Wtmp
#      #MASS::ginv(Wtmp)
#    },
#
#    #TODO: compare with R version?
#    computeW_orig = function(θ)
#    {
#      require(MASS)
#      dd <- d + d^2 + d^3
#      M <- Moments(θ)
#      Omega <- matrix( .C("Compute_Omega",
#        X=as.double(X), Y=as.double(Y), M=as.double(M),
#        pn=as.integer(n), pd=as.integer(d),
#        W=as.double(W), PACKAGE="morpheus")$W, nrow=dd, ncol=dd )
#      Omega
#      #MASS::ginv(Omega) #, tol=1e-4)
#    },
#
#    Moments = function(θ)
#    {
#      "Vector of moments, of size d+d^2+d^3"
#
#      p <- θ$p
#      β <- θ$β
#      λ <- sqrt(colSums(β^2))
#      b <- θ$b
#
#      # Tensorial products β^2 = β2 and β^3 = β3 must be computed from current β1
#      β2 <- apply(β, 2, function(col) col %o% col)
#      β3 <- apply(β, 2, function(col) col %o% col %o% col)
#
#      c(
#        β  %*% (p * .G(li,1,λ,b)),
#        β2 %*% (p * .G(li,2,λ,b)),
#        β3 %*% (p * .G(li,3,λ,b)))
#    },
#
