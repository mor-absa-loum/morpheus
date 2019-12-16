#' Wrapper function for OptimParams class
#'
#' @param K Number of populations.
#' @param link The link type, 'logit' or 'probit'.
#' @param X Data matrix of covariables
#' @param Y Output as a binary vector
#'
#' @return An object 'op' of class OptimParams, initialized so that \code{op$run(x0)}
#'   outputs the list of optimized parameters
#'   \itemize{
#'     \item p: proportions, size K
#'     \item β: regression matrix, size dxK
#'     \item b: intercepts, size K
#'   }
#'   θ0 is a vector containing respectively the K-1 first elements of p, then β by
#'   columns, and finally b: \code{θ0 = c(p[1:(K-1)],as.double(β),b)}.
#'
#' @seealso \code{multiRun} to estimate statistics based on β, and
#'   \code{generateSampleIO} for I/O random generation.
#'
#' @examples
#' # Optimize parameters from estimated μ
#' io = generateSampleIO(10000, 1/2, matrix(c(1,-2,3,1),ncol=2), c(0,0), "logit")
#' μ = computeMu(io$X, io$Y, list(K=2))
#' o <- optimParams(io$X, io$Y, 2, "logit")
#' θ0 <- list(p=1/2, β=μ, b=c(0,0))
#' par0 <- o$run(θ0)
#' # Compare with another starting point
#' θ1 <- list(p=1/2, β=2*μ, b=c(0,0))
#' par1 <- o$run(θ1)
#' o$f( o$linArgs(par0) )
#' o$f( o$linArgs(par1) )
#' @export
optimParams <- function(X, Y, K, link=c("logit","probit"))
{
  # Check arguments
  if (!is.matrix(X) || any(is.na(X)))
    stop("X: numeric matrix, no NAs")
  if (!is.numeric(Y) || any(is.na(Y)) || any(Y!=0 & Y!=1))
    stop("Y: binary vector with 0 and 1 only")
  link <- match.arg(link)
  if (!is.numeric(K) || K!=floor(K) || K < 2)
    stop("K: integer >= 2")

  # Build and return optimization algorithm object
  methods::new("OptimParams", "li"=link, "X"=X,
    "Y"=as.integer(Y), "K"=as.integer(K))
}

#' Encapsulated optimization for p (proportions), β and b (regression parameters)
#'
#' Optimize the parameters of a mixture of logistic regressions model, possibly using
#' \code{mu <- computeMu(...)} as a partial starting point.
#'
#' @field li Link function, 'logit' or 'probit'
#' @field X Data matrix of covariables
#' @field Y Output as a binary vector
#' @field K Number of populations
#' @field d Number of dimensions
#' @field W Weights matrix (iteratively refined)
#'
setRefClass(
  Class = "OptimParams",

  fields = list(
    # Inputs
    li = "character", #link function
    X = "matrix",
    Y = "numeric",
    Mhat = "numeric", #vector of empirical moments
    # Dimensions
    K = "integer",
    n = "integer",
    d = "integer",
    # Weights matrix (generalized least square)
    W = "matrix"
  ),

  methods = list(
    initialize = function(...)
    {
      "Check args and initialize K, d, W"

      callSuper(...)
      if (!hasArg("X") || !hasArg("Y") || !hasArg("K") || !hasArg("li"))
        stop("Missing arguments")

      # Precompute empirical moments
      M <- computeMoments(X, Y)
      M1 <- as.double(M[[1]])
      M2 <- as.double(M[[2]])
      M3 <- as.double(M[[3]])
      Mhat <<- c(M1, M2, M3)

      n <<- nrow(X)
      d <<- length(M1)
      W <<- diag(d+d^2+d^3) #initialize at W = Identity
    },

    expArgs = function(v)
    {
      "Expand individual arguments from vector v into a list"

      list(
        # p: dimension K-1, need to be completed
        "p" = c(v[1:(K-1)], 1-sum(v[1:(K-1)])),
        "β" = t(matrix(v[K:(K+d*K-1)], ncol=d)),
        "b" = v[(K+d*K):(K+(d+1)*K-1)])
    },

    linArgs = function(L)
    {
      "Linearize vectors+matrices from list L into a vector"

      # β linearized row by row, to match derivatives order
      c(L$p[1:(K-1)], as.double(t(L$β)), L$b)
    },

    computeW = function(θ)
    {
      require(MASS)
      dd <- d + d^2 + d^3
      M <- Moments(θ)
      Omega <- matrix( .C("Compute_Omega",
        X=as.double(X), Y=as.integer(Y), M=as.double(M),
        pn=as.integer(n), pd=as.integer(d),
        W=as.double(W), PACKAGE="morpheus")$W, nrow=dd, ncol=dd )
      MASS::ginv(Omega)
    },

    Moments = function(θ)
    {
      "Vector of moments, of size d+d^2+d^3"

      p <- θ$p
      β <- θ$β
      λ <- sqrt(colSums(β^2))
      b <- θ$b

      # Tensorial products β^2 = β2 and β^3 = β3 must be computed from current β1
      β2 <- apply(β, 2, function(col) col %o% col)
      β3 <- apply(β, 2, function(col) col %o% col %o% col)

      c(
        β  %*% (p * .G(li,1,λ,b)),
        β2 %*% (p * .G(li,2,λ,b)),
        β3 %*% (p * .G(li,3,λ,b)))
    },

    f = function(θ)
    {
      "Product t(hat_Mi - Mi) W (hat_Mi - Mi) with Mi(theta)"

      L <- expArgs(θ)
      A <- as.matrix(Mhat - Moments(L))
      t(A) %*% W %*% A
    },

    grad_f = function(θ)
    {
      "Gradient of f, dimension (K-1) + d*K + K = (d+2)*K - 1"

      L <- expArgs(θ)
      -2 * t(grad_M(L)) %*% W %*% as.matrix(Mhat - Moments(L))
    },

    grad_M = function(θ)
    {
      "Gradient of the vector of moments, size (dim=)d+d^2+d^3 x K-1+K+d*K"

      p <- θ$p
      β <- θ$β
      λ <- sqrt(colSums(β^2))
      μ <- sweep(β, 2, λ, '/')
      b <- θ$b

      res <- matrix(nrow=nrow(W), ncol=0)

      # Tensorial products β^2 = β2 and β^3 = β3 must be computed from current β1
      β2 <- apply(β, 2, function(col) col %o% col)
      β3 <- apply(β, 2, function(col) col %o% col %o% col)

      # Some precomputations
      G1 = .G(li,1,λ,b)
      G2 = .G(li,2,λ,b)
      G3 = .G(li,3,λ,b)
      G4 = .G(li,4,λ,b)
      G5 = .G(li,5,λ,b)

      # Gradient on p: K-1 columns, dim rows
      km1 = 1:(K-1)
      res <- cbind(res, rbind(
        sweep(as.matrix(β [,km1]), 2, G1[km1], '*') - G1[K] * β [,K],
        sweep(as.matrix(β2[,km1]), 2, G2[km1], '*') - G2[K] * β2[,K],
        sweep(as.matrix(β3[,km1]), 2, G3[km1], '*') - G3[K] * β3[,K] ))

      for (i in 1:d)
      {
        # i determines the derivated matrix dβ[2,3]

        dβ_left <- sweep(β, 2, p * G3 * β[i,], '*')
        dβ_right <- matrix(0, nrow=d, ncol=K)
        block <- i
        dβ_right[block,] <- dβ_right[block,] + 1
        dβ <- dβ_left + sweep(dβ_right, 2,  p * G1, '*')

        dβ2_left <- sweep(β2, 2, p * G4 * β[i,], '*')
        dβ2_right <- do.call( rbind, lapply(1:d, function(j) {
          sweep(dβ_right, 2, β[j,], '*')
        }) )
        block <- ((i-1)*d+1):(i*d)
        dβ2_right[block,] <- dβ2_right[block,] + β
        dβ2 <- dβ2_left + sweep(dβ2_right, 2, p * G2, '*')

        dβ3_left <- sweep(β3, 2, p * G5 * β[i,], '*')
        dβ3_right <- do.call( rbind, lapply(1:d, function(j) {
          sweep(dβ2_right, 2, β[j,], '*')
        }) )
        block <- ((i-1)*d*d+1):(i*d*d)
        dβ3_right[block,] <- dβ3_right[block,] + β2
        dβ3 <- dβ3_left + sweep(dβ3_right, 2, p * G3, '*')

        res <- cbind(res, rbind(dβ, dβ2, dβ3))
      }

      # Gradient on b
      res <- cbind(res, rbind(
        sweep(β,  2, p * G2, '*'),
        sweep(β2, 2, p * G3, '*'),
        sweep(β3, 2, p * G4, '*') ))

      res
    },

    run = function(θ0)
    {
      "Run optimization from θ0 with solver..."

      if (!is.list(θ0))
        stop("θ0: list")
      if (is.null(θ0$β))
        stop("At least θ0$β must be provided")
      if (!is.matrix(θ0$β) || any(is.na(θ0$β))
        || nrow(θ0$β) != d || ncol(θ0$β) != K)
      {
        stop("θ0$β: matrix, no NA, nrow = d, ncol = K")
      }
      if (is.null(θ0$p))
        θ0$p = rep(1/K, K-1)
      else if (!is.numeric(θ0$p) || length(θ0$p) != K-1
        || any(is.na(θ0$p)) || sum(θ0$p) > 1)
      {
        stop("θ0$p: length K-1, no NA, positive integers, sum to <= 1")
      }
      if (is.null(θ0$b))
        θ0$b = rep(0, K)
      else if (!is.numeric(θ0$b) || length(θ0$b) != K || any(is.na(θ0$b)))
        stop("θ0$b: length K, no NA")
      # TODO: stopping condition? N iterations? Delta <= epsilon ?
      loopMax <- 2
      for (loop in 1:loopMax)
      {
        op_res = constrOptim( linArgs(θ0), .self$f, .self$grad_f,
          ui=cbind(
            rbind( rep(-1,K-1), diag(K-1) ),
            matrix(0, nrow=K, ncol=(d+1)*K) ),
          ci=c(-1,rep(0,K-1)) )
        if (loop < loopMax) #avoid computing an extra W
          W <<- computeW(expArgs(op_res$par))
        #print(op_res$value) #debug
        #print(expArgs(op_res$par)) #debug
      }

      expArgs(op_res$par)
    }
  )
)

# Compute vectorial E[g^{(order)}(<β,x> + b)] with x~N(0,Id) (integral in R^d)
#                 = E[g^{(order)}(z)] with z~N(b,diag(λ))
# by numerically evaluating the integral.
#
# @param link Link, 'logit' or 'probit'
# @param order Order of derivative
# @param λ Norm of columns of β
# @param b Intercept
#
.G <- function(link, order, λ, b)
{
  # NOTE: weird "integral divergent" error on inputs:
  # link="probit"; order=2; λ=c(531.8099,586.8893,523.5816); b=c(-118.512674,-3.488020,2.109969)
  # Switch to pracma package for that (but it seems slow...)
  sapply( seq_along(λ), function(k) {
    res <- NULL
    tryCatch({
      # Fast code, may fail:
      res <- stats::integrate(
        function(z) .deriv[[link]][[order]](λ[k]*z+b[k]) * exp(-z^2/2) / sqrt(2*pi),
        lower=-Inf, upper=Inf )$value
    }, error = function(e) {
      # Robust slow code, no fails observed:
      sink("/dev/null") #pracma package has some useless printed outputs...
      res <- pracma::integral(
        function(z) .deriv[[link]][[order]](λ[k]*z+b[k]) * exp(-z^2/2) / sqrt(2*pi),
        xmin=-Inf, xmax=Inf, method="Kronrod")
      sink()
    })
    res
  })
}

# Derivatives list: g^(k)(x) for links 'logit' and 'probit'
#
.deriv <- list(
  "probit"=list(
    # 'probit' derivatives list;
    # NOTE: exact values for the integral E[g^(k)(λz+b)] could be computed
    function(x) exp(-x^2/2)/(sqrt(2*pi)),                     #g'
    function(x) exp(-x^2/2)/(sqrt(2*pi)) *  -x,               #g''
    function(x) exp(-x^2/2)/(sqrt(2*pi)) * ( x^2 - 1),        #g^(3)
    function(x) exp(-x^2/2)/(sqrt(2*pi)) * (-x^3 + 3*x),      #g^(4)
    function(x) exp(-x^2/2)/(sqrt(2*pi)) * ( x^4 - 6*x^2 + 3) #g^(5)
  ),
  "logit"=list(
    # Sigmoid derivatives list, obtained with http://www.derivative-calculator.net/
    # @seealso http://www.ece.uc.edu/~aminai/papers/minai_sigmoids_NN93.pdf
    function(x) {e=exp(x); .zin(e                                    /(e+1)^2)}, #g'
    function(x) {e=exp(x); .zin(e*(-e   + 1)                         /(e+1)^3)}, #g''
    function(x) {e=exp(x); .zin(e*( e^2 - 4*e    + 1)                /(e+1)^4)}, #g^(3)
    function(x) {e=exp(x); .zin(e*(-e^3 + 11*e^2 - 11*e   + 1)       /(e+1)^5)}, #g^(4)
    function(x) {e=exp(x); .zin(e*( e^4 - 26*e^3 + 66*e^2 - 26*e + 1)/(e+1)^6)}  #g^(5)
  )
)

# Utility for integration: "[return] zero if [argument is] NaN" (Inf / Inf divs)
#
# @param x Ratio of polynoms of exponentials, as in .S[[i]]
#
.zin <- function(x)
{
  x[is.nan(x)] <- 0.
  x
}
