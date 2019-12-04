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
#'   x0 is a vector containing respectively the K-1 first elements of p, then β by
#'   columns, and finally b: \code{x0 = c(p[1:(K-1)],as.double(β),b)}.
#'
#' @seealso \code{multiRun} to estimate statistics based on β, and
#'   \code{generateSampleIO} for I/O random generation.
#'
#' @examples
#' # Optimize parameters from estimated μ
#' io = generateSampleIO(10000, 1/2, matrix(c(1,-2,3,1),ncol=2), c(0,0), "logit")
#' μ = computeMu(io$X, io$Y, list(K=2))
#' o <- optimParams(io$X, io$Y, 2, "logit")
#' x0 <- list(p=1/2, β=μ, b=c(0,0))
#' par0 <- o$run(x0)
#' # Compare with another starting point
#' x1 <- list(p=1/2, β=2*μ, b=c(0,0))
#' par1 <- o$run(x1)
#' o$f( o$linArgs(par0) )
#' o$f( o$linArgs(par1) )
#' @export
optimParams = function(X, Y, K, link=c("logit","probit"))
{
	# Check arguments
  if (!is.matrix(X) || any(is.na(X)))
    stop("X: numeric matrix, no NAs")
  if (!is.numeric(Y) || any(is.na(Y)) || any(Y!=0 | Y!=1))
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
		M1 = "numeric",
		M2 = "numeric", #M2 easier to process as a vector
		M3 = "numeric", #same for M3
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
      M <- computeMoments(optargs$X,optargs$Y)
      M1 <<- as.double(M[[1]])
      M2 <<- as.double(M[[2]])
      M3 <<- as.double(M[[3]])

			n <<- nrow(X)
			d <<- length(M1)
      W <<- diag(d+d^2+d^3) #initialize at W = Identity
		},

		expArgs = function(x)
		{
			"Expand individual arguments from vector x"

			list(
				# p: dimension K-1, need to be completed
				"p" = c(x[1:(K-1)], 1-sum(x[1:(K-1)])),
				"β" = matrix(x[K:(K+d*K-1)], ncol=K),
				"b" = x[(K+d*K):(K+(d+1)*K-1)])
		},

		linArgs = function(o)
		{
			" Linearize vectors+matrices into a vector x"

			c(o$p[1:(K-1)], as.double(o$β), o$b)
		},

    getOmega = function(theta)
    {
      dim <- d + d^2 + d^3
      matrix( .C("Compute_Omega",
        X=as.double(X), Y=as.double(Y), pn=as.integer(n), pd=as.integer(d),
        p=as.double(theta$p), β=as.double(theta$β), b=as.double(theta$b),
        W=as.double(W), PACKAGE="morpheus")$W, nrow=dim, ncol=dim)
    },

    f = function(theta)
    {
			"Product t(Mi - hat_Mi) W (Mi - hat_Mi) with Mi(theta)"

      p <- theta$p
			β <- theta$β
			λ <- sqrt(colSums(β^2))
			b <- theta$b

			# Tensorial products β^2 = β2 and β^3 = β3 must be computed from current β1
			β2 <- apply(β, 2, function(col) col %o% col)
			β3 <- apply(β, 2, function(col) col %o% col %o% col)

			A <- matrix(c(
				β  %*% (p * .G(li,1,λ,b)) - M1,
				β2 %*% (p * .G(li,2,λ,b)) - M2,
				β3 %*% (p * .G(li,3,λ,b)) - M3), ncol=1)
      t(A) %*% W %*% A
    },

		grad_f = function(x)
		{
			"Gradient of f, dimension (K-1) + d*K + K = (d+2)*K - 1"

      # TODO: formula -2 t(grad M(theta)) . W . (Mhat - M(theta))
    }

    grad_M = function(theta)
    {
      # TODO: adapt code below for grad of d+d^2+d^3 vector of moments,
      # instead of grad (sum(Mhat-M(theta)^2)) --> should be easier

      P <- expArgs(x)
			p <- P$p
			β <- P$β
			λ <- sqrt(colSums(β^2))
			μ <- sweep(β, 2, λ, '/')
			b <- P$b

			# Tensorial products β^2 = β2 and β^3 = β3 must be computed from current β1
			β2 <- apply(β, 2, function(col) col %o% col)
			β3 <- apply(β, 2, function(col) col %o% col %o% col)

			# Some precomputations
			G1 = .G(li,1,λ,b)
			G2 = .G(li,2,λ,b)
			G3 = .G(li,3,λ,b)
			G4 = .G(li,4,λ,b)
			G5 = .G(li,5,λ,b)

			# (Mi - hat_Mi)^2 ' == (Mi - hat_Mi)' 2(Mi - hat_Mi) = Mi' Fi
			F1 = as.double( 2 * ( β  %*% (p * G1) - M1 ) )
			F2 = as.double( 2 * ( β2 %*% (p * G2) - M2 ) )
			F3 = as.double( 2 * ( β3 %*% (p * G3) - M3 ) )

			km1 = 1:(K-1)
			grad <- #gradient on p
			  t( sweep(as.matrix(β [,km1]), 2, G1[km1], '*') - G1[K] * β [,K] ) %*% F1 +
				t( sweep(as.matrix(β2[,km1]), 2, G2[km1], '*') - G2[K] * β2[,K] ) %*% F2 +
				t( sweep(as.matrix(β3[,km1]), 2, G3[km1], '*') - G3[K] * β3[,K] ) %*% F3

			grad_β <- matrix(nrow=d, ncol=K)
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

				grad_β[i,] <- t(dβ) %*% F1 + t(dβ2) %*% F2 + t(dβ3) %*% F3
			}
			grad <- c(grad, as.double(grad_β))

			grad = c(grad, #gradient on b
				t( sweep(β,  2, p * G2, '*') ) %*% F1 +
				t( sweep(β2, 2, p * G3, '*') ) %*% F2 +
				t( sweep(β3, 2, p * G4, '*') ) %*% F3 )

			grad
		},

    # TODO: rename x(0) into theta(0) --> θ
		run = function(x0)
		{
			"Run optimization from x0 with solver..."

	    if (!is.list(x0))
		    stop("x0: list")
      if (is.null(x0$β))
        stop("At least x0$β must be provided")
			if (!is.matrix(x0$β) || any(is.na(x0$β)) || ncol(x0$β) != K)
				stop("x0$β: matrix, no NA, ncol == K")
      if (is.null(x0$p))
        x0$p = rep(1/K, K-1)
      else if (length(x0$p) != K-1 || sum(x0$p) > 1)
        stop("x0$p should contain positive integers and sum to < 1")
      # Next test = heuristic to detect missing b (when matrix is called "beta")
      if (is.null(x0$b) || all(x0$b == x0$β))
        x0$b = rep(0, K)
      else if (any(is.na(x0$b)))
        stop("x0$b cannot have missing values")

			op_res = constrOptim( linArgs(x0), .self$f, .self$grad_f,
				ui=cbind(
					rbind( rep(-1,K-1), diag(K-1) ),
					matrix(0, nrow=K, ncol=(d+1)*K) ),
				ci=c(-1,rep(0,K-1)) )

      # We get a first non-trivial estimation of W: getOmega(theta)^{-1}
      # TODO: loop, this redefine f, so that we can call constrOptim again...
      # Stopping condition? N iterations? Delta <= ε ?

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
