#auxiliary to test diagonality
.computeMuCheckDiag = function(X, Y, K, jd_method, β_ref)
{
  d <- ncol(X)
  #TODO: redundant code, same as computeMu() main method. Comments are stripped
  M3 <- .Moments_M3(X,Y)
  M2_t <- array(dim=c(d,d,d))
  for (i in 1:d)
  {
    ρ <- c(rep(0,i-1),1,rep(0,d-i))
    M2_t[,,i] <- .T_I_I_w(M3,ρ)
  }
  jd <- jointDiag::ajd(M2_t, method=jd_method)
  V <- if (jd_method=="uwedge") jd$B else solve(jd$A)
  M2_t <- array(dim=c(d,d,K))
  for (i in 1:K)
    M2_t[,,i] <- .T_I_I_w(M3,V[,i])
  #END of computeMu() code

  max_error <- 1.0 #TODO: tune ?
  invβ <- MASS::ginv(β_ref)
  for (i in 1:K)
  {
    shouldBeDiag <- invβ %*% M2_t[,,i] %*% t(invβ)
    expect_lt(
      mean( abs(shouldBeDiag[upper.tri(shouldBeDiag) | lower.tri(shouldBeDiag)]) ),
      max_error)
  }
}

test_that("'jedi' and 'uwedge' joint-diagonalization methods return a correct matrix",
{
  n <- 10000
  d_K <- list( c(2,2), c(5,3), c(20,13) )

  for (dk_index in 1:length(d_K))
  {
    d <- d_K[[dk_index]][1]
    K <- d_K[[dk_index]][2]
    #NOTE: sometimes large errors if pr is not balanced enough (e.g. random);
    #      same note for β. However we could be more random than that...
    β_ref <- rbind(diag(K),matrix(0,nrow=d-K,ncol=K))
    io <- generateSampleIO(n, p=rep(1/K,K-1), β=β_ref, rep(0,K), link="logit")
#    .computeMuCheckDiag(io$X, io$Y, K, jd_method="uwedge", β_ref) #TODO: sometimes failing test
    #TODO: some issues with jedi method (singular system)
    #.computeMuCheckDiag(io$X, io$Y, K, jd_method="jedi", β_ref)
  }
})
