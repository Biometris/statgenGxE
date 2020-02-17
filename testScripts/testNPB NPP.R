test.npb <- function(Y, n.boot=1000, K = 2)  {
  I <- nrow(Y)
  J <- ncol(Y)
  M <- min(I-1, J-1)
  E <- sweep(Y, 1, rowMeans(Y))
  E <- sweep(E, 2, colMeans(Y))
  E <- E + mean(Y)
  E.svd <- svd(E)
  U <- E.svd$u[,1:M]
  V <- E.svd$v[,1:M]
  D <- diag(E.svd$d[1:M])
  lam <- E.svd$d[1:M]
  t.obs <- lam[K+1]^2/sum(lam[(K+1):M]^2)
  t.boot <- rep(NA, n.boot)
  if (K > 0) {
    U.K <- U[,1:K]
    D.K <- diag((lam[1:K]), nrow=K, ncol=K)
    V.K <- V[,1:K]
    theta.K <- U.K %*% D.K %*% t(V.K)
  } else {
    theta.K <- matrix(0, nrow= I, ncol=J)
  }
  U.B <- U[,(K+1):M]
  D.B <- diag((lam[(K+1):M]), nrow=(M-K), ncol=(M-K))
  V.B <- V[,(K+1):M]
  R.B <- U.B %*% D.B %*% t(V.B)
  for(bb in 1:n.boot) {
    R.b <- matrix(sample(R.B, I*J, replace=T), nrow=I, ncol=J)
    E.b <- theta.K + R.b
    E.bb <- sweep(E.b,  1, rowMeans(E.b))
    E.bb <- sweep(E.bb, 2, colMeans(E.b))
    E.bb <- E.bb + mean(E.b)
    lam.b <- svd(E.bb)$d
    t.boot[bb] <- lam.b[K+1]^2/sum(lam.b[(K+1):M]^2)
  }
  pvalue <- colMeans(t.boot > matrix(rep(t.obs, n.boot),
                                     nrow=n.boot, byrow=TRUE))
  cat("Test statistics of observed data:", t.obs, "\n")
  cat("p-value of AMMI", K , "using Non-parametric bootstrap:", pvalue)
}

test.ppb <- function(Y, n.boot=1000, K = 2)  {
  I <- nrow(Y)
  J <- ncol(Y)
  M <- min(I-1, J-1)
  E <- sweep(Y, 1, rowMeans(Y))
  E <- sweep(E, 2, colMeans(Y))
  E <- E + mean(Y)
  E.svd <- svd(E)
  U <- E.svd$u[,1:M]
  V <- E.svd$v[,1:M]
  D <- diag(E.svd$d[1:M])
  lam <- E.svd$d[1:M]
  t.obs <- lam[K+1]^2/sum(lam[(K+1):M]^2)
  t.boot <- rep(NA, n.boot)
  if (K > 0) {
    U.K <- U[,1:K]
    D.K <- diag((lam[1:K]), nrow=K, ncol=K)
    V.K <- V[,1:K]
    theta.K <- U.K %*% D.K %*% t(V.K)
  } else {
    theta.K <- matrix(0, nrow= I, ncol=J)
  }
  U.B <- U[,(K+1):M]
  D.B <- diag((lam[(K+1):M]), nrow=(M-K), ncol=(M-K))
  V.B <- V[,(K+1):M]
  R.B <- U.B %*% D.B %*% t(V.B)
  for(bb in 1:n.boot) {
    R.b <- matrix(sample(R.B, I*J, replace=F), nrow=I, ncol=J)
    E.b <- theta.K + R.b
    E.bb <- sweep(E.b,  1, rowMeans(E.b))
    E.bb <- sweep(E.bb, 2, colMeans(E.b))
    E.bb <- E.bb + mean(E.b)
    lam.b <- svd(E.bb)$d
    t.boot[bb] <- lam.b[K+1]^2/sum(lam.b[(K+1):M]^2)
  }
  pvalue <- colMeans(t.boot > matrix(rep(t.obs, n.boot),
                                     nrow=n.boot, byrow=TRUE))
  cat("Test statistics of observed data:", t.obs, "\n")
  cat("p-value of AMMI", K , "using Permutation based bootstrap:", pvalue)
}








test.npb(testDat, n.boot = 10000, K = 0)
test.npb(testDat, n.boot = 10000, K = 1)

test.ppb(testDat, n.boot = 10000, K = 0)
test.ppb(testDat, n.boot = 10000, K = 1)

maizeTot <- do.call(rbind, TDMaize)
testDat2 <- as.matrix(tapply(maizeTot$yld, list(maizeTot$genotype, maizeTot$trial), I))

test.npb(testDat2, n.boot = 10000, K = 0)
test.npb(testDat2, n.boot = 10000, K = 1)
test.npb(testDat2, n.boot = 10000, K = 2)
test.npb(testDat2, n.boot = 10000, K = 3)
test.npb(testDat2, n.boot = 10000, K = 4)
test.npb(testDat2, n.boot = 10000, K = 5)
test.npb(testDat2, n.boot = 10000, K = 6)

test.ppb(testDat2, n.boot = 10000, K = 0)
test.ppb(testDat2, n.boot = 10000, K = 1)
test.ppb(testDat2, n.boot = 10000, K = 2)
test.ppb(testDat2, n.boot = 10000, K = 3)
test.ppb(testDat2, n.boot = 10000, K = 4)
test.ppb(testDat2, n.boot = 10000, K = 5)
test.ppb(testDat2, n.boot = 10000, K = 6)

geAmMaizeNw <-  gxeAmmi(TDMaize, trait = "yld", nPC = 7)
summary(geAmMaizeNw)

testDat3[is.na(testDat3)] <- mean(testDat3, na.rm = TRUE)

test.npb(testDat3, n.boot = 10000, K = 0)
test.npb(testDat3, n.boot = 10000, K = 1)
test.npb(testDat3, n.boot = 10000, K = 2)
test.npb(testDat3, n.boot = 10000, K = 3)
test.npb(testDat3, n.boot = 10000, K = 4)
test.npb(testDat3, n.boot = 10000, K = 5)
test.npb(testDat3, n.boot = 10000, K = 6)
test.npb(testDat3, n.boot = 10000, K = 7)
test.npb(testDat3, n.boot = 10000, K = 8)
test.npb(testDat3, n.boot = 10000, K = 9)

test.ppb(testDat3, n.boot = 10000, K = 0)
test.ppb(testDat3, n.boot = 10000, K = 1)
test.ppb(testDat3, n.boot = 10000, K = 2)
test.ppb(testDat3, n.boot = 10000, K = 3)
test.ppb(testDat3, n.boot = 10000, K = 4)
test.ppb(testDat3, n.boot = 10000, K = 5)
test.ppb(testDat3, n.boot = 10000, K = 6)
test.ppb(testDat3, n.boot = 10000, K = 7)
test.ppb(testDat3, n.boot = 10000, K = 8)
test.ppb(testDat3, n.boot = 10000, K = 9)

geAmAusOrig <- gxeAmmiOrig(BLUEsAus, trait = "BLUEs_yield", nPC = 10)
geAmAusNw <- gxeAmmi(BLUEsAus, trait = "BLUEs_yield", nPC = 10)
