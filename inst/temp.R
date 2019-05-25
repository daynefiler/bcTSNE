bctsOG <- function(X, Z, k = max(2, round(NCOL(X)/10)), rescale = FALSE) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!is.matrix(Z)) stop("Z must be a matrix")
  if (nrow(X) != nrow(Z)) stop("X and Z must have the same rows")
  if (rescale) X <- scale(X)
  SVD <- svd(X, nu = k, nv = k)
  LM <- lm(SVD$u %*% diag(SVD$d[1:k]) ~ -1 + Z)
  S <- SVD$u %*% diag(SVD$d[1:k]) - Z %*% LM$coef
  list(S = S, U = t(SVD$v))
}

bctsGenSim <- function(n = 200, p = 3000, k = 20, seed = NULL, 
                   scaleX = TRUE, categorical = TRUE) {
  set.seed(seed)
  W <- matrix(rnorm(p*k), nrow = k)
  S <- matrix(rnorm(n*k), nrow = n)
  if (categorical) {
    z <- sample(1:3, replace = TRUE, size = n)
    Z <- model.matrix( ~ as.factor(z))[ , -1]
    b <- matrix(rnorm(k*NCOL(Z)), k)
    X <- (S - Z %*% t(b)) %*% W
  } else {
    Z <- matrix(rnorm(n))
    X <- sweep(S, 1, exp(Z)*1.5) %*% W
  }
  if (scaleX) X <- scale(X)
  ogX <- OG(X, Z, k = k)$S
  list(X = X, 
       z = z, 
       unadj = Rtsne(X, initial_dims = k)$Y, 
       og = Rtsne(ogX, pca = FALSE)$Y)
}

# s1 <- genSim(seed = 2731)
# s2 <- genSim(seed = 2731, categorical = FALSE)


bctsProcSimEmbed <- function(y, z, n = 10, categorical = TRUE) {
  getNMin <- function(x, n) Rfast::Order(x)[2:(n + 1)]
  yd <- as.matrix(dist(y))
  minInd <- apply(yd, 2, getNMin, n = n)
  minGrp <- matrix(z[minInd], ncol = ncol(minInd))
  d <- sweep(minGrp, 2, z)
  if (categorical) return(sum(d != 0))
  sum(abs(d))
}

# procEmbed(s$unadj, s$z, n = 5)
# procEmbed(s$og, s$z, n = 5)
