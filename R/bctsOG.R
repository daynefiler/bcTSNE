#' @name bctsOG
#' @title Perform OG
#' @param X numeric matrix
#' @param Z model matrix with covariate information
#' @param k integer of length 1
#' @details Will take out when SOG gets added to CRAN
#' @return numeric matrix (A + t(A))
#' @importFrom stats lm
#' @export

bctsOG <- function(X, Z, k = max(2, round(NCOL(X)/10)), rescale = FALSE) {
  
  ## Data checks
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.matrix(Z) || !is.numeric(Z)) stop("Z must be a numeric matrix")
  if (NROW(X) != NROW(Z)) stop("X and Z must have the same rows")
  if (k > NCOL(X)) stop("'k' must be <= ncol(X).")
  stopifnot(is.numeric(k))
  if (length(k) > 1) {
    k <- k[1]
    warning("length(k) > 1, only using first element.")
  }
  
  ## Perform algorithm
  if (rescale) X <- scale(X)
  SVD <- svd(X, nu = k, nv = k)
  LM <- lm(SVD$u %*% diag(SVD$d[1:k]) ~ -1 + Z)
  S <- SVD$u %*% diag(SVD$d[1:k]) - Z %*% LM$coef
  list(S = S, U = t(SVD$v))
  
}
