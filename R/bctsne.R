#' @name bctsne
#' @title Calculate BC t-SNE by orthogonal gradient descent
#' @param X numeric matrix, input matrix
#' @param Z numeric matrix, covariate matrix
#' @param outDim integer of length 1, the output dimension
#' @param perplexity numeric of length 1, the t-SNE perplexity
#' @details 
#' \code{X} should be preprocessed (e.g. PCA, centered and scaled). \code{Z}
#' is the full model matrix, including the intercept. 
#' 
#' TODO:
#' \itemize{
#'   \item Add cost function calculation to tsne C code
#'   \item Add parameters for max iterations, momentum values, etc. 
#'   \item Once cost function is calculated, return training data 
#'   \item Add vector to zero out betas
#'   \item Add verbose parameter and print statments in the C code
#' }
#' 
#' @return numeric matrix, t-SNE gradient
#' @useDynLib bcTSNE
#' @export

bctsne <- function(X, Z, outDim, perplexity) {

  ## Data checks
  stopifnot(is.numeric(X) & is.numeric(Z))
  stopifnot(is.numeric(outDim) & is.numeric(perplexity))
  if (length(outDim) > 1) {
    warning("length(outDim) > 1; only using the first element.")
    outDim <- outDim[1]
  }
  if (length(perplexity) > 1) {
    warning("length(perplexity) > 1; only using the first element.")
    perplexity <- perplexity[1]
  }
  
  ## Get dimensions
  N <- NROW(X)
  D <- NCOL(X)
  K <- NCOL(Z)
  
  ## More data checks  
  if (N != NROW(Z)) stop("'X' and 'Z' do not have the same number of rows.")
  if (K > D) stop("ncol(Z) > ncol(X) -- problem is underdetermined.")
  if (D > 100) warning("ncol(X) > 100; consider reducing X using PCA.")
  if (outDim == 0) stop("'outDim' must be >= 1.")
  if (outDim > 3) warning("'outDim' is greater than 3; consider values 1:3.")
  
  J <- as.integer(outDim)
  X <- as.vector(X)
  Z <- as.vector(Z)
  Y <- rep(0.0, N*J)
  res <- .C("tsne", X, N, D, Z, K, J, log(perplexity), Y)
  list(X = matrix(res[[1]], N, D), 
       Z = matrix(res[[4]], N, K), 
       perplexity = perplexity,
       Y = matrix(res[[8]], N, J))
  
}