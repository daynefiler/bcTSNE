#' @name bctsGenSim
#' @title Generate simulated data for the BC t-SNE algorithm 
#' @param n integer of length 1, number of observations
#' @param p integer of length 1, number of parameters
#' @param k integer of length 1, number of components to define output, 'X'
#' @param seed integer of length 1, starting seed
#' @param categorical logical of length 1, switch simulation type, see details
#' @details 
#' To add -- need to add data checks 
#' @examples 
#' ## Categorical example
#' s1 <- bctsGenSim(seed = 7231)
#' s1 <- bctsProcSim(s1)
#' p <- list(bty = "n", axes = FALSE, xlab = "", ylab = "", pch = 16)
#' m <- c("Unadjusted", "Pre-processed", "batch-corrected t-SNE")
#' par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
#' palette(c("#0072B2", "#E69F00", "gray60"))
#' mapply(plot, x = s1$embed, main = m, MoreArgs = c(col = list(s1$z), p))
#' title("Categorical group effect", outer = TRUE)
#' palette(value = "default")
#' 
#' ## Continuous example
#' s2 <- bctsGenSim(seed = 7231, categorical = FALSE)
#' s2 <- bctsProcSim(s2)
#' s2col <- rainbow(21)[findInterval(s2$z, quantile(s2$z, seq(0, 1, 0.05)))]
#' par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
#' mapply(plot, x = s2$embed, main = m, MoreArgs = c(col = list(s2col), p))
#' title("Continuous group effect", outer = TRUE)
#' 
#' @return 
#' List with the following:
#' \itemize{
#'   \item X the (n x p) simulated matrix
#'   \item z the vector of group membership
#'   \item seed and categorical parameters stored as attributes 
#' }
#' @importFrom stats lm rnorm model.matrix
#' @export

bctsGenSim <- function(n = 200, p = 3000, k = 20, seed = NULL, 
                       categorical = TRUE) {
  set.seed(seed)
  W <- matrix(rnorm(p*k), nrow = k)
  S <- matrix(rnorm(n*k), nrow = n)
  if (categorical) {
    z <- as.factor(sample(1:3, replace = TRUE, size = n))
    Z <- model.matrix( ~ z)[ , -1]
    b <- matrix(rnorm(k*NCOL(Z)), k)
    X <- (S - Z %*% t(b)) %*% W
  } else {
    z <- rnorm(n)
    X <- sweep(S, 1, exp(z)*1.5) %*% W
  }
  out <- list(X, z)
  attributes(out) <- list(names = c("X", "z"),
                          seed = seed, 
                          categorical = categorical)
  out
}

#' @name bctsProcEmbed
#' @title Analyze the embeddings of simulated data
#' @param Y numeric matrix with embeddings
#' @param z numeric vector of group membership
#' @param n interger of length 1, the number of nearest neighbors to analyze
#' @inheritParams bctsGenSim 
#' @details 
#' To add -- need to add data checks 
#' @return metric
#' @importFrom stats dist
#' @importFrom Rfast Order
#' @export

bctsProcEmbed <- function(Y, z, n = 10, categorical = TRUE) {
  z <- as.numeric(z)
  getNMin <- function(x, n) Order(x)[2:(n + 1)]
  Ydist <- as.matrix(dist(Y))
  minInd <- apply(Ydist, 2, getNMin, n = n)
  minGrp <- matrix(z[minInd], ncol = ncol(minInd))
  d <- sweep(minGrp, 2, z)
  if (categorical) return(sum(d != 0))
  sum(abs(d))
}

#' @name bctsProcSim 
#' @title Process simulated data 
#' @description Produce embeddings for the simulated data object
#' @param s list containing simulated data, returned by \link{\code{bctsGenSim}}
#' @param k integer of length 1, the number of initial dimensions to reduce X
#' @details
#' Need to add. -- need to add data checks to the function
#' @return 
#' Modified sim object with the three embeddings stored in s$embed
#' @importFrom stats model.matrix 
#' @importFrom Rtsne Rtsne
#' @export

bctsProcSim <- function(s, k = 50) {
  
  ## Create the model matrix, Z
  Z <- model.matrix( ~ s$z)
  
  ## Store the embeddings in a list
  s$embed <- list(uats = NULL, ppts = NULL, bcts = NULL)
  
  ## Perform unadjusted t-SNE with k initial dims
  s$embed$uats <- Rtsne::Rtsne(s$X, initial_dims = k)$Y
  
  ## Perform OG-preprocessed t-SNE
  s$Xtilde <- bctsOG(X = s$X, Z = Z, k = k)
  s$embed$ppts <- Rtsne::Rtsne(s$Xtilde$S, pca = FALSE)$Y
  
  ## Perform batch-corrected t-SNE
  Xsvd <- svd(s$X, nu = k, nv = k)
  s$Xred <- scale(Xsvd$u %*% diag(Xsvd$d[1:k]))
  s$embed$bcts <- bctsne(s$Xred, Z, 2, perplexity = 30)$Y
  
  ## Get the metric for how well the algorithm removed batch
  s$embedMetric <- sapply(s$embed, 
                          bctsProcEmbed, 
                          z = s$z, 
                          categorical = attr(s, "categorical"))
  
  s
  
}

