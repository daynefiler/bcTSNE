##----------------------------------------------------------------------------##
## Create some data for testing against python script
##----------------------------------------------------------------------------##

require(dlfUtils)
require(bcTSNE)

set.seed(1234)
X <- matrix(sample(30), 10, 3)
Y <- matrix(rnorm(20), ncol = 2)
Z <- model.matrix( ~ factor(sample(0:1, 10, replace = TRUE)))
P <- bcTSNE:::calcPvals(bcTSNE:::sqdist(X))*4
P[P < 1e-12] <- 1e-12

matToPyArray(X)
matToPyArray(Y)
matToPyArray(Z)




