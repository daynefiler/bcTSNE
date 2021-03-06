% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bctsne.R
\name{bctsne}
\alias{bctsne}
\title{Calculate BC t-SNE by orthogonal gradient descent}
\usage{
bctsne(X, Z, k = 50, outDim = 2, perplexity = 30, maxIter = 1000)
}
\arguments{
\item{X}{numeric matrix, input matrix}

\item{Z}{numeric matrix, covariate matrix}

\item{k}{integer of length 1, reduced dimension (number of eigenvectors)}

\item{outDim}{integer of length 1, the output dimension}

\item{perplexity}{numeric of length 1, the t-SNE perplexity}

\item{maxIter}{integer of length 1, the maximum iterations for the BC t-SNE
algorithm}
}
\value{
numeric matrix, t-SNE gradient
}
\description{
Calculate BC t-SNE by orthogonal gradient descent
}
\details{
\code{X} should be preprocessed (e.g. PCA, centered and scaled). \code{Z}
is the full model matrix, including the intercept. 

TODO:
\itemize{
  \item Add cost function calculation to tsne C code
  \item Add parameters for max iterations, momentum values, etc. 
  \item Once cost function is calculated, return training data 
  \item Add vector to zero out betas
  \item Add verbose parameter and print statments in the C code
}
}
