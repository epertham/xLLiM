\name{preprocess_data}
\alias{preprocess_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A proposition of function to process high dimensional data before running gllim, sllim or bllim
}
\usage{
preprocess_data(tapp,yapp,in_K,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tapp}{An \code{L x N} matrix of training responses with variables in rows and subjects in columns}
  \item{yapp}{An \code{D x N} matrix of training covariates with variables in rows and subjects in columns}
  \item{in_K}{Initial number of components or number of clusters}
  \item{...}{Other arguments of glmnet can be passed}
}
\value{
  \item{selected.variables}{Vector of the indexes of selected variables. Selection is made within clusters and merged hereafter.}
  \item{clusters}{Initialization clusters with k-means}
}
\description{
The goal of \code{preprocess_data()} is to get relevant clusters for G-, S-, or BLLiM initialization, coupled with a feature selection for high-dimensional datasets. This function is an alternative to the default initialization implemented in \code{gllim()}, \code{sllim()} and \code{bllim()}. 

In this function, clusters are initialized with K-means, and variable selection is performed with a LASSO (\code{glmnet}) within each clusters. Then selected features are merged to get a subset variables before running any prediction method of xLLiM. 
}
\references{
[1] E. Devijver, M. Gallopin, E. Perthame. Nonlinear network-based quantitative trait prediction from transcriptomic data. Submitted, 2017, available at \url{https://arxiv.org/abs/1701.07899}.
}
\author{
Emeline Perthame (emeline.perthame@pasteur.fr), Emilie Devijver (emilie.devijver@kuleuven.be), Melina Gallopin (melina.gallopin@u-psud.fr)
}
\seealso{
\code{\link[xLLiM]{xLLiM-package}}, \code{\link[glmnet]{glmnet-package}}, \code{\link{kmeans}}
}
\examples{
x <- 1
}
