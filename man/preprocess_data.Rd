\name{preprocess_data}
\alias{preprocess_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A proposition of function to process high dimensional data before running gllim, sllim or bllim
}
\description{
XX
}
\usage{
preprocess_data(tapp,yapp,in_K,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tapp}{An \code{L x N} matrix of training responses with variables in rows and subjects in columns}
  \item{yapp}{An \code{D x N} matrix of training covariates with variables in rows and subjects in columns}
  \item{in_K}{Initial number of components or number of clusters}
  \item{...}{Other arguments to run glmnet}
}
\value{
  \item{selected.variables}{Vector of the indexes of selected variables}
  \item{clusters}{XX}
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
