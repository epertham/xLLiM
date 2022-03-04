\name{Kfoldcv_xllim}
\alias{Kfoldcv_xllim}
\alias{createFolds_xllim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to ...
}
\description{
XX (expliquer pourquoi on a repris celle de caret pour justifier sinon ca se fait pas)
}
\usage{
Kfoldcv_xllim(yapp,tapp,func,verb=1,Kfold=10,B=10,...)
createFolds_xllim(y, k = 10, list = TRUE, returnTrain = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tapp}{An \code{L x N} matrix of training responses with variables in rows and subjects in columns}
  \item{yapp}{An \code{D x N} matrix of training covariates with variables in rows and subjects in columns}
  \item{func}{XX}
  \item{verb}{XX}
  \item{Kfold}{XX}
  \item{B}{XX}
  \item{y}{XX}
  \item{k}{XX}
  \item{list}{XX}
  \item{returnTrain}{XX}
  \item{...}{Other arguments inherited from \code{func} argument to run the cross validation}
}
\value{
Kfoldcv_xllim returns XX.
createFolds_xllim returns XX. 
}
\author{
Emeline Perthame (emeline.perthame@pasteur.fr), Emilie Devijver (emilie.devijver@kuleuven.be), Melina Gallopin (melina.gallopin@u-psud.fr)
}
\seealso{
\code{\link[xLLiM]{xLLiM-package}},\code{\link[caret]{createFolds}}
}
\examples{
x <- 1
}

