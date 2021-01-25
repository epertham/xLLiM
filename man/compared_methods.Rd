\name{bllim_cv}
\alias{bllim_cv}
\alias{mean_cv}
\alias{randomForest_cv}
\alias{lasso_cv}
\alias{mars_cv}
\alias{svm_cv}
\title{
Functions to learn a model on a training dataset and compute prediction on a testing dataset for multiple regression methods.  
}
\description{
These functions are wrappers for basic regression methods used for comparison in <arXiv:1701.07899>. 
}
\usage{
bllim_cv(trainx, trainy, testx, testy, K, verb = 0, alpha, nfolds, ...)
mean_cv(trainx, trainy, testx, testy)
randomForest_cv(trainx, trainy, testx, testy)
lasso_cv(trainx,trainy, testx, testy)
mars_cv(trainx,trainy, testx, testy)
svm_cv(trainx,trainy, testx, testy, kernel = "linear", type = "eps-regression")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainx}{An L x N matrix of training responses with variables in rows and subjects in columns}
  \item{trainy}{An D x N matrix of training covariates with variables in rows and subjects in columns}
  \item{testx}{An L x n matrix of testing responses with variables in rows and subjects in columns}
  \item{testy}{An D x n matrix of training covariates with variables in rows and subjects in columns}
  \item{K}{Initial number of components or number of clusters}
  \item{verb}{Verbosity: print out the progression of the algorithm. If \code{verb=0}, there is no print, if \code{verb=1}, the progression is printed out. Default is 0.}
  \item{alpha}{Parameter inherited from glmnet package, the elasticnet mixing parameter, see \code{\link[glmnet]{glmnet}}}
  \item{nfolds}{Parameter inherited from glmnet package, see \code{\link[glmnet]{cv.glmnet}}}
  \item{kernel}{Parameter inherited from e1071, see \code{\link[e1071]{svm}}}
  \item{type}{Parameter inherited from e1071, see \code{\link[e1071]{svm}}}
  \item{...}{Other arguments to tune bllim}
}
\value{
Returns a L x n matrix containing prediction for the provided testing matrix.  
}
\references{
[1] E. Devijver, M. Gallopin, E. Perthame. Nonlinear network-based quantitative trait prediction from transcriptomic data. Submitted, 2017, available at https://arxiv.org/abs/1701.07899.
}
\author{
Emeline Perthame (emeline.perthame@pasteur.fr), Emilie Devijver (emilie.devijver@kuleuven.be), Melina Gallopin (melina.gallopin@u-psud.fr)
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{bllim}}, \code{\link[glmnet]{glmnet-package}}, \code{\link[mda]{mda}}
}
\examples{
## Load training data
data(data.xllim)
## Load testing data
data(data.xllim.test)

trainx = t(data.xllim[1:2,]) # 2 responses in rows and 100 observations in columns
trainy = t(data.xllim[3:52,]) # 50 covariates in rows and 100 observations in columns
testx = NULL # testing responses
testy = t(data.xllim.test) # 50 testing covariates

pred <- lasso_cv(trainx, trainy, testx, testy)
}
