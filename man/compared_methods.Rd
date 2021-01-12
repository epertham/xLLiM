\name{bllim_cv}
\alias{lasso_cv}
\alias{mars_cv}
\alias{bllim_cv}
\title{
Functions to learn a model on a training dataset and compute prediction on a testing dataset for multiple regression methods.  
}
\description{
These functions are wrappers for basic regression methods used for comparison in <arXiv:1701.07899>. 
}
\usage{
bllim_cv(trainx,trainy,testx,testy)
lasso_cv(trainx,trainy,testx,testy)
mars_cv(trainx,trainy,testx,testy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainx}{An L x N matrix of training responses with variables in rows and subjects in columns}
  \item{trainy}{An D x N matrix of training covariates with variables in rows and subjects in columns}
  \item{testx}{An L x n matrix of testing responses with variables in rows and subjects in columns}
  \item{testy}{An D x n matrix of training covariates with variables in rows and subjects in columns}
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
\seealso{
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{bllim}}, \code{\link[glmnet]{glmnet-package}}, \code{\link[mda]{mda-package}}
}
\examples{
## Load training data
data(data.xllim)
## Load testing data
data(data.xllim.test)

trainx = data.xllim[1:2,] # 2 responses in rows and 100 observations in columns
trainy = data.xllim[3:52,] # 50 covariates in rows and 100 observations in columns
testx = NULL # testing responses
testy = data.xllim.test # 50 testing covariates

lasso_cv(trainx, trainy, testx, testy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
