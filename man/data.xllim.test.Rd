\name{data.xllim.test}
\alias{data.xllim.test}
\docType{data}
\title{Testing data to run examples of usage of \code{\link{gllim_inverse_map}} and \code{\link{sllim_inverse_map}} functions
}
\description{
\code{data.xllim.test} is a matrix of simulated testing data, generated under the same GLLiM model as \code{\link{data.xllim}}, from the true parameters available in object \code{\link{data.xllim.trueparameters}}. The goal is to train a GLLiM (resp. SLLiM and BLLiM) model on training data (see \code{\link{data.xllim}}) and to retrieve the unknown responses from data.xllim.test using  \code{\link{gllim_inverse_map}} (resp. \code{\link{sllim_inverse_map}}). 
}
\usage{data(data.xllim.test)}
\format{A matrix of simulated testing data with 50 rows (covariates) and 20 columns (observations). 
}
\examples{
  data(data.xllim.test)
  dim(data.xllim.test) # 50 20
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{data.xllim}}, \code{\link{gllim_inverse_map}}, \code{\link{sllim_inverse_map}}
}
