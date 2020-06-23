\name{data.xllim.trueparameters}
\alias{data.xllim.trueparameters}
\docType{data}
\title{True parameters used to simulate the datasets \code{\link{data.xllim}} and \code{\link{data.xllim.test}}
}
\description{
\code{data.xllim.trueparameters} is a list containing the true parameters of the GLLiM model used to generate the datasets \code{\link{data.xllim}} and \code{\link{data.xllim.test}}. We set the number of covariates to D=50, number of responses to L=2 and we simulated a GLLiM model with K=5 components. 
}
\usage{data(data.xllim.trueparameters)}
\format{A list with the following elements
\itemize{
\item{pi}{ A vector of length \code{K} of mixture weights i.e. prior probabilities for each  component}
\item{c}{ An \code{(L x K)} matrix of means of responses (X)}
\item{Gamma}{ An \code{(L x L x K)} array of \code{K} matrices of covariances of responses (X)}
\item{A}{ An \code{(D x L x K)} array of \code{K} matrices of linear transformation matrices}
\item{b}{ An \code{(D x K)} matrix in which affine transformation vectors are in columns}
\item{Sigma}{ An \code{(D x D x K)} array of covariances of \eqn{Y}}
}
\code{data.xllim.trueparameters} has the same that the values returned by \code{\link{gllim}} function.
}
\examples{
  data(data.xllim.trueparameters)
  ## data.xllim.trueparameters$pi # A vector with K=5 elements
  ## data.xllim.trueparameters$c # A matrix with dimension L=2 x K=5
  ## data.xllim.trueparameters$Gamma # An array with dimension L=2 x L=2 x K=5
  ## data.xllim.trueparameters$A # An array with dimension D=50 x L=2 x K=5
  ## data.xllim.trueparameters$b # A matrix with dimension D=50 x K=5
  ## data.xllim.trueparameters$Sigma # An array with dimension D=50 x D=50 x K=5
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{data.xllim}}, \code{\link{gllim_inverse_map}}, \code{\link{sllim_inverse_map}}
}
