\name{data.xllim}
\alias{data.xllim}
\docType{data}
\title{Simulated data to run examples of usage of \code{\link{gllim}} and \code{\link{sllim}} functions
}
\description{
Matrix of simulated data, generated under a GLLiM model, with K=5 clusters from the true parameters available in object \code{\link{data.xllim.trueparameters}}. The goal is to learn the non linear relation between the responses (Y) and the covariates (X) using \code{\link{gllim}}, \code{\link{bllim}} or \code{\link{sllim}}. Details are given hereafter. 
}
\details{
This dataset is generated under a GLLiM model with L=2, D=50 and N=100. 

First, the responses \eqn{Y} are generated according to a Gaussian Mixture model with K=5 clusters: 
\deqn{p(Y=y | Z=k)= N(y; c_k,\Gamma_k)}
where each \eqn{(c_k)_{k=1}^K} is a L-vector randomly sampled from a standardized Gaussian, \eqn{(\Gamma_k)_{k=1}^K} are LxL random correlation matrix and \eqn{Z} is a multinomial hidden variable which indicates the cluster membership of each observation:
\deqn{p(Z=k) =\pi_k}
where the probabilities \eqn{(\pi_k)_{k=1}^K} are sampled from a standard uniform distribution and normalized to sum to 1. 

Then, the covariates \eqn{X} are generated according to a Gaussian Mixture of regressions. It is recalled that GLLiM models the following inverse relation, which is used to generate \eqn{X}:
\deqn{X = \sum_{k=1}^{K=5} I_{Z=k}(A_kX+b_k+E_k)}
where \eqn{Y} is the vector of L responses and \eqn{X} is the vector of D covariates and \eqn{Z} is the hidden variable of cluster membership introduced above. Regression coefficients \eqn{A_k} and intercepts \eqn{b_k} are sampled from a standard Gaussian and the covariance matrix of the noise \eqn{\Sigma_k=Var(E_k)} is the identity.  

The goal is to retrieve \eqn{Y} from \eqn{X} using \code{\link{gllim}}, \code{\link{bllim}} or \code{\link{sllim}}.
}
\usage{data(data.xllim)}
\format{A matrix of simulated data with 52 rows and 100 columns (observations). The first 2 rows are responses (Y) and the last 50 rows are covariates (X). The goal is to retrieve \eqn{Y} from \eqn{X} using \code{\link{gllim}} or \code{\link{sllim}}.
}
\examples{
data(data.xllim)
dim(data.xllim) # 52 100
Y = data.xllim[1:2,] # responses # 2 100
X = data.xllim[3:52,] # covariates # 50 100
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{gllim}}, \code{\link{sllim}}, \code{\link{data.xllim.test}}
}
