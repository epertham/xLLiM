\name{bllim}
\alias{bllim}
\title{EM Algorithm for Block diagonal Gaussian Locally Linear Mapping
}
\usage{
bllim(tapp,yapp,in_K,in_r=NULL,ninit=20,maxiter=100,verb=0,in_theta=NULL,plot=TRUE)
}
\arguments{
\item{tapp}{An \code{L x N} matrix of training responses with variables in rows and subjects in columns}
\item{yapp}{An \code{D x N} matrix of training covariates with variables in rows and subjects in columns}
\item{in_K}{Initial number of components or number of clusters}
\item{in_r}{Initial assignments (default NULL). If NULL, the model is initialized with the best initialisation among 20, computed by a joint Gaussian mixture model on both response and covariates.}
\item{ninit}{Number of random initializations. Not used of \code{in_r} is specified. Default is 20 and the random initialization which maximizes the likelihood is retained.}
\item{maxiter}{Maximum number of iterations (default 100). The algorithm stops if the number of iterations exceeds \code{maxiter} or if the difference of likelihood between two iterations is smaller than a threshold fixed to \eqn{0.001 (max(LL)-min(LL))} where \eqn{LL} is the vector of log-likelihoods at the successive iterations.}
\item{verb}{Verbosity: print out the progression of the algorithm. If \code{verb=0}, there is no print, if \code{verb=1}, the progression is printed out. Default is 0.}
\item{in_theta}{Initial parameters (default NULL), same structure as the output of this function. The EM algorithm can be  initialized either with initial assignments or initial parameters values.}
\item{plot}{Displays plots to allow user to check that the slope heuristics can be applied confidently to select the conditional block structure of predictors, as in the \code{\link[capushe]{capushe-package}} package. Default is TRUE.}
}
\value{
Returns a list with the following elements:
\item{LLf}{Final log-likelihood}
\item{LL}{Log-likelihood value at each iteration of the EM algorithm}
\item{pi}{A vector of length \code{K} of mixture weights i.e. prior probabilities for each  component}
\item{c}{An \code{(L x K)} matrix of means of responses (Y)}
\item{Gamma}{An \code{(L x L x K)} array of \code{K} matrices of covariances of responses (Y) }
\item{A}{An \code{(D x L x K)} array of \code{K} matrices of linear transformation matrices}
\item{b}{An \code{(D x K)} matrix in which affine transformation vectors are in columns}
\item{Sigma}{An \code{(D x D x K)} array of covariances of \eqn{X}}
\item{r}{An \code{(N x K)} matrix of posterior probabilities}
\item{nbpar}{The number of parameters estimated in the model}
}
\description{EM Algorithm for Block diagonal Gaussian Locally Linear Mapping
}
\details{
The BLLiM model implemented in this function adresses the following non-linear mapping issue:
\deqn{ E(Y | X=x) = g(x),}
where \eqn{Y} is a L-vector of multivariate responses and \eqn{X} is a large D-vector of covariates' profiles such that \eqn{D \gg L}. As \code{\link{gllim}} and \code{\link{sllim}}, the \code{bllim} function aims at estimating the non linear regression function \eqn{g}.

First, the methods of this package are based on an inverse regression strategy. The inverse conditional relation \eqn{p(X | Y)} is specified in a way that the forward relation of interest \eqn{p(Y | X)} can be deduced in closed-from. Under some hypothesis on covariance structures, the large number \eqn{D} of covariates is handled by this inverse regression trick, which acts as a dimension reduction technique. The number of parameters to estimate is therefore drastically reduced. Second, we propose to approximate the non linear \eqn{g} regression function by a piecewise affine function. Therefore, a hidden discrete variable \eqn{Z} is introduced, in order to divide the space into \eqn{K} regions such that an affine model holds  between responses Y and variables X in each  region \eqn{k}:
\deqn{X = \sum_{k=1}^K I_{Z=k} (A_k Y + b_k + E_k)}
where \eqn{A_k} is a \eqn{D \times L} matrix of coeffcients for regression \eqn{k}, \eqn{b_k} is a D-vector of intercepts and \eqn{E_k} is a Gaussian noise with covariance matrix \eqn{\Sigma_k}. 

BLLiM is defined as the following hierarchical Gaussian mixture model for the inverse conditional density \eqn{(X | Y)}:
\deqn{p(X | Y=y,Z=k;\theta) = N(X; A_kx+b_k,\Sigma_k)}
\deqn{p(Y | Z=k; \theta) = N(Y; c_k,\Gamma_k)}
\deqn{p(Z=k)=\pi_k}
where \eqn{\Sigma_k} is a \eqn{D \times D} block diagonal covariance structure automatically learnt from data. \eqn{\theta} is the set of parameters \eqn{\theta=(\pi_k,c_k,\Gamma_k,A_k,b_k,\Sigma_k)_{k=1}^K}.
The forward conditional density of interest \eqn{p(Y | X)} is deduced from these equations and is also a Gaussian mixture of regression model.

For a given number of affine components (or clusters) K and a given block structure, the number of parameters to estimate is:
\deqn{(K-1)+ K(DL+D+L+ nbpar_{\Sigma}+L(L+1)/2)}
where \eqn{L} is the dimension of the response, \eqn{D} is the dimension of covariates and \eqn{nbpar_{\Sigma}} is the total number of parameters in the large covariance matrix \eqn{\Sigma_k} in each cluster. This number of parameters depends on the number and size of blocks in each matrices. 

Two hyperparameters must be estimated to run BLLiM: 
\itemize{ 
\item Number of mixtures components (or clusters) \eqn{K}: we propose to use BIC criterion or slope heuristics as implemented in \code{\link[capushe]{capushe-package}}
\item For a given number of clusters K, the block structure of large covariance matrices specific of each cluster: the size and the number of blocks of each \eqn{\Sigma_k} matrix is automatically learnt from data, using an extension of the shock procedure (see \code{\link[shock]{shock-package}}). This procedure is based on a successive thresholding of sample conditional covariance matrix within clusters, building a collection of block structure candidates. The final block structure is retained using slope heuristics. 
} 

BLLiM is not only a prediction model but also an interpretable tool. For example, it is useful for the analysis of transcriptomic data. Indeed, if covariates are genes and response is a phenotype, the model provides clusters of individuals based on the relation between gene expression data and the phenotype, and also leads to infer a gene regulatory network specific for each cluster of individuals. 
}
\references{
[1] E. Devijver, M. Gallopin, E. Perthame. Nonlinear network-based quantitative trait prediction from transcriptomic data. Submitted, 2017, available at \url{https://arxiv.org/abs/1701.07899}.
}
\examples{
data(data.xllim)

## Setting 5 components in the model
K = 5

## the model can be initialized by running an EM algorithm for Gaussian Mixtures (EMGM)
r = emgm(data.xllim, init=K); 
## and then the gllim model is estimated
responses = data.xllim[1:2,] # 2 responses in rows and 100 observations in columns
covariates = data.xllim[3:52,] # 50 covariates in rows and 100 observations in columns

## if initialization is not specified, the model is automatically initialized by EMGM
# mod = bllim(responses,covariates,in_K=K)

## Prediction can be performed using prediction function gllim_inverse_map
# pred = gllim_inverse_map(covariates,mod)$x_exp
}
\author{Emeline Perthame (emeline.perthame@pasteur.fr), Emilie Devijver (emilie.devijver@kuleuven.be), Melina Gallopin (melina.gallopin@u-psud.fr)
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{emgm}}, \code{\link{gllim_inverse_map}},\code{\link[capushe]{capushe-package}},\code{\link[shock]{shock-package}}
}

