\name{gllim}
\alias{gllim}
\title{EM Algorithm for Gaussian Locally Linear Mapping
}
\usage{
gllim(tapp,yapp,in_K,in_r=NULL,maxiter=100,Lw=0,cstr=NULL,verb=0,in_theta=NULL,...)
}
\arguments{
\item{tapp}{An \code{Lt x N} matrix of training responses with variables in rows and subjects in columns}
\item{yapp}{An \code{D x N} matrix of training covariates with variables in rows and subjects in columns}
\item{in_K}{Initial number of components}
\item{in_r}{Initial assignments (default NULL)}
\item{maxiter}{Maximum number of iterations (default 100). The algorithm stops if the number of iterations exceeds \code{maxiter} or if the difference of likelihood between two iterations is smaller than a threshold fixed to \eqn{0.001 (max(LL)-min(LL))} where \eqn{LL} is the vector of log-likelihoods at the successive iterations.}
\item{Lw}{Number of hidden components (default 0)}
\item{cstr}{Constraints on error covariance matrices. Must be a list as following \code{cstr=list(Sigma="i")} constraints \eqn{\Sigma_k} to be diagonal and isotropic, which is the default. See details section hereafter to see the other available options to constraint the covariance matrices. }
\item{verb}{Verbosity: print out the progression of the algorithm. If \code{verb=0}, there is no print, if \code{verb=1}, the progression is printed out. Default is 0.}
\item{in_theta}{The EM algorithm can be  initialized either with initial assignments or initial parameters values. In that case, the initial parameters (default NULL) must have the same structure as the output \code{theta} of this function.}
\item{...}{other arguments to be passed for internal use only}
}
\value{
Returns a list with the following elements:
\item{LLf}{Final log-likelihood}
\item{LL}{Log-likelihood value at each iteration of the EM algorithm}
\item{pi}{A vector of length \code{K} of mixture weights i.e. prior probabilities for each  component}
\item{c}{An \code{(L x K)} matrix of means of responses (Y) where \code{L=Lt+Lw}}
\item{Gamma}{An \code{(L x L x K)} array of \code{K} matrices of covariances of responses (Y) where \code{L=Lt+Lw}}
\item{A}{An \code{(D x L x K)} array of \code{K} matrices of linear transformation matrices where \code{L=Lt+Lw}}
\item{b}{An \code{(D x K)} matrix in which affine transformation vectors are in columns}
\item{Sigma}{An \code{(D x D x K)} array of covariances of \eqn{X}}
\item{r}{An \code{(N x K)} matrix of posterior probabilities}
\item{nbpar}{The number of parameters estimated in the model}
}
\description{EM Algorithm for Gaussian Locally Linear Mapping
}
\details{
The GLLiM model implemented in this function adresses the following non-linear mapping issue:
\deqn{ E(Y | X=x) = g(x),}
where \eqn{Y} is a L-vector of multivariate responses and \eqn{X} is a large D-vector of covariates' profiles such that \eqn{D \gg L}. The methods implemented in this package aims at estimating the non linear regression function \eqn{g}.

First, the methods of this package are based on an inverse regression strategy. The inverse conditional relation \eqn{p(X | Y)} is specified in a way that the forward relation of interest \eqn{p(Y | X)} can be deduced in closed-from. Under some hypothesis on covariance structures, the large number \eqn{D} of covariates is handled by this inverse regression trick, which acts as a dimension reduction technique. The number of parameters to estimate is therefore drastically reduced. Second, we propose to approximate the non linear \eqn{g} regression function by a piecewise affine function. Therefore, a hidden discrete variable \eqn{Z} is introduced, in order to divide the space into \eqn{K} regions such that an affine model holds  between responses Y and variables X in each  region \eqn{k}:
\deqn{X = \sum_{k=1}^K I_{Z=k} (A_k Y + b_k + E_k)}
where \eqn{A_k} is a \eqn{D \times L} matrix of coeffcients for regression \eqn{k}, \eqn{b_k} is a D-vector of intercepts and \eqn{E_k} is a Gaussian noise with covariance matrix \eqn{\Sigma_k}. 

GLLiM is defined as the following hierarchical Gaussian mixture model for the inverse conditional density \eqn{(X | Y)}:
\deqn{p(X | Y=y,Z=k;\theta) = N(X; A_kx+b_k,\Sigma_k)}
\deqn{p(Y | Z=k; \theta) = N(Y; c_k,\Gamma_k)}
\deqn{p(Z=k)=\pi_k}
where \eqn{\theta} is the set of parameters \eqn{\theta=(\pi_k,c_k,\Gamma_k,A_k,b_k,\Sigma_k)_{k=1}^K}.
The forward conditional density of interest \eqn{p(Y | X)} is deduced from these equations and is also a Gaussian mixture of regression model.

\code{\link{gllim}} allows the addition of  \eqn{L_w} latent variables in order to account for correlation among covariates or if it is supposed that responses are only partially observed. Adding latent factors is known to improve prediction accuracy, if \eqn{L_w} is not too large with regard to the number of covariates. When latent factors are added, the dimension of the response is \eqn{L=L_t+L_w} and \eqn{L=L_t} otherwise.

For GLLiM, the number of parameters to estimate is:
\deqn{(K-1)+ K(DL+D+L_t+ nbpar_{\Sigma}+nbpar_{\Gamma})}
where \eqn{L=L_w+L_t} and \eqn{nbpar_{\Sigma}} (resp. \eqn{nbpar_{\Gamma}}) is the number of parameters in each of the large (resp. small) covariance matrix \eqn{\Sigma_k} (resp. \eqn{\Gamma_k}). For example,
\itemize{
\item if the constraint on \eqn{\Sigma} is \code{cstr$Sigma="i"}, then \eqn{nbpar_{\Sigma}=1},which is the default constraint in the \code{gllim} function
\item if the constraint on \eqn{\Sigma} is \code{cstr$Sigma="d"}, then \eqn{nbpar_{\Sigma}=D}, 
\item if the constraint on \eqn{\Sigma} is \code{cstr$Sigma=""}, then \eqn{nbpar_{\Sigma}=D(D+1)/2},
\item if the constraint on \eqn{\Sigma} is \code{cstr$Sigma="*"}, then \eqn{nbpar_{\Sigma}=D(D+1)/(2K)}. 
}
The rule to compute the number of parameters of \eqn{\Gamma} is the same as \eqn{\Sigma}, replacing D by \eqn{L_t}. Currently the \eqn{\Gamma_k} matrices are not constrained and \eqn{nbpar_{\Gamma}=L_t(L_t+1)/2} because for indentifiability reasons the \eqn{L_w} part is set to the identity matrix.

The user must choose the number of mixtures components \eqn{K} and, if needed, the number of latent factors \eqn{L_w}. For small datasets (less than 100 observations), it is suggested to select both \eqn{(K,L_w)} by minimizing the BIC criterion. For larger datasets, it is suggested to save computational time, to set \eqn{L_w} using BIC while setting \eqn{K} to an arbitrary value large enough to catch non linear relations between responses and covariates and small enough to have several observations (at least 10) in each clusters. Indeed, for large datasets, the number of clusters should not have a strong impact on the results while it is sufficiently large.
}
\references{
[1] A. Deleforge, F. Forbes, and R. Horaud. High-dimensional regression with Gaussian mixtures and partially-latent response variables. Statistics and Computing, 25(5):893--911, 2015.

[2] E. Perthame, F. Forbes, and A. Deleforge. Inverse regression approach to robust nonlinear high-to-low dimensional mapping. Journal of Multivariate Analysis, 163(C):1--14, 2018. https://doi.org/10.1016/j.jmva.2017.09.009

Converted to R from the Matlab code of the GLLiM toolbox available on: \url{https://team.inria.fr/perception/gllim_toolbox/}
}
\examples{
data(data.xllim)

## Setting 5 components in the model
K =5

## the model can be initialized by running an EM algorithm for Gaussian Mixtures (EMGM)
r = emgm(data.xllim, init=K); 
## and then the gllim model is estimated
responses = data.xllim[1:2,] # 2 responses in rows and 100 observations in columns
covariates = data.xllim[3:52,] # 50 covariates in rows and 100 observations in columns
mod = gllim(responses,covariates,in_K=K,in_r=r);

## if initialization is not specified, the model is automatically initialized by EMGM
## mod = gllim(responses,covariates,in_K=K)

## Adding 1 latent factor 
## mod = gllim(responses,covariates,in_K=K,in_r=r,Lw=1)

## Some constraints on the covariance structure of \eqn{X} can be added
## mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="i")) 
# Isotropic covariances
# (same variance among covariates but different in each component)

## mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="d")) 
# Heteroskedastic covariances
# (variances are different among covariates and in each component)

## mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="")) 
# Unconstrained full matrix

## mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="*")) 
# Full matrix but equal between components
}
\author{Emeline Perthame (emeline.perthame@inria.fr), Florence Forbes (florence.forbes@inria.fr), Antoine Deleforge (antoine.deleforge@inria.fr)
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{emgm}}, \code{\link{gllim_inverse_map}}, \code{\link{sllim}}
}

