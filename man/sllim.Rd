\name{sllim}
\alias{sllim}
\title{EM Algorithm for Student Locally Linear Mapping
}
\usage{
sllim(tapp,yapp,in_K,in_r=NULL,maxiter=100,Lw=0,cstr=NULL,verb=0,in_theta=NULL,
 in_phi=NULL)
}
\arguments{
\item{tapp}{An \code{Lt x N} matrix of training responses with variables in rows and subjects in columns}
\item{yapp}{An \code{D x N} matrix of training covariates with variables in rows and subjects in columns}
\item{in_K}{Initial number of components}
\item{in_r}{Initial assignments (default NULL)}
\item{maxiter}{Maximum number of iterations (default 100). The algorithm stops if the number of iterations exceeds \code{maxiter} or if the difference of likelihood between two iterations is smaller than a threshold (fixed to \eqn{0.001(max(LL)-min(LL))} where \eqn{LL} is the vector of successive log-likelihood values at each iteration). }
\item{Lw}{Number of hidden components (default 0)}
\item{cstr}{Constraints on \eqn{X} covariance matrices. Must be a list as following \code{cstr=list(Sigma="i")} constraints \eqn{\Sigma} to be diagonal and isotropic, which is the default. See details section hereafter to see the other available options to constraint the covariance matrix.}
\item{verb}{Verbosity: print out the progression of the algorithm. If \code{verb=0}, there is no print, if \code{verb=1}, the progression is printed out. Default is 0.}
\item{in_theta}{Initial parameters (default NULL), same structure as the output of this function}
\item{in_phi}{Initial parameters (default NULL), same structure as the output of this function}
}
\value{
Returns a list with the following elements:
\item{LLf}{Final log-likelihood}
\item{LL}{Log-likelihood value at each iteration of the EM algorithm}
\item{theta}{A list containing the estimations of parameters as follows:}
\item{c}{An \code{L x K} matrix of means of responses (Y) where \code{L=Lt+Lw}}
\item{Gamma}{An \code{L x L x K} array of \code{K} matrices of covariances of responses (Y) where \code{L=Lt+Lw}}
\item{A}{An \code{D x L x K} array of \code{K} matrices of affine transformation matrices where \code{L=Lt+Lw}}
\item{b}{An \code{D x K} matrix in which affine transformation vectors are in columns}
\item{Sigma}{An \code{D x D x K} array of \eqn{X} covariances}
\item{nbpar}{The number of parameters estimated in the model}
\item{phi}{A list containing the estimations of parameters as follows:}
\item{r}{An \code{N x K} matrix of posterior probabilities}
\item{pi}{A vector of length \code{K} of mixture weights i.e. prior probabilities of all components}
\item{alpha}{A vector of length \code{K} of degree of freedom parameters (heaviness of the tail) for each Student component}
}
\description{EM Algorithm for Student Locally Linear Mapping
}
\details{
This function implements the robust counterpart of GLLiM model and should be applied when outliers are present in the data.

The SLLiM model implemented in this function addresses the following non-linear mapping issue:
\deqn{ E(Y | X=x) = g(x),}
where \eqn{Y} is a L-vector of multivariate responses and \eqn{X} is a large D-vector of covariates' profiles such that \eqn{D \gg L}. The methods implemented in this package aims at estimating the non linear regression function \eqn{g}.

First, the methods of this package are based on an inverse regression strategy. The inverse conditional relation \eqn{p(X | Y)} is specified in a way that the forward relation of interest \eqn{p(Y | X)} can be deduced in closed-from. Under some hypothesis on covariance structures, the large number \eqn{D} of covariates is handled by this inverse regression trick, which acts as a dimension reduction technique. The number of parameters to estimate is therefore drastically reduced. Second, we propose to approximate the non linear \eqn{g} regression function by a piecewise affine function. Therefore, an hidden discrete variable \eqn{Z} is introduced, in order to divide the space in \eqn{K} regions such that an affine model holds between responses Y and variables X, in each  region \eqn{k}:
\deqn{X = \sum_{k=1}^K I_{Z=k} (A_k Y + b_k + E_k)}
where \eqn{A_k} is a \eqn{D \times L} matrix of coefficients for regression \eqn{k}, \eqn{b_k} is a D-vector of intercepts and \eqn{E_k} is a noise with covariance matrix proportional to \eqn{\Sigma_k}.

SLLiM is defined as the following hierarchical generalized Student mixture model for the inverse conditional density \eqn{p(X | Y)}:
\deqn{p(X=x | Y=y,Z=k; \theta,\phi) = S(x; A_kx+b_k,\Sigma_k,\alpha_k^x,\gamma_k^x)}
\deqn{p(Y=y | Z=k; \theta,\phi) = S(y; c_k,\Gamma_k,\alpha_k,1)}
\deqn{p(Z=k | \phi)=\pi_k}
where \eqn{(\theta,\phi)} are the sets of parameters  \eqn{\theta=(c_k,\Gamma_k,A_k,b_k,\Sigma_k)_{k=1}^K} and \eqn{\phi=(\pi_k,\alpha_k)_{k=1}^K}. In the previous expression, \eqn{\alpha_k} and \eqn{(\alpha_k^x,\gamma_k^x)} determine the heaviness of the tail of the generalized Student distribution, which gives robustness to the model. Note that \eqn{\alpha_k^x=\alpha_k + L/2} and \eqn{\gamma_k^x=1 + 1/2 \delta(y,c_k,\Gamma_k)} where \eqn{\delta} is the Mahalanobis distance.
The forward conditional density of interest can be deduced from these equations and is also a Student mixture of regressions model.

Like \code{\link{gllim}}, \code{\link{sllim}} allows the addition of latent variables in order to account for correlation among covariates or if it is supposed that responses are only partially observed. Adding latent factors is known to improve prediction accuracy, if \code{Lw} is not too large with regard to the number of covariates. When latent factors are added, the dimension of the response is \code{L=Lt+Lw} and \code{L=Lt} otherwise.

For SLLiM, the number of parameters to estimate is:
\deqn{(K-1)+ K(1+DL+D+L_t+ nbpar_{\Sigma}+nbpar_{\Gamma})}
where \eqn{L=L_w+L_t} and \eqn{nbpar_{\Sigma}} (resp. \eqn{nbpar_{\Gamma}}) is the number of parameters in each of the large (resp. small) covariance matrix \eqn{\Sigma_k} (resp. \eqn{\Gamma_k}). For example,
\itemize{
\item if the constraint on \eqn{\Sigma_k} is \code{cstr$Sigma="i"}, then \eqn{nbpar_{\Sigma}=1},which is the default constraint in the \code{gllim} function
\item if the constraint on \eqn{\Sigma_k} is \code{cstr$Sigma="d"}, then \eqn{nbpar_{\Sigma}=D},
\item if the constraint on \eqn{\Sigma_k} is \code{cstr$Sigma=""}, then \eqn{nbpar_{\Sigma}=D(D+1)/2},
\item if the constraint on \eqn{\Sigma_k} is \code{cstr$Sigma="*"}, then \eqn{nbpar_{\Sigma}=D(D+1)/(2K)}.
}
The rule to compute the number of parameters of \eqn{\Gamma_k} is the same as \eqn{\Sigma_k}, replacing D by \eqn{L_t}. Currently the \eqn{\Gamma_k} matrices are not constrained and \eqn{nbpar_{\Gamma}=L_t(L_t+1)/2} because for indentifiability reasons the \eqn{L_w} part is set to the identity matrix.

The user must choose the number of mixtures components \eqn{K} and, if needed, the number of latent factors \eqn{L_w}. For small datasets (less than 100 observations), we suggest to select both \eqn{(K,L_w)} by minimizing the BIC criterion. For larger datasets, to save computation time, we suggest to set \eqn{L_w} using BIC while setting \eqn{K} to an arbitrary value large enough to catch non linear relations between responses and covariates and small enough to have several observations (at least 10) in each clusters. Indeed, for large datasets, the number of clusters should not have a strong impact on the results while it is sufficiently large.

}
\references{
[1] A. Deleforge, F. Forbes, and R. Horaud. High-dimensional regression with Gaussian mixtures and partially-latent response variables. Statistics and Computing, 25(5):893--911, 2015.

[2] E. Perthame, F. Forbes, and A. Deleforge. Inverse regression approach to robust nonlinear high-to-low dimensional mapping. Journal of Multivariate Analysis, 163(C):1--14, 2018. https://doi.org/10.1016/j.jmva.2017.09.009
}
\examples{
data(data.xllim)
responses = data.xllim[1:2,] # 2 responses in rows and 100 observations in columns
covariates = data.xllim[3:52,] # 50 covariates in rows and 100 observations in columns

## Setting 5 components in the model
K = 5

## the model can be initialized by running an EM algorithm for Gaussian Mixtures (EMGM)
r = emgm(rbind(responses, covariates), init=K); 
## and then the sllim model is estimated
mod = sllim(responses,covariates,in_K=K,in_r=r);

## if initialization is not specified, the model is automatically initialized by EMGM
## mod = sllim(responses,covariates,in_K=K)

## Adding 1 latent factor 
## mod = sllim(responses,covariates,in_K=K,in_r=r,Lw=1)

## Some constraints on the covariance structure of \eqn{X} can be added
## mod = sllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="i")) 
# Isotropic covariance matrices
# (same variance among covariates but different in each component)

## mod = sllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="d")) 
# Heteroskedastic covariance matrices
# (variances are different among covariates and in each component)

## mod = sllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="")) 
# Unconstrained full covariance matrices

## mod = sllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="*")) 
# Full covariance matrices but equal for all components
}
\author{Emeline Perthame (emeline.perthame@inria.fr), Florence Forbes (florence.forbes@inria.fr), Antoine Deleforge (antoine.deleforge@inria.fr)
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{emgm}}, \code{\link{sllim_inverse_map}}, \code{\link{gllim}}
}
