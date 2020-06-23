\name{emgm}
\alias{emgm}
\title{Perform EM algorithm for fitting a Gaussian mixture model (GMM)
}
\usage{
emgm(X, init, maxiter,verb)
}
\arguments{
\item{X}{An \code{(M x N)} matrix with variables in rows and observations in columns. \code{M} is \code{D+L} in the proposed approach}
\item{init}{This argument can be a number \eqn{K} of classes (integer), a matrix of posterior probabilities (\code{(N x K)} matrix) or a matrix of centers (\code{(M x K)} matrix)}
\item{maxiter}{Maximum number of iterations for estimation of the GMM}
\item{verb}{Print out the progression of the algorithm. If \code{verb=0}, there is no print, if \code{verb=1}, the progression is printed out. Default is 0.}
}
\value{
Returns a list with the following elements:
\item{label}{An \code{N} vector of class assignments provided by maximum a posteriori (MAP) on posterior probabilities to belong to each of the K components for each observation}
\item{model}{A list with the estimated parameters of the  GMM}
\item{model$mu}{An \code{(M x K)} matrix of estimations of means in each cluster of the joint GMM}
\item{model$Sigma}{An \code{(M x M x K)} array of estimations of covariance matrix in each cluster of the  GMM}
\item{model$weight}{An \code{K} vector of estimated prior probabilities of each cluster}
\item{llh}{A vector of values of the log-likelihood for each iteration of the algorithm}
\item{R}{An \code{N x K} matrix of estimations of posterior probabilities to belong to each of the K components for each observation}
}
\description{Perform EM algorithm for fitting a Gaussian mixture model (GMM). In the GLLiM context, this is done  jointly on both responses and covariates
}
\references{
[1] A. Deleforge, F. Forbes, and R. Horaud. High-dimensional regression with Gaussian mixtures and partially-latent response variables. Statistics and Computing,25(5):893--911, 2015.

[2] E. Perthame, F. Forbes, and A. Deleforge. Inverse regression approach to robust nonlinear high-to-low dimensional mapping. Journal of Multivariate Analysis, 163(C):1--14, 2018. https://doi.org/10.1016/j.jmva.2017.09.009

[3] Y. Qiao and N. Minematsu. Mixture of probabilistic linear regressions: A unified view of GMM-based mapping techiques. IEEE International Conference on Acoustics, Speech, and Signal Processing, 2009. 

Converted to R from the Matlab code of the GLLiM toolbox available on: \url{https://team.inria.fr/perception/gllim_toolbox/}
}

\examples{
# data(data.xllim)
# K=5
# r = emgm(data.xllim, init=K, verb=0);  
# r$R # estimation of posterior probabilities to belong to 
## each of the K components for each observation
}
\author{Emeline Perthame (emeline.perthame@inria.fr), Florence Forbes (florence.forbes@inria.fr), Antoine Deleforge (antoine.deleforge@inria.fr)
}
\seealso{\code{\link[xLLiM]{xLLiM-package}}, \code{\link{gllim}}, \code{\link{sllim}}
}

