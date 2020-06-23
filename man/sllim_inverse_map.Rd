\name{sllim_inverse_map}
\alias{sllim_inverse_map}
\title{Inverse Mapping from sllim parameters
}
\usage{
sllim_inverse_map(y,theta,verb=0)
}
\arguments{
\item{y}{An \code{D x N} matrix of input observations  with variables in rows and subjects on columns}
\item{theta}{An object returned by the \code{sllim} function}
\item{verb}{Verbosity: print out the progression of the algorithm. If \code{verb=0}, there is no print, if \code{verb=1}, the progression is printed out. Default is 0.}
}
\value{
Returns a list with the following elements:
\item{x_exp}{An \code{L x N} matrix of predicted responses by posterior mean. If \eqn{L_w} latent factors are added to the model, the first \eqn{Lt} rows (\eqn{1:Lt}) are predictions of responses and rows \eqn{(L_t+1):L} (recall that \eqn{L=L_t+L_w}) are estimations of latent factors.}
\item{alpha}{Weights of the posterior Gaussian mixture model}
}
\description{This function computes the prediction of a new response from the estimation of the SLLiM model, returned by the function \code{sllim}.}
\details{
This function computes the prediction of a new response from the estimation of a SLLiM model, returned by the function \code{sllim}.
Indeed, if the inverse conditional density \eqn{p(X | Y)} and the marginal density \eqn{p(Y)} are defined according to a SLLiM model (as described in \code{\link{xLLiM-package}} and \code{\link{sllim}}), the forward conditional density \eqn{p(Y | X)} can be deduced.

Under SLLiM model, it is recalled that the inverse conditional \eqn{p(X | Y)} is a mixture of Student regressions with parameters \eqn{(c_k,\Gamma_k,A_k,b_k,\Sigma_k)_{k=1}^K} and \eqn{(\pi_k,\alpha_k)_{k=1}^K}. Interestingly, the forward conditional \eqn{p(Y | X)} is also a mixture of Student regressions with parameters \eqn{(c_k^*,\Gamma_k^*,A_k^*,b_k^*,\Sigma_k^*)_{k=1}^K} and \eqn{(\pi_k,\alpha_k)_{k=1}^K}. These parameters have a closed-form expression depending only on \eqn{(c_k,\Gamma_k,A_k,b_k,\Sigma_k)_{k=1}^K} and \eqn{(\pi_k,\alpha_k)_{k=1}^K}. 

Finally, the forward density (of interest) has the following expression:
\deqn{p(Y | X=x) = \sum_k \frac{\pi_k S(x; c_k^*,\Gamma_k^*,\alpha_k,1)}{\sum_j \pi_j S(x; c_j^*,\Gamma_j^*,\alpha_j,1)} S(y; A_k^*x + b_k^*,\Sigma_k^*,\alpha_k^y,\gamma_k^y)}
where \eqn{(\alpha_k^y,\gamma_k^y)} determine the heaviness of the tail of the Generalized Student distribution.
Note that \eqn{\alpha_k^y= \alpha_k + D/2}  and \eqn{\gamma_k^y= 1 + 1/2 \delta(x,c_k^*,\Gamma_k^*)} where \eqn{\delta} is the Mahalanobis distance. A prediction of a new vector of responses is computed by:
\deqn{E (Y | X=x) = \sum_k \frac{\pi_k S(x; c_k^*,\Gamma_k^*,\alpha_k,1)}{\sum_j \pi_j S(x; c_j^*,\Gamma_j^*,\alpha_j,1)} (A_k^*x + b_k^*)}
where \eqn{x} is a new  vector of observed covariates.
}
\references{
[1] A. Deleforge, F. Forbes, and R. Horaud. High-dimensional regression with Gaussian mixtures and partially-latent response variables. Statistics and Computing, 25(5):893--911, 2015.

[2] E. Perthame, F. Forbes, and A. Deleforge. Inverse regression approach to robust nonlinear high-to-low dimensional mapping. Journal of Multivariate Analysis, 163(C):1--14, 2018. https://doi.org/10.1016/j.jmva.2017.09.009
}

\examples{
data(data.xllim)

## Setting 5 components in the model
K = 5

## the model can be initialized by running an EM algorithm for Gaussian Mixtures (EMGM)
r = emgm(data.xllim, init=K); 
## and then the sllim model is estimated
responses = data.xllim[1:2,] # 2 responses in rows and 100 observations in columns
covariates = data.xllim[3:52,] # 50 covariates in rows and 100 observations in columns
mod = sllim(responses,covariates,in_K=K,in_r=r);

# Prediction on a test dataset
data(data.xllim.test)
pred = sllim_inverse_map(data.xllim.test,mod)
## Predicted responses
print(pred$x_exp)

}
\author{Emeline Perthame (emeline.perthame@inria.fr), Florence Forbes (florence.forbes@inria.fr), Antoine Deleforge (antoine.deleforge@inria.fr)
}
\seealso{\code{\link{xLLiM-package}},\code{\link{sllim}} 
}

