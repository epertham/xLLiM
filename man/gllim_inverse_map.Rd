\name{gllim_inverse_map}
\alias{gllim_inverse_map}
\title{Inverse Mapping from gllim or bllim parameters
}
\usage{
gllim_inverse_map(y,theta,verb=0)
}
\arguments{
\item{y}{An \code{D x N} matrix of input observations with variables in rows and subjects on columns}
\item{theta}{An object returned by the \code{gllim} function corresponding to the learned GLLiM model}
\item{verb}{Verbosity: print out the progression of the algorithm. If \code{verb=0}, there is no print, if \code{verb=1}, the progression is printed out. Default is 0.}
}
\value{
Returns a list with the following elements:
\item{x_exp}{An \code{L x N} matrix of predicted responses by posterior mean. If \eqn{L_w} latent factors are added to the model, the first \eqn{Lt} rows (\eqn{1:Lt}) are predictions of responses and rows \eqn{(L_t+1):L} (recall that \eqn{L=L_t+L_w}) are estimations of latent factors.}
\item{alpha}{Weights of the posterior Gaussian mixture model}
}
\description{This function computes the prediction of a new response from the estimation of the GLLiM model, returned by the function \code{gllim}. Given an observed \eqn{X}, the prediction of the corresponding \eqn{Y} is obtained by setting \eqn{Y} to the mean of the distribution \eqn{p(Y | X)}.
}
\details{This function computes the prediction of a new response from the estimation of GLLiM or a BLLiM model, returned by functions \code{gllim} and \code{bllim}.
Indeed, if the inverse conditional density \eqn{p(X | Y)} and the marginal density \eqn{p(Y)} are defined according to a GLLiM model (or BLLiM) (as described on \code{\link{xLLiM-package}} and \code{\link{gllim}}), the forward conditional density \eqn{p(Y | X)} can be deduced.

Under GLLiM and BLLiM model, it is recalled that the inverse conditional \eqn{p(X | Y)} is a mixture of Gaussian regressions with parameters \eqn{(\pi_k,c_k,\Gamma_k,A_k,b_k,\Sigma_k)_{k=1}^K}. Interestingly, the forward conditional \eqn{p(Y | X)} is also a mixture of Gaussian regressions with parameters \eqn{(\pi_k,c_k^*,\Gamma_k^*,A_k^*,b_k^*,\Sigma_k^*)_{k=1}^K}. These parameters have a closed-form expression depending only on \eqn{(\pi_k,c_k,\Gamma_k,A_k,b_k,\Sigma_k)_{k=1}^K}. 

Finally, the forward density (of interest) has the following expression:
\deqn{p(Y | X=x) = \sum_{k=1}^K \frac{\pi_k N(x; c_k^*,\Gamma_k^*)}{\sum_j \pi_j N(x; c_j^*,\Gamma_j^*)} N(y; A_k^*x + b_k^*,\Sigma_k^*)}
and a prediction of a new vector of responses is computed as:
\deqn{E (Y | X=x) = \sum_{k=1}^K \frac{\pi_k N(x; c_k^*,\Gamma_k^*)}{\sum_j \pi_j N(x; c_j^*,\Gamma_j^*)} (A_k^*x + b_k^*)}
where \eqn{x} is a  new vector of observed covariates.

When applied on a BLLiM model (returned by function \code{bllim}), the prediction function \code{gllim_inverse_map} accounts for the block structure of covariance matrices of the model.  
}
\references{
[1] A. Deleforge, F. Forbes, and R. Horaud. High-dimensional regression with Gaussian mixtures and partially-latent response variables. Statistics and Computing, 25(5):893--911, 2015.

[2] E. Devijver, M. Gallopin, E. Perthame. Nonlinear network-based quantitative trait prediction from transcriptomic data. Submitted, 2017, available at \url{https://arxiv.org/abs/1701.07899}.

[3] E. Perthame, F. Forbes, and A. Deleforge. Inverse regression approach to robust nonlinear high-to-low dimensional mapping. Journal of Multivariate Analysis, 163(C):1--14, 2018. https://doi.org/10.1016/j.jmva.2017.09.009

Converted to R from the Matlab code of the GLLiM toolbox available on: \url{https://team.inria.fr/perception/gllim_toolbox/}
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
mod = gllim(responses,covariates,in_K=K,in_r=r);

## Charge testing data
data(data.xllim.test)
## Prediction on a test dataset
pred = gllim_inverse_map(data.xllim.test,mod)
## Predicted responses
print(pred$x_exp)

## Can also be applied on an object returned by bllim function 
## Learn the BLLiM model
# mod = bllim(responses,covariates,in_K=K,in_r=r);
## Prediction on a test dataset
# pred = gllim_inverse_map(data.xllim.test,mod)
## Predicted responses
# print(pred$x_exp)
}
\author{Emeline Perthame (emeline.perthame@inria.fr), Florence Forbes (florence.forbes@inria.fr), Antoine Deleforge (antoine.deleforge@inria.fr)
}
\seealso{\code{\link{xLLiM-package}},\code{\link{gllim}} 
}

