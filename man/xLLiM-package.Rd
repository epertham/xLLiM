\name{xLLiM-package}
\alias{xLLiM-package}
\docType{package}
\title{
High Dimensional Locally-Linear Mapping
}
\description{
Provides a tool for non linear mapping (non linear regression) using a mixture of regression model and an inverse regression strategy. The methods include the GLLiM model (see Deleforge et al (2015) <DOI:10.1007/s11222-014-9461-5>) based on Gaussian mixtures and a robust version of GLLiM, named SLLiM (see Perthame et al (2016) <https://hal.archives-ouvertes.fr/hal-01347455>) based on a mixture of Generalized Student distributions. The methods also include BLLiM (see Devijver et al (2017) <https://arxiv.org/abs/1701.07899>) which is an extension of GLLiM with a sparse block diagonal structure for large covariance matrices (particularly interesting for transcriptomic data).
}
\details{
\tabular{ll}{
Package: \tab xLLiM\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2017-05-23\cr
License: \tab GPL (>= 2)\cr
}
The methods implemented in this package adress the following non-linear mapping issue:
\deqn{ E(Y | X=x) = g(x),}
where \eqn{Y} is a L-vector of multivariate responses and \eqn{X} is a large D-vector of covariates' profiles such that \eqn{D \gg L}. The methods implemented in this package aims at estimating the non linear regression function \eqn{g}.

First, the methods of this package are based on an inverse regression strategy. The inverse conditional relation \eqn{p(X | Y)} is specified in a way that the forward relation of interest \eqn{p(Y | X)} can be deduced in closed-from. The large number \eqn{D} of covariates is handled by this inverse regression trick, which acts as a dimension reduction technique. The number of parameters to estimate is therefore drastically reduced.

Second, we propose to approximate the non linear \eqn{g} regression function by a piecewise affine function. Therefore, a hidden discrete variable \eqn{Z} is introduced, in order to separate the space in \eqn{K} regressions such that an affine model holds in each  region \eqn{k} between responses Y and variables X:
\deqn{X = \sum_{k=1}^K I_{Z=k} (A_k Y + b_k + E_k)}
where \eqn{A_k} is a \eqn{D \times L} matrix of coeffcients for regression \eqn{k}, \eqn{b_k} is a D-vector of intercepts and \eqn{E_k} is a random noise. 

All the models implemented in this package are based on mixture of regression models. The components of the mixture are Gaussian for GLLiM. SLLiM is a robust extension of GLLiM, based on Generalized Student mixtures. Indeed, Generalized Student distributions are heavy-tailed distributions which improve the robustness of the model compared to their Gaussian counterparts. BLLiM is an extension of GLLiM designed to provide an interpretable prediction tool for the analysis of transcriptomic data. It assumes a block diagonal dependence structure between covariates (genes) conditionally to the response. The block structure is automatically chosen among a collection of models using the slope heuristics.

For both GLLiM and SLLiM, this package provides the possibility to add \eqn{L_w} latent variables, when the responses are partially observed. In this situation, the vector \eqn{Y=(T,W)} is split into an observed \eqn{L_t}-vector \eqn{T} and an unobserved \eqn{L_w}-vector \eqn{W}. The total size of the response is therefore \eqn{L=L_t+L_w} where \eqn{L_w} is chosen by the user. See [1] for details but this amounts to consider factors and allows to add structure in the large dimensional covariance matrices. The user must choose the number of mixtures components \eqn{K} and, if needed, the number of latent factors \eqn{L_w}. For small datasets (less than 100 observations), we suggest to select both \eqn{(K,L_w)} by minimizing the BIC criterion. For larger datasets, we suggest to set \eqn{L_w} using BIC while setting \eqn{K} to an arbitrary value large enough to catch non linear relations between responses and covariates and small enough to have several observations (at least 10) in each clusters. Indeed, for large datasets, the number of clusters should not have a strong impact on the results provided it is sufficiently large.

We propose to assess the prediction accuracy of a new response \eqn{x_{test}} by computing the NRMSE (Normalized Root Mean Square Error) which is the RMSE normalized by the RMSE of prediction by the mean of training responses:
\deqn{NRMSE = \frac{|| \hat{y} - x_{test}||_2}{|| \bar{y} - x_{test} ||_2}}
where \eqn{\hat{y}} is the predicted response, \eqn{x_{test}} is the true testing response and \eqn{\bar{y}} is the mean of training responses. 

The functions available in this package are used in this order: 
\itemize{
\item Step 1 (optional): Initialization of the algorithm using  a Multivariate Gaussian mixture model and an EM algorithm implemented in the \code{\link{emgm}} function. Responses and covariates must be concatenated as described in the documentation of \code{\link{emgm}} which corresponds to a joint Gaussian Mixture Model (see Qiao et al, 2009).
\item Step 2: Estimation of a regression model using one of the available models (\code{\link{gllim}}, \code{\link{sllim}} or \code{\link{bllim}}). User must specify the following arguments
\itemize{
\item for GLLiM or SLLiM: constraint on the large covariance matrices of covariates named \eqn{\Sigma_k}. These matrices can be supposed diagonal and homoskedastic (isotropic) by setting \code{cstr=list(Sigma="i")} which is the default. Other constraints are diagonal and heteroskedastic \code{(Sigma="d")}, full matrix \code{(Sigma="")} or full but equal for each class \code{(Sigma="*")}. Except for the last constraint, in all previous constraints the matrices have their own parameterization.
\item number of components \eqn{K} in the  model. 
\item for GLLiM or SLLiM: if needed, number of latent factors \eqn{L_w}
}
\item Step 3: Prediction of responses for a testing dataset using the \code{\link{gllim_inverse_map}} or \code{\link{sllim_inverse_map}} functions. 
}
}

\author{Emeline Perthame (emeline.perthame@inria.fr), Florence Forbes (florence.forbes@inria.fr), Antoine Deleforge (antoine.deleforge@inria.fr)

}
\seealso{\code{\link[shock]{shock-package}},\code{\link[capushe]{capushe-package}}
}
\references{
[1] A. Deleforge, F. Forbes, and R. Horaud. High-dimensional regression with Gaussian mixtures and partially-latent response variables. Statistics and Computing, 25(5):893--911, 2015.

[2] E. Devijver, M. Gallopin, E. Perthame. Nonlinear network-based quantitative trait prediction from transcriptomic data. Submitted, 2017, available at \url{https://arxiv.org/abs/1701.07899}.

[3] E. Perthame, F. Forbes, and A. Deleforge. Inverse regression approach to robust nonlinear high-to-low dimensional mapping. Journal of Multivariate Analysis, 163(C):1--14, 2018. https://doi.org/10.1016/j.jmva.2017.09.009

[4] X. Qiao and N. Minematsu. Mixture of probabilistic linear regressions: A unified view of GMM-based mapping techiques. IEEE International Conference on Acoustics, Speech, and Signal Processing, 2009.

The \code{\link{gllim}} and \code{\link{gllim_inverse_map}} functions have been converted to R from the original Matlab code of the GLLiM toolbox available on: \url{https://team.inria.fr/perception/gllim_toolbox/}
}
\examples{
### Not run

## Example of inverse regression with GLLiM model
# data(data.xllim)  
# dim(data.xllim) #  size 52 y 100
# responses = data.xllim[1:2,] # 2 responses in rows and 100 observations in columns
# covariates = data.xllim[3:52,] # 50 covariates in rows and 100 observations in columns

## Set 5 components in the model
# K = 5

## Step 1: initialization of the posterior probabilities (class assignments) 
## via standard EM for a joint Gaussian Mixture Model
# r = emgm(rbind(responses, covariates), init=K); 

## Step 2: estimation of the model
## Default Lw=0 and cstr$Sigma="i"
# mod = gllim(responses,covariates,in_K=K,in_r=r)

## Skip Step 1 and go to Step 2: automatic initialization and estimation of the model
# mod = gllim(responses,covariates,in_K=K)

## Alternative: Add Lw=1 latent factor to the model
# mod = gllim(responses,covariates,in_K=K,in_r=r,Lw=1)

## Different constraints on the large covariance matrices can be added: 
## see details in the documentation of the GLLiM function
## description
# mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="i")) #default
# mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="d"))
# mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma=""))
# mod = gllim(responses,covariates,in_K=K,in_r=r,cstr=list(Sigma="*"))
## End of example of inverse regression with GLLiM model

## Step 3: Prediction on a test dataset
# data(data.xllim.test) size 50 y 20
# pred = gllim_inverse_map(data.xllim.test,mod)
## Predicted responses using the mean of \eqn{p(y | x)}.
# pred$x_exp

## Example of leave-ntest-out (1 fold cross-validation) procedure 
# n = ncol(covariates)
# ntest=10
# id.test = sample(1:n,ntest)
# train.responses = responses[,-id.test]
# train.covariates = covariates[,-id.test]
# test.responses = responses[,id.test]
# test.covariates = covariates[,id.test]

## Learn the model on training data
# mod = gllim(train.responses, train.covariates,in_K=K)

## Predict responses on testing data 
# pred = gllim_inverse_map(test.covariates,mod)$x_exp

## nrmse : normalized root mean square error to measure prediction performance 
## the normalization term is the rmse of the prediction by the mean of training responses
## an nrmse larger than 1 means that the procedure performs worse than prediction by the mean
# norm_term = sqrt(rowMeans(sweep(test.responses,1,rowMeans(train.responses),"-")^2))
## Returns 1 value for each response variable 
# nrmse = sqrt(rowMeans((test.responses-pred)^2))/norm_term 

}
