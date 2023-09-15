# Overview

Provides a tool for non linear mapping (non linear regression) using a mixture of regression model and an inverse regression strategy. The methods include the GLLiM model (see Deleforge et al (2015) \doi{10.1007/s11222-014-9461-5}) based on Gaussian mixtures and a robust version of GLLiM, named SLLiM (see Perthame et al (2016) \doi{10.1016/j.jmva.2017.09.009}) based on a mixture of Generalized Student distributions. The methods also include BLLiM (see Devijver et al (2017) <arXiv:1701.07899>) which is an extension of GLLiM with a sparse block diagonal structure for large covariance matrices (particularly interesting for transcriptomic data).

# Installation

```
# To get xLLiM from CRAN
install.packages("xLLiM")
library(xLLiM)
```

# Or the development version from GitHub

```
# install.packages("devtools")
devtools::install_github("epertham/xLLiM", ref = "master")
library(xLLiM)
```
