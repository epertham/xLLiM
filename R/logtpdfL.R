#% Arellano-Valle and Bolfarine generalized t distribution for dimension L
logtpdfL = function(logDetGamma, mahalx, alpha, L){
    y = lgamma(alpha+L/2) - lgamma(alpha) - (L/2)*log(2*pi) - logDetGamma - (alpha+L/2)*(log(mahalx/2 + 1)); #% 1x1 normalization constant
return(y)}