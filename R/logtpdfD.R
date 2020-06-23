#% Arellano-Valle and Bolfarine generalized t distribution for dimension D
logtpdfD = function(logDetSigma, mahaly, alpha, gammakn, D){
    y = lgamma(alpha+D/2) - lgamma(alpha) - (D/2)*log(2*pi) - (D/2)*log(gammakn) - logDetSigma - (alpha+D/2)*(log(mahaly/(2*gammakn) + 1)); # % 1x1 normalization constant
return(y)}