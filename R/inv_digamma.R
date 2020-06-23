inv_digamma = function(y,niter=5){
# % INV_DIGAMMA    Inverse of the digamma function.
# %
# % inv_digamma(y) returns x such that digamma(x) = y.

#% never need more than 5 iterations

#% Newton iteration to solve digamma(x)-y = 0
x = exp(y)+1/2;
i = which(y <= -2.22);
x[i] = -1/(y[i] - digamma(1));

for (iter in 1:niter){ x = x - (digamma(x)-y)/trigamma(x);}

return(x)
}