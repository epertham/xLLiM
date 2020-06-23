logsumexp = function(x,dimension=c(1,2)){
# % Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
# %   By default dim = 2 (columns).
	if (is.null(dimension)) dimension=1 ;
# % subtract the largest in each column
	y = apply(x,1,max)
	x = x-y#sweep(x,1,y,"-")
	s = y+log(rowSums(exp(x)))
	i = is.infinite(y)
	if (sum(i) > 0) s[i]=y[i]
	return(s)
}