loggausspdf = function(X, mu, Sigma){
	d = nrow(X) ; n = ncol(X) ;
	if (ncol(mu) == n)
	{X=X-mu;}#sweep(X,2,mu,"-")}
	else {X=sweep(X,1,mu,"-")}
	#X = #scale(X,center=mu,scale=FALSE) # d x n  ### X ou t(X)
	p = is.positive.definite(Sigma)
	if (! p) {
		print("SNPD !");
		y= rep(-Inf,n)
	} else {
	U = chol(Sigma) # d x d
	Q = solve(t(U),X) # d x n
	q = colSums(Q^2) # 1xn quadratic term (M distance)
	c = d*log(2*pi)+2*sum(log(diag(U))) #1x1 normalization constant
	y = -(c+q)/2;} # 1xn
	return(y)
}