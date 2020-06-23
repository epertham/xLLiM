mahalanobis_distance = function(x,mu,Sigma){
    #%Compute the mahalanobis distance
    
    n = ncol(x);
    if (ncol(mu) == n)
	{x=x-mu;}#sweep(X,2,mu,"-")}
	else {x=sweep(x,1,mu,"-")}
    # [U,num] = cholcov(sigma); 
    # if num > 0
        # fprintf(1,'toto');
    # end
    p = is.positive.definite(Sigma)
    if (p) {U = chol(Sigma);    
    			Q = solve(t(U)) %*% x; #% dxn
   			y = colSums(Q^2);} #% 1xn quadratic term (M distance)
	else { print('SNPD! ');
        y=rep(Inf,n);} # % 1xn %%% EP=-Inf a l'origine ; plutot valeur Inf??
    return(y)
}