emgm = function(X, init, maxiter=100,verb=0){
# % Perform EM algorithm for fitting the Gaussian mixture model.
# %   X: d x n data matrix
# %   init: k (1 x 1) or posteriors (n x k) or center (d x k)
# % Written in Matlab by Michael Chen (sth4nth@gmail.com)
# % Converted to R by Emeline Perthame (perthame.emeline@gmail.com)
# %% initialization
 
initialization = function(X, init){
  d = nrow(X) ; n = ncol(X);
  if (is.list(init))   #% initialize with a model
  {R  = expectation(X,init);}
  else {
    if (length(init) == 1) #% random initialization
    {k = init;
    idx = sample(1:n,k,replace=FALSE);
    m = X[,idx,drop=FALSE];
    label = max.col(t(sweep(t(m)%*%X,1,colSums(m^2)/2,"-")));
    u = sort(unique(label));
    count=0;
    while (k != length(u) && count<20){
    count=count+1;
    k=length(u);
    idx = sample(1:n,k,replace=FALSE);
    ###
    m = X[,idx];
    m = as.matrix(m);
    label = max.col(t(sweep(t(m)%*%X,1,colSums(m^2)/2,"-")));
    u = sort(unique(label));
    }
    k=length(u);
    R = as.matrix(sparseMatrix(i=1:n,j=label,x=rep(1,n),dims=c(n,k)))         
    }
    else {  
     if (nrow(init) == n) 
      {R = init;}
     else {
      if (nrow(init) == d) 
      {k = ncol(init);
       m = init;
       m = as.matrix(m);
       label = max.col(t(sweep(t(m)%*%X,2,colSums(m^2)/2,"-")));
       R = as.matrix(sparseMatrix(i=1:n,j=label,x=rep(1,n),dims=c(n,k)))  ;
       }
	else {stop('ERROR: init is not valid.');}
	}
	}
}
return(R)
}

expectation = function(X, model){
  mu = model$mu;
  Sigma = model$Sigma;
  w = model$weight;

  n = ncol(X);
  k = ncol(mu);
  logRho = matrix(0,n,k);

  for (i in 1:k){
    logRho[,i] = loggausspdf(X,mu[,i,drop=FALSE],Sigma[,,i]);
  }
  logRho = sweep(logRho,2,log(w),"+") 
  TT = logsumexp(logRho,2);
  llh = sum(TT)/n;
  logR= sweep(logRho,1,TT,"-")
  R = exp(logR);
  return(list(R=R, llh=llh))
}

maximization = function(X, R){
  d = nrow(X) ; n = ncol(X)
  k = ncol(R) ;

  nk = colSums(R);
  w = nk/n;
  mu = sweep(X%*%R,2,1/nk,"*") ### attention risque d'erreur ici

  Sigma = array(0,dim=c(d,d,k));
  sqrtR = sqrt(R);
  for (i in 1:k){
    Xo = sweep(X,1,mu[,i],"-") 
    Xo = sweep(Xo,2,sqrtR[,i],"*")
    Sigma[,,i] = tcrossprod(Xo)/nk[i];
    #% add a prior for numerical stability
    Sigma[,,i] = Sigma[,,i]+diag(d)*(1e-08);
  }
return(list(mu=mu,Sigma=Sigma,weight=w))
}


  if(verb>=1) print('     EM for Gaussian mixture: running ... ');

  R = initialization(X,init);
  label = max.col(R) 
  R = R[,sort(unique(label))];

  tol = 1e-14;
  llh = rep(-Inf, maxiter)
  converged = FALSE;
  t = 0;
while (!converged & t < maxiter){
    t = t+1;
    if(verb>=1) print(paste('     Step ',t,sep=""));
    model = maximization(X,as.matrix(R));
    tmp = expectation(X,model);
    R=tmp$R; llh[t]=tmp$llh
   
    label = max.col(R)
    u = unique(label);  
    if (ncol(R) != length(u))
    { R = R[,u]; } else {converged = ((llh[t+1]-llh[t]) < tol*abs(llh[t+1]));} #% remove empty components    
}

if(verb>=1) {
    if (converged)
        {print(paste('Converged in ',t,' steps.',sep=""));}
    else
        {print(paste('Did not converge in ',maxiter,' steps.',sep=""));}
			}

return(list(label= label, model= model, llh= llh[1:t], R=R))
}