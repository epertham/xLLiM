sllim_inverse_map = function(y,theta,verb=0){
# %%%%%%%%%%%%%%%%% Inverse Mapping from Gllim Parameters %%%%%%%%%%%%%%%%%%%
# %%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
# %% Converted to R: Emeline Perthame (2016) - perthame.emeline@gmail.com %%%
# % Description: Map N observations y using the inverse conditional
# % expectation E[x|y;theta] of the sllim model with parameters theta.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% Input %%%%
# % - y (DxN)               % Input observations to map
# % - theta  (struct)       % Gllim model parameters
# %   - theta.c (LxK)       % Gaussian means of X's prior
# %   - theta.Gamma (LxLxK) % Gaussian covariances of X's prior
# %   - theta.pi (1xK)      % Gaussian weights of X's prior
# %   - theta.A (DxLxK)     % Affine transformation matrices
# %   - theta.b (DxK)       % Affine transformation vectors
# %   - theta.Sigma (DxDxK) % Error covariances
# % - verb {0,1,2}          % Verbosity (def 1)
# %-  phi  (struct)         % Estimated parameters
# %   - phi.pi (1xK)        % t weights of X
# %   - phi.alpha (1xK)     % Arellano-Valle and Bolfarine's Generalized t
# %   distrib - alpha parameter
# %   - phi.gamma (1xK)   % Arellano-Valle and Bolfarine's Generalized t
# %   distrib - gamma parameter
# %%%% Output %%%%
# % - x_exp (LxN)           % Posterior mean estimates E[xn|yn;theta]
# % - alpha (NxK)           % Weights of the posterior GMMs
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	phi = theta$phi
	theta = theta$theta

	D = nrow(y) ; N = ncol(y) ;
	L = nrow(theta$c) ; K = ncol(theta$c) ;
	
# % ======================Inverse density parameters=========================
	if(verb>=1) print('Compute K projections to X space and weights')
	proj=array(NaN,dim=c(L,N,K));     
	logalpha=matrix(0,N,K); 
	for (k in 1:K){
    if(verb>=2) print(paste("k=",k,sep=""));
    if(verb>=2) print('AbcG ');
    if (L==1) {Ak=theta$A[,,k,drop=FALSE];} else {Ak=theta$A[,,k];} # % DxL
    bk=theta$b[,k]; # % Dx1
    Sigmak=theta$Sigma[,,k]; # %DxD  ## OK 
    if (L==1) {ck=theta$c[,k,drop=FALSE];} else {ck=theta$c[,k];} # % Lx1
    Gammak=theta$Gamma[,,k]; # % LxL ## OK 
    
    if(verb>=2) print('cks ');
    cks=Ak%*%ck+bk; ## OK 
    #####theta_star.c(:,k) = cks;
    
    if(verb>=2) print('Gks ');
    Gammaks=Sigmak+Ak%*%tcrossprod(Gammak,Ak); ## OK 
    ######theta_star.Gamma(:,:,k) = Gammaks;
    
    if(verb>=2) print('iSks ');
    iSk = solve(Sigmak) #diag(1/diag(Sigmak))
    invSigmaks2=diag(L)+Gammak%*%crossprod(Ak, iSk)%*%Ak; 
    ######theta_star.Sigma(:,:,k) = invSigmaks2;

    if(verb>=2) print('Aks ');
    Aks=solve(invSigmaks2,Gammak)%*%crossprod(Ak, iSk)      
    ######theta_star.A(:,:,k) = Aks;

    if(verb>=2) print('bks ');
    bks= solve(invSigmaks2) %*% (ck-Gammak%*%crossprod(Ak, iSk%*%bk))
    ######theta_star.b(:,k) = bks;
    
    if(verb>=2) print('projections '); 
    proj[,,k]=sweep(Aks%*%y,1,bks,"+"); 

    if(verb>=2) print('logalpha '); 
    
    p = is.positive.definite(Gammaks)
	if (!p) {sldR = -Inf;}
	else {R = chol(Gammaks) ; sldR=sum(log(diag(R)));}
    logalpha[,k]=log(phi$pi[k])+ logtpdfL(sum(log(diag(R))),mahalanobis_distance(y,cks,Gammaks),phi$alpha[k], D);
 
    if(verb>=2) print(''); 
    }
    
	den=logsumexp(logalpha,2);
	logalpha= sweep(logalpha,1,den,"-") 
	alpha=exp(logalpha); 

	x_exp = lapply(1:L,function(l) proj[l,,]*alpha)
	x_exp = lapply(x_exp,rowSums)
	x_exp = do.call(rbind, x_exp)
	
	return(list(x_exp=x_exp,alpha=alpha))
}