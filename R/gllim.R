gllim = function(tapp,yapp,in_K,in_r=NULL,maxiter=100,Lw=0,cstr=NULL,verb=0,in_theta=NULL,...){
# %%%%%%%% General EM Algorithm for Gaussian Locally Linear Mapping %%%%%%%%%
# %%% Author: Antoine Deleforge (April 2013) - antoine.deleforge@inria.fr %%%
# % Description: Compute maximum likelihood parameters theta and posterior
# % probabilities r=p(z_n=k|x_n,y_n;theta) of a gllim model with constraints
# % cstr using N associated observations t and y.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% Input %%%%
# %- t (LtxN)               % Training latent variables
# %- y (DxN)                % Training observed variables
# %- in_K (int)             % Initial number of components
# % <Optional>
# %- Lw (int)               % Dimensionality of hidden components (default 0)
# %- maxiter (int)          % Maximum number of iterations (default 100)
# %- in_theta (struct)      % Initial parameters (default NULL)
# %                         | same structure as output theta
# %- in_r (NxK)             % Initial assignments (default NULL)
# %- cstr (struct)          % Constraints on parameters theta (default NULL,'')
# %   - cstr$ct             % fixed value (LtxK) or ''=uncons.
# %   - cstr$cw             % fixed value (LwxK) or ''=fixed to 0
# %   - cstr$Gammat         % fixed value (LtxLtxK) or ''=uncons.
# %                         | or {'','d','i'}{'','*','v'} [1]
# %   - cstr$Gammaw         % fixed value (LwxLwxK) or ''=fixed to I
# %   - cstr$pi             % fixed value (1xK) or ''=uncons. or '*'=equal
# %   - cstr$A             % fixed value (DxL) or ''=uncons.
# %   - cstr$b             % fixed value (DxK) or ''=uncons.
# %   - cstr$Sigma         % fixed value (DxDxK) or ''=uncons.
# %                         | or {'','d','i'}{'','*'} [1]
# %- verb {0,1,2}           % Verbosity (default 1)
# %%%% Output %%%%
# %- theta  (struct)        % Estimated parameters (L=Lt+Lw)
# %   - theta.c (LxK)       % Gaussian means of X
# %   - theta.Gamma (LxLxK) % Gaussian covariances of X
# %   - theta.pi (1xK)      % Gaussian weights of X
# %   - theta.A (DxLxK)     % Affine transformation matrices
# %   - theta.b (DxK)       % Affine transformation vectors
# %   - theta.Sigma (DxDxK) % Error covariances
# %- r (NxK)                % Posterior probabilities p(z_n=k|x_n,y_n;theta) 
# %%% [1] 'd'=diag., 'i'=iso., '*'=equal for all k, 'v'=equal det. for all k
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % ======================Input Parameters Retrieval=========================
# A faire plus tard mais pas forcement indispensable maintenant? 
# [Lw, maxiter, in_theta, in_r, cstr, verb] = ...
#     process_options(varargin,'Lw',0,'maxiter',100,'in_theta',NULL,...
#                              'in_r',NULL,'cstr',struct(),'verb',1);
                      
# % ==========================Default Constraints============================
if(! "ct" %in% names(cstr)) cstr$ct=NULL;
if(! "cw" %in% names(cstr)) cstr$cw=NULL;
if(! "Gammat" %in% names(cstr)) cstr$Gammat=NULL;
if(! "Gammaw" %in% names(cstr)) cstr$Gammaw=NULL;
if(! "pi" %in% names(cstr)) cstr$pi=NULL;
if(! "A" %in% names(cstr)) cstr$A=NULL;
if(! "b" %in% names(cstr)) cstr$b=NULL;
if(! "Sigma" %in% names(cstr)) cstr$Sigma="i";

if (ncol(tapp) != ncol(yapp)) {stop("Observations must be in columns and variables in rows")}

covParBloc_EM =function(yCen,rk_bar, model){
  covY = tcrossprod(yCen)/rk_bar
  Dim = dim(covY)[1]
  covEmpB = matrix(0,nrow = Dim,ncol = Dim)
  for (b in 1:max(model)){
    ind = which(model == b)
    covEmpB[ind,ind] = covY[ind,ind]
  }
  return(covEst = covEmpB )
}

ExpectationZ = function(tapp,yapp,th,verb){
    if(verb>=1) print('  EZ');
    if(verb>=3) print(' k='); 
    D= nrow(yapp)
    N = ncol(yapp) 
    K=length(th$pi);
    Lt = nrow(tapp)
    L = nrow(th$c)
    Lw=L-Lt;
    logr=matrix(NaN,N,K);
    for (k in 1:K){
        if(verb>=3) print(k);    
        muyk=th$b[,k,drop=FALSE]; #% Dx1
        covyk= th$Sigma[,,k]; #% DxD
        if(Lt>0)
            {if (L==1) {Atk=th$A[,1:Lt,k,drop=FALSE];} else {Atk=th$A[,1:Lt,k]} #% DxLt    
            muyk= sweep(Atk%*%tapp,1,muyk,"+");#% DxN
        }
        if(Lw>0)
            {Awk=matrix(th$A[,(Lt+1):L,k,drop=FALSE],ncol=Lw,nrow=D); #% DxLw
            Gammawk=th$Gamma[(Lt+1):L,(Lt+1):L,k]; #% LwxLw
            cwk=th$c[(Lt+1):L,k]; #% Lwx1
            covyk=covyk+Awk%*%Gammawk%*%t(Awk); #% DxD
            muyk=sweep(muyk,1,Awk%*%cwk,"+"); #% DxN
        }
        logr[,k] = log(th$pi[k]);#*rep(1,N); #N x K
        logr[,k] = logr[,k] + loggausspdf(yapp,muyk,covyk); 
        if (Lt>0)
            logr[,k] = logr[,k]+ loggausspdf(tapp,th$c[1:Lt,k,drop=FALSE],th$Gamma[1:Lt,1:Lt,k]);        
    }
    lognormr=logsumexp(logr,2);
    LL=sum(lognormr);
    r=exp(sweep(logr,1, lognormr,"-"));
    
    # % remove empty clusters
    ec=rep(TRUE,K); #% false if component k is empty.
    for (k in 1:K){
        if(sum(r[,k])==0 | !is.finite(sum(r[,k])))
            {ec[k]=FALSE;
            if(verb>=1) {print(paste('     WARNING: CLASS ',k,' HAS BEEN REMOVED'));}  
		}      
    }
    if (sum(ec)==0){
      print('REINIT! ');
      r = emgm(rbind(tapp,yapp), K, 2, verb)$R;
      ec=rep(TRUE,ncol(r));
      } else {
        r=r[,ec,drop=FALSE];
      }
return(list(r=r,LL=LL,ec=ec))
}

ExpectationW=function(tapp,yapp,th,verb){
    if(verb>=1) print('  EW');
    if(verb>=3) print(' k='); 
    D = nrow(yapp) ; N=ncol(yapp)
    K=length(th$pi);
    Lt = nrow(tapp);
    L = nrow(th$c)
    Lw=L-Lt;
    if(Lw==0)
        {muw=NULL;
        Sw=NULL;}
    Sw=array(0,dim=c(Lw,Lw,K));
    muw=array(0,dim=c(Lw,N,K));
    for (k in 1:K){
        if(verb>=3) print(k)        
        Atk=th$A[,1:Lt,k]; #%DxLt 
        Sigmak=th$Sigma[,,k]; #%DxD
        if (Lw==0) 
        {Awk = NULL ; Gammawk=NULL ;cwk =NULL;invSwk=NULL} else {Awk=th$A[,(Lt+1):L,k]; Gammawk=th$Gamma[(Lt+1):L,(Lt+1):L,k];cwk=th$c[(Lt+1):L,k];invSwk=diag(Lw)+tcrossprod(Gammawk,Awk) %*% solve(Sigmak)%*%Awk;}    #%DxLw # gerer le cas ou Lw=0 Matlab le fait tout seul
        if (!is.null(tapp))
            {Atkt=Atk%*%tapp;}
        else
            {Atkt=0;}
        if (Lw==0) {muw=NULL;Sw=NULL;} else {
        	#invSwk\bsxfun(@plus,Gammawk*Awk'/Sigmak*bsxfun(@minus,y-Atkt,th.b(:,k)),cwk)
        muw[,,k]= solve(invSwk,sweep(Gammawk %*% t(Awk) %*% solve(Sigmak) %*% sweep(yapp-Atkt,1,th$b[,k],"-"),1,cwk,"+")); #%LwxN
        Sw[,,k]=solve(invSwk,Gammawk);}
    }
return(list(muw=muw,Sw=Sw))
}


Maximization = function(tapp,yapp,r,muw,Sw,cstr,verb,model){
    if(verb>=1) print('  M'); 
    if(verb>=3) print(' k='); 
    K = ncol(r);
    D = nrow(yapp);N=ncol(yapp)
	Lt = nrow(tapp)
	Lw = ifelse(is.null(muw),0,nrow(muw))
    L=Lt+Lw;
    th = list()
    th$c=matrix(NaN,nrow=L,ncol=K)
    th$Gamma=array(0,dim=c(L,L,K));
    if(Lw>0)
        {th$c[(Lt+1):L,]=cstr$cw; #% LwxK
        th$Gamma[(Lt+1):L,(Lt+1):L,]=cstr$Gammaw;} #% LwxLwxK}
    th$pi=rep(NaN,K);    
    th$A=array(NaN,dim=c(D,L,K));
    th$b=matrix(NaN,nrow=D,ncol=K);
    th$Sigma= array(NaN,dim=c(D,D,K));  
    
    rk_bar=rep(0,K);
    for (k in 1:K){
        if(verb>=3) print(k);    
      #  % Posteriors' sums
        rk=r[,k]; #% 1xN         
        rk_bar[k]=sum(rk); #% 1x1
        
        if(Lt>0)
           {
           	if(verb>=3) {print('c');}
            #% Compute optimal mean ctk  
            if(is.null(cstr$ct))
                {th$c[1:Lt,k]=rowSums(sweep(tapp,2,rk,"*"))/rk_bar[k];}# % Ltx1
            else {th$c[1:Lt,k]=cstr$ct[,k];}
            #% Compute optimal covariance matrix Gammatk
            if(verb>=3) {print('Gt');}
            diffGamma <- sweep(sweep(tapp,1,th$c[1:Lt,k],"-"),2,sqrt(rk),"*");    #% LtxN
            if( is.null(cstr$Gammat) || (length(cstr$Gammat)==1 & cstr$Gammat=='*')) # | ou ||?
               # %%%% Full Gammat
                {th$Gamma[1:Lt,1:Lt,k]=tcrossprod(diffGamma)/rk_bar[k]; #% DxD
                }                    
            else
            	{
            		if( !is.character(cstr$Gammat))
                #%%%% Fixed Gammat
                {th$Gamma[1:Lt,1:Lt,k]=cstr$Gammat[,,k];  }          
            		else
            			{
            			if(cstr$Gammat[1]=='d' | cstr$Gammat[1]=='i')
	                		#% Diagonal terms   
		                {gamma2=rowSums(diffGamma^2)/rk_bar[k]; #%Ltx1
		                if(cstr$Gammat[1]=='d')
		                    #%%% Diagonal Gamma
		                    {th$Gamma[1:Lt,1:Lt,k]=diag(gamma2);} #% LtxLt  
		                else
		                    #%%% Isotropic Gamma
		                    {th$Gamma[1:Lt,1:Lt,k]=mean(gamma2)*diag(Lt);} #% LtxLt
	                		}
            			else
		            		{if(cstr$Gammat[1]=='v')
		                #%%%% Full Gamma
		                {th$Gamma[1:Lt,1:Lt,k]=tcrossprod(diffGamma)/rk_bar[k];} #% LtxLt
		            		else {# cstr$Gammat,
		                stop('  ERROR: invalid constraint on Gamma.'); }
	                		}
            		}
			}				
           }        
				
        # % Compute optimal weight pik
        th$pi[k]=rk_bar[k]/N; #% 1x1

        if(Lw>0)
            {x=rbind(tapp,muw[,,k]); #% LxN
            Skx=rbind(cbind(matrix(0,Lt,Lt),matrix(0,Lt,Lw)),cbind(matrix(0,Lw,Lt),Sw[,,k])); }#% LxL    
        else
            {x=tapp; #% LtxN
            Skx=matrix(0,Lt,Lt);} #%LtxLt
        
        if(verb>=3) {print('A');}
        if(is.null(cstr$b))
            {# % Compute weighted means of y and x
            yk_bar=rowSums(sweep(yapp,2,rk,"*"))/rk_bar[k]; #% Dx1
            if(L>0)
                {xk_bar= rowSums(sweep(x,2,rk,"*"))/rk_bar[k];} #% Lx1
            else
                {xk_bar=NULL;}
            }
        else
            {yk_bar=cstr$b[,k];
            xk_bar=rep(0,L);
            th$b[,k]=cstr$b[,k]; 
            } 
        #% Compute weighted, mean centered y and x
        weights=sqrt(rk); #% 1xN        
        y_stark=sweep(yapp,1,yk_bar,"-"); #% DxN #col or row? 
        y_stark= sweep(y_stark,2,weights,"*"); #% DxN  #col or row?     
        if(L>0)
           { x_stark=sweep(x,1,xk_bar,"-"); #% LxN  
            x_stark= sweep(x_stark,2,weights,"*"); #% LxN
            }            
        else
            {x_stark=NULL;}
        
        #% Robustly compute optimal transformation matrix Ak
        #warning off MATLAB:nearlySingularMatrix;
        if(!all(Skx==0))
            {if(N>=L & det(Skx+tcrossprod(x_stark)) >10^(-8))
                {th$A[,,k]=tcrossprod(y_stark,x_stark) %*% qr.solve(Skx+tcrossprod(x_stark));} #% DxL
            else
                {th$A[,,k]=tcrossprod(y_stark,x_stark) %*% ginv(Skx+tcrossprod(x_stark));} #%DxL
            }
        else
        		{if(!all(x_stark==0))
	            {if(N>=L & det(tcrossprod(x_stark))>10^(-8))
	               {th$A[,,k]=tcrossprod(y_stark,x_stark) %*% qr.solve(tcrossprod(x_stark));} #% DxL
	            else
		            {if(N<L && det(crossprod(x_stark))>10^(-8)) 
		               {th$A[,,k]=y_stark %*% solve(crossprod(x_stark)) %*% t(x_stark);} #% DxL
		            else
		                {if(verb>=3) print('p') 
		                th$A[,,k]=y_stark %*% ginv(x_stark);}  #% DxL
		            }}
       		 else
            {#% Correspond to null variance in cluster k or L=0:
            if(verb>=1 & L>0) print('null var\n');
            th$A[,,k]=0; # % DxL
            }
			}  

        if(verb>=3)print('b'); 
       # % Intermediate variable wk=y-Ak*x
        if(L>0)
            {wk=yapp-th$A[,,k]%*%x;} #% DxN  
        else
            {wk=yapp;}

        #% Compute optimal transformation vector bk
        if(is.null(cstr$b))
            th$b[,k]=rowSums(sweep(wk,2,rk,"*"))/rk_bar[k]; #% Dx1 

        if(verb>=3) print('S');
        #% Compute optimal covariance matrix Sigmak
        if(Lw>0)
           { Awk=th$A[,(Lt+1):L,k];
            Swk=Sw[,,k];                
            ASAwk=Awk%*%tcrossprod(Swk,Awk);}
        else
            ASAwk=0;

        diffSigma=sweep(sweep(wk,1,th$b[,k],"-"),2,sqrt(rk),"*"); #%DxN
        
        if (cstr$Sigma %in% c("","*")) 
            {#%%%% Full Sigma  
            th$Sigma[,,k]=tcrossprod(diffSigma)/rk_bar[k]; #% DxD
            th$Sigma[,,k]=th$Sigma[,,k]+ASAwk;  }                  
        else 
        {
	        	if(!is.character(cstr$Sigma))
	            #%%%% Fixed Sigma
	            {th$Sigma=cstr$Sigma;}
	        else {
		        		if(cstr$Sigma[1]=='d' || cstr$Sigma[1]=='i')
		            #% Diagonal terms   
		            {sigma2=rowSums(diffSigma^2)/rk_bar[k]; #%Dx1
			            if(cstr$Sigma[1]=='d')
			                {#%%% Diagonal Sigma
			                th$Sigma[,,k]=diag(sigma2,ncol=D,nrow=D); #% DxD
			                	if (is.null(dim(ASAwk))) {th$Sigma[,,k]=th$Sigma[,,k] + diag(ASAwk,ncol=D,nrow=D)}
			                		else {th$Sigma[,,k]=th$Sigma[,,k]+diag(diag(ASAwk));}    
			                }            
			            else
			                {#%%% Isotropic Sigma
			                th$Sigma[,,k]=mean(sigma2)*diag(D); #% DxD
			                		if (is.null(dim(ASAwk))) {th$Sigma[,,k]=th$Sigma[,,k]+sum(diag(ASAwk,ncol=D,nrow=D))/D*diag(D);}
			                		else {th$Sigma[,,k]=th$Sigma[,,k]+sum(diag(ASAwk))/D*diag(D);}
			                }  
		             }                       
	          else {	
	            if (cstr$Sigma=='bSHOCK') {
	              th$Sigma[,,k] = covParBloc_EM(yCen=diffSigma,rk_bar=rk_bar[k], model=model[[k]]);}
	            else {stop('  ERROR: invalid constraint on Sigma.');}}
	        }
        }
			       
        #% Avoid numerical problems on covariances:
        if(verb>=3) print('n');
        if(! is.finite(sum(th$Gamma[1:Lt,1:Lt,k]))) {th$Gamma[1:Lt,1:Lt,k]=0;}
        th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k]+1e-8*diag(Lt);
        if(! is.finite(sum(th$Sigma[,,k]))) {th$Sigma[,,k]=0;}
        th$Sigma[,,k]=th$Sigma[,,k]+1e-8*diag(D);
        if(verb>=3) print(',');
    } 
    
    if(verb>=3) print('end');

    if (cstr$Sigma=="*")
        {#%%% Equality constraint on Sigma
        th$Sigma=sweep(th$Sigma ,3,rk_bar,"*"); 
        th$Sigma=array(apply(th$Sigma,c(1,2),mean),dim=c(D,D,K)) 
    		}
        
    if( !is.null(cstr$Gammat) && cstr$Gammat=='v')
        {#%%% Equal volume constraint on Gamma
        detG=rep(0,K);
        for (k in 1:K){
        	if (D==1) {detG[k]=th$Gamma[1:Lt,1:Lt,k]}
             else {detG[k]=det(th$Gamma[1:Lt,1:Lt,k]);} #% 1x1
        th$Gamma[1:Lt,1:Lt,k] = th$Gamma[1:Lt,1:Lt,k] / detG[k]
        }
        th$Gamma[1:Lt,1:Lt,]=sum(detG^(1/Lt)*th$pi)*th$Gamma[1:Lt,1:Lt,];
    		}
    
    if(is.character(cstr$Gammat) && !is.null(cstr$Gammat) && cstr$Gammat[length(cstr$Gammat)]=='*')
        {#%%% Equality constraint on Gammat
        for (k in 1:K){
        th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k]%*%diag(rk_bar);    
        th$Gamma[1:Lt,1:Lt,k]=matrix(1,Lt,Lt) * sum(th$Gamma[1:Lt,1:Lt,k])/N;  
        }  
    		}
    		
    if( ! is.character(cstr$pi) || is.null(cstr$pi))
        {if(! is.null(cstr$pi)) {th$pi=cstr$pi;}} else {
    	if (!is.null(cstr$pi) && cstr$pi[1]=='*') 
    	{th$pi=1/K*rep(1,K);} else {stop('  ERROR: invalid constraint on pi.');} 
        }              
return(th)
}

remove_empty_clusters= function(th,cstr,ec){
  if(sum(ec) != length(ec))
  {if( !is.null(cstr$ct) && !is.character(cstr$ct))
    cstr$ct=cstr$ct[,ec];  
  if(!is.null(cstr$cw) && !is.character(cstr$cw))
    cstr$cw=cstr$cw[,ec];    
  if(!is.null(cstr$Gammat) && !is.character(cstr$Gammat))
    cstr$Gammat=cstr$Gammat[,,ec];
  if(!is.null(cstr$Gammaw) && !is.character(cstr$Gammaw))
    cstr$Gammaw=cstr$Gammaw[,,ec];    
  if(!is.null(cstr$pi) && !is.character(cstr$pi))
    cstr$pi=cstr$pi[,ec];       
  if(!is.null(cstr$A) && !is.character(cstr$A))
    cstr$A=cstr$A[,,ec];
  if(!is.null(cstr$b) && !is.character(cstr$b))
    cstr$b=cstr$b[,ec];     
  if(!is.null(cstr$Sigma) && !is.character(cstr$Sigma))
    cstr$Sigma=cstr$Sigma[,,ec];        
  
  th$c=th$c[,ec];
  th$Gamma=th$Gamma[,,ec];
  th$pi=th$pi[ec];
  th$A=th$A[,,ec];
  th$b=th$b[,ec];
  th$Sigma=th$Sigma[,,ec]; 
  }
  return(list(th=th,cstr=cstr))
}

# % ==========================EM initialization==============================
Lt=nrow(tapp)
L=Lt+Lw;
D = nrow(yapp) ; N = ncol(yapp);

if(verb>=1) {print('EM Initializations');}
if(!is.null(in_theta)) {
    theta=in_theta;
    K=length(theta$pi);
    if(is.null(cstr$cw))
        {cstr$cw=matrix(0,L,K);} # % Default value for cw 
    if(is.null(cstr$Gammaw))
       { cstr$Gammaw=array(diag(Lw),dim=c(Lw,Lw,K));}  #% Default value for Gammaw
    tmp = ExpectationZ(tapp,yapp,theta,verb);
    r = tmp$r ;
    ec = tmp$ec ;
    tmp = remove_empty_clusters(theta,cstr,ec);
    theta = tmp$th ;
    cstr = tmp$cstr ;
	tmp = ExpectationW(tapp,yapp,theta,verb);
	muw = tmp$muw
	Sw = tmp$Sw  
    if(verb>=1) print(""); 
    } else {if(is.null(in_r)){ 
			 r = emgm(rbind(tapp,yapp), in_K, 1000, verb=verb)$R; 
    } else {r=in_r$R;}
    
    if(Lw==0) {Sw=NULL; muw=NULL;} else {
        # % Start by running an M-step without hidden variables (partial
        # % theta), deduce Awk by local weighted PCA on residuals (complete
        # % theta), deduce r, muw and Sw from E-steps on complete theta.
        theta = Maximization(tapp,yapp,r,NULL,NULL,cstr,verb);
        #print(colMeans(theta$A)) OK no problem here : error is fater
        K=length(theta$pi);
        if(is.null(cstr$cw))
            {cstr$cw=matrix(0,Lw,K);} 
        theta$c=rbind(theta$c,cstr$cw[,1:K]); 
        Gammaf=array(0,dim=c(L,L,K));
        Gammaf[1:Lt,1:Lt,]=theta$Gamma; 
        if(is.null(cstr$Gammaw))
            {cstr$Gammaw=array(diag(Lw),dim=c(Lw,Lw,K));} 
        Gammaf[(Lt+1):L,(Lt+1):L,]=cstr$Gammaw[,,1:K]; #%LwxLwxK    
        theta$Gamma=Gammaf;    
        # % Initialize Awk with local weighted PCAs on residuals:
        Aw=array(0,dim=c(D,Lw,K));
        for (k in 1:K)
            {rk_bar=sum(r[,k]);
            bk=theta$b[,k];
            w=sweep(yapp,1,bk,"-"); #%DxN  
            if(Lt>0)
                {Ak=theta$A[,,k];
                w=w-Ak%*%tapp;}
            w=sweep(w,2,sqrt(r[,k]/rk_bar),"*"); #%DxN
            C=tcrossprod(w); #% Residual weighted covariance matrix
            tmp = eigen(C) ##svd?
           	U = tmp$vectors[,1:Lw]
          	Lambda = tmp$values[1:Lw] #% Weighted residual PCA U:DxLw
            #% The residual variance is the discarded eigenvalues' mean           
            sigma2k=(sum(diag(C))-sum(Lambda))/(D-Lw); #scalar
            #print(sigma2k) #OK here
            theta$Sigma[,,k]=sigma2k * diag(D);
            Aw[,,k]=U%*%sqrt(diag(Lambda,ncol=length(Lambda),nrow=length(Lambda))-sigma2k*diag(Lw));} 
        theta$A=abind(theta$A,Aw,along=2); #%DxLxK
        tmp = ExpectationZ(tapp,yapp,theta,verb);
        r =tmp$r ;
        ec=tmp$ec;
		tmp = remove_empty_clusters(theta,cstr,ec);
		theta = tmp$th ; 
		cstr = tmp$cstr ;
        tmp = ExpectationW(tapp,yapp,theta,verb);
        muw = tmp$muw ;
        Sw = tmp$Sw ;
        if(verb>=1) print("");      
	}}

# %===============================EM Iterations==============================


if(verb>=1) print('       Running EM');
LL = rep(-Inf,maxiter);
iter = 0;
converged= FALSE;
while ( !converged & iter<maxiter){
    iter = iter + 1;

    if(verb>=1) print(paste('       Iteration ',iter,sep=""));   
    # % =====================MAXIMIZATION STEP===========================
    theta = Maximization(tapp,yapp,r,muw,Sw,cstr,verb,...);    
    
    # % =====================EXPECTATION STEPS=========================== 
    tmp = ExpectationZ(tapp,yapp,theta,verb);
    r = tmp$r ;
    LL[iter] =tmp$LL;
    if (verb>=1) {print(LL[iter]);}
    ec=tmp$ec
	  tmp = remove_empty_clusters(theta,cstr,ec);
	  theta = tmp$th
	  cstr = tmp$cstr
    
    tmp = ExpectationW(tapp,yapp,theta,verb);  
    muw=tmp$muw
    Sw=tmp$Sw
    
    if(iter>=3) {
      deltaLL_total=max(LL[1:iter])-min(LL[1:iter]);
      deltaLL=LL[iter]-LL[iter-1];
      converged=(deltaLL <= (0.001*deltaLL_total));
      }       
    
    if(verb>=1) print("");
	}
# %%% Final log-likelihood %%%%
LLf=LL[iter];

# % =============================Final plots===============================
if(verb>=1) print(paste('Converged in ',iter,' iterations',sep=""));

theta$r = r
theta$LLf=LLf
theta$LL = LL[1:iter]

if (cstr$Sigma == "i") {nbparSigma = 1}
if (cstr$Sigma == "d") {nbparSigma = D}  
if (cstr$Sigma == "") {nbparSigma = D*(D+1)/2}
if (cstr$Sigma == "*") {nbparSigma = D*(D+1)/(2*in_K)} 
if (cstr$Sigma == "bSHOCK") {nbparSigma = Inf} 
  
if (!is.null(cstr$Gammat)){
  if (cstr$Gammat == "i") {nbparGamma = 1}
  if (cstr$Gammat == "d") {nbparGamma = Lt}  
  if (cstr$Gammat == "") {nbparGamma = Lt*(Lt+1)/2}
  if (cstr$Gammat == "*") {nbparGamma = Lt*(Lt+1)/(2*in_K)}    
}  
if (is.null(cstr$Gammat)){
  nbparGamma = Lt*(Lt+1)/2  
}

theta$nbpar = (in_K-1) + in_K*(D*L + D + Lt + nbparSigma + nbparGamma)

return(theta)
}