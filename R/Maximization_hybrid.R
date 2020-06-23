Maximization_hybrid = function(tapp,yapp,r,u,phi,muw,Sw,mahalt,mahaly,cstr,verb){
    if(verb>=1) print('  M'); 
    if(verb>=3) print(' k='); 
    K = ncol(r);
    D = nrow(yapp);N=ncol(yapp)
	Lt = nrow(tapp)
	Lw = ifelse(is.null(muw),0,nrow(muw))
    L=Lt+Lw;
    
    th = list()
    ph = list()
    th$c=matrix(NaN,nrow=L,ncol=K)
    th$Gamma=array(0,dim=c(L,L,K));
    if(Lw>0)
        {th$c[(Lt+1):L,]=cstr$cw; #% LwxK
        th$Gamma[(Lt+1):L,(Lt+1):L,]=cstr$Gammaw;} #% LwxLwxK}
    ph$pi=rep(NaN,K);    
    th$A=array(NaN,dim=c(D,L,K));
    th$b=matrix(NaN,nrow=D,ncol=K);
    th$Sigma= array(NaN,dim=c(D,D,K));  

    rk_bar=rep(0,K);
    for (k in 1:K){
        if(verb>=3) print(k);    
      #  % Posteriors' sums
        rk=r[,k]; #% 1xN         
        rk_bar[k]=sum(rk); #% 1x1

		uk=u[,k];  #% 1xN
		rk_tilde = rk * uk;
		rk_bar_tilde = sum(rk_tilde);

        if(Lt>0)
           {
           	if(verb>=3) {print('c');}
            #% Compute optimal mean ctk  
            if(is.null(cstr$ct))
                {th$c[1:Lt,k]=rowSums(sweep(tapp,2, rk_tilde,"*"))/rk_bar_tilde;}# % Ltx1 
            else {th$c[1:Lt,k]=cstr$ct[,k];}
			#% Compute optimal covariance matrix Gammatk
            if(verb>=3) {print('Gt');}
            diffGamma= sweep(sweep(tapp,1,th$c[1:Lt,k],"-"),2,sqrt(rk_tilde),"*");    #% LtxN 
            if( is.null(cstr$Gammat) || (length(cstr$Gammat)==1 & cstr$Gammat=='*')) # | ou ||?
               # %%%% Full Gammat
                {th$Gamma[1:Lt,1:Lt,k]=tcrossprod(diffGamma)/rk_bar[k]; #% DxD
                }#th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k];                    
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
           
		#% Compute optimal weight pik
		ph$pi[k]=rk_bar[k]/N; #% 1x1

        if(Lw>0)
            {x=rbind(tapp,muw[,,k]); #% LxN      
		Skx=rbind(cbind(matrix(0,Lt,Lt),matrix(0,Lt,Lw)),cbind(matrix(0,Lw,Lt),Sw[,,k])); }#% LxL    
        else
            {x=tapp; #% LtxN
            Skx=matrix(0,Lt,Lt);} #%LtxLt

        if(verb>=3) {print('A');}
        if(is.null(cstr$b))
            {# % Compute weighted means of y and x
            yk_bar=rowSums(sweep(yapp,2,rk_tilde,"*"))/rk_bar_tilde; #% Dx1
            if(L>0)
                xk_bar= rowSums(sweep(x,2, rk_tilde,"*"))/rk_bar_tilde #% Lx1
            else
                {xk_bar=NULL;}
            }
        else
            {yk_bar=cstr$b[,k];
            xk_bar=rep(0,L);
            th$b[,k]=cstr$b[,k]; 
            } 
		#% Compute weighted, mean centered y and x
		weights=sqrt(rk_tilde)/sqrt(rk_bar[k]); #% 1xN  
        y_stark=sweep(yapp,1,yk_bar,"-"); #% DxN #col or row? 
        y_stark= sweep(y_stark,2,weights,"*"); #% DxN  #col or row?     
        if(L>0)
           { x_stark=sweep(x,1,xk_bar,"-"); #% LxN  
            x_stark= sweep(x_stark,2,weights,"*"); #% LxN
            }            
        else
            {x_stark=NULL;}
        
       # % Robustly compute optimal transformation matrix Ak
       # warning off MATLAB:nearlySingularMatrix;
if(!all(Skx==0)) 
            {if(N>=L & det(Skx+tcrossprod(x_stark))>10^(-8))
                {th$A[,,k]=tcrossprod(y_stark,x_stark)%*%solve(Skx+tcrossprod(x_stark));} #% DxL
            else
                {th$A[,,k]=tcrossprod(y_stark,x_stark)%*%ginv(Skx+tcrossprod(x_stark));} #%DxL
            }
        else
        		{if(!all(x_stark==0))
	            {if(N>=L & det(tcrossprod(x_stark))>10^(-8))
	               {th$A[,,k]=tcrossprod(y_stark,x_stark)%*% solve(tcrossprod(x_stark));} #% DxL
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
            {wk=yapp-th$A[,,k]%*%x;} #% DxN #attention au reshape? 
        else
            {wk=yapp;}

        #% Compute optimal transformation vector bk
        if(is.null(cstr$b))
            th$b[,k]=rowSums(sweep(wk,2,rk_tilde,"*"))/rk_bar_tilde; #% Dx1 #col ou row?
 
        if(verb>=3) print('S');
        #% Compute optimal covariance matrix Sigmak
        if(Lw>0)
           { Awk=th$A[,(Lt+1):L,k];
            Swk=Sw[,,k];                
            ASAwk=Awk%*%tcrossprod(Swk,Awk);}
        else
            ASAwk=0;

        diffSigma=sweep(sweep(wk,1,th$b[,k],"-"),2,sqrt(rk_tilde),"*"); #%DxN
                    
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
			else {	cstr$Sigma ;
			            stop('  ERROR: invalid constraint on Sigma.');}
					}
				}
				
		#% Avoid numerical problems on covariances:
         if(verb>=3) print('n');
        if(! is.finite(sum(th$Gamma[1:Lt,1:Lt,k]))) {th$Gamma[1:Lt,1:Lt,k]=0;}
        th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k]+1e-8*diag(Lt);
        if(! is.finite(sum(th$Sigma[,,k]))) {th$Sigma[,,k]=0;}
        th$Sigma[,,k]=th$Sigma[,,k]+1e-8*diag(D);
        if(verb>=3) print(',');
        
        # % Compute phi.alpha %
 if(!is.null(mahalt) && !is.null(mahaly)) { ph$alpha[k] = inv_digamma((digamma(phi$alpha[k] + (D+Lt)/2) - (1/rk_bar[k]) * sum( rk * log(1 + (1/2) * (mahaly[,k] + mahalt[,k]))))); ## potentielle erreur?passage difficile
if(verb>=3) print(paste("K",k,"-> alpha=",ph$alpha[k]));} 
else {ph$alpha = phi$alpha;}

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
        th$Gamma[1:Lt,1:Lt,]=sum(detG^(1/Lt)*ph$pi)*th$Gamma[1:Lt,1:Lt,];
    		}
 
    if(is.character(cstr$Gammat) && !is.null(cstr$Gammat) && cstr$Gammat[length(cstr$Gammat)]=='*')
        {#%%% Equality constraint on Gammat
        for (k in 1:K){
        th$Gamma[1:Lt,1:Lt,k]=th$Gamma[1:Lt,1:Lt,k]%*%diag(rk_bar);    
        th$Gamma[1:Lt,1:Lt,k]=matrix(1,Lt,Lt) * sum(th$Gamma[1:Lt,1:Lt,k])/N;  
        }  
    		}

    if( ! is.character(cstr$pi) || is.null(cstr$pi))
        {if(! is.null(cstr$pi)) {ph$pi=cstr$pi;}} else {
    	if (!is.null(cstr$pi) && cstr$pi[1]=='*') 
    	{ph$pi=1/K*rep(1,K);} else {stop('  ERROR: invalid constraint on pi.');} 
        }    
return(list(ph=ph,th=th))
}
