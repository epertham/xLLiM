bllim = function(tapp,yapp,in_K,in_r=NULL,ninit=20,maxiter=100,verb=0,in_theta=NULL,plot=TRUE){
  # %%%%%%%% General EM Algorithm for Block diagonal gaussian Locally Linear Mapping %%%%%%%%%
  # %%% Author: ED, MG, EP (April 2017) - emeline.perthame@pasteur.fr %%%
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
  
  if (ncol(tapp) != ncol(yapp)) {stop("Observations must be in columns and variables in rows")}
  
  
  # % ==========================EM initialization==============================
  L <- nrow(tapp)
  D <- nrow(yapp) ; N = ncol(yapp);
  
  if (is.null(in_r)){
    ## Step 1 A)
    if (verb) {print("Initialization ... ")}
    ## we perform 10 initialisations of the model
    ## with diagonal Sigma and keep the best initialisation to 
    r = emgm(as.matrix(rbind(tapp, yapp)), init=in_K, 1000, verb=0)
    LLinit = sapply(1:ninit,function(it){
      while (min(table(r$label)) < 2){
        r = emgm(as.matrix(rbind(tapp, yapp)), init=in_K, 1000, verb=0)
      }
      modInit = gllim(as.matrix(tapp),as.matrix(yapp),in_K=in_K,in_r=r,cstr=list(Sigma="d"),verb=0,maxiter = 5,in_theta=in_theta)
      return(list(r=r,LLinit=modInit$LLf))
    })
    ## we select the model maximizing the logliklihood among the 10 models
    ## this model will be used in step 2
    r = LLinit[1,][[which.max(unlist(LLinit[2,]))]]
  } else {r=in_r}
  
  ## Step 1 B)
  ## using the best initialisation selected in step 1 "r"
  ## we rerun the model estimation "modInit"
  ## "modInit" is a glimm model with diagonale matrices cstr=list(Sigma="d")
  modInit = gllim(as.matrix(tapp),as.matrix(yapp),in_K=in_K,in_r=r,cstr=list(Sigma="d"),verb=0,maxiter = maxiter,in_theta=in_theta)
  ## save modInit for further 
  FinalModInit <- modInit
  ## STEP 2 INITIALIZATION Sigma_k
  ## construct the collection of block diagonal matrices for Sigma_k
  ## to get a first estimation of Sigma_k for k in 1,...,K
  ## based on the classification of observations obtained (modInit)
  
  ## STEP 2 A) we extract the classification of individuals (drosophiles) based on the modInit                            
  affecIndModInit = sapply(1:N,function(i)which.max(modInit$r[i,]))
  
  ## STEP 2 B) we compute full sigma_k (the covariance associated to the noise) 
  ## for each group of individuals                             
  listKfullSigma = lapply(1:in_K,function(dum){
    tmp = cov(t(yapp[,affecIndModInit == dum, drop = FALSE])) -
      matrix(modInit$A[,,dum], ncol = L)%*%modInit$Gamma[,,dum]%*%t(matrix(modInit$A[,,dum],ncol = L))
    return(tmp )
  })
  
  ## List of K matrices with full matrix 
  
  ## optional plot of estimated matrix Sigma
  ## lapply(listKfullSigma,image.plot)
  ## lapply(listKfullSigma,hist)
  if (verb) {print("Building model collection ... ")}
  ## STEP 2 C) we threshold each Sigma_k from Sigma_k full to diagonal 
  ## partsSparseSig is a list of size K of matrix (N_struc x D) that contain in each row a partition 
  ## of variables obtained by threshold of the matrix Sigma_k
  partsSparseSig = lapply(listKfullSigma,function(x) do.call(rbind,thresholdAbsSPath2(cov2cor(x))$partitionList))
  ## We clean each matrix to suprress extra complex models
  partsSparseSig = lapply(partsSparseSig,function(x) as.matrix(x[(nrow(x)-min(sapply(partsSparseSig,nrow))+1):nrow(x),]))
  
  ## STEP 3 : Compute Glimm-shock model
  ## it initizalizes the count of the number of model 
  it = 0
  ## LL will contain the log-likelihood for each partition of Sigma_k 
  ## (ED) for each model or for each partition of Sigma_k ? Pour moi le dernier veut dire qu'on a un vecteur de taille K
  LL = rep(0,nrow(partsSparseSig[[1]]))
  ## nbpar will contain the total number of parameters for each partition of Sigma_k
  nbpar = rep(0,nrow(partsSparseSig[[1]]))
  ## save the number of obs in each cluster, used for cleaning after
  effectifs = matrix(0,nrow=in_K,ncol=nrow(partsSparseSig[[1]]))
  ## save the size of the largest clique for each partition
  Pg  = matrix(0,nrow=in_K,ncol=nrow(partsSparseSig[[1]]))
  ## Collect the corresponding estimated glimm model with shock
  modelSparSig <- list()
  ## Collect the corresponding partition used to estimate glimm model with shock
  strucSig <- list()
  
  if (verb) {print("Running BLLiM ... ")}
  ## this loop estimate a gllim model with each possible partition of (Sigma_k), from the more simple (each matrix is diagonal) to the more complex (each matrix is full)
  ## indPart index the partition of variables describing the structure of sigma 
  if (verb) {pb <- progress_bar$new(total = nrow(partsSparseSig[[1]]))}
  for (indPart in 1:nrow(partsSparseSig[[1]])){ 
    # if (verb) print(paste0("Model ",indPart))
    if (verb) {pb$tick()}
    ## strucSig is a list of size K (for each groups of individuals)
    ## containing the partition of variables associated to cluster k 
    ## given by thresholding of sigma_k  
    strucSig[[indPart]] = lapply(1:in_K,function(dum) partsSparseSig[[dum]][indPart,])
    
    ## estimation of parameters taking into account
    ## the variable to predict as.matrix(tapp)
    ## the covariables as.matrix(yapp)
    ## K number of clusters of individuals
    ## r estimation of parameters from gaussian mixture model for initialization
    ## model describes the partition of variables for each clusters
    modelSparSig[[indPart]] = gllim(as.matrix(tapp),as.matrix(yapp),in_K=in_K,in_r=list(R=modInit$r),Lw=0,maxiter=maxiter,verb=0,cstr=list(Sigma="bSHOCK"),model = strucSig[[indPart]])
    ## size of blocks 
    pg = lapply(strucSig[[indPart]],table)
    ##  number of parameters in each block     
    dimSigma = do.call(sum,lapply(pg,function(x)x*(x+1)/2))
    ## K = dim(modShock$Sigma)[3]
    ## par(mar = rep(2, 4),mfrow=c(1,5))  
    ## for (k in 1:K)image.plot(modShock$Sigma[,,k]) 
    LL[indPart] = modelSparSig[[indPart]]$LLf
    nbpar[indPart] = (in_K-1) + in_K*(D*L + D + L + L*(L+1)/2) + dimSigma
    effectifs[,indPart] = modelSparSig[[indPart]]$pi*N
    ## size of the largest clique in each class
    Pg[,indPart] = sapply(pg,max)
  }
  
  
  
  ## put the interesting objet in the big lists 
  BigModelSparSig <- modelSparSig ## modShock ## partsSparseSig is more interesting than modShock?? # Used for cleaning partsSparseSig # Used after selection by capushe to run the selected model
  BigStrucSig <- strucSig
  ### BigModShock[[k]] <- modShock ### devenu inutile non? 
  BigLL<- LL # Used for capushe after
  BigNbpar <- nbpar # Used for capushe after
  BigEffectifs <- effectifs # Used for cleaning SparseSigColl
  BigPg <- Pg # Used for cleaning SparseSigColl 
  
  
  ## STEP 4 : select the right level of sparsity for sigma_k
  ## i.e. supressing the structure of sigma_k matrix that would 
  ## be unrealistic to estimate, given the size of the cluster
  
  ## if the size of the cluster k (number of individuals)
  ## is too small to estimate sigma_k, we suppress the model
  supressModel = list()
  for(kk in 1:in_K){
    pg = Pg[kk,]
    supressModel[[kk]] <- which(pg > effectifs[kk,])
  }
  ## supressModel contains all models to supress
  supressModel <- unique(unlist(supressModel))
  
  ## once we have the vector supressModel, we actually clean the model collection "modelSparSig"
  ## the collection of partitions "strucSig" to keep only the right model
  if(length(supressModel)>1){
    modelSparSigClean <- modelSparSig[-supressModel]
    strucSigClean <- strucSig[-supressModel]
    nbparClean <- nbpar[-supressModel]
    LLClean <- LL[-supressModel]
  } else {
    modelSparSigClean <- modelSparSig
    strucSigClean <- strucSig
    nbparClean <- nbpar
    LLClean <- LL  
  }
  
  BigModelSparSigClean <- modelSparSigClean 
  BigStrucSigClean <- strucSigClean
  BigLLClean <- LLClean # Used for capushe after
  BigNbparClean <- nbparClean # Used for capushe after
  
  ## STEP 5 : model selection based on the slope heuristic
  
  ## We create TabForCapClean, an input for capushe
  ## to indice the model, we simply use the model dimension
  ## to avoid mistake when selected the right model
  
  TabForCapClean <- cbind(nbparClean,nbparClean,nbparClean,-LLClean)
  
  ## select the id of model
  if (length(nbparClean)>10) {
    resCap <- capushe(TabForCapClean,n=N,psi.rlm = "lm")
    id_model_final = which(nbparClean==as.integer(resCap@DDSE@model))
    FinalStrucSig = strucSigClean[id_model_final]
    FinalModelSparSig = modelSparSigClean[id_model_final]
    FinalLL = LLClean[id_model_final]
    FinalNbpar = nbparClean[id_model_final]
    
    if (plot) {
      x <- resCap
      scoef=x@Djump@ModelHat$Kopt/x@Djump@ModelHat$kappa[x@Djump@ModelHat$JumpMax+1]
      
      
      leng=length(x@Djump@graph$model)
      mleng=length(x@Djump@ModelHat$model_hat)
      
      Intervalslope=scoef*x@DDSE@interval$interval
      Absmax=max((scoef+1)/scoef*x@Djump@ModelHat$Kopt,5/4*Intervalslope[2])
      
      
      Ordmin=max(which((x@Djump@ModelHat$kappa<=Absmax)==TRUE))
      plot(x=x@Djump@ModelHat$Kopt,y=x@Djump@graph$complexity[x@Djump@graph$Modopt],xlim=c(0,Absmax),ylim=c(x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[Ordmin]],x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[1]]),ylab="Model dimension",xlab=expression(paste("Values of the penalty constant ",kappa)),yaxt="n",pch=4,col="blue",main="",lwd=3,xaxt="n")
      complex=x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[1:Ordmin]]
      for (i in 1:(mleng-1)){
        lines(x=c(x@Djump@ModelHat$kappa[i],x@Djump@ModelHat$kappa[i+1]),y=c(complex[i],complex[i]))
      }
      if (x@Djump@ModelHat$kappa[i+1]<Absmax){
        lines(x=c(x@Djump@ModelHat$kappa[mleng],Absmax),y=c(complex[mleng],complex[mleng]))
      }
      lines(x=c(x@Djump@ModelHat$Kopt/scoef,x@Djump@ModelHat$Kopt/scoef),y=c(complex[x@Djump@ModelHat$JumpMax],complex[x@Djump@ModelHat$JumpMax+1]),col="blue",lty=2)
      ordon=paste(as.character(complex),"(",as.character(x@Djump@graph$model[x@Djump@ModelHat$model_hat[1:Ordmin]]),")")
      par(cex.axis=0.6)
      
      complex2=complex
      pas=(max(complex2)-min(complex2))/39*5.6/par("din")[2]
      leng2=length(complex2)
      Coord=c()
      i=leng2
      j=leng2-1
      while (j>0){
        while ((j>1)*((complex2[j]-complex2[i])<pas)){
          Coord=c(Coord,-j)
          j=j-1
        }
        if ((j==1)*((complex2[j]-complex2[i])<pas)){
          Coord=c(Coord,-j)
          j=j-1
        }
        i=j
        j=j-1
      }
      
      if (length(Coord)>0){
        complex2=complex2[Coord]
        ordon=ordon[Coord]
      }
      axis(2,complex2,las=2)
      
      
      par(cex.axis=1)
      axis(2,labels="Model dimension",outer=TRUE,at=(x@Djump@graph$complexity[leng]+x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[Ordmin]])/2,lty=0,cex=3)
      
      axis(1,labels=expression(hat(kappa)^{dj}),at=x@Djump@ModelHat$Kopt/scoef)
      axis(1,labels=  expression(kappa[opt]),at=x@Djump@ModelHat$Kopt)
      
      
      absci=seq(0,Absmax,by=x@Djump@ModelHat$Kopt/scoef/3)[c(-4)]
      Empty=c(which(abs(absci-x@Djump@ModelHat$Kopt)<absci[2]/4*(12.093749/par("din")[1])^2*scoef),which(abs(absci-x@Djump@ModelHat$Kopt/scoef)<absci[2]/4*(12.093749/par("din")[1])^2*scoef))
      absci=absci[-Empty]
      axis(1,absci,signif(absci,2))
      
      plength=length(x@DDSE@graph$model)-1
      p=plength-x@DDSE@interval$point_using+2
      
      Model=x@DDSE@ModelHat$model_hat[x@DDSE@ModelHat$point_breaking[x@DDSE@ModelHat$imax]]
      
      Couleur=c(rep("blue",(p-1)),rep("red",(plength-p+2)))
      plot(x=x@DDSE@graph$pen,y=-x@DDSE@graph$contrast,main="",xlab="Model dimension",ylab="log-likelihood",col=Couleur,pch=19,lwd=0.1)
      
      labelpen=x@DDSE@graph$pen
      labelmodel=x@DDSE@graph$model
      pas=(max(labelpen)-min(labelpen))/50*11.5/par("din")[1]
      leng=length(labelpen)
      Coord=c()
      i=1
      j=2
      while (j<leng){
        while ((j<leng)*((labelpen[j]-labelpen[i])<pas)){
          Coord=c(Coord,-j)
          j=j+1
        }
        i=j
        j=j+1
      }
      if (length(Coord)>0){
        labelpen=labelpen[Coord]
        labelmodel=labelmodel[Coord]
      }
      #axis(1,labelpen,labelmodel,las=2)
      
      Cof=x@DDSE@graph$reg$coefficients
      abline(a=Cof[1],b=Cof[2],col="red",lwd=2)
      abscisse=seq(plength+1,2)
    }
    
  } else {
    bicfinal  = -2*LLClean + log(N)*nbparClean
    if (plot) {plot(bicfinal ~ nbparClean,pch=16,col="blue",main="Plot of BIC criterion",xlab="# of parameters",ylab="BIC")}
    # id_model_final = which.max(LLClean)
    id_model_final = which.min(bicfinal)
    FinalStrucSig = strucSigClean[id_model_final]
    FinalModelSparSig = modelSparSigClean[id_model_final]
    FinalLL = LLClean[id_model_final]
    FinalNbpar = nbparClean[id_model_final] 
  }
  
  ##################################################
  ## STEP 3 : perform prediction 
  ##################################################
  ModFinal = FinalModelSparSig[[1]]
  ModFinal$nbpar = FinalNbpar
  
  return(ModFinal)
}