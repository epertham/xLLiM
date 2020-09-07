Kfoldcv_xllim = function(yapp,tapp,func,verb=1,Kfold=10,B=10,...){
  ## func must be a function with func(cov.train,resp.train,cov.test,resp.test,...) and returns a prediction value 
  n = nrow(yapp)
  l = ncol(tapp)
  
  pred = list()
  
  for (b in 1:B){
    if (verb) print(paste0("Repetition ",b,"/",B))
    pred[[b]] = matrix(0,nrow=n,ncol=l)
    
    fold <- createFolds_xllim(1:n,k=Kfold)
    
    for (i in 1:length(fold)){
      if (verb) print(paste0("Fold ",i,"/",length(fold)))
      
      yapp.train = yapp[-fold[[i]],]
      yapp.test = yapp[fold[[i]],]
      tapp.train = tapp[-fold[[i]],]
      tapp.test = tapp[fold[[i]],]
      
      pred[[b]][fold[[i]],] = func(tapp.train,yapp.train,tapp.test,yapp.test,...)
    }
  }
  return(pred)
}