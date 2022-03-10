################### Dumb comparison to the mean #############
mean_cv = function(trainx,trainy,testx,testy){
  pred <- colMeans(trainx)
  pred <- matrix(pred,ncol=ncol(trainx),nrow=nrow(testy),byrow=TRUE)
  return(pred)
}

######################### randomForest ##########################
randomForest_cv = function(trainx,trainy,testx,testy){
  pred = matrix(0,ncol=ncol(trainx),nrow=nrow(testx))
  for (k in 1:ncol(trainx)){
    mod = randomForest(x=trainy,y=trainx[,k])
    pred[,k] = predict(mod,testy)
  } 
  return(pred)
}

################### LASSO #############
lasso_cv <- function(trainx,trainy,testx,testy){
  cv <- cv.glmnet(as.matrix(trainy),as.matrix(trainx),family="mgaussian")
  mod <- glmnet(as.matrix(trainy),as.matrix(trainx),family="mgaussian",lambda=cv$lambda.min)
  pred <- predict(mod,as.matrix(testy))
  return(pred[,,1])
}

####################### spline regression #######################
mars_cv = function(trainx,trainy,testx,testy){
  mod = mars(trainy,trainx)
  testy = data.frame(testy)
  pred = predict(mod,testy)
  return(pred)
}

############################## svm #############################
svm_cv = function(trainx,trainy,testx,testy,kernel="linear",type="eps-regression"){
  pred = matrix(0,ncol=ncol(trainx),nrow=nrow(testx))
  for (k in 1:ncol(trainx)){
    tmp = data.frame(trainy,x=trainx[,k])
    mod = svm(x ~ .,data=tmp,kernel=kernel,type=type)
    testy = data.frame(testy)
    pred[,k] = predict(mod,testy)
  } 
  return(pred)
}

################### BLLiM #############
bllim_cv <- function(trainx,trainy,testx,testy,K,verb=0,alpha, nfolds,...){
  prep_data <- preprocess_data(trainx,trainy,in_K=K,alpha = alpha, nfolds = nfolds)
  mod <- bllim(t(trainx), t(trainy[,prep_data$selected.variables,drop=FALSE]), in_K=K,maxiter=100, in_r=list(R=prep_data$clusters),plot=FALSE,verb=FALSE)
  pred <- gllim_inverse_map(t(testy[,prep_data$selected.variables,drop=FALSE]),mod)$x_exp
  return(t(pred))
}


####################### spls regression #######################
mixOmics_cv = function(trainx,trainy,testx,testy){
  X <- trainy # omics data 
  Y <- trainx # pheno data
  # set range of test values for number of variables to use from trainy dataframe
  list.keepX <- c(seq(20, 50, 5))
  # set range of test values for number of variables to use from Y dataframe
  list.keepY <- c(ncol(Y)) 
  # tune parameters 
  tune.spls.res <- tune.spls(X, Y, ncomp = 2:6,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'regression', measure = 'cor') 
  optimal.keepX <- tune.spls.res$choice.keepX # extract optimal number of variables for X dataframe
  optimal.keepY <- tune.spls.res$choice.keepY # extract optimal number of variables for Y datafram
  optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components
  
  # use all tuned values from above
  final.spls.res <- spls(X, Y, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "regression") # explanitory approach being used
  return(predict(final.spls.res  , newdata=testy)$predict[,,optimal.ncomp])
}


####################### gllim #######################

gllim_cv = function(trainx,trainy,testx,testy,K,Lw=0){
  
  D = nrow(trainy)
  n = length(trainx)
  Lt = ncol(trainx)
    
  ## Run the final model, with the selected values of K and Lw 
  mod = gllim(trainx,trainy,in_K=K,Lw=Lw,cstr=list(Sigma="d"),verb=0)
  for (i in 1:10){
    tmp = gllim(trainx,trainy,in_K=K,Lw=Lw,cstr=list(Sigma="d"),verb=0)
    if (tmp$LLf > mod$LLf) {mod = tmp}
  }
  
  ## Compute prediction accuracy through root mean square error
  pred = gllim_inverse_map(testy,mod,verb=0)$x_exp ## compute prediction 
  
  return(pred[1:Lt,])
}
