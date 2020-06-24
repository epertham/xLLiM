################### Dumb comparison to the mean #############
mean_cv = function(trainx,trainy,testx,testy){
  pred <- colMeans(trainx)
  pred <- matrix(pred,ncol=ncol(trainx),nrow=nrow(testy),byrow=TRUE)
  return(pred)
}

######################### randomForest ##########################
# randomForest_EP = function(trainx,trainy,testx,testy){
#   norm_term = 1 #mean((testx)^2);
#   mod = randomForest(x=trainy,y=trainx)
#   pred = predict(mod,testy)
#   #return(sqrt((pred-testx)^2/norm_term))   
#   return(pred)
# }
randomForest_cv = function(trainx,trainy,testx,testy){
  pred = matrix(0,ncol=ncol(trainx),nrow=nrow(testx))
  for (k in 1:ncol(trainx)){
    mod = randomForest(x=trainy,y=trainx[,k])
    pred[,k] = predict(mod,testy)
  } 
  return(pred)
}

# library(ranger)
# ranger_EPm = function(trainx,trainy,testx,testy){
#   pred = matrix(0,ncol=ncol(trainx),nrow=nrow(testx))
#   for (k in 1:ncol(trainx)){
#     mod = ranger(y ~ .,data=data.frame(y=trainx[,k],trainy))
#     pred[,k] = predict(mod,data.frame(testy))$predictions
#   } 
#   return(pred)
# }

################### LASSO #############
lasso_cv <- function(trainx,trainy,testx,testy){
  pred <-  matrix(0,ncol=ncol(trainx),nrow=nrow(testx))
  cv <- cv.glmnet(as.matrix(trainy),as.matrix(trainx),family="mgaussian")
  mod <- glmnet(as.matrix(trainy),as.matrix(trainx),family="mgaussian",lambda=cv$lambda.min)
  pred <- predict(mod,as.matrix(testy))
  return(pred)
}

####################### spline regression #######################
mars_cv = function(trainx,trainy,testx,testy){
  mod = mars(trainy,trainx)
  testy = data.frame(testy)
  pred = predict(mod,testy)
  return(pred)
}

############################## svm #############################

# svm_EP = function(trainx,trainy,testx,testy,kernel="linear",type="eps-regression"){
#   norm_term = 1 #mean((testx)^2);
#   tmp = data.frame(trainy,x=trainx)
#   mod = svm(x ~ .,data=tmp,kernel=kernel,type=type)
#   testy = data.frame(testy)
#   pred = predict(mod,testy)
#   #return(sqrt((pred-testx)^2/norm_term))   
#   return(pred)
# }

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
bllim_cv <- function(tapp.train,yapp.train,tapp.test,yapp.test,K,verb=0,alpha, nfolds,...){
  prep_data <- preprocess_data(tapp.train,yapp.train,in_K=K,alpha = alpha, nfolds = nfolds)
  mod <- bllim(t(tapp.train), t(yapp.train[,prep_data$selected.variables]), in_K=K,maxiter=100, in_r=list(R=prep_data$clusters),plot=FALSE,verb=FALSE)
  pred <- gllim_inverse_map(t(yapp.test[,prep_data$selected.variables,drop=FALSE]),mod)$x_exp
  return(t(pred))
}