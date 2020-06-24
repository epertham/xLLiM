preprocess_data = function(tapp,yapp,in_K,alpha,nfolds){
  init.kmeans = kmeans(cbind(tapp,yapp),in_K)
  ind = c()
  for (k in 1:in_K){
    cv <- cv.glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                    as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian", 
                    alpha=alpha,nfolds = nfolds)
    mod <- glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                  as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian",
                  lambda=cv$lambda.1se, alpha=alpha)
    for (l in 1:dim(tapp)[2]){
      ind = c(ind, which(mod$beta[[l]] !=0))
    }
  }
  ind = unique(ind)
  return(list(selected.variables=ind,clusters=model.matrix(~ -1 + factor(init.kmeans$cluster))))
}