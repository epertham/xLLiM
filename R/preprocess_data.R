preprocess_data = function(tapp,yapp,in_K,...){
  L <- ncol(tapp)
  init.kmeans = kmeans(cbind(tapp,yapp),in_K)
  ind = c()
  
  probs <- kmeans_probs(cbind(tapp,yapp), init.kmeans$centers)
  for (k in 1:in_K){
    
    cv <- cv.glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                    as.matrix(tapp[init.kmeans$cluster== k,]), family="mgaussian",...)
    
    lambda <- min(cv$lambda.1se, max(cv$lambda[cv$nzero >= L]))
    mod <- glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                  as.matrix(tapp[init.kmeans$cluster== k,]), family="mgaussian",
                  lambda=lambda, ...)
    indk <- c()
    for (l in 1:dim(tapp)[2]){
      indk = c(indk, which(mod$beta[[l]] !=0))
    }
    indk <- unique(indk)
    
    if (length(indk) == 0){
      mod <- glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                    as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian",
                    lambda=cv$lambda, ...)
      mod <- glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                    as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian",
                    lambda=max(cv$lambda[mod$dfmat[1,] >= 1]), ...)
    }
    
    for (l in 1:dim(tapp)[2]){
      ind = c(ind, which(mod$beta[[l]] !=0))
    }
  }
  ind <- unique(ind)
  if (in_K == 1) {
    clusters <- data.frame(cluster1 = init.kmeans$cluster)
  } else {
    clusters <- probs
  }
  colnames(clusters) <- paste0("cluster",1:ncol(clusters)) 
  return(list(selected.variables=sort(ind),clusters=clusters))
}
