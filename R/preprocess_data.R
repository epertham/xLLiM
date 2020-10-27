preprocess_data = function(tapp,yapp,in_K,...){
  init.kmeans = kmeans(cbind(tapp,yapp),in_K)
  ind = c()
  for (k in 1:in_K){
    print(paste0("k=",k))
    
    cv <- cv.glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                    as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian",...)
    mod <- glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                  as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian",
                  lambda=cv$lambda.1se, ...)
    indk <- c()
    for (l in 1:dim(tapp)[2]){
      indk = c(indk, which(mod$beta[[l]] !=0))
    }
    indk <- unique(indk)
    # print(paste0("indk1=",indk))
    # print(paste0("length indk1=",length(indk)))
    
    if (length(indk) == 0){
      print(paste0("nzero",c(min(cv$nzero),max(cv$nzero))))
      mod <- glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                    as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian",
                    lambda=cv$lambda, ...)
      mod <- glmnet(as.matrix(yapp[init.kmeans$cluster== k,]),
                    as.matrix(tapp[init.kmeans$cluster== k,]),family="mgaussian",
                    lambda=max(cv$lambda[mod$dfmat[1,] >= 1]), ...)
      # indk <- c()
      # for (l in 1:dim(tapp)[2]){
      #   indk = c(indk, which(mod$beta[[l]] !=0))
      # }
      # indk <- unique(indk)
      # print(paste0("indk2=",indk))
      # print(paste0("length indk2=",length(indk)))
    }
    
    for (l in 1:dim(tapp)[2]){
      ind = c(ind, which(mod$beta[[l]] !=0))
    }
  }
  ind = unique(ind)
  return(list(selected.variables=ind,clusters=model.matrix(~ -1 + factor(init.kmeans$cluster))))
}
