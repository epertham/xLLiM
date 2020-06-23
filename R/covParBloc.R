covParBloc = function(yCen, model){
  covY = cov(yCen)
  Dim = dim(covY)[1]
  covEmpB = matrix(0,nrow = Dim,ncol = Dim)
  for (b in 1:max(model)){
    ind = which(model == b)
    covEmpB[ind,ind] = covY[ind,ind]
  }
  return(covEst = covEmpB )
}