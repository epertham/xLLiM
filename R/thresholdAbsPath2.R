thresholdAbsSPath2 = function (Sabs) {
  labelsPath <- list()
  valThres <- Sabs[upper.tri(Sabs)]
  orderValue <- c(0,valThres[order(valThres)])
  labelsPath <- lapply(orderValue,function(lambdaR){
    E <- Sabs
    E[Sabs > lambdaR] <- 1
    E[Sabs < lambdaR] <- 0
    E[Sabs == lambdaR] <- 0
    goutput <- graph.adjacency(E, mode = "undirected", weighted = NULL)
    return(clusters(goutput)$membership)
  })
  return(list(partitionList = unique(labelsPath), lambdaPath = unique(orderValue)))
}