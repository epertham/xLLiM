kmeans_probs <- function(mat,m){
  tmp <- apply(m,1,function(x){rowSums(sweep(mat,2,x,"-")^2)})
  tmp <- tmp/rowSums(tmp)
  return(tmp)
}