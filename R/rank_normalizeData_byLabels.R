rank_normalizeData_byLabels <- function(X,labels,full=FALSE){
  
  X = as.matrix(X)
  ls = unique(labels)
  for (l in ls){
    il = which(labels == l)
    if (full){
      X[il,] = sweep(X[il,],2,colMeans(X[il,]),FUN="-")%*%inverse_sqrt(cov(X[il,]))
    } else{
      for (i in seq_len(ncol(X))){
        if (sd(X[il,i]) == 0){
          X[il,i] = rank(X[il,i])
          X[il,i] = X[il,i]-mean(X[il,i])
        } else{
          X[il,i] = rank(X[il,i])
          X[il,i] = (X[il,i]-mean(X[il,i]))/sd(X[il,i]) # mean 0, sd 1
        }
      }
    }
  }
  
  
  return(X)
}
