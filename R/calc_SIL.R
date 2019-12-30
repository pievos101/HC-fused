calc.SIL <- function(mat, size, method="ward.D"){

  DIST <- mat 
  SI   <- numeric(length(2:(size)))
  count<-1
  hc <- fastcluster::hclust(DIST, method=method)
  for(xx in 2:(size)){
   si <- silhouette(cutree(hc,k=xx), DIST)
   SI[count] <- mean(si[,3])
   count <- count + 1 
  }
  names(SI) <- as.character(2:(size))
  return(SI)

}
