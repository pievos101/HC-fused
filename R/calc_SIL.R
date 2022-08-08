calc.SIL <- function(mat, size, fix.k=NaN, method="ward.D"){

  if(method=="kmeans"){

   cat("Using kmeans SIL calculation ...\n")

   DIST <- mat 
   SI   <- numeric(length(2:(size)))
   count<-1
   for(xx in 2:(size)){
    cl <- kmeans(DIST, xx)$cluster
    si <- silhouette(cl, DIST)
    SI[count] <- mean(si[,3])
    count <- count + 1 
   }
   names(SI) <- as.character(2:(size))
   return(SI)

  }else{ # hierarchical clustering

   if(is.na(fix.k)){
    DIST  <- mat 
    SI    <- numeric(length(2:(size)))
    count <-1
    hc <- fastcluster::hclust(DIST, method=method)
    for(xx in 2:(size)){
      si <- silhouette(cutree(hc,k=xx), DIST)
      SI[count] <- mean(si[,3])
      count <- count + 1 
    }
    names(SI) <- as.character(2:(size))
    return(SI)
    # fixed.k
   }else{
    DIST  <- mat 
    hc <- fastcluster::hclust(DIST, method=method)
    si <- silhouette(cutree(hc,k=fix.k), DIST)
    SI <- mean(si[,3])
    names(SI) <- as.character(fix.k)
    return(SI)
   }

  }


}
