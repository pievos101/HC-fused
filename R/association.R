association <- function(cl){
  
  clustering = cl
  BIN <- matrix(0,length(clustering), length(clustering))
  
  cl  <- unique(clustering)
  
  for(xx in 1:length(cl)){
    
    ids <- which(clustering==cl[xx])
    BIN[ids,ids] <- 1
    
  }
  
  return(BIN)
  
}
