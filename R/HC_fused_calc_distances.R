HC_fused_calc_distances <- function(obj, MAT=list()){
  
 pairs     <- combn(length(obj),2)
  
 DISTANCES <- matrix(NaN, length(MAT), choose(length(obj),2))    
  
 for(xx in 1:length(MAT)){

  mat <- MAT[[xx]]
  mat_distances <- apply(pairs,2,function(x){
    
    o1 <- obj[[x[1]]]
    o2 <- obj[[x[2]]]
    
    o1_zero <- sum(mat[o1,]==0)
    o1_one  <- sum(mat[o1,]==1)
    
    o2_zero <- sum(mat[o2,]==0)
    o2_one  <- sum(mat[o2,]==1)
    
    res     <- mean(mat[o1,o2]) #+ mean(mat1[o1,o1]) + mean(mat1[o2,o2])  
    
    return(res)
  })

 DISTANCES[xx,] <- mat_distances 

 }#end of for loop
  

  #matAND UND MATRIX
  matAND <- Reduce("&", MAT)
  
  matAND[matAND] <- 1
  matAND_distances <- apply(pairs,2,function(x){
    
    o1 <- obj[[x[1]]]
    o2 <- obj[[x[2]]]
    
    o1_zero <- sum(matAND[o1,]==0)
    o1_one  <- sum(matAND[o1,]==1)
    
    o2_zero <- sum(matAND[o2,]==0)
    o2_one  <- sum(matAND[o2,]==1)
    
    res     <- mean(matAND[o1,o2]) #+ mean(mat12[o1,o1]) + mean(mat12[o2,o2])
    
    return(res)
    
  })
  
  return(rbind(DISTANCES,matAND_distances))    
}

