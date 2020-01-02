sim1 <- function(binary=TRUE, my_var=0.1, mode="HC"){
  
  #mat1
  g1 <- matrix(rnorm(2*100,-10,my_var),2 ,100)
  g2 <- matrix(rnorm(2*100,-10,my_var),2 ,100)
  g3 <- matrix(rnorm(4*100,0,my_var), 4 ,100)
  g4 <- matrix(rnorm(2*100,10,my_var), 2 ,100)
  g5 <- matrix(rnorm(2*100,10,my_var), 2 ,100)

  mat  <- rbind(g1,g2,g3,g4,g5)
  mat1 <- mat
  sil  <- calc.SIL(dist(mat),7)
  
  id <- which.max(sil)
  k  <- as.numeric(names(sil)[id])
  
  if(mode=="HC"){
  hc  <- hclust(dist(mat), method="ward.D")
  cl  <- cutree(hc, k)
  }
  if(mode=="convex"){
  carp_fit <- CARP(mat, status=FALSE)
  hc  <- as.hclust(carp_fit)
  cl  <- cutree(hc, k)
  }

  bin1  <- calc.BINARY(cl)
  
  #mat2
  g1 <- matrix(rnorm(4*100,0,my_var),4 ,1000)
  g2 <- matrix(rnorm(4*100,0,my_var),4 ,1000)
  g3 <- matrix(rnorm(2*100,10,my_var),2 ,1000)
  g4 <- matrix(rnorm(2*100,10,my_var),2 ,1000)
  
  
  mat  <- rbind(g1,g2,g3,g4)
  mat2 <- mat
  
  sil <- calc.SIL(dist(mat),7)
  
  id <- which.max(sil)
  k  <- as.numeric(names(sil)[id])
  
  if(mode=="HC"){
  hc  <- hclust(dist(mat), method="ward.D")
  cl  <- cutree(hc, k)
  }
  if(mode=="convex"){
  carp_fit <- CARP(mat, status=FALSE)
  hc  <- as.hclust(carp_fit)
  cl  <- cutree(hc, k)
  }

  
  bin2 <- calc.BINARY(cl)
  
  if(binary){
   return(list(mat1=bin1, mat2=bin2))
  }else{
   return(list(mat1=mat1, mat2=mat2))
  }
  
}
sim2 <- function(binary=TRUE, my_var=0.1, mode="HC"){
  
  #mat1
  g1 <- matrix(rnorm(1*100,-10,my_var),1 ,100)
  g2 <- matrix(rnorm(1*100,-10,my_var),1 ,100)
  g3 <- matrix(rnorm(6*100,0,my_var),6 ,100)
  g4 <- matrix(rnorm(1*100,10,my_var),1 ,100)
  g5 <- matrix(rnorm(1*100,10,my_var),1 ,100)
  
  mat  <- rbind(g1,g2,g3,g4,g5)
  mat1 <- mat
  sil  <- calc.SIL(dist(mat),7)
  
  id <- which.max(sil)
  k  <- as.numeric(names(sil)[id])
  
  if(mode=="HC"){
  hc  <- hclust(dist(mat), method="ward.D")
  cl  <- cutree(hc, k)
  }
  if(mode=="convex"){
  carp_fit <- CARP(mat, status=FALSE)
  hc  <- as.hclust(carp_fit)
  cl  <- cutree(hc, k)
  }

  
  bin1  <- calc.BINARY(cl)
  
  #mat2
  g1 <- matrix(rnorm(2*100,-10,my_var),2 ,1000)
  g2 <- matrix(rnorm(6*100,0,my_var),6 ,1000)
  g3 <- matrix(rnorm(1*100,10,my_var),1 ,1000)
  g4 <- matrix(rnorm(1*100,30,my_var),1 ,1000)
  
  mat  <- rbind(g1,g2,g3,g4)
  mat2 <- mat
  
  sil <- calc.SIL(dist(mat),7)
  
  id <- which.max(sil)
  k  <- as.numeric(names(sil)[id])
  
  if(mode=="HC"){
  hc  <- hclust(dist(mat), method="ward.D")
  cl  <- cutree(hc, k)
  }
  if(mode=="convex"){
  carp_fit <- CARP(mat, status=FALSE)
  hc  <- as.hclust(carp_fit)
  cl  <- cutree(hc, k)
  }

  
  bin2 <- calc.BINARY(cl)
  
  if(binary){
    return(list(mat1=bin1, mat2=bin2))
  }else{
    return(list(mat1=mat1, mat2=mat2))
  }
  
  
}

sim3 <- function(binary=TRUE, my_var=0.1, mode="HC"){
  
  #mat1
  g1 <- matrix(rnorm(1*100,-10,my_var),1 ,100)
  g2 <- matrix(rnorm(1*100,-10,my_var),1 ,100)
  g3 <- matrix(rnorm(6*100,0,my_var),6 ,100)
  g4 <- matrix(rnorm(1*100,10,my_var),1 ,100)
  g5 <- matrix(rnorm(1*100,10,my_var),1 ,100)
  
  mat  <- rbind(g1,g2,g3,g4,g5)
  mat1 <- mat
  sil  <- calc.SIL(dist(mat),7)
  
  id <- which.max(sil)
  k  <- as.numeric(names(sil)[id])
  
  if(mode=="HC"){
  hc  <- hclust(dist(mat), method="ward.D")
  cl  <- cutree(hc, k)
  }
  if(mode=="convex"){
  carp_fit <- CARP(mat, status=FALSE)
  hc  <- as.hclust(carp_fit)
  cl  <- cutree(hc, k)
  }

  
  bin1  <- calc.BINARY(cl)
  
  #mat2
  g1 <- matrix(rnorm(2*100,-10,my_var),2 ,1000)
  g2 <- matrix(rnorm(6*100,0,my_var),6 ,1000)
  g3 <- matrix(rnorm(1*100,30,my_var),1 ,1000)
  g4 <- matrix(rnorm(1*100,40,my_var),1 ,1000)
  g5 <- matrix(rnorm(1*100,50,my_var),1 ,1000)
  
  mat  <- rbind(g1,g2,g3,g4,g5)
  mat2 <- mat
  
  sil <- calc.SIL(dist(mat),7)
  
  id <- which.max(sil)
  k  <- as.numeric(names(sil)[id])
  
  if(mode=="HC"){
  hc  <- hclust(dist(mat), method="ward.D")
  cl  <- cutree(hc, k)
  }
  if(mode=="convex"){
  carp_fit <- CARP(mat, status=FALSE)
  hc  <- as.hclust(carp_fit)
  cl  <- cutree(hc, k)
  }

  
  bin2 <- calc.BINARY(cl)
  
  if(binary){
    return(list(mat1=bin1, mat2=bin2))
  }else{
    return(list(mat1=mat1, mat2=mat2))
  }
  
  
}
