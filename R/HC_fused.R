
# HC_fused  ###############################################################
#HC_fused <- function(obj, MAT=list(), n.iter=100, print=FALSE){

HC_fused <- function(MAT=list(), n.iter=100, use_opt_code=TRUE, print=FALSE){
  
  if(use_opt_code==TRUE){

    P <- matrix(unlist(HC_fused_cpp_opt6(MAT, n.iter)), 
                nrow=dim(omics_binary[[1]])[1], byrow = TRUE)
    S <- NULL

    return(list(NETWORK=P, SOURCE=S))
  }


  obj   <- as.list(1:dim(MAT[[1]])[1])

  n.patients <- dim(MAT[[1]])[1]
  
  #init 
  NETWORK <- matrix(0,n.patients,n.patients)
  SOURCE  <- vector("list",n.patients)
  
  #init 
  sss <- sapply(1:(length(MAT)),function(x){paste("mat",x, sep="")})
  sss <- c(sss,"matAND")

  for(xx in 1:length(SOURCE)){
    SOURCE[[xx]] <- sss
  }
  
  #print(SOURCE)  

  for(xx in 1:n.iter){
  
    cat(xx," of ",n.iter,"done! \n")  
    
    obj_dyn <- obj
    K       <- vector("list",length(obj))
    K[[1]]  <- obj
    # build the tree
    count   <- 2
  
   while(length(obj_dyn)>1){
    
    distances  <- HC_fused_calc_distances(obj_dyn, MAT)
    #distances <- HC_calc_distances(obj_dyn, mat1, mat2)
    
    #print(distances)
    
    #MAP matrix
    pairs        <- combn(length(obj_dyn),2)
    #MAP_matrix_1 <- rbind(rep(1,dim(distances)[2]),rep(2,dim(distances)[2]),rep(3,dim(distances)[2]))
    MAP_matrix_1 <- t(sapply(1:dim(distances)[1],function(x){rep(x,dim(distances)[2])}))

    MAP_matrix_2 <- apply(pairs,2,function(x){paste(x[1],"-",x[2], sep="")})
    colnames(distances) <- MAP_matrix_2

    
    #print(distances)   

    MAP_matrix_2    <- matrix(0,dim(distances)[2],dim(distances)[1])
    MAP_matrix_2[,] <- colnames(distances)
    MAP_matrix_2    <- t(MAP_matrix_2)
    #MAP_matrix_2 <- rbind(MAP_matrix_2,MAP_matrix_2,MAP_matrix_2)
    
    #print(MAP_matrix_2)   
    #stop()
 

    #print(MAP_matrix_1)
    #print(MAP_matrix_2)
    
    if(print){
     print(distances)
    }
      
    #get the min-distance/max-similarities 
    ids   <- which(distances==max(distances))
    #ids2  <- which(distances==1,arr.ind = TRUE) #this worked also
    ids2  <- which(distances==max(distances), arr.ind = TRUE)
    
    #print(ids)
    #print(ids2)
    
    is3   <- which(ids2[,1]==dim(distances)[1])    
    
    #print(is3)
    
    if(length(is3)!=0){
      
      if(length(is3)>1){
        id_min       <- ids[sample(is3,1)] 
      }else{
        id_min       <- ids[is3] 
      }
      
      #print(id_min)    
      
    }else{
      
      if(length(ids)>1){  
        id_min       <- sample(ids,1)
      }else{
        id_min       <- ids
      }
    }
    
    map_info_mat  <- MAP_matrix_1[id_min]
    map_info_pair <- as.numeric(strsplit(MAP_matrix_2[id_min],"-")[[1]])
    
    #print(MAP_matrix_1)
    #print("---------------")
    
    
    #print(map_info_mat)
    #print("---------------")
    
    if(print){
     print(map_info_pair)
     print("---------------")
    }
    
    #merge 
    ele <- c(obj_dyn[[map_info_pair[1]]],obj_dyn[[map_info_pair[2]]])
    obj_dyn[[map_info_pair[1]]] <- ele
    
    #save the source
    for(xx in ele){
     
       SOURCE[[xx]] <- c(SOURCE[[xx]],sss[map_info_mat])
    
    } 


    obj_dyn[[map_info_pair[2]]] <- NULL
   
    if(print){ 
     print(obj_dyn)
    }
      
    #Increment when samples are in the same cluster 
    for(xx in 1:length(obj_dyn)){
      x <- obj_dyn[[xx]]
      if(length(x)>1){
        NETWORK[x,x]<-NETWORK[x,x] + 1
        }
    }
    
    K[[count]] <- obj_dyn 
    count <- count + 1
    
   } 
   
} 
  
  #Get the contribution  
  NUM_SOURCE  <- sapply(SOURCE,table)

  return(list(NETWORK=NETWORK, SOURCE=NUM_SOURCE))
  
}

