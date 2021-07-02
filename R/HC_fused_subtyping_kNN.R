
HC_fused_subtyping_kNN <- function(omics=list(), this_method="ward.D", HC.iter=20, k=5){

max.k=10

kNN_LIST <- vector("list", length(k))

# Calculate Fused Matrix for each k
for (xx in 1:length(k)){

	cat(paste("FUSE for k=",k[xx],"\n", sep=""))
	P 	         <- HC_fused_subtyping_kNN_sub(omics=omics, this_method="ward.D", HC.iter=HC.iter, k=k[xx]) 
	kNN_LIST[[xx]]  <- HC_fused_kNNGraph2(P, k=k[xx]) #Graph2 because Input are Distances
}

if(length(kNN_LIST)>1){
	P <- matrix(unlist(HC_fused_cpp_opt6(kNN_LIST, HC.iter)), nrow=dim(kNN_LIST[[1]])[1], byrow = TRUE)
	# Normalize 
	P <- 1 - (P/max(P))
}else{
	P <- P
}

S <- NULL

#C        <- estimateNumberOfClustersGivenGraph(1-P)[[2]]
#cl_fused <- spectralClustering(1-P,C);

# Cluster the P matrix 
sil_fused  <- calc.SIL(as.dist(P), max.k, method=this_method)
k_fused    <- as.numeric(names(which.max(sil_fused)))

if(this_method=="kmeans"){

 cat("Using kmeans for final clustering ...\n")  
 cl_fused   <- kmeans(as.dist(P), k_fused)$cluster


}else{ # hierarchical clustering

 cat("Using hierarchical clustering for final clustering ...\n")
 hc_fused   <- hclust(as.dist(P), method=this_method)
 cl_fused   <- cutree(hc_fused, k_fused)

}

# S
if(length(S)>0){
 ss       <- apply(S,2,sum)
 for(xx in 1:length(ss)){S[,xx] <- S[,xx]/ss[xx]}

 mm <- paste("omic",1:(dim(S)[1]-1),sep="")

 rownames(S) <- c(mm,"omic_AND")
}

#return(list(cluster=cl_fused, P=P))
return(list(cluster=cl_fused, P=P, S=S, SIL=max(sil_fused)))

}# end of function


HC_fused_subtyping_kNN_sub <- function(omics=list(), this_method="ward.D", max.k=10,
					HC.iter=20, k=5){

#this_method="ward.D" 

parallel=FALSE


omics_binary <- vector("list", length(omics))


# Create kNN Graph
for (xx in 1:length(omics)){
        
	omics_binary[[xx]] <- HC_fused_kNNGraph(omics[[xx]], k) 
}


# Now, fuse these binary matrices 

if(parallel==TRUE){

	cores=detectCores()
	cl <- makeCluster(cores[1]-1) #not to overload your computer
	registerDoParallel(cl)

	res <- foreach(i=1:HC.iter) %dopar% {
		#sink("log.txt", append=TRUE)
		#cat(paste("Starting iteration",i,"\n"))
		HC_fused(omics_binary, n.iter=1)
	}

	P   <- lapply(res,function(x){return(x$NETWORK)})
	P   <- Reduce("+", P)

       S   <- lapply(res,function(x){return(x$SOURCE)})
	S   <- Reduce("+", S)
  
	stopCluster(cl)

}else{

       P <- matrix(unlist(HC_fused_cpp_opt6(omics_binary, HC.iter)), nrow=dim(omics_binary[[xx]])[1], byrow = TRUE)
       S <- NULL
       
}

#print(S)

# Normalize 
P        <- 1 - (P/max(P))

return(P)

}
