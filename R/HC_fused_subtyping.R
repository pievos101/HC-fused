HC_fused_subtyping <- function(omics=list(), max.k=10, this_method="ward.D", 
					HC.iter=20, use_opt_code=TRUE){

parallel=FALSE


omics_binary <- vector("list", length(omics))


 for (xx in 1:length(omics)){

	omic  <- omics[[xx]]	

	sil   <- calc.SIL(dist(omic), max.k, method=this_method)
	id    <- which.max(sil)
 	k     <- as.numeric(names(sil)[id])

	hc    <- hclust(dist(omic), method=this_method)
	cl    <- cutree(hc, k)

	mat   <- calc.BINARY(cl)	
        
	omics_binary[[xx]] <- mat 
       
 }


# Now, fuse thes binary matrices 

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

       if(use_opt_code==FALSE){ 
        res <- HC_fused(omics_binary, n.iter = HC.iter)
        P   <- res$NETWORK
	 S   <- res$SOURCE
       }else{
        P <- matrix(unlist(HC_fused_cpp_opt6(omics_binary, HC.iter)), nrow=dim(omics_binary[[xx]])[1], byrow = TRUE)
        S <- NULL
       }
}

#print(S)

# Normalize 
P        <- 1 - (P/max(P))

# Cluster the P matrix 
sil_fused  <- calc.SIL(as.dist(P),max.k, method=this_method)
k_fused    <- as.numeric(names(which.max(sil_fused)))

cat(paste("Hierarchical clustering with ", this_method, "\n", sep=""))

hc_fused   <- hclust(as.dist(P), method=this_method)
cl_fused   <- cutree(hc_fused, k_fused)

# S
if(length(S)>0){
 ss       <- apply(S,2,sum)
 for(xx in 1:length(ss)){S[,xx] <- S[,xx]/ss[xx]}

 mm <- paste("omic",1:(dim(S)[1]-1),sep="")

 rownames(S) <- c(mm,"omic_AND")
}

return(list(cluster=cl_fused, P=P, S=S, SIL=max(sil_fused)))

}# end of function
