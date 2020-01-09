HC_fused_subtyping <- function(omics=list(), max.k=10, this_method="ward.D", HC.iter=10, parallel=FALSE){

omics_binary <- vector("list", length(omics))

# First, convert the omics data sets into binary views

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

        res <- HC_fused(omics_binary, n.iter = HC.iter)
        P   <- res$NETWORK
	S   <- res$SOURCE

}

#print(S)

# Normalize 
P        <- 1 - (P/max(P))

ss       <- apply(S,2,sum)
for(xx in 1:length(ss)){S[,xx] <- S[,xx]/ss[xx]}


# Cluster the P matrix 
sil_fused  <- calc.SIL(as.dist(P),max.k, method=this_method)
k_fused    <- as.numeric(names(which.max(sil_fused)))

hc_fused   <- hclust(as.dist(P), method=this_method)
cl_fused   <- cutree(hc_fused, k_fused)

#ASSIGN NAMES

mm <- paste("omic",1:(dim(S)[1]-1),sep="")

rownames(S) <- c(mm,"omic_AND")

return(list(cluster=cl_fused, P=P, S=S))

}# end of function
