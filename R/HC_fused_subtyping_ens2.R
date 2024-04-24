
HC_fused_subtyping_ens2 <- function(omics=list(), this_method="ward.D", HC.iter=20, max.k=10, fix.k=NaN){

ens_LIST <- vector("list", length(this_method))

# Calculate Fused Matrix for each method
for (xx in 1:length(this_method)){

	cat(paste("FUSE for method=",this_method[xx],"\n", sep=""))

	res              <- HC_fused_subtyping(omics=list(omics[[xx]]), this_method=this_method[xx], HC.iter=HC.iter, max.k=max.k, fix.k=fix.k, use_opt_code=TRUE) 
	ens_LIST[[xx]]   <- calc.BINARY(res$cluster)	
}

if(length(ens_LIST)>1){
	# Fuse
	P <- matrix(unlist(HC_fused_cpp_opt6(ens_LIST, HC.iter)), nrow=dim(ens_LIST[[1]])[1], byrow = TRUE)
	# Normalize 
	P <- 1 - (P/max(P))
}else{
	P <- P
}

S <- NULL

# Cluster fused matrix and combine solutions
CL_LIST     <- vector("list", length(this_method))
SIL_x = rep(NaN, length(this_method))

for(xx in 1:length(this_method)){
	if(is.na(fix.k)){
		sil_fused  <- calc.SIL(as.dist(P), max.k, method=this_method[xx])
		k_fused    <- as.numeric(names(which.max(sil_fused)))
	}else{
		sil_fused  <- calc.SIL(as.dist(P), max.k, fix.k=fix.k, method=this_method[xx])
		k_fused    <- fix.k
		SIL_x[xx]  <- sil_fused
	}
 hc_fused   <- hclust(as.dist(P), method=this_method[xx])
 cl_fused   <- cutree(hc_fused, k_fused)
 CL_LIST[[xx]]  <- cl_fused
}


if(is.na(fix.k)){
 cl_fused <- combine_clusters(CL_LIST)
}else{
 v_id = which.max(SIL_x) #length(CL_LIST) # sample(1:length(CL_LIST), 1)	
 sil_fused = max(SIL_x)
 cl_fused = CL_LIST[[v_id]]
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


