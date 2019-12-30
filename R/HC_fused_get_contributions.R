HC_fused_get_contributions <- function(S, cluster){

 
 cl.ids    <- unique(cluster)
 n.cluster <- length(cl.ids)


 CONTRIB   <- matrix(NaN, dim(S)[1], n.cluster) 


 for(xx in 1:n.cluster){

 	ids <- which(cluster==cl.ids[xx])
        aa  <- apply(S[,ids, drop=FALSE],1,sum)
        aa  <- aa/sum(aa) 
	CONTRIB[,xx] <- aa
 
 }

nn <- paste("cluster",1:n.cluster,sep="")

mm <- paste("omic",1:(dim(S)[1]-1),sep="")

colnames(CONTRIB) <- nn

rownames(CONTRIB) <- c(mm,"omic_AND")

return(CONTRIB)
}
