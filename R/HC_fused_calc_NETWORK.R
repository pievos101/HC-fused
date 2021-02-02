HC_fused_calc_NETWORK <- function(omics=list(), max.k=10, this_method="ward.D"){

omics_binary <- vector("list", length(omics))

# Convert the omics data sets into binary views

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

return(omics_binary)

}
