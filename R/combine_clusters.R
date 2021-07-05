combine_clusters <- function(CL_LIST){

n.samples   <- length(CL_LIST[[1]])

cl_final    <- rep(NaN, n.samples)

n.solutions <- length(CL_LIST)


k <- 1
for (xx in 1:n.samples){

 ids <- which(CL_LIST[[1]]==CL_LIST[[1]][xx])

 for (yy in 2:n.solutions){	

	m   <- which(CL_LIST[[yy]]==CL_LIST[[yy]][xx])
	ids <- intersect(ids,m) 

 }
  if(all(is.na(cl_final[ids]))){
		cl_final[ids] <- k
		k <- k + 1
  }
  	
}

return(cl_final)

}