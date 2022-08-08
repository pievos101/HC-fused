# Ensemble clustering with genetic algorithm

# Fitness function
check_ensemble2 <- function(x, methods=FALSE, omics_in=FALSE, fix.k=NaN){
	ens  <- round(x)
	ens  <- methods[ens]
	res  <- HC_fused_subtyping_ens2(omics=omics_in, this_method=ens)
	#print(res$SIL)
	return(res$SIL)
} # end of fitness function


HC_fused_subtyping_ga2 <- function(INPUT=list(), HC.iter=30, max.k=10, fix.k=NaN){

require(GA)

# Available methods
methods = c("single", "complete", "average", "mcquitty", "ward.D",
"ward.D2", "centroid", "median")

# Perform the genetic algorithm
gann <- ga(
	type = "real-valued", 
	fitness = check_ensemble2, methods, INPUT, fix.k=fix.k, 
	lower = rep(1,length(INPUT)), upper = rep(8,length(INPUT)), 
	#seed = 1234, 
	elitism = 20, maxiter = 20, popSize = 20, 
	run = 20, parallel=FALSE)


return(gann)

}# end of function


