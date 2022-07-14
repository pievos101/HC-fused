# Ensemble clustering with genetic algorithm

# Fitness function
check_ensemble <- function(x, methods=FALSE, omics_in=FALSE){
	ens1 <- round(x[1])
	ens2 <- round(x[2])
	ens1 <- methods[ens1]
	ens2 <- methods[ens2]
	res  <- HC_fused_subtyping_ens(omics=omics_in, this_method=c(ens1, ens2))
	print(res$SIL)
	return(res$SIL)
} # end of fitness function


HC_fused_subtyping_ga <- function(INPUT=list(), HC.iter=30, max.k=10){

require(GA)

# Available methods
methods = c("single", "complete", "average", "mcquitty", "ward.D",
"ward.D2", "centroid", "median")

# Perform the genetic algorithm
gann <- ga(
	type = "real-valued", 
	fitness = check_ensemble, methods, INPUT, lower = c(1,1), upper = c(8,8), 
	#seed = 1234, 
	elitism = 20, maxiter = 20, popSize = 20, 
	run = 20, parallel=FALSE)


return(gann)

}# end of function


