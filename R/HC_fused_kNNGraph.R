# Calculation of a kNN-Graph

HC_fused_kNNGraph <- function(M, k){

DDist    <- as.matrix(dist(M)) # dist can be optimized

RANK     <- apply(DDist,1,rank, ties.method="random")

knearest <- apply(RANK,2, function(x){which(x<=k)})

# init binary matrix
BMatrix <- matrix(0, dim(knearest)[2], dim(knearest)[2])

for(xx in 1:dim(BMatrix)[1]){

	BMatrix[knearest[,xx], knearest[,xx]] <- 1

}

return(BMatrix)

}

HC_fused_kNNGraph2 <- function(DDist, k){

#DDist    <- as.matrix(dist(M)) # dist can be optimized

RANK     <- apply(DDist,1,rank, ties.method="random")

knearest <- apply(RANK,2, function(x){which(x<=k)})

# init binary matrix
BMatrix <- matrix(0, dim(knearest)[2], dim(knearest)[2])

for(xx in 1:dim(BMatrix)[1]){

	BMatrix[knearest[,xx], knearest[,xx]] <- 1

}

return(BMatrix)

}
