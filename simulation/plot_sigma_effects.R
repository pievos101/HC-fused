 
library(HCfused)
source("~/GitHub/HC-fused/sim.R")

#var c(0.1,0.5,1,5,10,15,20)

#SIM1 ###################################################

par(mfrow=c(2,2))

var=20
lim <- 70

res  <- sim1(binary=FALSE,my_var=var)

plot(res$mat1[,1],res$mat1[,2], xlim=c(-lim,lim), ylim=c(-lim,lim), pch=19, col=c(rep("light blue",4),rep("blue",4),rep("light green",4)), xlab="feature values", ylab="feature value", main="view 1")


for (xx in 1:100){

pairs <- sample(1:100,2)

points(res$mat1[,pairs[1]],res$mat1[,pairs[2]],pch=19, col=c(rep("light blue",4),rep("blue",4),rep("light green",4)))


}


## add the actual points 
mu <- apply(res$mat1,1,mean)
text(mu,mu, label=1:12, cex=0.7)

#legend("topleft",as.character(1:12),
#	pch=1:12, bty = "n",cex = 0.80)
legend("bottomright",c("cluster 1", "cluster 2","cluster 3"),
	col=c("light blue","blue","light green"),	
	pch=19, bty = "n",cex = 0.80)


plot(res$mat2[,1],res$mat2[,2], xlim=c(-lim,lim), ylim=c(-lim,lim), pch=19, col=c(rep("light blue",4),rep("light blue",4),rep("light green",4)), xlab="feature values", ylab="feature value", main="view 2")


for (xx in 1:100){

pairs <- sample(1:100,2)

points(res$mat2[,pairs[1]],res$mat2[,pairs[2]], pch=19, col=c(rep("light blue",4),rep("light blue",4),rep("light green",4)))


}

## add the actual points 
mu <- apply(res$mat2,1,mean)
text(mu,mu, label=1:12, cex=0.7)

#legend("topleft",as.character(1:12),
#	pch=1:12, bty = "n",cex = 0.80)
legend("bottomright",c("cluster 1", "cluster 2"),
	col=c("light blue","light green"),	
	pch=19, bty = "n",cex = 0.80)


#SIM2 ####################################################
##########################################################
##########################################################
#var c(0.1,0.5,1,5,10,15,20)

par(mfrow=c(2,2))

var=20
lim <- 70


res  <- sim2(binary=FALSE,my_var=var)

plot(res$mat1[,1],res$mat1[,2], xlim=c(-lim,lim), ylim=c(-lim,lim), pch=19, col=c(rep("light blue",2),rep("blue",6),rep("light green",2)), xlab="feature values", ylab="feature value", main="view 1")


for (xx in 1:100){

pairs <- sample(1:100,2)

points(res$mat1[,pairs[1]],res$mat1[,pairs[2]], pch=19, col=c(rep("light blue",2),rep("blue",6),rep("light green",2)))


}


## add the actual points 
mu <- apply(res$mat1,1,mean)
text(mu,mu, label=1:10, cex=0.7)


#legend("topleft",as.character(1:10),
#	pch=1:10, bty = "n", cex=0.8)
legend("bottomright",c("cluster 1", "cluster 2","cluster 3"),
	col=c("light blue","blue","light green"),	
	pch=19, bty = "n",cex = 0.80)


plot(res$mat2[,1],res$mat2[,2], xlim=c(-lim,lim), ylim=c(-lim,lim), pch=19, col=c(rep("light blue",2),rep("blue",6),rep("light grey",2)), xlab="feature values", ylab="feature value", main="view 2")


for (xx in 1:100){

pairs <- sample(1:100,2)

points(res$mat2[,pairs[1]],res$mat2[,pairs[2]], pch=19, col=c(rep("light blue",2),rep("blue",6),rep("light grey",2)))


}

## add the actual points 
mu <- apply(res$mat2,1,mean)
text(mu,mu, label=1:10, cex=0.7)

#legend("topleft",as.character(1:10),
#	pch=1:10, bty = "n", cex=0.8)
legend("bottomright",c("cluster 1", "cluster 2"),
	col=c("light blue","blue"),	
	pch=19, bty = "n",cex = 0.80)

