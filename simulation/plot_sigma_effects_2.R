 
library(HCfused)
source("~/GitHub/HC-fused/sim.R")

#var c(0.1,0.5,1,5,10,15,20)

###################
#SIM1 ###################################################
####################

var=15
lim <- 60

res  <- sim1(binary=FALSE,my_var=var)

#VIEW1
plot(res$mat1[,1],res$mat1[,2], xlim=c(-lim,lim), ylim=c(-lim,lim),type='n', xlab="feature values", ylab="feature value", main="view 1")
#text(res$mat1[,1],res$mat1[,2],label=1:12)

for(xx in 1:10){
res  <- sim1(binary=FALSE,my_var=var)

for(xx in 1:5){
pairs <- sample(1:100,2)
HCres          <- HC_fused_subtyping(list(res$mat1,res$mat2),
			max.k=10,
			HC.iter=50)
print(HCres$cluster)
coll=c("light blue","blue","light green","grey")

points(res$mat1[,pairs[1]],res$mat1[,pairs[2]], col=coll[HCres$cluster], cex=1.5)
text(res$mat1[,pairs[1]],res$mat1[,pairs[2]],label=1:12, cex=0.5)

}
}


#VIEW2
plot(res$mat2[,1],res$mat2[,2], xlim=c(-lim,lim), ylim=c(-lim,lim),type='n', xlab="feature values", ylab="feature value", main="view 2")
#text(res$mat1[,1],res$mat1[,2],label=1:12)

for(xx in 1:10){
res  <- sim1(binary=FALSE,my_var=var)

for(xx in 1:5){
pairs <- sample(1:100,2)
HCres          <- HC_fused_subtyping(list(res$mat1,res$mat2),
			max.k=10,
			HC.iter=50)
print(HCres$cluster)
coll=c("light blue","blue","light green","grey")

points(res$mat2[,pairs[1]],res$mat2[,pairs[2]], col=coll[HCres$cluster], cex=1.5)
text(res$mat2[,pairs[1]],res$mat2[,pairs[2]],label=1:12, cex=0.5)

}
}

#############################
## SIM2
#############################

var <- 5
lim <- 70

res  <- sim2(binary=FALSE,my_var=var)

#VIEW1
plot(res$mat1[,1],res$mat1[,2], xlim=c(-lim,lim), ylim=c(-lim,lim),type='n', xlab="feature values", ylab="feature value", main="view 1")
#text(res$mat1[,1],res$mat1[,2],label=1:12)

for(xx in 1:10){
res  <- sim2(binary=FALSE,my_var=var)

for(xx in 1:5){
pairs <- sample(1:100,2)
HCres          <- HC_fused_subtyping(list(res$mat1,res$mat2),
			max.k=9,
			HC.iter=50)
print(HCres$cluster)
coll=c("light blue","blue","light green","grey")

points(res$mat1[,pairs[1]],res$mat1[,pairs[2]], col=coll[HCres$cluster], cex=1.5)
text(res$mat1[,pairs[1]],res$mat1[,pairs[2]],label=1:10, cex=0.5)

}
}


#VIEW2
plot(res$mat2[,1],res$mat2[,2], xlim=c(-lim,lim), ylim=c(-lim,lim),type='n', xlab="feature values", ylab="feature value", main="view 2")
#text(res$mat1[,1],res$mat1[,2],label=1:12)

for(xx in 1:10){
res  <- sim2(binary=FALSE,my_var=var)

for(xx in 1:5){
pairs <- sample(1:100,2)
HCres          <- HC_fused_subtyping(list(res$mat1,res$mat2),
			max.k=9,
			HC.iter=50)
print(HCres$cluster)
coll=c("light blue","blue","light green","grey")

points(res$mat2[,pairs[1]],res$mat2[,pairs[2]], col=coll[HCres$cluster], cex=1.5)
text(res$mat2[,pairs[1]],res$mat2[,pairs[2]],label=1:10, cex=0.5)

}
}





