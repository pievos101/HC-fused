perm_features <- function(mat, n){

ids     <- 1:(dim(mat)[2])
n.pat   <- dim(mat)[1]
ids_pat <- 1:n.pat

per_ids_feat <- sample(ids,n,replace=FALSE)
per_ids_pat  <- sample(ids_pat,n.pat,replace=FALSE)

mat[,per_ids_feat] <- mat[per_ids_pat,per_ids_feat]

return(mat)
}


###
library(aricode)
library(SNFtool)
library(cluster)
library(iClusterPlus)
library(clustRviz)
library(PINSPlus)
library(Spectrum)
library(HCfused)

source("~/GitHub/HC-fused/NEMO.R")

source("~/GitHub/HC-fused/sim.R")

mode="HC"
this_method = "ward.D"

n.iter <- 100
my_var <- 1

VVV <- c(1,5,10,50)

OURS   <- vector("list",length(VVV))
SNF    <- vector("list",length(VVV))
CON    <- vector("list",length(VVV))
BAY    <- vector("list",length(VVV))
#SPEC   <- vector("list",length(VVV))
PINS   <- vector("list",length(VVV))
X_NEMO   <- vector("list",length(VVV))


OURS_sil   <- vector("list",length(VVV))
SNF_sil    <- vector("list",length(VVV))
CON_sil    <- vector("list",length(VVV))
BAY_sil    <- vector("list",length(VVV))
#SPEC_sil   <- vector("list",length(VVV))
PINS_sil   <- vector("list",length(VVV))
NEMO_sil   <- vector("list",length(VVV))

truelab <- c(rep(1,4),rep(2,4),rep(3,4))

count <- 0
for(my_perm in VVV){

count <- count + 1 

for(xx in 1:n.iter){
  
  print(my_perm)
  print(xx)
  #HC-fused
  res  <- sim1(FALSE, my_var, mode=mode)
  mat1 <- res$mat1
  mat2 <- res$mat2
  
  #permute
  mat1 <- perm_features(mat1,my_perm)
  mat2 <- perm_features(mat2,my_perm)

  #obj   <- as.list(1:dim(mat2)[1])
  #res   <- HC_fused(obj, mat1, mat2, n.iter=10)  
  #m <- (res$NETWORK/max(res$NETWORK))

  HCres  <- HC_fused_subtyping(omics = list(mat1 , mat2),
  max.k = 9, this_method = "ward.D", HC.iter = 50)
 
  P <- HCres$P  

  hc_final   <- hclust(as.dist(P), method=this_method)
  sil        <- calc.SIL(as.dist(P),9, method=this_method)
  k          <- as.numeric(names(which.max(sil)))
  cl         <- cutree(hc_final,k)
  #print(cl)
  OURS[[count]]     <- c(OURS[[count]],clustComp(cl,truelab)$ARI) #rand.index(cl,truelab)
  OURS_sil[[count]] <- c(OURS_sil[[count]],sil[2]) 
  
  #SNF
  K=7
  alpha=0.5
  TTT=20
  C=3
  
  res   <- sim1(FALSE, my_var)
  Data1 <- res$mat1
  Data2 <- res$mat2
  
  #permute
  Data1 <- perm_features(Data1,my_perm)
  Data2 <- perm_features(Data2,my_perm)

  Data1 = standardNormalization(Data1)
  Data2 = standardNormalization(Data2)
  
  Dist1 = (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2)
  Dist2 = (dist2(as.matrix(Data2),as.matrix(Data2)))^(1/2)
  
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  W = SNF(list(W1,W2), K, TTT)
  
  sil      <- calc.SIL(as.dist(1-W),9)
  #k        <- as.numeric(names(which.max(sil)))
  
  estimationResult <- estimateNumberOfClustersGivenGraph(W, 2:9);
  k <- estimationResult[[1]] 
  
  group = spectralClustering(W,k)
  #print(group)
  
  SNF[[count]]     <-  c(SNF[[count]],clustComp(group,truelab)$ARI)
  SNF_sil[[count]] <-  c(SNF_sil[[count]],sil[2])
 #rand.index(group,truelab)
  
  #concatenate
  res   <- sim1(FALSE, my_var)
  Data1 <- res$mat1
  Data2 <- res$mat2
  
  #permute
  Data1 <- perm_features(Data1,my_perm)
  Data2 <- perm_features(Data2,my_perm)


  Data   <- cbind(Data1,Data2)
  
  hc_con <- hclust(dist(Data))
  sil      <- calc.SIL(dist(Data),9)
  k        <- as.numeric(names(which.max(sil)))
  cl_con   <- cutree(hc_con,k)
  #print(cl)
  CON[[count]]      <- c(CON[[count]],clustComp(cl_con,truelab)$ARI) #rand.index(cl,truelab)
  CON_sil[[count]]  <- c(CON_sil[[count]],sil[2]) 
  
  #iClusterPlus
  #BIC <- rep(NaN,2)
  #for(kk in 1:2){
  # fit.single=iClusterPlus(dt1=Data1,dt2=Data2,type=c("gaussian","gaussian"),
  #                        K=kk,maxiter=10)
  # BIC[kk] <- fit.single$BIC
  #}
  #fit.single=iClusterPlus(dt1=Data1,dt2=Data2,type=c("gaussian","gaussian"),
  #                        K=which.min(BIC),maxiter=10)
  #cl_bayes           <- fit.single$clusters
  #BAY[[count]]       <- c(BAY[[count]],clustComp(cl_bayes,truelab)$ARI) #rand.index(cl,truelab)
  #BAY_sil[[count]]  <- c(BAY_sil[[count]],sil[2]) 
  
 
  ## PINSPLUS

  result  <- try(SubtypingOmicsData(dataList = list(Data1,Data2), iterMin=20,verbose=FALSE),silent=TRUE)
  
  if(class(result) == "try-error"){
  cl_pins <- rep(1,dim(Data1)[1])   
  }else{
  #print(result)
  cl_pins <- result$cluster2
  }
  PINS[[count]]       <- c(PINS[[count]],clustComp(cl_pins,truelab)$ARI)

  #NEMO
  omics_list = list(as.data.frame(t(Data1)),as.data.frame(t(Data2))) 
  cl_nemo    = nemo.clustering(omics_list,num.neighbors=7)
  nag        = nemo.affinity.graph(omics_list, k = 7)
  sil_nemo   <- calc.SIL(as.dist(1-nag),9)
  X_NEMO[[count]]    <- c(X_NEMO[[count]],clustComp(cl_nemo,truelab)$ARI)
  NEMO_sil[[count]]  <- c(NEMO_sil[[count]],sil_nemo[2])
  ####

  #Spectrum#########################
  #res_spec      <- Spectrum(list(as.data.frame(t(Data1)),as.data.frame(t(Data2))), showres=FALSE, silent=TRUE)
  #cl_spec       <- res_spec$assignments
  #SPEC[[count]] <- c(SPEC[[count]],clustComp(cl_spec,truelab)$ARI) #rand.index(cl,truelab)

}

}

#par(mfrow=c(2,4))
#boxplot(OURS, ylim=c(0,1), ylab="ARI", names=as.character(VVV), xlab="sd", main="HC-fused",col="grey")
#boxplot(SNF, ylim=c(0,1), ylab="ARI", names=as.character(VVV), xlab="sd", main="SNF",col="grey")
#boxplot(CON, ylim=c(0,1), ylab="ARI", names=as.character(VVV), xlab="sd", main="Concatenate",col="grey")
#boxplot(BAY, ylim=c(0,1), ylab="ARI", names=as.character(VVV), xlab="sd", main="iClusterPlus",col="grey")

#boxplot(OURS_sil, ylim=c(0,1), ylab="SIL (k=3)", names=as.character(VVV), xlab="sd", main="HC-fused",col="grey")
#boxplot(SNF_sil, ylim=c(0,1), ylab="SIL (k=3)", names=as.character(VVV), xlab="sd", main="SNF",col="grey")
#boxplot(CON_sil, ylim=c(0,1), ylab="SIL (k=3)", names=as.character(VVV), xlab="sd", main="Concatenate",col="grey")
#boxplot(BAY_sil, ylim=c(0,1), ylab="SIL (k=3)", names=as.character(VVV), xlab="sd", main="iClusterPlus",col="grey")

par(mfrow=c(2,3))

plot(sapply(OURS,mean), ylim=c(0,1), ylab="ARI", xaxt="n", xlab="permuted features", main="HC-fused",col="black", type="b", pch=19)
points(sapply(OURS,sd), col="dark grey", type="b", pch=19,lty=2)
axis(1,1:length(VVV),as.character(VVV))
plot(sapply(SNF,mean), ylim=c(0,1), ylab="ARI", xaxt="n", xlab="permuted features", main="SNF",col="blue", type="b", pch=9)
points(sapply(SNF,sd), col="dark grey", type="b", pch=19,lty=2)
axis(1,1:length(VVV),as.character(VVV))
#points(sapply(BAY,mean), ylim=c(0,1), ylab="ARI", xaxt="n", xlab="variance", main="",col="green", type="b", pch=10)
#axis(1,1:length(VVV),as.character(VVV))
plot(sapply(CON,mean), ylim=c(0,1), ylab="ARI", xaxt="n", xlab="permuted features", main="HC-concatenate",col="orange", type="b", pch=11)
points(sapply(CON,sd), col="dark grey", type="b", pch=19,lty=2)
axis(1,1:length(VVV),as.character(VVV))
plot(sapply(PINS,mean), ylim=c(0,1), ylab="ARI", xaxt="n", xlab="permuted features", main="PINSPLUS",col="red", type="b", pch=12)
points(sapply(PINS,sd), col="dark grey", type="b", pch=19,lty=2)
axis(1,1:length(VVV),as.character(VVV))
plot(sapply(X_NEMO,mean), ylim=c(0,1), ylab="ARI", xaxt="n", xlab="permuted features", main="NEMO",col="green", type="b", pch=5)
points(sapply(X_NEMO,sd), col="dark grey", type="b", pch=19,lty=2)
axis(1,1:length(VVV),as.character(VVV))
#plot(sapply(SPEC,mean), ylim=c(0,1), ylab="ARI", xaxt="n", xlab="variance", main="SPECTRUM",col="black", type="b", pch=19)
#axis(1,1:length(VVV),as.character(VVV))


#legend
#legend("bottomleft", c("HC-fused", "SNF","HC-concatenate","PINSPlus"),
#col = c("black", "blue", "orange","red"), text.col = "black", 
#lty = c(1, 1, 1), pch = c(19, 9, 11, 12),
#merge = TRUE, bg = FALSE,bty = "n")


plot(sapply(OURS_sil,mean), ylim=c(0,1), ylab="SIL (k=3)", xaxt="n", xlab="permuted features", main="",col="black", type="b", pch=19)
axis(1,1:length(VVV),as.character(VVV))
points(sapply(SNF_sil,mean), ylim=c(0,1), ylab="SIL (k=3)", xaxt="n", xlab="permuted features", main="",col="blue", type="b", pch=9)
#axis(1,1:length(VVV),as.character(VVV))
#points(sapply(BAY_sil,mean), ylim=c(0,1), ylab="SIL (k=3)", xaxt="n", xlab="variance", main="",col="black", type="b", pch=10)
#axis(1,1:length(VVV),as.character(VVV))
points(sapply(CON_sil,mean), ylim=c(0,1), ylab="SIL (k=3)", xaxt="n", xlab="permuted features", main="",col="orange", type="b", pch=11)
#axis(1,1:length(VVV),as.character(VVV))
points(sapply(NEMO_sil,mean), ylim=c(0,1), ylab="SIL (k=3)", xaxt="n", xlab="permuted features", main="",col="green", type="b", pch=5)
#axis(1,1:length(VVV),as.character(VVV))

#legend
legend("topright", c("HC-fused", "SNF","HC-concatenate","NEMO"),
col = c("black", "blue", "orange","green"), text.col = "black", 
lty = c(1, 1, 1,1), pch = c(19, 9, 11,5),
merge = TRUE, bg = FALSE,bty = "n")

