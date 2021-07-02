#https://gdc.cancer.gov/about-data/publications/coadread_2012

#library(dendextend)
#library(foreach)
###
library(aricode)
library(SNFtool)
library(cluster)
#library(iClusterPlus)
#library(clustRviz)
library(PINSPlus)
#library(Spectrum)

library(survival)
library(SNFtool)
library(HCfused)

#library(hcfusedpkg)

source("~/GitHub/HC-fused/NEMO.R")

do.LOG <- FALSE
do.PCA <- FALSE

cat("Reading in TCGA data ... \n")

#aml, gbm, lung, sarcoma, colon, liver, ovarian, breast, kidney, melanoma
cancertype <- "aml"
LOC <- paste("~/TCGA_data/NAR Data/",cancertype,"/", sep="")

#mRNA
mRNAX  <- t(read.table(paste(LOC,"exp", sep="")))

#Methy 
MethyX <- t(read.table(paste(LOC,"methy", sep="")))

#miRNA 
miRNAX <- t(read.table(paste(LOC,"mirna", sep="")))

#CLIN
survivalX <- read.table(paste(LOC,"survival", sep=""), header=TRUE)

# define the patients 
patientsX <- intersect(intersect(rownames(mRNAX),rownames(MethyX)),rownames(miRNAX))

n.iter=20

P_SNF      <- rep(NaN,n.iter)
P_PINS     <- rep(NaN,n.iter)
P_FUSED    <- rep(NaN,n.iter)
#P_SPECTRUM <- rep(NaN,n.iter)
P_NEMO     <- rep(NaN,n.iter)

this_method = "ward.D"

for (xx in 1:n.iter){

cat(xx, "of", n.iter,"\n")

patients <- sample(patientsX,100)

mRNA  <- mRNAX[patients,]
Methy <- MethyX[patients,]
miRNA <- miRNAX[patients,]

if(do.LOG==TRUE){

ids <- which(mRNA==0)
if(length(ids)!=0){
mRNA[ids] <- min(mRNA[-ids])
mRNA <- log(mRNA) 
}

ids <- which(Methy==0)
if(length(ids)!=0){
Methy[ids] <- min(Methy[-ids])
Methy <- log(Methy) 
}

ids <- which(miRNA==0)
if(length(ids)!=0){
miRNA[ids] <- min(miRNA[-ids])
miRNA <- log(miRNA) 
}


}# End of do LOG

#normalization
mRNA  = standardNormalization(mRNA)
Methy = standardNormalization(Methy)
miRNA = standardNormalization(miRNA)

#PCA
if(do.PCA==TRUE){
zero     <- apply(mRNA,2,sum)
zero.ids <- which(zero==0)
if(length(zero.ids)==0){
pca      <- prcomp(mRNA[,], center = TRUE, scale. = TRUE)
}else{
pca      <- prcomp(mRNA[,-zero.ids], center = TRUE, scale. = TRUE)
}
summ     <- summary(pca)
csum     <- cumsum(summ$importance[2,])
id       <- which(csum>0.90)[1]
#mRNA     <- pca$x[,1:id]
mRNA     <- pca$x[,]

zero     <- apply(Methy,2,sum)
zero.ids <- which(zero==0)
if(length(zero.ids)==0){
pca      <- prcomp(Methy[,], center = TRUE, scale. = TRUE)
}else{
pca      <- prcomp(Methy[,-zero.ids], center = TRUE, scale. = TRUE)
}
summ     <- summary(pca)
csum     <- cumsum(summ$importance[2,])
id       <- which(csum>0.90)[1]
#Methy    <- pca$x[,1:id]
Methy    <- pca$x[,]


zero     <- apply(miRNA,2,sum)
zero.ids <- which(zero==0)
if(length(zero.ids)==0){
pca      <- prcomp(miRNA[,], center = TRUE, scale. = TRUE)
}else{
pca      <- prcomp(miRNA[,-zero.ids], center = TRUE, scale. = TRUE)
}
summ     <- summary(pca)
csum     <- cumsum(summ$importance[2,])
id       <- which(csum>0.90)[1]
#miRNA    <- pca$x[,1:id]
miRNA    <- pca$x[,]

}# End of IF PCA

if(is.element(cancertype,c("breast","lung","gbm"))){
 patients <- tolower(patients)
 patients <- strsplit(patients,".", fixed=TRUE)
 patients <- sapply(patients, function(x){paste(x[1:3],collapse=".")})
 ids <- match(patients,as.character(survivalX[,1]))
###########
}else{
#for all other
 ids <- match(gsub(".", "-", patients, fixed=TRUE),as.character(survivalX[,1]))
}


survival <- survivalX[ids,]

#print(table(survival[,3]))

#### CLUSTERING 

## HCfused - HC.iter=20 !!!
res        <- HC_fused_subtyping(list(mRNA,Methy,miRNA), max.k=10, 
                          HC.iter=50, use_opt_code=TRUE, use_kNNGraph=TRUE, k=5)

cl_fused   <- res$cluster

## SNF ######################################################
  ## First, set all the parameters:
  K = 20;##number of neighbors, usually (10~30)
  alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
  NIT = 10; ###Number of Iterations, usually (10~20)

  #datGE = standardNormalization(mRNA)
  #datME = standardNormalization(Methy)
  #datMI = standardNormalization(miRNA)
  
  PSMgeneE = dist2(as.matrix(mRNA),as.matrix(mRNA));
  PSMmethy = dist2(as.matrix(Methy),as.matrix(Methy));
  PSMmir   = dist2(as.matrix(miRNA),as.matrix(miRNA));
  
  W1 = affinityMatrix(PSMgeneE, K, alpha)
  W2 = affinityMatrix(PSMmethy, K, alpha)
  W3 = affinityMatrix(PSMmir, K, alpha)
  W = SNF(list(W1,W2,W3), K, NIT)

#Groups with SNF
  C = estimateNumberOfClustersGivenGraph(W)[[1]]
  
  groupSNF = spectralClustering(W,C);


## PINSPLUS
 #result  <- SubtypingOmicsData(dataList = list(mRNA,Methy,miRNA), iterMin=20)
 result  <- SubtypingOmicsData(dataList = list(mRNA,Methy,miRNA), iterMin=20, verbose=FALSE)
  
 
 cl_pins <- result$cluster2
  
#NEMO
  omics_list = list(as.data.frame(t(mRNA)),as.data.frame(t(Methy)),as.data.frame(t(miRNA))) 
  cl_nemo    = nemo.clustering(omics_list,num.neighbors=20)


#Spectrum#########################

#one   <- as.data.frame(t(x_mRNA))
#two   <- as.data.frame(t(x_Methy))
#three <- as.data.frame(t(x_miRNA))

#spec_list <- list(one=one, one=two, three=three)

# res_spec <- Spectrum::Spectrum(spec_list, showres=FALSE, silent=TRUE)
  
# cl_spec       <- res_spec$assignments

#rm(res_spec)
#gc()
###
  

### Survival

groups   <- factor(groupSNF) 
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_snf  <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_snf)
P_SNF[xx] = round(summary(coxFit)$sctest[3],digits = 40);


groups   <- factor(cl_pins) 
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_pins <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_pins)
P_PINS[xx] = round(summary(coxFit)$sctest[3],digits = 40);


groups <- factor(cl_fused)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_nemo)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_nemo <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_NEMO[xx] = round(summary(coxFit)$sctest[3],digits = 40);

#groups <- factor(cl_spec)
#names(groups) = rownames(survival)
#coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
#cox_spectrum      <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_spectrum)
#P_SPECTRUM[xx] = round(summary(coxFit)$sctest[3],digits = 40);

#print(P_SNF[count])
#print(P_PINS[count])
#print(P_FUSED[count])
#print(P_SPECTRUM[count])

print(cbind(P_SNF,P_PINS,P_NEMO,P_FUSED))


}#end of loop

RESULT  <- cbind(P_SNF, P_PINS, P_NEMO, P_FUSED)
RESULT_log <- -log10(RESULT)
colnames(RESULT_log) <- c("SNF","PINSplus", "NEMO", "HC-FUSED")

boxplot(RESULT_log, col="grey", ylab="-log10(logrank p-value)", outline=FALSE)
abline(h=-log10(0.05), col="red")

stop("Allet jut!")