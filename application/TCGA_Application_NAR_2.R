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

source("~/GitHub/HC-fused/application/NEMO.R")
source("~/GitHub/HC-fused/application/TCGA_clinical_enrichment.R")

do.LOG <- FALSE
do.PCA <- FALSE

do.SNF <- TRUE
do.PINSPLUS <- TRUE
do.NEMO <- TRUE
do.HCfused <- TRUE
do.PAREA1 <- TRUE
do.PAREA2 <- TRUE

cat("Reading in TCGA data ... \n")

#aml, gbm, lung, sarcoma, colon, liver, ovarian, breast, kidney, melanoma
cancertype <- "kidney"
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

n.iter=30

P_SNF      <- rep(NaN,n.iter)
P_PINS     <- rep(NaN,n.iter)
P_FUSED    <- rep(NaN,n.iter)
#P_SPECTRUM <- rep(NaN,n.iter)
P_NEMO     <- rep(NaN,n.iter)
P_PAREA1   <- rep(NaN,n.iter)
P_PAREA2   <- rep(NaN,n.iter)

K_SNF      <- rep(NaN,n.iter)
K_PINS     <- rep(NaN,n.iter)
K_FUSED    <- rep(NaN,n.iter)
#P_SPECTRUM <- rep(NaN,n.iter)
K_NEMO     <- rep(NaN,n.iter)
K_PAREA1   <- rep(NaN,n.iter)
K_PAREA2   <- rep(NaN,n.iter)


C_SNF      <- rep(NaN,n.iter)
C_PINS     <- rep(NaN,n.iter)
C_FUSED    <- rep(NaN,n.iter)
#P_SPECTRUM <- rep(NaN,n.iter)
C_NEMO     <- rep(NaN,n.iter)
C_PAREA1   <- rep(NaN,n.iter)
C_PAREA2   <- rep(NaN,n.iter)



CLIN_SNF    <- vector("list", n.iter)
CLIN_PINS   <- vector("list", n.iter)
CLIN_FUSED  <- vector("list", n.iter)
CLIN_NEMO   <- vector("list", n.iter)
CLIN_PAREA1 <- vector("list", n.iter)
CLIN_PAREA2 <- vector("list", n.iter)


this_method = "ward.D"

for (xx in 1:n.iter){

cat(xx, "of", n.iter,"\n")

patients <- sample(patientsX,100)

# DELETE again !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#patients <- patientsX

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
#res        <- HC_fused_subtyping(list(mRNA,Methy,miRNA), max.k=10, 
#                          HC.iter=30, use_opt_code=TRUE)

#cl_fused   <- res$cluster

## SNF ######################################################
if(do.SNF){ 
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
  W  = SNF(list(W1,W2,W3), K, NIT)

#Groups with SNF
  C = estimateNumberOfClustersGivenGraph(W)[[1]] #1 is eigen gap
  
  groupSNF = spectralClustering(W,C);
  names(groupSNF)   <- survival$PatientID

  # Clinical enrichment
  CLIN_SNF[[xx]] <-  check.clinical.enrichment(groupSNF, 
                        subtype.name=cancertype)

}

## PINSPLUS
 #result  <- SubtypingOmicsData(dataList = list(mRNA,Methy,miRNA), iterMin=20)
if(do.PINSPLUS){ 
 result  <- SubtypingOmicsData(dataList = list(mRNA,Methy,miRNA), iterMin=20,
             verbose=FALSE, kMax=10, 
             clusteringMethod="kmeans")
  
 
 cl_pins <- result$cluster2
 names(cl_pins)   <- survival$PatientID

 # Clinical enrichment
 CLIN_PINS[[xx]] <-  check.clinical.enrichment(cl_pins, 
                        subtype.name=cancertype)
}  
#NEMO
if(do.NEMO){
  omics_list = list(as.data.frame(t(mRNA)),as.data.frame(t(Methy)),as.data.frame(t(miRNA))) 
  cl_nemo    = nemo.clustering(omics_list,num.neighbors=20)
  names(cl_nemo)   <- survival$PatientID

 # Clinical enrichment
 CLIN_NEMO[[xx]] <-  check.clinical.enrichment(cl_nemo, 
                        subtype.name=cancertype)

}

if(do.HCfused){
  HC.iter=30
  res                 <- HC_fused_subtyping(list(mRNA,Methy,miRNA), max.k=10, 
                          HC.iter=HC.iter, use_opt_code = TRUE)
  cl_fused            <- res$cluster

  names(cl_fused)   <- survival$PatientID


 # Clinical enrichment
 CLIN_FUSED[[xx]] <-  check.clinical.enrichment(cl_fused, 
                        subtype.name=cancertype)

}

if(do.PAREA1){

## HCfused - original
HC.iter=30

# Available methods
methods = c("single", "complete", "average", "mcquitty", "ward.D",
"ward.D2", "centroid", "median")

#print("##################")
#print("GENTIC ALGORITHM")
#print("##################")

# Perform the genetic algorithm
res1       <- HC_fused_subtyping_ga(list(mRNA, Methy, miRNA))

#print(round(res1@solution))
sel       <- methods[round(res1@solution)]

#print("####################")
#print("Selected Methods:")
#print(sel)
#print("####################")


res2       <- HC_fused_subtyping_ens(list(mRNA,Methy,miRNA), 
              max.k=5, 
              this_method=c(sel[1],sel[2]),
              HC.iter=HC.iter)

cl_parea1  <- res2$cluster
names(cl_parea1)   <- survival$PatientID


 # Clinical enrichment
CLIN_PAREA1[[xx]] <-  check.clinical.enrichment(cl_parea1, 
                        subtype.name=cancertype)

}


if(do.PAREA2){

## HCfused - original
HC.iter=30

# Available methods
methods = c("single", "complete", "average", "mcquitty", "ward.D",
"ward.D2", "centroid", "median")

# Perform the genetic algorithm
res1       <- HC_fused_subtyping_ga2(list(mRNA, Methy, miRNA))

sel       <- methods[round(res1@solution)]

res2       <- HC_fused_subtyping_ens2(list(mRNA,Methy,miRNA), 
              max.k=5, 
              this_method=sel,#c(sel[1],sel[2]),
              HC.iter=HC.iter)

cl_parea2  <- res2$cluster
names(cl_parea2)   <- survival$PatientID

# Clinical enrichment
CLIN_PAREA2[[xx]] <-  check.clinical.enrichment(cl_parea2, 
                        subtype.name=cancertype)

}

#Spectrum #########################

#one   <- as.data.frame(t(x_mRNA))
#two   <- as.data.frame(t(x_Methy))
#three <- as.data.frame(t(x_miRNA))

#spec_list <- list(one=one, one=two, three=three)

# res_spec <- Spectrum::Spectrum(spec_list, showres=FALSE, silent=TRUE)
  
# cl_spec       <- res_spec$assignments

#rm(res_spec)
#gc()
###

#surv <- Surv(survival, censor)
#sum.surv <- summary(coxph(surv ~ group))
#c_index <- sum.surv$concordance  


### Survival
if(do.SNF){
groups   <- factor(groupSNF) 
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_snf  <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_snf)
P_SNF[xx] = round(summary(coxFit)$sctest[3],digits = 40);
K_SNF[xx] = length(unique(groups))
C_SNF[xx] = summary(coxFit)$concordance 
}

if(do.PINSPLUS){
groups   <- factor(cl_pins) 
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_pins <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_pins)
P_PINS[xx] = round(summary(coxFit)$sctest[3],digits = 40);
K_PINS[xx] = length(unique(groups))
C_PINS[xx] = summary(coxFit)$concordance 
}

if(do.NEMO){
groups <- factor(cl_nemo)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_nemo <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_NEMO[xx] = round(summary(coxFit)$sctest[3],digits = 40);
K_NEMO[xx] = length(unique(groups))
C_NEMO[xx] = summary(coxFit)$concordance 
}

if(do.HCfused){
groups <- factor(cl_fused)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_nemo <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED[xx] = round(summary(coxFit)$sctest[3],digits = 40);
K_FUSED[xx] = length(unique(groups))
C_FUSED[xx] = summary(coxFit)$concordance 
}

if(do.PAREA1){
groups <- factor(cl_parea1)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_nemo <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_PAREA1[xx] = round(summary(coxFit)$sctest[3],digits = 40);
K_PAREA1[xx] = length(unique(groups))
C_PAREA1[xx] = summary(coxFit)$concordance 
}

if(do.PAREA2){
groups <- factor(cl_parea2)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_nemo <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_PAREA2[xx] = round(summary(coxFit)$sctest[3],digits = 40);
K_PAREA2[xx] = length(unique(groups))
C_PAREA2[xx] = summary(coxFit)$concordance 
}


#RESX <- cbind(P_SNF,P_PINS,P_NEMO,P_FUSED, P_PAREA1, P_PAREA2)
#RESX <- cbind(K_SNF,K_PINS,K_NEMO,K_FUSED,K_PAREA1,K_PAREA2)
RESX <- cbind(C_SNF,C_PINS,C_NEMO,C_FUSED,C_PAREA1,C_PAREA2)
colnames(RESX) <- c("SNF","PINS","NEMO","HCFUSED","PAREA1","PAREA2")

print(CLIN_FUSED)
# What the hack does manipulate the seed?
rm(.Random.seed, envir=globalenv())

}#end of loop


## save clinical enrichment
CLIN_ENRICH = Reduce('rbind',CLIN_FUSED)
write.table(CLIN_ENRICH, paste("HCfused_CLIN_ENRICH2_",cancertype,".txt", sep=""))
CLIN_ENRICH = Reduce('rbind',CLIN_SNF)
write.table(CLIN_ENRICH, paste("SNF_CLIN_ENRICH2_",cancertype,".txt", sep=""))
CLIN_ENRICH = Reduce('rbind',CLIN_PINS)
write.table(CLIN_ENRICH, paste("PINS_CLIN_ENRICH2_",cancertype,".txt", sep=""))
CLIN_ENRICH = Reduce('rbind',CLIN_NEMO)
write.table(CLIN_ENRICH, paste("NEMO_CLIN_ENRICH2_",cancertype,".txt", sep=""))
CLIN_ENRICH = Reduce('rbind',CLIN_PAREA1)
write.table(CLIN_ENRICH, paste("PAREA1_CLIN_ENRICH2_",cancertype,".txt", sep=""))
CLIN_ENRICH = Reduce('rbind',CLIN_PAREA2)
write.table(CLIN_ENRICH, paste("PAREA2_CLIN_ENRICH2_",cancertype,".txt", sep=""))



stop("Allet jut!")

RESULT  <- cbind(P_SNF, P_PINS, P_NEMO)
RESULT_log <- -log10(RESULT)
colnames(RESULT_log) <- c("SNF","PINSplus", "NEMO")

boxplot(RESULT_log, col="grey", ylab="-log10(logrank p-value)", outline=FALSE)
abline(h=-log10(0.05), col="red")

stop("Allet jut!")

library(ggplot2)
library(reshape)
## Paper plots

par(mfrow=c(2,2))
GBM_OTHER <- read.table("GBM_OTHER.txt")
colnames(GBM_OTHER) <- c("SNF","PINSplus","NEMO")
GBM_ENS <- read.table("GBM_ENSEMBLE.txt")
KIDNEY_OTHER <- read.table("KIDNEY_OTHER.txt")
colnames(KIDNEY_OTHER) <- c("SNF","PINSplus","NEMO")
KIDNEY_ENS <- read.table("KIDNEY_ENSEMBLE.txt")
LIVER_OTHER <- read.table("LIVER_OTHER.txt")
colnames(LIVER_OTHER) <- c("SNF","PINSplus","NEMO")
LIVER_ENS <- read.table("LIVER_ENSEMBLE.txt")
SARCOMA_OTHER <- read.table("SARCOMA_OTHER.txt")
colnames(SARCOMA_OTHER) <- c("SNF","PINSplus","NEMO")
SARCOMA_ENS <- read.table("SARCOMA_ENSEMBLE.txt")


boxplot(-log10(cbind(GBM_OTHER,GBM_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",3),rep("cadetblue",4)),
  main="GBM")
abline(h=-log10(0.05), col="red")

boxplot(-log10(cbind(KIDNEY_OTHER,KIDNEY_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",3),rep("cadetblue",4)),
  main="KIRC")
abline(h=-log10(0.05), col="red", main="KIRC")

boxplot(-log10(cbind(LIVER_OTHER,LIVER_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",3),rep("cadetblue",4)),
  main="LIHC")
abline(h=-log10(0.05), col="red", main="LIHC")

boxplot(-log10(cbind(SARCOMA_OTHER,SARCOMA_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",3),rep("cadetblue",4)),
  main="SARC")
abline(h=-log10(0.05), col="red", main="SARC")
