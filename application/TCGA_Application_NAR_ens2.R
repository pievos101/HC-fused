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

ENS1 <- FALSE
ENS2 <- FALSE
ENS3 <- FALSE
ENS4 <- TRUE


cat("Reading in TCGA data ... \n")

# aml is slow
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

P_FUSED      <- rep(NaN,n.iter)
P_FUSED_3    <- rep(NaN,n.iter)
P_FUSED_5    <- rep(NaN,n.iter)
P_FUSED_10   <- rep(NaN,n.iter)
P_FUSED_combined   <- rep(NaN,n.iter)
P_FUSED_combined_1   <- rep(NaN,n.iter)
P_FUSED_combined_2   <- rep(NaN,n.iter)
P_FUSED_combined_3   <- rep(NaN,n.iter)
P_FUSED_combined_4   <- rep(NaN,n.iter)
P_FUSED_combined_5   <- rep(NaN,n.iter)
P_FUSED_combined_6   <- rep(NaN,n.iter)
P_FUSED_combined_7   <- rep(NaN,n.iter)

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

## HCfused - original
HC.iter=30
if(ENS1){
res                 <- HC_fused_subtyping_ens(list(mRNA,Methy,miRNA), max.k=10, 
                          this_method=c("ward.D","ward.D2"),#Ward
                          HC.iter=HC.iter)

cl_fused            <- res$cluster
}

if(ENS2){
res_combined        <- HC_fused_subtyping_ens(list(mRNA,Methy,miRNA),
                          this_method=c("complete","single"),#Linkage
                          HC.iter=HC.iter) 
                          
cl_fused_combined_1 <- res_combined$cluster
}

if(ENS3){
res_combined        <- HC_fused_subtyping_ens(list(mRNA,Methy,miRNA),
                          this_method=c("centroid","median"),#PGMC
                          HC.iter=HC.iter)

cl_fused_combined_2 <- res_combined$cluster
}

if(ENS4){
res_combined        <- HC_fused_subtyping_ens(list(mRNA,Methy,miRNA),
                          this_method=c("mcquitty","average"),#PGMA
                          HC.iter=HC.iter)

cl_fused_combined_3 <- res_combined$cluster
}

#################################################################

if(ENS1){
groups <- factor(cl_fused)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED[xx] = round(summary(coxFit)$sctest[3],digits = 40);
}

if(ENS2){
groups <- factor(cl_fused_combined_1)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_1[xx] = round(summary(coxFit)$sctest[3],digits = 40);
}

if(ENS3){
groups <- factor(cl_fused_combined_2)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_2[xx] = round(summary(coxFit)$sctest[3],digits = 40);
}

if(ENS4){
groups <- factor(cl_fused_combined_3)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_3[xx] = round(summary(coxFit)$sctest[3],digits = 40);
}

RESULT <- cbind(P_FUSED, P_FUSED_combined_1,P_FUSED_combined_2, P_FUSED_combined_3)

colnames(RESULT) <- c("ward.D-ward.D2","complete-single","centroid-median","mcquitty-average")

print(RESULT)

}#end of loop

RESULT_log <- -log10(RESULT)

#colnames(RESULT_log) <- c("HC_FUSED","HC_FUSED_kNN_1","HC_FUSED_kNN_2",
#  "HC_FUSED_kNN_3","HC_FUSED_kNN_4")

boxplot(RESULT_log, col="grey", ylab="-log10(logrank p-value)", las=1,
  outline=FALSE, cex.axis=0.9)
abline(h=-log10(0.05), col="red")

stop("Allet jut!")


boxplot(-log10(cbind(ALL,ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",8),rep("cadetblue",4)))
abline(h=-log10(0.05), col="red")


## Paper plots

par(mfrow=c(2,2))
GBM_ALL <- read.table("GBM_ALL.txt")
GBM_ENS <- read.table("GBM_ENSEMBLE.txt")
KIDNEY_ALL <- read.table("KIDNEY_ALL.txt")
KIDNEY_ENS <- read.table("KIDNEY_ENSEMBLE.txt")
LIVER_ALL <- read.table("LIVER_ALL.txt")
LIVER_ENS <- read.table("LIVER_ENSEMBLE.txt")
SARCOMA_ALL <- read.table("SARCOMA_ALL.txt")
SARCOMA_ENS <- read.table("SARCOMA_ENSEMBLE.txt")


boxplot(-log10(cbind(GBM_ALL,GBM_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",8),rep("cadetblue",4)),
  main="GBM")
abline(h=-log10(0.05), col="red")

boxplot(-log10(cbind(KIDNEY_ALL,KIDNEY_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",8),rep("cadetblue",4)),
  main="KIRC")
abline(h=-log10(0.05), col="red", main="KIRC")

boxplot(-log10(cbind(LIVER_ALL,LIVER_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",8),rep("cadetblue",4)),
  main="LIHC")
abline(h=-log10(0.05), col="red", main="LIHC")

boxplot(-log10(cbind(SARCOMA_ALL,SARCOMA_ENS)), ylab="-log10(logrank p-value)", 
  las=2,outline=FALSE, cex.axis=0.6, col=c(rep("grey",8),rep("cadetblue",4)),
  main="SARC")
abline(h=-log10(0.05), col="red", main="SARC")


# Speed

TIME <- c(45.84, 49.35, 86.23, 129, 131.2, 137.7, 143, 162.6)
barplot(TIME, ylab="Elapsed Time (Seconds)", 
  name=c("NEMO","SNF","HCfused","UPGMC-WPGMC","ward.D-ward.D2",
          "complete-single","UPGMA-WPGMA","PINSPlus"), las=2, cex.names=0.6, 
  cex.axis=.8)