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

# aml is slow
#aml, gbm, lung, sarcoma, colon, liver, ovarian, breast, kidney, melanoma
cancertype <- "colon"
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

#print(table(survival[,3]))

#### CLUSTERING 

## HCfused - 
#res_3        <- HC_fused_subtyping_kNN(list(mRNA,Methy,miRNA),
                          #this_method="kmeans",
#                          HC.iter=20, k=3)

#cl_fused_3  <- res_3$cluster

#res_5        <- HC_fused_subtyping_kNN(list(mRNA,Methy,miRNA),
                          #this_method="kmeans",
#                          HC.iter=20, k=5)

#cl_fused_5  <- res_5$cluster

#res_10        <- HC_fused_subtyping_kNN(list(mRNA,Methy,miRNA),
                          #this_method="kmeans",
#                          HC.iter=20, k=10)

#cl_fused_10  <- res_10$cluster

## HCfused - original
HC.iter=30
res                 <- HC_fused_subtyping(list(mRNA,Methy,miRNA), max.k=10, 
                          this_method="ward.D",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused            <- res$cluster

res_combined        <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
                          this_method="ward.D2",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused_combined_1 <- res_combined$cluster


res_combined        <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
                          this_method="single",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused_combined_2 <- res_combined$cluster

#
res_combined        <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
                          this_method="complete",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused_combined_3 <- res_combined$cluster

res_combined        <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
                          this_method="average",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused_combined_4 <- res_combined$cluster

res_combined        <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
                          this_method="mcquitty",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused_combined_5 <- res_combined$cluster

res_combined        <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
                          this_method="median",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused_combined_6 <- res_combined$cluster

res_combined        <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
                          this_method="centroid",
                          HC.iter=HC.iter, use_opt_code=TRUE)

cl_fused_combined_7 <- res_combined$cluster


#################################################################

groups <- factor(cl_fused)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_fused_combined_1)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_1[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_fused_combined_2)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_2[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_fused_combined_3)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_3[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_fused_combined_4)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_4[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_fused_combined_5)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_5[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_fused_combined_6)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_6[xx] = round(summary(coxFit)$sctest[3],digits = 40);

groups <- factor(cl_fused_combined_7)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED_combined_7[xx] = round(summary(coxFit)$sctest[3],digits = 40);


RESULT <- cbind(P_FUSED, P_FUSED_combined_1,P_FUSED_combined_2,
              P_FUSED_combined_3,P_FUSED_combined_4,P_FUSED_combined_5,
              P_FUSED_combined_6,P_FUSED_combined_7)

colnames(RESULT) <- c("ward.D","ward.D2","single","complete",
                      "average","mcquitty","median","centroid")

print(RESULT)

}#end of loop

RESULT_log <- -log10(RESULT)

#colnames(RESULT_log) <- c("HC_FUSED","HC_FUSED_kNN_1","HC_FUSED_kNN_2",
#  "HC_FUSED_kNN_3","HC_FUSED_kNN_4")

boxplot(RESULT_log, col="grey", ylab="-log10(logrank p-value)", las=2,
  outline=FALSE)
abline(h=-log10(0.05), col="red")

stop("Allet jut!")