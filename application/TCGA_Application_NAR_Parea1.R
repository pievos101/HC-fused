#https://gdc.cancer.gov/about-data/publications/coadread_2012

#library(dendextend)
#library(foreach)
###

library(SNFtool)
library(aricode)
library(cluster)
library(survival)
library(HCfused)
#require(parallel)
source("~/GitHub/HC-fused/application/TCGA_clinical_enrichment.R")

do.sampling <- TRUE

do.LOG <- FALSE
do.PCA <- FALSE

cat("Reading in TCGA data ... \n")

# aml is slow
# melanoma = SKCM
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
CLIN_FUSED   <- vector("list", n.iter)


if(!do.sampling){
  n.iter=1
}

for (xx in 1:n.iter){

cat(xx, "of", n.iter,"\n")

if(do.sampling){
  patients <- sample(patientsX,100)
}else{
  patients <- patientsX
}

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

# Available methods
methods = c("single", "complete", "average", "mcquitty", "ward.D",
"ward.D2", "centroid", "median")

print("##################")
print("GENTIC ALGORITHM")
print("##################")

# Perform the genetic algorithm
res1       <- HC_fused_subtyping_ga(list(mRNA, Methy, miRNA))

print(round(res1@solution))
sel       <- methods[round(res1@solution)]

print("####################")
print("Selected Methods:")
print(sel)
print("####################")

res2       <- HC_fused_subtyping_ens(list(mRNA,Methy,miRNA), 
              max.k=10, 
              this_method=c(sel[1],sel[2]),
              HC.iter=HC.iter)

cl_fused  <- res2$cluster
names(cl_fused)   <- survival$PatientID

#################################################################

groups <- factor(cl_fused)
names(groups) = rownames(survival)
coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups, data = survival, ties="exact")
cox_fused <- round(summary(coxFit)$sctest[3],digits = 40);
#print(cox_fused)
P_FUSED[xx] = round(summary(coxFit)$sctest[3],digits = 40);

print(P_FUSED)

# Clinical enrichment
CLIN_FUSED[[xx]] <-  check.clinical.enrichment(cl_fused, 
                        subtype.name=cancertype)

print(CLIN_FUSED)


}#end of loop

RESULT_log <- -log10(P_FUSED)

boxplot(RESULT_log, col="grey", ylab="-log10(logrank p-value)", las=1,
  outline=FALSE, cex.axis=0.9)
abline(h=-log10(0.05), col="red")

# Clinical enrichment
CLIN_ENRICH = Reduce('rbind',CLIN_FUSED)
write.table(CLIN_ENRICH, paste("Parea1_CLIN_ENRICH_",cancertype,".txt", sep=""))

stop("Allet jut!")

## Plots
# Loading
library("survminer") 
fit<- survfit(Surv(time=Survival, event=Death) ~ groups, data = survival)
# Drawing survival curves
ggsurvplot(fit, pval = TRUE, risk.table = TRUE)
