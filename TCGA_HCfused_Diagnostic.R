
library(HCfused)
library(SNFtool)

#mRNA
mRNAX  <- t(read.table("~/TCGA_data/NAR Data/kidney/exp"))

#Methy 
MethyX <- t(read.table("~/TCGA_data/NAR Data/kidney/methy"))

#miRNA 
miRNAX <- t(read.table("~/TCGA_data/NAR Data/kidney/mirna"))

# define the patients 
patientsX <- intersect(intersect(rownames(mRNAX),rownames(MethyX)),rownames(miRNAX))

patients <- sample(patientsX,100)

mRNA  <- mRNAX[patients,]
Methy <- MethyX[patients,]
miRNA <- miRNAX[patients,]

#normalization
mRNA  = standardNormalization(mRNA)
Methy = standardNormalization(Methy)
miRNA = standardNormalization(miRNA)

#for breast, lung and gbm  !!
#patients <- tolower(patients)
#patients <- strsplit(patients,".", fixed=TRUE)
#patients <- sapply(patients, function(x){paste(x[1:3],collapse=".")})
#ids <- match(patients,as.character(survivalX[,1]))
###########

ITER        <- c(2,3,5,10,20,30,40,50,100,150,200)
HC_CLUST    <- list()
SIL         <- list()

for(xx in 1:length(ITER)){
## HCfused
res             <- HC_fused_subtyping(list(mRNA,Methy,miRNA),
			max.k=20,
			HC.iter=ITER[xx])
HC_CLUST[[xx]]  <- res$cluster
SIL[[xx]]  <- res$SIL
}

n.cl       <- sapply(HC_CLUST,function(x){length(unique(x))})
# PLOTS 
par(mfrow=c(3,3))
plot(n.cl, type="b", pch=19, xaxt="n", ylab="Number of clusters", xlab="Iterations")
axis(1,1:length(ITER), as.character(ITER), las=2)
plot(unlist(SIL), type="b", pch=19, xaxt="n", xlab="Iterations", ylab="SIL")
axis(1,1:length(ITER), as.character(ITER), las=2)
plot(unlist(SIL)/n.cl, type="b", pch=19, xaxt="n", xlab="Iterations", ylab="SIL/n.cluster")
axis(1,1:length(ITER), as.character(ITER), las=2)



