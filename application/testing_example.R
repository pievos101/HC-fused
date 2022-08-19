install.packages("devtools")
library(devtools)

small_networks <- as.list(read.csv("networks_small.csv", header = TRUE))
large_networks <- as.list(read.csv("networks_large.csv", header = TRUE))

build()
install()

library(HCFUSEDoptimized)
# This package contains the 6th optimized C++ version of HC_fused

#105 patients
start_time <- Sys.time()
res_cpp_networks_small = matrix(unlist(HC_fused_cpp_opt6(small_networks, 10)),nrow=105,byrow = TRUE)
end_time <- Sys.time()
end_time - start_time


#849 patients
start_time <- Sys.time()
res_cpp_networks_large = matrix(unlist(HC_fused_cpp_opt6(large_networks, 10)),nrow=849,byrow = TRUE)
end_time <- Sys.time()
end_time - start_time

##### The following is needed for the running of the R code

install_github("pievos101/HC-fused")

# Loading the libraries
library(HCfused)
library(fastcluster)

# mRNA
mRNA <- read.table("BREAST_Gene_Expression.txt")
mRNA <- t(mRNA)
dim(mRNA)

# Methy
Methy <- read.table("BREAST_Methy_Expression.txt")
Methy <- t(Methy)
dim(Methy)

# Get the network (binary) structured data (n.patients x n.patients)
omics        <- list(mRNA, Methy)
networks     <- HC_fused_calc_NETWORK(omics)

start_time <- Sys.time()
res_r_networks = HC_fused_new(networks,10)
end_time <- Sys.time()
end_time - start_time

