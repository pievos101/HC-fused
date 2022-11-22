# Parea
Parea <- function(omics=list(), this_method="ward.D", HC.iter=20, max.k=10, fix.k=NaN, type=1){

    if(type==1){

        ens_LIST <- vector("list", length(this_method))

        # Calculate Fused Matrix for each method
        for (xx in 1:length(this_method)){

            cat(paste("FUSE for method=",this_method[xx],"\n", sep=""))

            res              <- HC_fused_subtyping(omics=omics, this_method=this_method[xx], HC.iter=HC.iter, max.k=max.k, fix.k=fix.k, use_opt_code=TRUE) 
            ens_LIST[[xx]]   <- calc.BINARY(res$cluster)	
        }

        if(length(ens_LIST)>1){
            # Fuse
            P <- matrix(unlist(HC_fused_cpp_opt6(ens_LIST, HC.iter)), nrow=dim(ens_LIST[[1]])[1], byrow = TRUE)
            # Normalize 
            P <- 1 - (P/max(P))
        }else{
            P <- P
        }

        S <- NULL

        # Cluster fused matrix and combine solutions
        CL_LIST     <- vector("list", length(this_method))

        for(xx in 1:length(this_method)){
        if(is.na(fix.k)){
            sil_fused  <- calc.SIL(as.dist(P), max.k, method=this_method[xx])
            k_fused    <- as.numeric(names(which.max(sil_fused)))
        }else{
            sil_fused  <- calc.SIL(as.dist(P), max.k, fix.k=fix.k, method=this_method[xx])
            k_fused    <- fix.k
        }
        hc_fused   <- hclust(as.dist(P), method=this_method[xx])
        cl_fused   <- cutree(hc_fused, k_fused)
        CL_LIST[[xx]]  <- cl_fused
        }

        cl_fused <- combine_clusters(CL_LIST)


        # S
        if(length(S)>0){
        ss       <- apply(S,2,sum)
        for(xx in 1:length(ss)){S[,xx] <- S[,xx]/ss[xx]}

        mm <- paste("omic",1:(dim(S)[1]-1),sep="")

        rownames(S) <- c(mm,"omic_AND")
        }

        #return(list(cluster=cl_fused, P=P))
        return(list(cluster=cl_fused, P=P, S=S, SIL=max(sil_fused)))
    }# end of type 1

    if(type==2){

        ens_LIST <- vector("list", length(this_method))

        # Calculate Fused Matrix for each method
        for (xx in 1:length(this_method)){

            cat(paste("FUSE for method=",this_method[xx],"\n", sep=""))

            res              <- HC_fused_subtyping(omics=list(omics[[xx]]), this_method=this_method[xx], HC.iter=HC.iter, max.k=max.k, fix.k=fix.k, use_opt_code=TRUE) 
            ens_LIST[[xx]]   <- calc.BINARY(res$cluster)	
        }

        if(length(ens_LIST)>1){
            # Fuse
            P <- matrix(unlist(HC_fused_cpp_opt6(ens_LIST, HC.iter)), nrow=dim(ens_LIST[[1]])[1], byrow = TRUE)
            # Normalize 
            P <- 1 - (P/max(P))
        }else{
            P <- P
        }

        S <- NULL

        # Cluster fused matrix and combine solutions
        CL_LIST     <- vector("list", length(this_method))

        for(xx in 1:length(this_method)){
            if(is.na(fix.k)){
                sil_fused  <- calc.SIL(as.dist(P), max.k, method=this_method[xx])
                k_fused    <- as.numeric(names(which.max(sil_fused)))
            }else{
                sil_fused  <- calc.SIL(as.dist(P), max.k, fix.k=fix.k, method=this_method[xx])
                k_fused    <- fix.k
            }
        hc_fused   <- hclust(as.dist(P), method=this_method[xx])
        cl_fused   <- cutree(hc_fused, k_fused)
        CL_LIST[[xx]]  <- cl_fused
        }

        cl_fused <- combine_clusters(CL_LIST)


        # S
        if(length(S)>0){
        ss       <- apply(S,2,sum)
        for(xx in 1:length(ss)){S[,xx] <- S[,xx]/ss[xx]}

        mm <- paste("omic",1:(dim(S)[1]-1),sep="")

        rownames(S) <- c(mm,"omic_AND")
        }

        #return(list(cluster=cl_fused, P=P))
        return(list(cluster=cl_fused, P=P, S=S, SIL=max(sil_fused)))

    }# end of type 2

}# end of parea

DeepParea <- function(omics=list(), neurons=10, pooling=2, this_method="ward.D", HC.iter=20, max.k=10, fix.k=NaN){

    count <- 1

    # Create new imput
    omics_shuffle <- vector("list", neurons)

    for(xx in 1:neurons){

        id <- sample(1:length(omics),1)
        omics_shuffle[[xx]] <- omics[[id]]
        #feat_ids <- sample(1:dim(omics[[id]])[2],dim(omics[[id]])[2], replace = TRUE)
        #omics_shuffle[[xx]] <- omics_shuffle[[xx]][,feat_ids] 
    }
   
    ens_LIST <- vector("list", neurons)

    # Calculate Binary Matrix for each method and view
    for (xx in 1:neurons){

        #cat(paste("Calc Binary Matrix for method=", this_method[xx],"\n", sep=""))

        res              <- HC_fused_subtyping(omics=list(omics_shuffle[[xx]]), this_method=this_method[count], HC.iter=HC.iter, max.k=max.k, fix.k=fix.k, use_opt_code=TRUE)
        count <- count + 1  
        ens_LIST[[xx]]   <- calc.BINARY(res$cluster)	
    }

    # Start the pooled fusion
    ## create the grouping

    while(TRUE){

    groups  <- split(1:neurons, sort(rep(1:(neurons/pooling), pooling)))
    #print(groups)
    CL_LIST <- vector("list", length(groups))

    for(xx in 1:length(groups)){
    
        #print(groups)
        #print(ens_LIST)
        ens_LIST_group = ens_LIST[groups[[xx]]]
        #print(ens_LIST_group)

        # Fuse
        P <- matrix(unlist(HC_fused_cpp_opt6(ens_LIST_group, HC.iter)), nrow=dim(ens_LIST_group[[1]])[1], byrow = TRUE)
        # Normalize 
        P <- 1 - (P/max(P))
        # Cluster the fused matrix 
        if(is.na(fix.k)){
            sil_fused  <- calc.SIL(as.dist(P), max.k, method=this_method[count])
            k_fused    <- as.numeric(names(which.max(sil_fused)))
        }else{
            sil_fused  <- calc.SIL(as.dist(P), max.k, fix.k=fix.k, method=this_method[count])
            k_fused    <- fix.k
        }
        hc_fused   <- hclust(as.dist(P), method=this_method[count])
        count <- count + 1
        cl_fused   <- cutree(hc_fused, k_fused)
        CL_LIST[[xx]]  <- cl_fused
        #print(count)
    
    }

    #print(CL_LIST)
    #print(length(CL_LIST))
    if(length(CL_LIST)<=3){
        cl_fused <- combine_clusters(CL_LIST)
        S <- NULL
        #print(count)
    return(list(cluster=cl_fused, P=P, S=S, SIL=max(sil_fused)))
    }else{
        # Calculate Binary Matrix for each method and view
        ens_LIST <- vector("list", length(CL_LIST))
        for (yy in 1:length(CL_LIST)){
            ens_LIST[[yy]]   <- calc.BINARY(CL_LIST[[yy]])	
        }
        neurons = length(ens_LIST)
        #print(neurons)
    }
    
    }# while loop

}# end of DeepParea
#Call
#res <- DeepParea(list(view1, view2), this_method=rep("ward.D",20))