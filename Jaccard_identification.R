###########################################################################################
#  Jaccard_identification()   function                    2024/10/25  Atsushi Kimura  
#
#   usage:
#       input:        Reference_hclust.cnum,    Original_hclust.cnum  
#       output:       Renumbered_hclust.cnum
#       Renumbered_hclust.cnum <- Jaccard_identification(Reference_hclut.cnum, Original_hclust.cnum)
#
#   example:
#       #Draw Japan Map
#       Ref_hclust.cnum <- cutree(Reference_hclust.result, k=12)
#       Org_hclust.cnum <- cutree(Original_hclust.result, k=12) 
#       Renumbered_hclust.cnum <- Jaccard_identification(Ref_hclust.cnum, Org_hclust.cnum)
#       library(NipponMap)
#       library(RColorBrewer)
#       NipponMap::JapanPrefMap(brewer.pal(12,"Set3")[Ref_hclust.cnum])
#       NipponMap::JapanPrefMap(brewer.pal(12,"Set3")[Org_hclust.cnum])
#       NipponMap::JapanPrefMap(brewer.pal(12,"Set3")[Renumbered_hclust.cnum])
#
################
disassemble_vec <- function(CN){
  MAT_CN <- matrix(data=rep(x=0, times=max(CN)*length(CN)), nrow = max(CN), ncol = length(CN))
  for(i in 1:length(CN)){
    MAT_CN[CN[i],i] <- as.numeric(1)
  }
  return(MAT_CN)
}

################
assemble_vec <- function(MAT_CN){
  ASEM_CN <- rep(x=0, times=ncol(MAT_CN))
  for(i in 1:ncol(MAT_CN)){
    count <- as.numeric(1)
    while(MAT_CN[count,i] == 0){
      count <- count + 1
    }
    ASEM_CN[i] <- count
  }
  return(ASEM_CN)
}
################
diff_MAT <- function(MAT_RefCN, MAT_CN){
  D_MAT_RefCN_CN <- matrix(0,nrow = nrow(MAT_RefCN), ncol = nrow(MAT_RefCN))
  for(j in 1:nrow(MAT_CN)){
    for(i in 1:nrow(MAT_RefCN)){
      D_MAT_RefCN_CN[i,j] <- sum((MAT_RefCN[i,]-MAT_CN[j,])^2)
    }
  }
  return(D_MAT_RefCN_CN)
}  

####################
match_MAT <- function(MAT_RefCN, MAT_CN){
  E_MAT_RefCN_CN <- MAT_RefCN %*% t(MAT_CN)
  return(E_MAT_RefCN_CN) 
} 
####################
Jaccard_MAT <- function(E_MAT_RefCN_CN, D_MAT_RefCN_CN){
  Jaccard_MAT <- E_MAT_RefCN_CN / (E_MAT_RefCN_CN + D_MAT_RefCN_CN)
  return(Jaccard_MAT)
}  

###################
which_max_for_matrix <- function(MAT){
  location_vec <- c(1,1)
  nrow_MAT <- nrow(MAT)
  ncol_MAT <- ncol(MAT)
  location_MAT <- which(MAT == max(MAT))[1]
  j <- 1
  while(location_MAT - nrow_MAT*j > 0){
    j<- j+1
  }
  i <- location_MAT - (j-1)*nrow_MAT
  location_vec <- c(i,j)
  return(location_vec)
}
#######################
cluster_min_defferencial_value <- function(RefCN, CN){
  MAT_RefCN <- matrix(data=rep(x=0, times=max(RefCN)*length(RefCN)), nrow = max(RefCN), ncol = length(RefCN))
  MAT_RefCN <- disassemble_vec(RefCN)
  MAT_CN <- matrix(data=rep(x=0, times=max(RefCN)*length(RefCN)), nrow = max(RefCN), ncol = length(RefCN))
  MAT_CN <- disassemble_vec(CN)
  MAT_SUM <- matrix(data=rep(x=0, times=max(RefCN)*length(RefCN)), nrow = max(RefCN), ncol = length(RefCN))
  MAT_SUM <- (MAT_RefCN - MAT_CN)^2
  SUM_DIFF_VAL <- 1/2 * sum(MAT_SUM)
  return(SUM_DIFF_VAL)
}

######################### MAIN ##############################
Jaccard_identification <- function(RefCN, CN){
  J_RefCM <- disassemble_vec(RefCN)
  J_CM    <- disassemble_vec(CN)
  J_Dij <- diff_MAT(J_RefCM, J_CM)
  J_Eij <- match_MAT(J_RefCM, J_CM)
  Jaccard_Dij <- Jaccard_MAT(J_Eij, J_Dij)
  J_MAT <- Jaccard_Dij
  CMm <- J_CM
  k <- ncol(Jaccard_Dij) -1
  for(i in 1:k){
    max_Jaccard_s <- which_max_for_matrix(J_MAT)[1]
    max_Jaccard_t <- which_max_for_matrix(J_MAT)[2]
    CMm_temp_s <- CMm[max_Jaccard_s,]
    CMm_temp_t <- CMm[max_Jaccard_t,]
    CMm[max_Jaccard_s,] <- CMm_temp_t
    CMm[max_Jaccard_t,] <- CMm_temp_s
    J_temp_s <- J_MAT[,max_Jaccard_s]
    J_temp_t <- J_MAT[,max_Jaccard_t]
    J_MAT[,max_Jaccard_s] <- J_temp_t
    J_MAT[,max_Jaccard_t] <- J_temp_s
    J_MAT[max_Jaccard_s,] <- -1
    J_MAT[,max_Jaccard_s] <- -1
  }
  J_CNmin <- assemble_vec(CMm)
  names(J_CNmin) <- names(RefCN)
  return(J_CNmin)
}  
######################### Jaccard_identification() function END  ######################
