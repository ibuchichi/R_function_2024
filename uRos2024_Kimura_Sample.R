####################################################################################################
# GitHub ibuchichi/R_function_2024 for uRos2024 Greece    sample script   2024/10/26
#    <uRos2024_Kimura.Rproj>
#                                           
####################################################################################################

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
######################### Jaccard_identification() function MAIN END  ######################

####################################################################################################
#    hclust_LIV()    function                       2024/10/25  Atsushi Kimura    
#
#      usage:
#          input:        std_data
#          output:       LIV_matrix
#
#      example:
#          # Draw LIV pareto chart
#          library(qcc)
#
#          LIV_matrix <- hclust_LIV(std_data)
#
#          for(i in 1: nrow(LIV_matrix)){
#            Pareto_vector <- LIV_matrix[i,]
#            names(Pareto_vector) <- colnames(LIV_matrix)
#            png_Pareto_file_name <- paste("LIV_",i ,"_pareto_qcc.png", sep="")
#            png(png_Pareto_file_name, width=3800, height=770)
#            ylab_words <- paste("LIV about ",i ,"th agglomaration")
#            pareto.chart(Pareto_vector, ylab = ylab_words)
#            dev.off()
#          }
#
#####################################################################################################
Elements <- function(Hclust.result, Data_std){
  Elements <- matrix(0, nrow(Data_std)-1, nrow(Data_std))
  for(i in 1: (nrow(Data_std)-1)){   
    if( Hclust.result$merge[i,1] < 0){ 
      Elements[i,(Hclust.result$merge[i,1])*(-1)] <- 1
    }else{                             
      Elements[i,] <- Elements[i,] | Elements[Hclust.result$merge[i,1],]
    }
    if( Hclust.result$merge[i,2] < 0){ 
      Elements[i,(Hclust.result$merge[i,2])*(-1)] <- 1
    }else{                             
      Elements[i,] <- Elements[i,] | Elements[Hclust.result$merge[i,2],]
    }
  }
  colnames(Elements) <- rownames(Data_std)
  return(Elements)
}

#################
hcluster_LIV <- function(Hclust.result, Elements_matrix, Data_std){
  LIV_matrix <- matrix(0, nrow(Data_std)-1, ncol(Data_std))
  for(i in 1: (nrow(Data_std)-1)){   
    Cp <- numeric(ncol(Data_std))   
    if( Hclust.result$merge[i,1] < 0){ 
      Cp <- Data_std[(-1)*Hclust.result$merge[i,1],]
    }else{                             
      for(j in 1: (nrow(Data_std))){
        if(Elements_matrix[Hclust.result$merge[i,1],j] == 1){
          Cp <- Data_std[j,] + Cp
        }      
      }
      Cp <- Cp / sum(Elements_matrix[Hclust.result$merge[i,1],]) 
    }
    Cq <- numeric(ncol(Data_std))   
    if( Hclust.result$merge[i,2] < 0){ 
      Cq <- Data_std[(-1)*Hclust.result$merge[i,2],]
    }else{                             
      for(j in 1: (nrow(Data_std))){
        if(Elements_matrix[Hclust.result$merge[i,2],j] == 1){
          Cq <- Data_std[j,] + Cq
        }      
      }
      Cq <- Cq / sum(Elements_matrix[Hclust.result$merge[i,2],]) 
    }
    vect_dist <- sum((Cp -Cq)^2)
    for(k in 1:ncol(Data_std)){
      LIV_matrix[i,k] <- LIV_matrix[i,k] + (Cp[k]-Cq[k])^2 / vect_dist  
    }
  }
  colnames(LIV_matrix) <- colnames(Data_std)
  LIV_matrix <- LIV_matrix * 100 
  return(LIV_matrix)
}

############################   hclust_LIV() function    MAIN   ##########################
hclust_LIV <- function(Data_std){
  # hclust()
  Data_std.dist <- sqrt(1/2)*dist(Data_std)
  methname <- "ward.D2"
  hclust_result <- hclust(Data_std.dist, method=methname)
  
  # cluster elements
  Elements_matrix <- Elements(hclust_result, Data_std)
  write.csv(Elements_matrix, file="Elements_result", quote = FALSE)
  
  # LIV calc.
  LIV_matrix <- hcluster_LIV(hclust_result, Elements_matrix, Data_std)
  colnames(LIV_matrix) <- colnames(Data_std)    
  write.csv(LIV_matrix, file="LIV_matrix.csv", quote = FALSE)
  return(LIV_matrix)
}

#######################   hclust_LIV()  function  MAIN END      #############





library(dplyr)
###########################################################################################################
# Loading FIES-from2017to2019.csv for Reference Data
#
#ssdse_c_Ref <- read.csv("FIES-from2017to2019.csv", fileEncoding="Shift-JIS")
ssdse_c_Ref <- read.csv("FIES-from2017to2019.csv")
rownames(ssdse_c_Ref) <- ssdse_c_Ref$Prefecture
ssdse_small_Ref <- ssdse_c_Ref |> select(-Prefecture) 

#############################################################
# Rationing and Standardisation of Reference Data
#
ssdse_small_Ref.p <- (ssdse_small_Ref / ssdse_small_Ref$Food )* 100
ssdse_small_Ref.std <- ssdse_small_Ref.p |> select(-Food) |> scale()

###############################
# Hierarchical Cluster Analysis (Ward Method)
#
ssdse_small_Ref.dist <- sqrt(1/2)*dist(ssdse_small_Ref.std)
methname <- "ward.D2"
ssdse_small_Ref.hclust <- hclust(ssdse_small_Ref.dist, method=methname)
ssdse12_Ref.cnum <- cutree(ssdse_small_Ref.hclust, k=12)

###########################
# Loading FIES(Family Income and Expenditure Survey Japan) Data Set to be Analysed 
#
#ssdse_c_year <- "FIES-from2007to2009"
#ssdse_c_year <- "FIES-from2010to2012"
#ssdse_c_year <- "FIES-from2011to2013"
#ssdse_c_year <- "FIES-from2012to2014"
#ssdse_c_year <- "FIES-from2015to2017"
#ssdse_c_year <- "FIES-from2016to2018"
#ssdse_c_year <- "FIES-from2017to2019"
ssdse_c_year <- "FIES-from2020to2022"
#ssdse_c_year <- "FIES-from2021to2023"

ssdse_data_path <- paste(ssdse_c_year, ".csv", sep="")
#ssdse_c <- read.csv(ssdse_data_path, fileEncoding="Shift-JIS")
ssdse_c <- read.csv(ssdse_data_path)

rownames(ssdse_c) <- ssdse_c$Prefecture
ssdse_small <- ssdse_c |> select(-Prefecture) 


#############################################################
# Rationing and Standardisation of Analysis Data
#
ssdse_small.p <- (ssdse_small / ssdse_small$Food )* 100
std_data <- ssdse_small.p |> select(-Food) |> scale()
ssdse_c_pareto_small_name <- colnames(std_data)

STD_file_name <- paste(ssdse_c_year , "_STD.csv",sep="")
write.csv(std_data, file=STD_file_name, quote = FALSE)


###############################################
# Hierarchical Cluster Analysis (Ward Method)
#
std_data.dist <- sqrt(1/2)*dist(std_data)
methname <- "ward.D2"
hclust.result <- hclust(std_data.dist, method=methname)


############### Reorder Process  hclust.result → hclust.result_reorder
#   from  <OUJ_estrela202103forOUJVideo.Rproj>
hclust.result_reorder <- hclust.result

dd <- as.dendrogram(hclust.result)
dd.reorder <- reorder(dd, 47:1, agglo.FUN = mean)
hc.dd.reorder <- as.hclust(dd.reorder)
hclust.result_reorder$order <- hc.dd.reorder$order


#########################
# Drawing the Dendrogram
#
Map_file_name <- paste(ssdse_c_year , "_dendrogram_reorder.png",sep="")
png(Map_file_name, width=920, height=770)
plot(hclust.result_reorder)
dev.off()

cluster_num <- 12

Map_file_name <- paste(ssdse_c_year , "_dendrogram_reorder_", cluster_num, ".png",sep="")
png(Map_file_name, width=920, height=770)
plot(hclust.result_reorder)
rect.hclust(hclust.result_reorder, k= cluster_num)
dev.off()


####################################################
# Caluculation of LIV values (LIV_matrix)
####################################################
LIV_matrix <- hclust_LIV(std_data)
colnames(LIV_matrix) <- colnames(std_data)
write_file_name_LIV <- paste("hcluster_LIV_", ssdse_c_year, ".csv", sep="") 
write.csv(LIV_matrix, file=write_file_name_LIV, quote = FALSE)

####################################################
# Creation and Display of Pareto Chart
####################################################
library(qcc)

for(i in 1: nrow(LIV_matrix)){
  Pareto_vector <- LIV_matrix[i,]
  names(Pareto_vector) <- colnames(LIV_matrix)
  png_Pareto_file_name <- paste(ssdse_c_year, "-LIV_",i ,"_pareto_qcc.png", sep="")
  png(png_Pareto_file_name, width=3800, height=770)
  ylab_words <- paste("LIV about ",i ,"th agglomaration")
  pareto.chart(Pareto_vector, ylab = ylab_words)
  dev.off()
}


############################################################
# Processing to create continuity in the cluster colours
#              of the agglomerating Japanese divisional map.
#
#  input:
#      ssdseN.cnum   : Cluster vector brfore agglomeration
#      ssdseN1.cnum  : Cluster vector  after agglomeration
#      ssdseNrenum.cnum  :  Renumbered Cluster vectore 
#                                     before agglomeration
#  output:
#      ssdseN1renum.cnum :  Renumbered Cluster vectore
#                                      after agglomeration
#
Cluster_Renumber_for_continuity <- function(ssdseN.cnum, ssdseN1.cnum, ssdseNrenum.cnum){
  #
  pref_name_back <- names(ssdseN.cnum[ssdseN.cnum != ssdseN1.cnum][1])
  back_cluster_num <- ssdseNrenum.cnum[pref_name_back]
  pref_name_forward <- names(ssdseN.cnum[ssdseN.cnum == ssdseN1.cnum[pref_name_back]][1])
  forward_cluster_num <- ssdseNrenum.cnum[pref_name_forward]
  
  #
  back_cluster_element_num <- length(ssdseN.cnum[ssdseN.cnum == ssdseN.cnum[pref_name_back]])
  forward_cluster_element_num <- length(ssdseN.cnum[ssdseN.cnum == ssdseN.cnum[pref_name_forward]])
  
  #
  if(forward_cluster_element_num >= back_cluster_element_num){
    ssdseNrenum.cnum[ssdseNrenum.cnum == back_cluster_num] <- forward_cluster_num
  }else{
    ssdseNrenum.cnum[ssdseNrenum.cnum == forward_cluster_num] <- back_cluster_num
  }
  
  ssdseN1renum.cnum <- ssdseNrenum.cnum
  
  return(ssdseN1renum.cnum)
}
################################################################################


ssdse12.cnum <- cutree(hclust.result, k=12)
ssdse12renum.cnum <- Jaccard_identification(ssdse12_Ref.cnum, ssdse12.cnum)

ssdse11.cnum <- cutree(hclust.result, k=11)
ssdse11renum.cnum <- Cluster_Renumber_for_continuity(ssdse12.cnum, ssdse11.cnum, ssdse12renum.cnum)

ssdse10.cnum <- cutree(hclust.result, k=10)
ssdse10renum.cnum <- Cluster_Renumber_for_continuity(ssdse11.cnum, ssdse10.cnum, ssdse11renum.cnum)

ssdse9.cnum <- cutree(hclust.result, k=9)
ssdse9renum.cnum <- Cluster_Renumber_for_continuity(ssdse10.cnum, ssdse9.cnum, ssdse10renum.cnum)

ssdse8.cnum <- cutree(hclust.result, k=8)
ssdse8renum.cnum <- Cluster_Renumber_for_continuity(ssdse9.cnum, ssdse8.cnum, ssdse9renum.cnum)

ssdse7.cnum <- cutree(hclust.result, k=7)
ssdse7renum.cnum <- Cluster_Renumber_for_continuity(ssdse8.cnum, ssdse7.cnum, ssdse8renum.cnum)

ssdse6.cnum <- cutree(hclust.result, k=6)
ssdse6renum.cnum <- Cluster_Renumber_for_continuity(ssdse7.cnum, ssdse6.cnum, ssdse7renum.cnum)

ssdse5.cnum <- cutree(hclust.result, k=5)
ssdse5renum.cnum <- Cluster_Renumber_for_continuity(ssdse6.cnum, ssdse5.cnum, ssdse6renum.cnum)

ssdse4.cnum <- cutree(hclust.result, k=4)
ssdse4renum.cnum <- Cluster_Renumber_for_continuity(ssdse5.cnum, ssdse4.cnum, ssdse5renum.cnum)

ssdse3.cnum <- cutree(hclust.result, k=3)
ssdse3renum.cnum <- Cluster_Renumber_for_continuity(ssdse4.cnum, ssdse3.cnum, ssdse4renum.cnum)

ssdse2.cnum <- cutree(hclust.result, k=2)
ssdse2renum.cnum <- Cluster_Renumber_for_continuity(ssdse3.cnum, ssdse2.cnum, ssdse3renum.cnum)

ssdse1.cnum <- cutree(hclust.result, k=1)
ssdse1renum.cnum <- Cluster_Renumber_for_continuity(ssdse2.cnum, ssdse1.cnum, ssdse2renum.cnum)


##############################################################################
####　　　　 Drawing of Japan Divisional Maps 　　　　　######################
##############################################################################
library(sf)

samplek <- st_read("Japan_H26_05.shp")

prefS <- c("WAKAYAMA","HIROSHIMA","OKAYAMA","SHIMANE","KOCHI","OITA","KUMAMOTO","SAGA","FUKUOKA","MIE",
           "KAGOSHIMA","NARA","YAMAGATA","AKITA","AOMORI","TOCHIGI","TOYAMA","FUKUSHIMA","SAITAMA","MIYAGI",
           "MIYAZAKI","IWATE","HYOGO","EHIME","NIGATA","TOKYO","GUNMA","ISHIKAWA","KANAGAWA","HOKKAIDO",
           "AICHI","SHIZUOKA","GIFU","OSAKA","CHIBA","YAMANASHI","FUKUI","IBARAKI","NAGASAKI","NAGANO",
           "KAGAWA","TOTTORI","KYOTO","SHIGA","OKINAWA","TOKUSHIMA","YAMAGUCHI")  

#########################
# Set Display Japan Map Area
xlim <- c(126, 146) # longitude
ylim <- c(26, 46)   # latitude

library(RColorBrewer)
##########################
#Draw Japan Map by ssdse12_Ref.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse12_Ref.cnum[prefS])])

Map_file_name <- paste("Reference FIES-from2017to2019" , "cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse12_Ref.cnum[prefS])])
dev.off()


##########################
#Draw Japan Map by ssdse12renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse12renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "12renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse12renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse11renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse11renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "11renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse11renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse10renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse10renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "10renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse10renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse9renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse9renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "9renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse9renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse8renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse8renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "8renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse8renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse7renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse7renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "7renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse7renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse6renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse6renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "6renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse6renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse5renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse5renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "5renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse5renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse4renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse4renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "4renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse4renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse3renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse3renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "3renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse3renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse2renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse2renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "2renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse2renum.cnum[prefS])])
dev.off()

##########################
#Draw Japan Map by ssdse1renum.cnum
par(mfrow=c(1L,1L))
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse1renum.cnum[prefS])])

Map_file_name <- paste(ssdse_c_year, "1renum_cluster_Map.png",sep="-")
png(Map_file_name, width=920, height=770)
plot(samplek, xlim=xlim, ylim=ylim, col=brewer.pal(12,"Set3")[as.integer(ssdse1renum.cnum[prefS])])
dev.off()

