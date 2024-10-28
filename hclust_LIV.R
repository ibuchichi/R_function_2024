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
