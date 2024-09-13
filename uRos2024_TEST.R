####################################################################################################
# Test R script for uRos2024AthensGreece27-29Nov     SSDSE-C-2024e.csv
#    2024/9/13 　　　　　　　　　　                                           
####################################################################################################


###########################
# Read SSDSE-C-202Xe.csv
#
ssdse_c_year <- "SSDSE-C-2023e"

ssdse_data_path <- paste(ssdse_c_year, ".csv", sep="")

### Read processing
items <- read.csv(ssdse_data_path, header=FALSE, nrows=1, fileEncoding="Shift-JIS") 
ssdse_c <- read.csv(ssdse_data_path, header=FALSE, skip=2, col.names=items, fileEncoding="Shift-JIS")

survey_years <- colnames(ssdse_c)[2]
colnames(ssdse_c)[2] <- "Prefecture" 
rownames(ssdse_c) <- ssdse_c$Prefecture

##########################################################
# Data preprocessing 1
#
ssdse_small <- ssdse_c |> filter(Prefecture != "全国") |> select(matches(c("LB00" , "LB\\d{6}")))


ssdse_c_pareto <- read.csv(ssdse_data_path, header=FALSE, col.names=items, fileEncoding="Shift-JIS")
colnames(ssdse_c_pareto)[2] <- "Prefecture"
ssdse_c_pareto_small <- ssdse_c_pareto |> filter(Prefecture != "全国") |> select(matches("LB\\d{6}"))
ssdse_c_pareto_small_name <- ssdse_c_pareto_small[2,]    

　
#############################################################
# Data preprocessing 2
#
ssdse_small.p <- (ssdse_small / ssdse_small$LB00 )* 100
ssdse_small.std <- ssdse_small.p |> select(-LB00) |> scale()

ssdse_small_foodName.std <- ssdse_small.std
colnames(ssdse_small_foodName.std) <- ssdse_c_pareto_small_name

STD_file_name <- paste(ssdse_c_year , "_STD.csv",sep="")
write.csv(ssdse_small.std, file=STD_file_name, quote = FALSE)

STD_foodName_file_name <- paste(ssdse_c_year , "_STD_foodName.csv",sep="")
write.csv(ssdse_small_foodName.std, file=STD_foodName_file_name, quote = FALSE)

###############################
# ward法階層的クラスター分析を実施して日本地図アニメーション用のssdse12.cnumを作成する
#
ssdse_small.dist <- sqrt(1/2)*dist(ssdse_small.std)
methname <- "ward.D2"
ssdse_small.hclust <- hclust(ssdse_small.dist, method=methname)
ssdse12.cnum <- cutree(ssdse_small.hclust, k=12)
############################################################################################################


############### Reorder processing ssdse_small.hclust → ssdse_small.hclust_reorder
ssdse_small.hclust_reorder <- ssdse_small.hclust

dd <- as.dendrogram(ssdse_small.hclust)
dd.reorder <- reorder(dd, 47:1, agglo.FUN = mean)
hc.dd.reorder <- as.hclust(dd.reorder)
ssdse_small.hclust_reorder$order <- hc.dd.reorder$order
###############################################################################

####################
# Draw dendrogram
Map_file_name <- paste(ssdse_c_year , "_dendrogram_reorder.png",sep="")
png(Map_file_name, width=920, height=770)
plot(ssdse_small.hclust_reorder)
dev.off()


