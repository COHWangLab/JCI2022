.libPaths(c("/home/jmeiling/R/x86_64-pc-linux-gnu-library/3.5","/opt/R/3.5.1/lib64/R/library"))


################################################
### loading libraries
################################################
#install.packages('doMC')
#install.packages('Cairo')
#install.packages('factoextra')
#install.packages('survminer')
#BiocManager::install("impute")
#BiocManager::install("ggpubr")
library("doMC")
library("Cairo")
library("gtools")
library("gdata")
library("caTools")
library("gplots")
library("factoextra")
library("ggplot2")
library("stringr")
library("survival")
library("survminer")
library("reshape2")
library("ComplexHeatmap")
library("circlize")
library("impute")
library("ggpubr")
################################################
### generate input
################################################
work_dir <- "."
prefix <- "MCL"
setwd(work_dir)


################################################
### run cluster
################################################

## before run , please delete all temporary files generated before.

input_data <- paste0(prefix,"_input_matrix28.txt") ## input table name, format is each column is one sample each row is one alteration, colname of alteration should be "gene"


source('src/GDAC_TopgenesforCluster/Topgenes_v1.R')
result <- main("-s./src/GDAC_TopgenesforCluster/",paste0("-minput_data/",input_data),"-uALL",paste0("-o",prefix)) ## read data delete genes more than 2 na, impute missing values
source("src/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
result <- main("-ssrc/GDAC_NmfConsensusClustering/",paste0("-m",prefix,".expclu.gct"),paste0("-o",prefix),'-u2','-v10')
source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
result <- main(paste0("-u",prefix,".expclu.gct"),"-mPearson",paste0("-c",prefix,".cophenetic.coefficient.txt"),paste0("-w",prefix,".membership.txt"),paste0("-voutput_dir/",prefix),paste0("-pinput_data/",input_data))

## in the membership file choose the model with 4 clusters "membership.2"
## 1 is cluster 2, 2 is cluster4 , 3 is cluster 3, 4 is cluster 1
## then use the feature output to plot out the heatmap. 

if(FALSE){

data <- read.table(paste0("./input_data/",input_data),sep="\t",quote = "",stringsAsFactors = F,header = T,row.names = 1)## should be matrix 
info <- read.table("./survival.txt",sep="\t",quote = "",stringsAsFactors = F,header = T,row.names = 1)

Status_col <- c("Treatment_naive"="#006400","Relapse"="#FF0000","Normal"="blue","NA"="white")
IGHV_col <- c("Unmutated"="#977C00","Mutated"="#64A600","NA"="white")
MIPI_col <- c("High_risk"="#F00078","Intermediate_risk"="#FF79BC", "Low_risk"="#FFD9EC","NA"="white")
clinal_col <- c("Indolent"="#9AFF02","Classic"="#E1E100","Normal"="blue","NA"="white") ## heatmap annotation color 
## if have more annotations please set the annotation table and color here
anno <- info[colnames(data),]
anno <- anno[,c("Clinical_group","MIPI","Status","IGHV_status")]


survival <- function(x,surv,stat){
  
  tmp <- na.omit(info[,c(x,surv,stat)])
  colnames(tmp)[2:3] <- c("Time_to_First_Tx","Treatment")
  surv_object <- Surv(time = tmp$Time_to_First_Tx, event = tmp$Treatment)
  fit_not_ajusted <- do.call(survfit,list(surv_object~tmp[,x], data=tmp))
  p_tmp <- summary(fit_not_ajusted)
  
  p <- survdiff(surv_object ~ tmp[,x] , data = tmp)
  p <- 1 - pchisq(p$chisq, length(p$n) - 1)
  #return(p)
  
  
  ##By default survdiff is used to calculate regular log-rank test (with weights == 1)
  
  a <- ggsurvplot(fit_not_ajusted,
                  conf.int=FALSE, 
                  pval=TRUE, 
                  risk.table=TRUE, 
                  legend.title = x,  
                  palette="Set1", 
                  title=paste0(x,"_",surv), 
                  risk.table.height=.25,
                  surv.plot.height=.75
  )
  pdf(paste0(x,"_",surv,"_simple_KM.pdf"))
  print(a)
  dev.off()
}


cluster <- read.table(paste0(prefix,".membership.txt"),header = T,sep="\t",stringsAsFactors = F)
info <- cbind(info,cluster[rownames(info),1:4]) 
colnames(info) <- c(colnames(info)[1:8],"cluster_2","cluster_3","Cluster_4","Cluster_5")
com_factor <- c("Status","Clinical_group","IGHV_status","MIPI" ,"cluster_2","cluster_3","Cluster_4","Cluster_5")
pfs_cluster <- sapply(com_factor,function(x){survival(x,surv = "PFS_month",stat = "progress")})
os_cluster <- sapply(com_factor,function(x){survival(x,surv = "OS_month",stat = "Death")})
data_pfs <- data.frame(cluster=names(pfs_cluster),surv_q <- pfs_cluster)
data_os <- data.frame(cluster=names(os_cluster),surv_q <- os_cluster)



### heatmap for each cluster and marker genes

marker <- read.table(paste0("./output_dir/",prefix,".selectmarker.txt"),header = T,sep="\t",stringsAsFactors = F)
marker <- marker[order(marker$q),]
marker <- marker[order(marker$difference,decreasing = T),]
marker <- marker[!duplicated(marker$Hybridization.REF),]
rownames(marker) <- marker$Hybridization.REF
rownames(data) <- toupper(rownames(data))

row_ma <- data.frame(row.names = rownames(data),
                     Cluster = marker[rownames(data),"subclass"],
                     q_value = -log10(marker[rownames(data),"q"]),stringsAsFactors = F)
row_ma <- row_ma[order(row_ma$Cluster,-row_ma$q_value),]





rowanno_col <- list(cluster= c('1' = "red",'2' = "green", '3' = 'blue'))

col_ma <- info[,c("Clinical_group", "IGHV_status", "MIPI", "Status", "Cluster_4")]
col_ma <- col_ma[order(col_ma$Cluster_4),]



tmp <- as.matrix(data[rownames(row_ma),rownames(col_ma)])



cluster_col <- list(cluster = c("1" = "red","2" = "green","3" = "blue","4"="black"),
                    Status = Status_col,
                    IGHV_status = IGHV_col,
                    MIPI_Risk = MIPI_col,
                    Clinical_group = clinal_col)

topanno <- HeatmapAnnotation(Cluster = as.character(col_ma$Cluster_4),
                             Clinical_group = col_ma$Clinical_group,
                             IGHV_status = col_ma$IGHV_status,
                             MIPI_Risk =col_ma$MIPI,
                             Status = col_ma$status,
                             col= cluster_col,
                             na_col = "white")

ht_list = Heatmap(tmp, name = "Genetic Alterations", row_dend_side = "left",na_col = "white",show_column_names = F,show_row_names = T,
                  cluster_rows = F,cluster_columns = F, 
                  col = c("0"="white","1"="black","2"="blue"),
                  top_annotation = topanno, row_title = NULL, show_row_dend = FALSE) 
row =    HeatmapAnnotation(q_value= anno_barplot(row_ma$q_value),width = unit(4, "cm"),which = "row")




pdf(paste0(prefix,"_cluster_heatmap.pdf"),height = 20,width = 15)
draw(row+ht_list, ht_gap = unit(1, "mm"))
dev.off()
}



