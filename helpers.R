library(umap)
library(ggplot2)
library(ggrepel)
Refdata <- readRDS("data/MCL_cluster_reference.rds")

read_input <- function(inFile){
  input_data <- read.table(inFile$datapath,header=T,sep="\t",quote="",stringsAsFactors = F,row.names = 1)
  input_merge <- cbind(Refdata$ref, input_data[rownames(Refdata$ref),] )
  
  return(input_merge)
}

project_sample <- function(dataInput){
  S_new <- t(dataInput) %*% Refdata$W
  umap_new <- data.frame(umap(S_new)$layout,
                         Refdata$cluster[colnames(dataInput),"final_cluster"])
  
  colnames(umap_new) <- c("UMAP1", "UMAP2", "Cluster")
  umap_new[is.na(umap_new$Cluster),"Cluster"] <- "Input_sample"
  umap_new$sample <- rownames(umap_new)
  umap_new[umap_new$Cluster != "Input_sample","sample"] <- NA
  return(umap_new)
  
  
}

plot_umap <- function(umap_new,RefDot, InDot ){
  
  ggplot(umap_new) + 
    geom_point(aes(x=UMAP1, y=UMAP2,colour= Cluster),alpha=0.5, size=RefDot) +
    geom_point(data=subset(umap_new,Cluster == "Input_sample"),
               aes(x=UMAP1, y=UMAP2,colour=Cluster),alpha=0.8, size=InDot ,shape = 18) +
    geom_text_repel(aes(x=UMAP1, y=UMAP2,label=sample),min.segment.length = 0,max.overlaps = 100)+
    ggtitle("UMAP Plot") + 
    theme_bw() +
    theme(aspect.ratio=1, plot.title = element_text(face = "bold", size=12, hjust = 0.5), legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size=5)))
  
  
  
}



