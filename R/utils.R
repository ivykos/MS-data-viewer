library(hdf5r)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(patchwork)
library(stringr)

transfer_clusters_merged <- function(seuratObj, data_set, tissue_csv){
  #Again, make sure that the standard seurat guided clustering workflow has been followed.
  
  data.subset <-  subset(x=seuratObj, subset = orig.ident==data_set)
  
  #Retrieve the cluster labels and project ID from the Seurat object and put them in tables
  write.table(data.subset@active.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  idents <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  write.table(data.subset@meta.data$orig.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  orig <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  
  #Get the coordinates for plotting the tissue image
  positions <- read.csv(tissue_csv, header = F)
  
  #Make sure only the spots within the tissue are included for this analysis
  positions <- positions[positions$V2 == 1,]
  
  cells <- idents$V1
  cells <- str_remove(cells, "_1") #Remove _1 suffix from cell names
  idents$Cells <- cells
  
  pos.ordered <- positions[order(positions$V1),]
  
  
  #Plot it
  
  plt2 <- ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, color=as.factor(idents$V2))) + 
    geom_point() + theme_linedraw()
  
  plt2
  
}

transfer_clusters <- function(obj, csv){
  positions <- data.frame(read.csv(csv, header = F))
  positions <- positions[positions$V2 == 1,]
  idents <- data.frame(obj@active.ident)
  idents <- rownames_to_column(idents)
  pos.ordered <- positions[order(positions$V1),]
  idents.ordered <- idents[order(idents$rowname),]
  
  plt2 <- ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, 
                                  color=as.factor(idents.ordered$obj.active.ident))) + 
    geom_point() + theme_linedraw()
  
  plt2
}
get_expression <- function(seuratObj, proj, feature, tissue_csv){
  
  #Get tissue positions again 
  write.table(seuratObj@active.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  idents <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  write.table(seuratObj@meta.data$orig.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  orig <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  
  positions <- read.csv(tissue_csv, header = F)
  positions <- positions[positions$V2 == 1,]
  
  idents["Origin"] = orig$V2
  table <- idents[idents$Origin == proj,]
  cluster.ordered <- table[order(table$V1),]
  pos.ordered <-positions[order(positions$V1),]
  
  #Extract the counts matrix from the Seurat object
  tmp <- as.matrix(GetAssayData(object = seuratObj, slot = "counts"))
  feature_count <- tmp[feature,]
  rotated <- as.data.frame(t(feature_count))
  rotated <- as.data.frame(t(rotated))
  table["Expr"] <- rotated$V1
  
  ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
    geom_point(aes(color = table$Expr), size = 2) + 
    scale_color_viridis(option = "inferno") + theme_bw() +ggtitle(as.character(paste(proj, feature)))
}

