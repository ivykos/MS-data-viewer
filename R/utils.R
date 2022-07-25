library(hdf5r)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(patchwork)
library(stringr)


# Function for overlaying cluster labels on Visium tissue slide
transfer_clusters <- function(obj, csv){
  positions <- data.frame(read.csv(csv, header = F))
  positions <- positions[positions$V2 == 1,]
  idents <- data.frame(obj@active.ident)
  idents <- rownames_to_column(idents)
  pos.ordered <- positions[order(positions$V1),]
  idents.ordered <- idents[order(idents$rowname),]
  identity <- as.factor(idents.ordered$obj.active.ident)
  
  plt2 <- ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, 
                                  color=identity)) + 
    geom_point() + theme_linedraw() + xlab("X") + ylab("Y") 
  
  plt2
}

# Function for determining gene expression spatially
get_expression <- function(obj, feature, csv){
  list_of_genes <- rownames(obj)
  if (feature %in% list_of_genes == FALSE){
    stop("Gene not present")
  }
  positions <- data.frame(read.csv(csv, header = F))
  positions <- positions[positions$V2 == 1,]
  idents <- data.frame(obj@active.ident)
  idents <- rownames_to_column(idents)
  pos.ordered <- positions[order(positions$V1),]
  idents.ordered <- idents[order(idents$rowname),]
  
  tmp <- as.matrix(GetAssayData(object = obj, slot = "counts"))
  feature_count <- tmp[feature,]
  rotated <- as.data.frame(t(feature_count))
  rotated <- as.data.frame(t(rotated))
  
  expression <- rotated$V1
  ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
    geom_point(aes(color = expression), size = 2) + 
    scale_color_viridis(option = "inferno") + theme_bw() + xlab("X") + ylab("Y") +
    ggtitle(as.character(feature))
  
}