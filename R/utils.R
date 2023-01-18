library(hdf5r)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(patchwork)
library(stringr)
library(tidyverse)
library(rlang)
library(data.table)
#library(ggtips)


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

#Function for plotting predicted Cell2Location proportions
cell2loc <- function(obj, predictions, csv, celltype){
  x <- parse_expr(celltype)
  cell <- read.csv(predictions)
  cell <- cell[cell$sample %like% as.character(obj),]
  tissue <- data.frame(read.csv(csv, header=F))
  tissue <- tissue[tissue$V2 == 1,]
  cells_ordered <- cell[order(cell$barcode),]
  pos_ordered <- tissue[order(tissue$V1),]
  type <- as.numeric(unlist(select(cells_ordered, !! x)))
  
  plt <- ggplot(pos_ordered, aes(x=pos_ordered$V3, y=pos_ordered$V4)) +
    geom_point(aes(color=type), size=2) + scale_color_viridis(option="viridis") +
    theme_bw() + xlab("X") + ylab("Y") + labs(color = "Proportion") +ggtitle(celltype)
  
  plt
}


tissue_patch <- function(obj, tiss, obj2, tiss2){
  p1 <- transfer_clusters(obj, tiss)
  p2 <- transfer_clusters(obj2, tiss2)
  p1 / p2
}

umap_patch <- function(obj, obj2){
  p1 <- DimPlot(obj, reduction = "umap")
  p2 <- DimPlot(obj2, reduction = "umap")
  p1 / p2
}

gene_patch <- function(obj, obj2, tiss, tiss2, feat){
  p1 <- get_expression(obj, feat, tiss)
  p2 <- get_expression(obj2, feat, tiss2)
  p1 / p2
}

vln_patch <- function(obj, obj2, feat){
  p1 <- VlnPlot(obj, features = feat, group.by = "seurat_clusters")
  p2 <- VlnPlot(obj2,features = feat, group.by = "seurat_clusters")
  p1 / p2
}

c2l_patch <- function(sample, sample2, csv, csv2, celltype, pred){
  p1 <- cell2loc(sample, pred, csv, celltype)
  p2 <- cell2loc(sample2, pred, csv2, celltype)
  p1 / p2
}