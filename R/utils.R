library(hdf5r)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(patchwork)
library(stringr)
library(rlang)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggridges)

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


cell2loc_violin <- function(obj, sample, predictions, csv, celltype){
  x <- parse_expr(celltype)
  cell <- read.csv(predictions)
  cell <- cell[cell$sample %like% sample,]
  tissue <- data.frame(read.csv(csv, header=F))
  tissue <- tissue[tissue$V2 == 1,]
  cells_ordered <- cell[order(cell$barcode),]
  pos_ordered <- tissue[order(tissue$V1),]
  cell_ident <- as.data.frame(obj@active.ident)
  colnames(cell_ident)[1] <- "x"
  cell_ident$names <- rownames(cell_ident)
  cell_ident_ordered <- cell_ident[order(cell_ident$x),]
  cells_ordered$ident <- cell_ident_ordered$x
  
  #class(celltype)
  #celltype <- parse(text=celltype)
  p<-ggplot(cells_ordered, aes(x=ident, y=eval(parse(text = celltype)), fill=ident)) +
    geom_violin() +ylab("Proportion")
  p
}


bulk_plot <- function(cell, region, gene_list){
  
  #Read data
  ctrl_string <- paste(cell, region, "CTRL.csv", sep = "_")
  prim_string <- paste(cell, region, "Primary.csv", sep = "_")
  sec_string <- paste(cell, region, "Secondary.csv", sep = "_")
  
  ctrl <- read.csv(ctrl_string, row.names = 1)
  prim <- read.csv(prim_string, row.names = 1)
  sec <- read.csv(sec_string, row.names = 1)
  
  ctrl_genes <- filter(ctrl, rownames(ctrl) %in% gene_list)
  prim_genes <- filter(prim, rownames(prim) %in% gene_list)
  sec_genes <- filter(sec, rownames(sec) %in% gene_list)
  
  table <- as.data.frame(log(t(ctrl_genes)))
  table2 <- as.data.frame(log(t(prim_genes)))
  table3 <- as.data.frame(log(t(sec_genes)))
  table3 <- cbind(table3, Stat=rep("Secondary_dmyel",length(rownames(table3))))
  table2 <- cbind(table2, Stat=rep("Primary_dmyel",length(rownames(table2))))
  table <- cbind(table, Stat=rep("No_dmyel",length(rownames(table))))
  full <- rbind(table3, table2)
  full <- rbind(full, table)
  
  reshaped <- melt(full,id.vars = "Stat", measure.vars = gene_list)
  
  ggplot(reshaped, aes(x=variable, y=value, fill=Stat)) +
    geom_boxplot() +
    theme_bw() +
    labs(y="", x = "", fill="") +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=10)) +
    # theme(axis.text.x = element_blank()) + # no x-axis labels
    # facet_wrap(~Demyel, scales = ‘free’, ncol = 1) +
    scale_fill_manual(values = c("forestgreen","magenta","lightskyblue")) +ylab("log10")
  
}
