setwd("./")
library(Seurat)
library(harmony)
library(liger)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(SeuratObject)
library(SeuratWrappers)
library(tidyverse)
library(DoubletFinder)

# Read Seurat objects for different time points
H0<-readRDS("rds/0h.rds")
H12<-readRDS("rds/12h.rds")
H24<-readRDS("rds/24h.rds")
H48<-readRDS("rds/48h.rds")
H72<-readRDS("rds/72h.rds")
H96<-readRDS("rds/96h.rds")
H120<-readRDS("rds/120h.rds")
H144<-readRDS("rds/144h.rds")
H168<-readRDS("rds/168h.rds")

# Generate violin plots for quality control metrics (nFeature_RNA and nCount_RNA)
p1 = VlnPlot(H0, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p2 = VlnPlot(H12, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p3 = VlnPlot(H24, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p4 = VlnPlot(H48, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p5 = VlnPlot(H72, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p6 = VlnPlot(H96, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p7 = VlnPlot(H120, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p8 = VlnPlot(H144, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 
p9 = VlnPlot(H168, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1) 

samples <- list(
  H0 = "H0@assays$RNA@counts",
  H12 = 'H12@assays$RNA@counts',
  H24 = 'H24@assays$RNA@counts',
  H48 = 'H48@assays$RNA@counts',
  H72 = 'H72@assays$RNA@counts',
  H96 = 'H96@assays$RNA@counts',
  H120 = 'H120@assays$RNA@counts',
  H144 = 'H144@assays$RNA@counts',
  H168 = 'H168@assays$RNA@counts'
)

samples <- list(
  H0 = H0,
  H12 = H12,
  H24 = H24,
  H48 = H48,
  H72 = H72,
  H96 = H96,
  H120 = H120,
  H144 = H144,
  H168 = H168
)

# Batch processing for each sample: filtering and basic preprocessing
for (name in names(samples)) {
  seurat_obj <- samples[[name]]
  
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2000)
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 300 & nCount_RNA < 10000)
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), dims = 1:30)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = c(seq(0, 1.5, by = 0.1)))
  
  file_name <- paste0(name,  "filtered.rds")
  saveRDS(seurat_obj, file = file_name)
}

rm(list = ls())

### DoubletFinder for doublet removal ###
library(DoubletFinder)

H0<-readRDS("H0filtered.rds")
H12<-readRDS("H12filtered.rds")
H24<-readRDS("H24filtered.rds")
H48<-readRDS("H48filtered.rds")
H72<-readRDS("H72filtered.rds")
H96<-readRDS("H96filtered.rds")
H120<-readRDS("H120filtered.rds")
H144<-readRDS("H144filtered.rds")
H168<-readRDS("H168filtered.rds")

samples <- list(
  H0 = H0,
  H12 = H12,
  H24 = H24,
  H48 = H48,
  H72 = H72,
  H96 = H96,
  H120 = H120,
  H144 = H144,
  H168 = H168
)

library(Seurat)
library(DoubletFinder)

for (name in names(samples)) {
  seurat_obj <- samples[[name]]
  
  sweep.res.list <- paramSweep(seurat_obj,PCs = 1:10,sct = F)
  sweep.stats <-summarizeSweep(sweep.res.list,GT = F)
  bcmvn<-find.pK(sweep.stats)
  pk_bcmvn<-bcmvn$pK[which.max(bcmvn$BCmetric)]%>%as.character()%in%as.numeric()
  
  DoubletRate <- 0.15
  je.prop<-modelHomotypic(seurat_obj$seurat_clusters)
  nExp_poi <-round(DoubletRate*ncol(seurat_obj))
  nExp_poi.adj <-round(nExp_poi*(1-je.prop))
  
  seurat_obj <- doubletFinder(seurat_obj, 
                              PCs = 1:10,
                              pN = 0.25, 
                              pK = pk_bcmvn, 
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE,sct = F)
  
  saveRDS(seurat_obj, file = paste0( name, "_dbfilter.rds"))
}

rm(list = ls())

H0<-readRDS("H0_dbfilter.rds")
H0<-H0[,H0$DF.classifications_0.25_FALSE_1830%in%c("Singlet")]

H12<-readRDS("H12_dbfilter.rds")
H12<-H12[,H12$DF.classifications_0.25_FALSE_3520%in%c("Singlet")]

H24<-readRDS("H24_dbfilter.rds")
H24<-H24[,H24$DF.classifications_0.25_FALSE_4979%in%c("Singlet")]

H48<-readRDS("H48_dbfilter.rds")
H48<-H48[,H48$DF.classifications_0.25_FALSE_3000%in%c("Singlet")]

H72<-readRDS("H72_dbfilter.rds")
H72<-H72[,H72$DF.classifications_0.25_FALSE_3396%in%c("Singlet")]

H96<-readRDS("H96_dbfilter.rds")
H96<-H96[,H96$DF.classifications_0.25_FALSE_2558%in%c("Singlet")]

H120<-readRDS("H120_dbfilter.rds")
H120<-H120[,H120$DF.classifications_0.25_FALSE_4221%in%c("Singlet")]

H144<-readRDS("H144_dbfilter.rds")
H144<-H144[,H144$DF.classifications_0.25_FALSE_6202%in%c("Singlet")]

H168<-readRDS("H168_dbfilter.rds")
H168<-H168[,H168$DF.classifications_0.25_FALSE_4572%in%c("Singlet")]

H0_matrix <- H0@assays$RNA@counts
H12_matrix <- H12@assays$RNA@counts
H24_matrix <- H24@assays$RNA@counts
H48_matrix <- H48@assays$RNA@counts
H72_matrix <- H72@assays$RNA@counts
H96_matrix <- H96@assays$RNA@counts
H120_matrix <- H120@assays$RNA@counts
H144_matrix <- H144@assays$RNA@counts
H168_matrix <- H168@assays$RNA@counts

data<-list(H0_matrix,H12_matrix,H24_matrix,H48_matrix,
           H72_matrix,H96_matrix,H120_matrix,H144_matrix,
           H168_matrix)
name<-c("H0","H12","H24","H48","H72","H96","H120","H144","H168")
names(data)<-name

scRNA.list<-list()
for (i in name) {
  scRNA.list[[i]] <- data[[i]]
  scRNA.list[[i]]<-CreateSeuratObject(scRNA.list[[i]],project = i)}
scRNA <- Reduce(merge,scRNA.list)

scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by = "orig.ident", do.center = FALSE)
scRNA <- RunPCA(scRNA)
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:30,verbose = T)
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:30,verbose = T)
scRNA <- FindClusters(scRNA, resolution = c(seq(0,2, by = 0.1)), verbose = T)

desired_order <- c("H0","H12","H24","H48","H72","H96","H120","H144","H168") 
scRNA$orig.ident <- factor(scRNA$orig.ident, levels = desired_order)

saveRDS(scRNA,file = "jellyfish_regeneration.rds")

nice_colors<-c('#6BAED6',
               '#FDAE61',
               '#FB9A99',
               '#A1D99B',
               '#FFDD57',
               '#9ECAE1',
               '#B3CDE3',
               '#9AD9C7',
               '#87CEEB',
               '#C5D9E8',
               '#A7C7E7',
               '#BDE4A7',
               '#FDD0A2',
               '#C5D9E8',
               '#7BCCC4',
               '#80B1D3')

VlnPlot(scRNA,features = "nFeature_RNA",cols = nice_colors,pt.size = 0)+ggtitle("")

DimPlot(scRNA, cols = nice_colors, group.by = "RNA_snn_res.0.1", 
        label = TRUE, label.size = 6, pt.size = 1.5, raster.dpi = c(1024, 1024)) + 
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"),
    plot.background = element_rect(
      color = "black",
      size = 1
    )
  ) + 
  labs(title = NULL)

### monocle3 pseudotime analysis ###
library(monocle3)
library(sf)
run_monocle3 <- function(sce){
  data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
  cell_metadata <- sce@meta.data
  data <- data[,row.names(cell_metadata)]
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  
  cds <- monocle3::new_cell_data_set(data,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)
  
  cds <- monocle3::preprocess_cds(cds, num_dim = 50)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "PCA")
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  
  cds.embed <- cds@int_colData$reducedDims$PCA
  int.embed <- Embeddings(sce, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  
  cds <- monocle3::cluster_cells(cds,reduction_method = 'UMAP') 
  cds <- monocle3::learn_graph(cds)
}

cds<-run_monocle3(sce)
cds = order_cells(cds)

library(RColorBrewer)
blues <- brewer.pal(9, "PiYG")
blue_gradient <- colorRampPalette(blues)
color_vector <- blue_gradient(n = 256)
light_purple_to_grey_palette <-c("#EF6548","#FC8D59", "#FDBB84", "#FDD49E", "#A6BDDB", "#74A9CF", "#3690C0")
color_vector <- rev(light_purple_to_grey_palette)
color_vector <- colorRampPalette(c("#023858", "#3690C0", "#A6BDDB", "#E0F3F8", "#F7FCF0", "#D9F0A3", "#78C679", "#238443"))(256)
color_vector <- colorRampPalette(c("#023858", "#3690C0", "#A6BDDB",  "#74A9CF", "#FC8D59", "#FDBB84", "#EF6548"))(256)
color_vector <- colorRampPalette(c("#F7FBFF", "#D2E3F3", "#A6BDDB", "#6C8EBF", "#4A6EA8", "#3B4CC0", "#2E2F92"))(256)

p1= monocle3::plot_cells(cds, 
                         color_cells_by="SCT_snn_res.0.1",
                         trajectory_graph_segment_size = 2,
                         label_cell_groups = F,
                         label_leaves = F,
                         label_roots = F,
                         label_branch_points = F,
                         cell_size = 1,
                         min_expr = 0,
                         alpha = 1)
p1

p2 = DimPlot(cells, label = T,label.size = 5)
p1+p2