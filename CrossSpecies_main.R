# ==============================================================================
# Single-Cell RNA-seq Analysis Pipeline for Aurelia
# Main Tasks:
# 1. Subset cell clusters
# 2. Cell stemness analysis (CytoTRACE2)
# 3. Lineage/trajectory inference (Slingshot, Monocle3)
# 4. Cross-species gene mapping & integration
# 5. Visualization of UMAP, trajectories, and gene expression
# ==============================================================================
# Load required libraries
library(Seurat)
library(tidyverse)
library(CytoTRACE2)
library(SeuratDisk)
library(slop)
library(monocle3)
library(writexl)
library(rtracklayer)
library(slingshot)
ZS_nosplit <- read_rds('~/Aurelia/ZS_nosplit.rds')
ZS_nosplit_sub <- subset(ZS_nosplit, RNA_snn_res.0.6 %in%
                           c(0,1,4,7,8,10,17,18,21,22,28))
set.seed(1234)
ZS_nosplit_down <- ZS_nosplit[, sample(1:ncol(ZS_nosplit), 10000)]
ZS_nosplit_down[['RNA']] <- as(ZS_nosplit_down[['RNA']], Class = 'Assay')
write_rds(ZS_nosplit_down, '~/Aurelia/ZS_nosplit_down.rds')
# ==============================================================================
#Lineage Inference with Slingshot
# ==============================================================================
pdf('~/Aurelia/lineage3.pdf', width = 8, height = 8)
# Run Slingshot trajectory analysis
ZS_nosplit_slingshot <- slop::RunSlingshot(
  srt = subset(ZS_nosplit_down, RNA_snn_res.0.6 %in%
                 c(0,1,4,7,8,10,17,18,21,22,28)),
  group.by = "RNA_snn_res.0.6",
  start = 1,
  reduction = "umap"
)
dev.off()

# Visualize lineages on UMAP
slop::CellDimPlot(
  ZS_nosplit_slingshot,
  group.by = "RNA_snn_res.0.6",
  reduction = "umap",
  lineages = paste0("Lineage", 1:3),
  lineages_span = 0.1,
  lineages_trim = c(0.05, 0.95)
)
# Plot pseudotime for each lineage
umap_limits <- list(
  x = range(ZS_nosplit@reductions$umap@cell.embeddings[,1]),
  y = range(ZS_nosplit@reductions$umap@cell.embeddings[,2])
)

p1 <- FeatureDimPlot(ZS_nosplit_slingshot, features = "Lineage1",
                     reduction = "umap", pt.size = 0.5) +
  xlim(umap_limits$x) + ylim(umap_limits$y)

p2 <- FeatureDimPlot(ZS_nosplit_slingshot, features = "Lineage2",
                     reduction = "umap", pt.size = 0.5) +
  xlim(umap_limits$x) + ylim(umap_limits$y)

p3 <- FeatureDimPlot(ZS_nosplit_slingshot, features = "Lineage3",
                     reduction = "umap", pt.size = 0.5) +
  xlim(umap_limits$x) + ylim(umap_limits$y)

# Combine and save lineage plots
pdf('~/Aurelia/reversion/lineage_all.pdf', width = 16, height = 5)
p1 + p2 + p3 +p4
dev.off()
# ==============================================================================
# Re-process Subsetted Data (Normalization, PCA, UMAP)
# ==============================================================================
ZS_nosplit_sub <- JoinLayers(ZS_nosplit_sub)
sceobj_exp <- LayerData(ZS_nosplit_sub, assay = "RNA", layer = 'counts')

# Create new Seurat object
sceobj_new <- CreateSeuratObject(counts = sceobj_exp,
                                 meta.data = ZS_nosplit_sub@meta.data)

#Split assay by sample for integration
sceobj_new[["RNA"]] <- split(sceobj_new[["RNA"]], f = sceobj_new$orig.ident)
# ==============================================================================
#: Export Cluster2 for Cross-Species Analysis
# ==============================================================================
C2 <- read_rds('~/Aurelia/reversion/cluster2.rds')
C2[['RNA']] <- as(C2[['RNA']], Class = 'Assay')
# Modify gene names to avoid conflicts
C2_exp <- C2@assays$RNA@counts
rownames(C2_exp) <- paste0(stringr::str_replace_all(rownames(C2_exp), 'gene', 'rna'), '.1')
# Create annotated Seurat object
C2 <- CreateSeuratObject(counts = C2_exp, meta.data = C2@meta.data)
C2 <- NormalizeData(C2)
C2 <- ScaleData(C2)
C2$celltype <- 'C2'

# Export as h5Seurat and h5ad (for Python/scanpy)
SeuratDisk::SaveH5Seurat(C2, filename = "~/Aurelia/reversion/C2_ann.h5Seurat")
SeuratDisk::Convert("~/Aurelia/reversion/C2_ann.h5Seurat", dest = "h5ad")
write_rds(C2, '~/Aurelia/reversion/C2.rds')
# ==============================================================================
# Cross-Species Gene Mapping Functions
# ==============================================================================
trans_species_obj <- function(query_obj, subject_obj, map_gene) {
  subject_exp <- subject_obj@assays$RNA$counts
  subject_exp_duplicate_new <- subject_exp
  for (i in 1:5) {
    temp_exp <- subject_exp
    rownames(temp_exp) <- paste0(rownames(subject_exp), "-", i)
    subject_exp_duplicate_new <- rbind(temp_exp, subject_exp_duplicate_new)
  }
  map_gene <- map_gene[map_gene$subject_id %in% rownames(subject_exp_duplicate_new), ]
  subject_exp_duplicate_new <- subject_exp_duplicate_new[map_gene$subject_id, ]
  rownames(subject_exp_duplicate_new) <- map_gene$query_id
  subject_obj_new <- CreateSeuratObject(counts = subject_exp_duplicate_new, meta.data = subject_obj@meta.data)
  merge_obj <- merge(query_obj, subject_obj_new)
  return(merge_obj)
}
trans_gene_species <- function(species1_vs_species2, species2_vs_species1) {
  colnames(species1_vs_species2) <- c(
    "query_id", "subject_id", "percent_identity", "alignment_length",
    "mismatch_count", "gap_open_count", "q_start", "q_end",
    "s_start", "s_end", "e_value", "bit_score"
  )
  colnames(species2_vs_species1) <- colnames(species1_vs_species2)
  
  best_hits_species1 <- species1_vs_species2 %>%
    group_by(query_id) %>%
    top_n(4, bit_score) %>%
    ungroup()
  
  best_hits_species2 <- species2_vs_species1 %>%
    group_by(query_id) %>%
    top_n(4, bit_score) %>%
    ungroup()
  
  rbh <- inner_join(best_hits_species1, best_hits_species2,
                    by = c("query_id" = "subject_id", "subject_id" = "query_id")) %>%
    select(query_id, subject_id)
  return(rbh)
}

# ==============================================================================
#  Cross-Species Integration
# ==============================================================================
# Xenopus
ROCs_sce <- SeuratDisk::LoadH5Seurat("~/Aurelia/regen/Xenopus/ROCs_sce.h5Seurat")
Au_to_Xp <- read_table('~/Aurelia/regen/maps/XeAu/Au_to_Xe.txt', col_names = F)
Xp_to_Au <- read_table('~/Aurelia/regen/maps/XeAu/Xe_to_Au.txt', col_names = F)
map_gene_Au_Xp <- trans_gene_species(Au_to_Xp, Xp_to_Au)
ROCs_sce$celltype <- 'Xenopus_ROCs'
ROCs_sce$species <- 'Xenopus'
# Hydra
Hy_iCell <- SeuratDisk::LoadH5Seurat("~/Aurelia/regen/Hydra/Hy_iCell.h5Seurat")
Au_to_Hy <- read_table('~/Aurelia/regen/maps/HyAu/Au_to_Hy.txt', col_names = F)
Hy_to_Au <- read_table('~/Aurelia/regen/maps/HyAu/Hy_to_Au.txt', col_names = F)
map_gene_Au_Hy <- trans_gene_species(Au_to_Hy, Hy_to_Au)
Hy_iCell$celltype <- 'Hydra_iCell'
Hy_iCell$species <- 'Hydra'
# Schmidtea mediterranea
Sc_neoblast <- SeuratDisk::LoadH5Seurat("~/Aurelia/regen/Schmidtea_mediterranea/Sc_neoblast.h5Seurat")
Au_to_Sc <- read_table('~/Aurelia/regen/maps/ScAu/Au_to_Sc.txt', col_names = F)
Sc_to_Au <- read_table('~/Aurelia/regen/maps/ScAu/Sc_to_Au.txt', col_names = F)
map_gene_Au_Sc <- trans_gene_species(Au_to_Sc, Sc_to_Au)
Sc_neoblast$celltype <- 'Sc_neoblast'
Sc_neoblast$species <- 'Schmidtea_mediterranea'
# Danio rerio
Danio_rerio_sce_Me <- SeuratDisk::LoadH5Seurat("~/Aurelia/regen/Danio_rerio/Danio_re_Mesenchymal.h5Seurat")
Au_to_Da <- read_table('~/Aurelia/regen/maps/DaAu/Au_to_Da.txt', col_names = F)
Da_to_Au <- read_table('~/Aurelia/regen/maps/DaAu/Da_to_Au.txt', col_names = F)
map_gene_Au_Da <- trans_gene_species(Au_to_Da, Da_to_Au)
Danio_rerio_sce_Me$celltype <- 'Danio_rerio_Mesenchymal'
Danio_rerio_sce_Me$species <- 'Danio_rerio'
# Acoel
Ac_neoblast <- SeuratDisk::LoadH5Seurat("~/Aurelia/regen/Acoel/Ac_neoblast1.h5Seurat")
Au_to_Ac <- read_table('~/Aurelia/regen/maps/AcAu/Au_to_Ac.txt', col_names = F)
Ac_to_Au <- read_table('~/Aurelia/regen/maps/AcAu/Ac_to_Au.txt', col_names = F)
map_gene_Au_Ac <- trans_gene_species(Au_to_Ac, Ac_to_Au)
Ac_neoblast$celltype <- 'Ac_neoblast'
Ac_neoblast$species <- 'Acoel'
# ==============================================================================
#  Merge All Species and Integration
# ==============================================================================
query_exp <- C2@assays$RNA$counts 
query_exp_duplicate <- query_exp
query_exp_duplicate_new <- query_exp
for (i in 1:5) {
  row.names(query_exp_duplicate) <- paste0(row.names(query_exp), "_", i)
  query_exp_duplicate_new <- rbind(query_exp_duplicate,query_exp_duplicate_new)
}
identical(colnames(query_exp_duplicate_new),
          row.names(C2@meta.data))
query_obj <- CreateSeuratObject(counts = query_exp_duplicate_new,
                                meta.data = C2@meta.data)
trans_Au_Ac_obj <- trans_species_obj(query_obj = query_obj,
                                     subject_obj = Ac_neoblast,
                                     map_gene = map_gene_Au_Ac)
trans_Au_Da_obj <- trans_species_obj(query_obj = query_obj,
                                     subject_obj = Danio_rerio_sce_Me,
                                     map_gene = map_gene_Au_Da)
trans_Au_Hy_obj <- trans_species_obj(query_obj = query_obj,
                                     subject_obj = Hy_iCell,
                                     map_gene = map_gene_Au_Hy)
trans_Au_Xe_obj <- trans_species_obj(query_obj = query_obj,
                                     subject_obj = ROCs_sce,
                                     map_gene = map_gene_Au_Xp)
trans_Au_Sc_obj <- trans_species_obj(query_obj = query_obj,
                                     subject_obj = Sc_neoblast,
                                     map_gene = map_gene_Au_Sc)
trans_Au_La_obj <- trans_species_obj(query_obj = query_obj,
                                     subject_obj = Sc_neoblast,
                                     map_gene = map_gene_Au_Sc)

C2$species <- 'AU'
merge_all_species <- merge(C2,
                           list(trans_Au_Ac_obj,
                                trans_Au_Da_obj,
                                trans_Au_Hy_obj,
                                trans_Au_Xe_obj,
                                trans_Au_Sc_obj))
merge_all_species <- JoinLayers(merge_all_species)
merge_all_species <- NormalizeData(merge_all_species) %>% ScaleData()
Idents(merge_all_species) <- 'species'
merge_all_species_exp <- AverageExpression(merge_all_species,assays = 'RNA')$RNA
write.csv(merge_all_species_exp,'~/Aurelia/regen/merge_all_species_exp.csv')
write_rds(merge_all_species,file = '~/Aurelia/regen/merge_all_species.rds')
# ==============================================================================
#  Integration (Harmony, CCA, RPCA)
# ==============================================================================
sce <- read_rds('~/Aurelia/regen/merge_all_species.rds')
sce <- JoinLayers(sce)
sce <- CreateSeuratObject(sce@assays$RNA, meta.data = sce@meta.data)
sce[["RNA"]] <- split(sce[["RNA"]], f = sce$species)

sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)

# Harmony
sce <- IntegrateLayers(object = sce, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = "integrated.harmony")
sce <- RunUMAP(sce, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

# CCA
sce_cca <- IntegrateLayers(object = sce, method = CCAIntegration,
                           orig.reduction = "pca", new.reduction = "integrated.cca")
sce_cca <- RunUMAP(sce_cca, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

# RPCA
sce_rpca <- IntegrateLayers(object = sce, method = RPCAIntegration,
                            orig.reduction = "pca", new.reduction = "integrated.rpca")
sce_rpca <- RunUMAP(sce_rpca, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

# Clustering
sce <- FindNeighbors(sce, reduction = "integrated.harmony", dims = 1:30)
sce <- FindClusters(sce, resolution = c(0.2, 0.5, 1))

# Plot UMAP
pdf('~/Aurelia/regen/merge_all_species_umap.pdf', width = 10, height = 10)
DimPlot(sce, group.by = 'species')
DimPlot(sce_cca, reduction = 'umap.cca', group.by = 'species')
DimPlot(sce_rpca, reduction = 'umap.rpca', group.by = 'species')
dev.off()

write_rds(sce_cca, '~/Aurelia/regen/merge_all_species.rds')

# Export results
exp <- AverageExpression(sce)
write_rds(exp$RNA, '~/Aurelia/regen/exp.rds')
writexl::write_xlsx(
  list(map_gene_Au_Sc = map_gene_Au_Sc,
       map_gene_Au_Ac = map_gene_Au_Ac,
       map_gene_Au_Da = map_gene_Au_Da,
       map_gene_Au_Hy = map_gene_Au_Hy,
       map_gene_Au_Xp = map_gene_Au_Xp),
  path = '~/Aurelia/regen/all_species_mapping_gene.xlsx'
)
# ==============================================================================