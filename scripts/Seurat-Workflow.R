library(dplyr)
library(Seurat)
library(patchwork)
library(loupeR)
library(clustree)

dataDir <- gene-barcode matrices

set.seed(50)
spleen.data <- Read10X(data.dir = dataDir)
spleen <- CreateSeuratObject(counts = spleen.data, project = "spleen_mmc", min.cells = 3, min.features = 200) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
    SCTransform(vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    RunTSNE(dims = 1:30) %>%
    FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) 

# Stability of cluster resolution 

stability_df <- clustree(spleen, prefix = "SCT_snn_res.", node_colour = "sc3_stability")
stability_df$data %>% group_by(SCT_snn_res.) %>% summarise(mean(sc3_stability))

png("FigureM1_clustree.png", width= 3500, height = 2500, res = 300)
stability_df
dev.off()

### Create Loupe File for Loupe Browser

create_loupe_from_seurat(spleen)

##### Splitting data to Parent and Subtype Resolutions

library(reshape2)
spleen.0.4 <- CreateSeuratObject(counts = spleen.data, project = "spleen_mmc", min.cells = 3, min.features = 200) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
    SCTransform(vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30, reduction.name = "umap") %>%
    FindClusters(resolution = c(0.4)) 

##### Finding all markers per split file

spleen.markers.0.4 <- FindAllMarkers(spleen.0.4, logfc.threshold = 0, return.thresh = 0.05) 
SCluster <- group_by(spleen.markers.0.4, cluster)
write.csv(dcast(SCluster, formula = p_val_adj + gene ~ cluster,value.var = "avg_log2FC"), "Subtype_Res.0.4_dge.csv")

##### Renaming clusters in parent and subtype resolution (Placeholder for now)
    
new.cluster.ids <- c("S0", "S1", "S2", "S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17")
names(new.cluster.ids) <- levels(spleen.0.4)
spleen.0.4 <- RenameIdents(spleen.0.4, new.cluster.ids)

##### UMAP Dimplots

library(ggplot2)
DimPlot(spleen.0.4, reduction = "umap", label = TRUE, label.size = 4) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 12), legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 6)))

###### DoHeatmaps

spleen.markers.0.4 %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 2) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

tiff("FigureM2b_Subtype_Res.0.4_doheatmap.tiff", width = 20000, height = 20000, compression = "none", res = 600, pointsize = 12)
DoHeatmap(spleen.0.4, features = top10$gene, size = 3.5, angle = 90, slot = "scale.data")
dev.off()

##### Feature Counts

VlnPlot(spleen.0.4, features = "nFeature_RNA") + ylab("Number of gene features") + ggtitle("Subtype Resolution") + xlab("Cluster")

###### Finding doublets

install.packages('BiocManager')
BiocManager::install('glmGamPoi')

spleen.doublets <- CreateSeuratObject(counts = spleen.data, project = "spleen_mmc", min.cells = 3, min.features = 200) %>%
    SCTransform() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30, reduction.name = "umap") 

pK <- paramSweep(spleen.doublets, PCs = 1:10, sct = TRUE) %>% summarizeSweep(GT = FALSE) %>% find.pK()
maxpK <- 0.2
nExp <- round(0.008 * ncol(spleen.doublets))

spleen_doublets <- doubletFinder(
  spleen.doublets, PCs = 1:30, pN = 0.25, pK = maxpK,
  nExp = nExp, sct = True
)

obj <- spleen_doublets

## 1) Get the DoubletFinder column
df_col <- grep("^DF\\.classifications_", colnames(obj@meta.data), value = TRUE)[1]

## 2) B-only program score (keep just a B set)
bcell_genes <- c("ENSDARG00000078975","ENSDARG00000093272","cxcr4a","ebf1a","cd74a","cd74b","pax5")
bcell_genes <- intersect(rownames(obj), bcell_genes)
obj <- AddModuleScore(obj, features = list(bcell_genes), name = "BcellScore")
bscore <- obj$BcellScore1

b_any_cells <- colnames(obj)[obj$B_to_ANY %in% TRUE]
length(b_any_cells)  # how many you found

obj_b_any <- subset(obj, cells = b_any_cells)

# Define groups: B_to_ANY vs B singlets (strong B-score + DF singlet)
bscore <- obj$BcellScore1
b_thr  <- median(bscore, na.rm = TRUE) + 0.15
B_singlets_idx <- which(bscore >= b_thr & obj@meta.data[[df_col]] == "Singlet")

Idents(obj) <- factor(ifelse(colnames(obj) %in% b_any_cells, "B_to_ANY",
                      ifelse(colnames(obj) %in% colnames(obj)[B_singlets_idx], "B_singlet", "Other")),
                      levels = c("B_to_ANY","B_singlet","Other"))

# working (below)

################################################################## Build “top-50 per cluster” gene sets from the reference

all_markers <- FindAllMarkers(
  spleen.0.4,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Make cluster IDs character and take top 50 by effect size (change to p_val_adj if preferred)
all_markers$cluster <- as.character(all_markers$cluster)

all_markers_filt <- all_markers %>%
  filter(p_val_adj < 0.05)

top50_by_cluster <- all_markers_filt %>%
   group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup()

# One gene set per cluster
cluster_gene_sets <- split(top50_by_cluster$gene, top50_by_cluster$cluster)

# Keep only genes present in obj (ID-space match)
cluster_gene_sets <- lapply(cluster_gene_sets, function(gs) intersect(gs, rownames(obj)))

# Drop any empty sets (rare, but safe)
cluster_gene_sets <- cluster_gene_sets[vapply(cluster_gene_sets, length, 1L) > 0]

cluster_ids <- names(cluster_gene_sets)  # e.g., "0","1","2",...
length(cluster_ids)

obj <- AddModuleScore(obj, features = cluster_gene_sets, name = "RefClust50")

# Build a tidy mapping from score columns to reference cluster IDs
score_cols <- grep("^RefClust50", colnames(obj@meta.data), value = TRUE)
stopifnot(length(score_cols) == length(cluster_ids))
cluster_score_map <- setNames(score_cols, cluster_ids)  # named vector: cluster_id -> score column
head(cluster_score_map)

