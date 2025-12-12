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

## ========= doublets_pairs_lineage_agnostic_trueIDs_conf.R =========
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## ------------ PARAMETERS ------------
TOP_N      <- 100          # markers per cluster
P_ADJ_MAX  <- 0.05        # significance filter for markers
PREFIX     <- "RefClust1_AGN"  # unique prefix for THIS scoring run
GAP_THR    <- 0.30        # min (top1 - top2) to avoid ambiguous calls
TOP1_FLOOR <- 0.50        # absolute floor for top1 (will be max with Q75)
TOP2_FLOOR <- 0.20        # require non-trivial second program
RATIO_MAX  <- 2.5         # top1/top2 should not be too huge (avoid singlet-like)
ENT_K      <- 5           # how many top scores to use for entropy

## ------------ INPUT ------------
## - Reference Seurat object: spleen.0.4  (Idents: "SCT_snn_res.0.4")
## - Working Seurat object:  obj          (contains DF classifications in meta.data)

## ------------ 0) REFERENCE MARKERS ------------
Idents(spleen.0.4) <- "SCT_snn_res.0.4"

#all_markers <- FindAllMarkers(
 # object = spleen.0.4,
 # only.pos = TRUE,
 # min.pct = 0.25,
 # logfc.threshold = 0.25
#)
all_markers$cluster <- as.character(all_markers$cluster)

topN_by_cluster <- all_markers %>%
  filter(p_val_adj < P_ADJ_MAX) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

## ------------ 1) BUILD GENE SETS & SCORE (unique prefix) ------------
cluster_gene_sets <- split(topN_by_cluster$gene, topN_by_cluster$cluster)
cluster_gene_sets <- lapply(cluster_gene_sets, function(gs) intersect(gs, rownames(obj)))

lens <- vapply(cluster_gene_sets, length, 1L)
if (any(lens == 0L)) {
  message("Dropping empty gene sets: ", paste(names(cluster_gene_sets)[lens == 0], collapse=", "))
}
cluster_gene_sets <- cluster_gene_sets[lens > 0]
cluster_ids <- names(cluster_gene_sets)  # TRUE IDs, e.g., "0","1","2",...
stopifnot(length(cluster_ids) > 0)

# remove old columns from same PREFIX (optional)
old_cols <- grep(paste0("^", PREFIX), colnames(obj@meta.data), value = TRUE)
if (length(old_cols) > 0) {
  obj@meta.data <- obj@meta.data[, setdiff(colnames(obj@meta.data), old_cols), drop = FALSE]
}

# score and capture only columns added now
pre_ncol <- ncol(obj@meta.data)
set.seed(1)
obj <- AddModuleScore(
  object = obj,
  features = cluster_gene_sets,
  name = PREFIX,
  nbin = 24, ctrl = 100, seed.use = 1
)
all_new <- colnames(obj@meta.data)[(pre_ncol + 1):ncol(obj@meta.data)]
stopifnot(length(all_new) >= length(cluster_gene_sets))
new_cols <- all_new[seq_along(cluster_gene_sets)]  # one score col per gene set, ordered

# map score column (this run) -> TRUE cluster ID (no parsing)
col2id <- setNames(cluster_ids, new_cols)  # e.g., "RefClust50_AGN1" -> "0"

## ------------ 2) LINEAGE-AGNOSTIC PAIRS (top1/top2 on DF doublets) ------------
df_col <- grep("^DF\\.classifications_", colnames(obj@meta.data), value = TRUE)[1]
stopifnot(!is.na(df_col))

# use ONLY this run's columns
ref_cols <- new_cols
stopifnot(length(ref_cols) > 1)

ref_scores_full <- obj@meta.data[, ref_cols, drop = FALSE]

# collapse by TRUE ID in case of duplicates (rare)
by_id <- split(colnames(ref_scores_full), col2id[colnames(ref_scores_full)])
ref_scores_cons <- do.call(
  cbind,
  lapply(by_id, function(cols) rowMeans(ref_scores_full[, cols, drop = FALSE], na.rm = TRUE))
)
colnames(ref_scores_cons) <- names(by_id)  # TRUE IDs as column names

# subset to DF doublets
doublet_cells <- rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Doublet"]
rs <- as.matrix(ref_scores_cons[doublet_cells, , drop = FALSE])
stopifnot(nrow(rs) > 0, ncol(rs) > 1)

# top1/top2 and gap
ord_idx <- t(apply(rs, 1, function(v) order(v, decreasing = TRUE)[1:2]))
c1_idx  <- ord_idx[, 1]; c2_idx <- ord_idx[, 2]
c1_id    <- colnames(rs)[c1_idx]                   # TRUE IDs like "0","1",...
c2_id    <- colnames(rs)[c2_idx]
c1_score <- rs[cbind(rownames(rs), c1_id)]
c2_score <- rs[cbind(rownames(rs), c2_id)]
gap      <- c1_score - c2_score

pairs_df <- data.frame(
  cell_barcode = rownames(rs),
  c1_id        = c1_id,
  c1_score     = as.numeric(c1_score),
  c2_id        = c2_id,
  c2_score     = as.numeric(c2_score),
  gap          = as.numeric(gap),
  stringsAsFactors = FALSE
)
pairs_df$pair <- mapply(function(a,b) paste(sort(c(a,b)), collapse = "+"),
                        pairs_df$c1_id, pairs_df$c2_id)

## ------------ 3) CONFIDENCE FILTERS (gap/top1/top2/ratio/entropy) ------------
# thresholds
all_scores_vec <- as.numeric(as.matrix(obj@meta.data[, ref_cols, drop = FALSE]))
t1_q75 <- quantile(all_scores_vec, 0.75, na.rm = TRUE)
top1_thr <- max(TOP1_FLOOR, t1_q75)

pairs_df$ratio_t1_t2 <- pairs_df$c1_score / pmax(pairs_df$c2_score, 1e-9)
pairs_df$pass_gap    <- pairs_df$gap      >= GAP_THR
pairs_df$pass_t1     <- pairs_df$c1_score >= top1_thr
pairs_df$pass_t2     <- pairs_df$c2_score >= TOP2_FLOOR
pairs_df$pass_ratio  <- pairs_df$ratio_t1_t2 <= RATIO_MAX

# entropy of top K scores (captures “two-peak-ness” vs singlet/diffuse)
entropy <- function(x) {
  x <- pmax(x, 0); s <- sum(x)
  if (s <= 0) return(NA_real_)
  p <- x / s; p <- p[p > 0]
  -sum(p * log(p))
}
rs_row <- ref_scores_cons[pairs_df$cell_barcode, , drop = FALSE]
pairs_df$entropy_topK <- sapply(rownames(rs_row), function(cb){
  v <- sort(as.numeric(rs_row[cb, ]), decreasing = TRUE)[1:min(ENT_K, ncol(rs_row))]
  entropy(v)
})
ent_low  <- quantile(pairs_df$entropy_topK, 0.25, na.rm = TRUE)
ent_high <- quantile(pairs_df$entropy_topK, 0.90, na.rm = TRUE)
pairs_df$pass_entropy <- pairs_df$entropy_topK >= ent_low & pairs_df$entropy_topK <= ent_high

pairs_df$Pair_Confidence <- ifelse(
  pairs_df$pass_gap & pairs_df$pass_t1 & pairs_df$pass_t2 & pairs_df$pass_ratio & pairs_df$pass_entropy,
  "strong",
  ifelse(pairs_df$pass_gap & pairs_df$pass_t2, "moderate", "weak")
)

## ------------ 4) WRITE RESULTS INTO obj@meta.data ------------
pairs_df <- pairs_df[pairs_df$cell_barcode %in% colnames(obj), , drop = FALSE]
rownames(pairs_df) <- pairs_df$cell_barcode

for (cc in c("Pair_Top1ID","Pair_Top1Score","Pair_Top2ID","Pair_Top2Score","Pair_Gap",
             "Pair_Unordered","Pair_Confidence")) {
  if (!cc %in% colnames(obj@meta.data)) obj@meta.data[[cc]] <- NA
}
obj@meta.data[rownames(pairs_df), "Pair_Top1ID"]       <- pairs_df$c1_id
obj@meta.data[rownames(pairs_df), "Pair_Top1Score"]    <- pairs_df$c1_score
obj@meta.data[rownames(pairs_df), "Pair_Top2ID"]       <- pairs_df$c2_id
obj@meta.data[rownames(pairs_df), "Pair_Top2Score"]    <- pairs_df$c2_score
obj@meta.data[rownames(pairs_df), "Pair_Gap"]          <- pairs_df$gap
obj@meta.data[rownames(pairs_df), "Pair_Unordered"]    <- pairs_df$pair
obj@meta.data[rownames(pairs_df), "Pair_Confidence"]   <- pairs_df$Pair_Confidence

## ------------ 5) SUMMARIES & SAVES ------------
df_col <- grep("^DF\\.classifications_", colnames(obj@meta.data), value = TRUE)[1]
tot_df_doublets <- sum(obj@meta.data[[df_col]] == "Doublet", na.rm = TRUE)
with_pair       <- sum(!is.na(obj@meta.data$Pair_Unordered))
message("DF doublets: ", tot_df_doublets, " | with top2 pairs: ", with_pair)

pair_counts_all <- sort(table(obj@meta.data$Pair_Unordered), decreasing = TRUE)
pair_frac_all   <- round(100 * prop.table(pair_counts_all), 2)

strong_pairs <- pairs_df$pair[pairs_df$Pair_Confidence == "strong"]
strong_counts <- sort(table(strong_pairs), decreasing = TRUE)
strong_frac   <- round(100 * prop.table(strong_counts), 2)

# saves (OS-safe filenames)
write.csv(pairs_df, "AllDoublets_Top2Pairs_perCell_trueIDs_withConfidence.csv", row.names = FALSE)
write.csv(data.frame(pair=names(pair_counts_all), count=as.integer(pair_counts_all),
                     frac_percent=as.numeric(pair_frac_all[names(pair_counts_all)])),
          "AllDoublets_PairCounts_trueIDs_ALL.csv", row.names = FALSE)
write.csv(data.frame(pair=names(strong_counts), count=as.integer(strong_counts),
                     frac_percent=as.numeric(strong_frac[names(strong_counts)])),
          "AllDoublets_PairCounts_trueIDs_STRONG.csv", row.names = FALSE)

## ------------ 6) OPTIONAL UMAP VIEWS ------------
if (!is.null(obj@reductions$umap)) {
  dbl_cells <- rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Doublet"]
  print(
    DimPlot(obj, cells = dbl_cells, group.by = "Pair_Top1ID", label = FALSE, pt.size = 0.7) +
      ggtitle("Doublets by Top1 cluster (true IDs; lineage-agnostic)")
  )
  topK <- names(head(strong_counts, n = min(6, length(strong_counts))))
  if (length(topK) > 0) {
    subset_cells <- rownames(obj@meta.data)[obj@meta.data$Pair_Unordered %in% topK &
                                             obj@meta.data$Pair_Confidence == "strong"]
    print(
      DimPlot(subset(obj, cells = subset_cells),
              group.by = "Pair_Unordered",
              split.by = "Pair_Unordered",
              ncol = 3, pt.size = 0.8) +
        ggtitle("Top doublet pairs (STRONG only)")
    )
  }
}

## ========= end =========


########################################################################################################## Graph
## ===== stacked_bar_doublet_pairs_by_confidence.R =====
## Builds a pair summary from obj@meta.data and plots a stacked bar chart
## Requirements in obj@meta.data: Pair_Unordered (e.g., "0+6"), Pair_Confidence (weak/moderate/strong)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## 0) Sanity checks -----------------------------------------------------------
req_cols <- c("Pair_Unordered", "Pair_Confidence")
stopifnot(all(req_cols %in% colnames(obj@meta.data)))

## 1) Build per-cell pairs table (collapse mirrored pairs; drop self-pairs) ---
pairs_df <- obj@meta.data[, req_cols, drop = FALSE] %>%
  tibble::rownames_to_column("cell_barcode") %>%
  rename(pair = Pair_Unordered, Pair_Confidence = Pair_Confidence) %>%
  filter(!is.na(pair), !is.na(Pair_Confidence)) %>%
  separate(pair, into = c("i","j"), sep = "\\+") %>%
  mutate(i_num = as.numeric(i), j_num = as.numeric(j)) %>%
  # keep numeric cluster IDs & remove diagonal (i == j)
  filter(!is.na(i_num), !is.na(j_num), i_num != j_num) %>%
  # collapse mirrored pairs: enforce i<j
  rowwise() %>%
  mutate(pair2 = paste(sort(c(i_num, j_num)), collapse = "+")) %>%
  ungroup()

## 2) Count by pair × confidence, and order pairs by TOTAL count --------------
pair_flat <- pairs_df %>%
  count(pair2, Pair_Confidence, name = "count")

pair_totals <- pair_flat %>%
  group_by(pair2) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total))

# ensure stacked order weak < moderate < strong (bottom→top)
pair_flat <- pair_flat %>%
  mutate(Pair_Confidence = factor(Pair_Confidence,
                                  levels = c("weak","moderate","strong"))) %>%
  # order x-axis by total descending
  mutate(pair2 = factor(pair2, levels = pair_totals$pair2))

## 3) Choose y-axis tick density (more numbers) --------------------------------
ymax <- max(pair_totals$total, 0)
ystep <- if (ymax <= 10) 1 else if (ymax <= 20) 2 else if (ymax <= 50) 5 else 10
ybreaks <- seq(0, ceiling(ymax / ystep) * ystep, by = ystep)

## 4) Segment labels (counts) + top total labels ------------------------------
# data for total labels on top of bars
top_labels <- pair_totals %>%
  mutate(pair2 = factor(pair2, levels = levels(pair_flat$pair2)))

p <- ggplot(pair_flat, aes(x = pair2, y = count, fill = Pair_Confidence)) +
  geom_bar(stat = "identity") +
  # print per-segment counts inside stacks
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  # total label at the top of each bar
  geom_text(data = top_labels,
            aes(x = pair2, y = total, label = total),
            vjust = -0.5, fontface = "bold", size = 3.2, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(weak="#d9d9d9", moderate="#a6bddb", strong="#1c9099"),
    name   = "Confidence"
  ) +
  scale_y_continuous(
    name   = "Count",
    breaks = ybreaks,
    expand = expansion(mult = c(0.02, 0.10))  # a bit of headroom for top labels
  ) +
  labs(
    title = "Doublet partner pairs by confidence (collapsed A+B ≡ B+A)",
    x = "Pair (cluster IDs)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    plot.title    = element_text(face = "bold")
  )

print(p)

## 5) (Optional) save to file --------------------------------------------------
# ggsave("Doublet_Pairs_byConfidence_stacked.png", p, width = 12, height = 6, dpi = 300)
# ggsave("Doublet_Pairs_byConfidence_stacked.pdf", p, width = 12, height = 6)
## ===== end =====


