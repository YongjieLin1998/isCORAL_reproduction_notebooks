# Set ISCORAL_PROJECT_ROOT, ISCORAL_SAMPLE_ID, ISCORAL_INPUT_RDATA,
# ISCORAL_SCAN_RESULTS, ISCORAL_PACKAGE_DIR, and
# ISCORAL_EMBEDDING_OUTPUT_DIR as needed.
BASE_DIR  <- normalizePath(Sys.getenv("ISCORAL_PROJECT_ROOT", unset = "."), mustWork = TRUE)
SAMPLE_ID <- Sys.getenv("ISCORAL_SAMPLE_ID", unset = "<SAMPLE_ID>")
if (identical(SAMPLE_ID, "<SAMPLE_ID>")) {
  stop("Set ISCORAL_SAMPLE_ID or replace <SAMPLE_ID> before running this template.")
}

PACKAGE_DIR <- Sys.getenv(
  "ISCORAL_PACKAGE_DIR",
  unset = file.path(BASE_DIR, "packages", "isCORAL_v6.1.4")
)
INPUT_RDATA <- Sys.getenv(
  "ISCORAL_INPUT_RDATA",
  unset = file.path(BASE_DIR, "1_data", "raw", paste0(SAMPLE_ID, ".RData"))
)
SCAN_RES_FILE <- Sys.getenv(
  "ISCORAL_SCAN_RESULTS",
  unset = file.path(BASE_DIR, "3_results_new", "<SCAN_RESULTS_FILE>.rds")
)
BINSIZE_KEY <- Sys.getenv("ISCORAL_BINSIZE_KEY", unset = "binsize_200")

devtools::load_all(PACKAGE_DIR)

global_start_time <- Sys.time()
message(sprintf("\n🚀 [Pipeline Start] Task Start Time: %s", format(global_start_time, "%Y-%m-%d %H:%M:%S")))

# ------------------------------------------------------------------------------
# 1. Environment Setup and Library Loading
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(irlba)
  library(glmnet)
  library(matrixStats)
  library(doParallel)
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggrepel)
  library(cowplot)
  library(scales)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggdist)
  library(scCustomize)
  library(extrafont)
})

theme_set(theme_classic(base_family = "Arial")) 
update_geom_defaults("text", list(family = "Arial"))
update_geom_defaults("label", list(family = "Arial"))

export_dir <- Sys.getenv(
  "ISCORAL_EMBEDDING_OUTPUT_DIR",
  unset = file.path(BASE_DIR, "5_embeddings", SAMPLE_ID)
)
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

CELL_TYPE_COL <- Sys.getenv("ISCORAL_CELL_TYPE_COL", unset = "cluster_integrated")
BARCODE_COL   <- Sys.getenv("ISCORAL_BARCODE_COL", unset = "larry_sapphire_merged")

my_colors <- c("HSC" = "#E69F00", "GMP" = "#56B4E9", "Mk" = "#009E73", "Other" = "#000000")
highlight_colors <- my_colors[c("HSC", "GMP", "Mk", "Other")]

# ------------------------------------------------------------------------------
# 2. Raw Data Loading and Graph Probability Merging
# ------------------------------------------------------------------------------
t_load_start <- Sys.time()

temp_env <- new.env()
if (!file.exists(INPUT_RDATA)) {
  stop(sprintf("Input RData file not found: %s", INPUT_RDATA))
}
load(INPUT_RDATA, envir = temp_env)
obj_names <- ls(temp_env)
is_seurat <- sapply(obj_names, function(x) inherits(temp_env[[x]], "Seurat"))
if(sum(is_seurat) == 0) stop("🚨 Seurat object not found in RData!")
sample <- temp_env[[obj_names[is_seurat][1]]]
rm(temp_env, obj_names, is_seurat)

required_cols <- c(CELL_TYPE_COL, BARCODE_COL)
missing_cols <- setdiff(required_cols, colnames(sample@meta.data))
if(length(missing_cols) > 0) {
  warning(sprintf("🚨 Warning: Missing metadata columns [%s]", paste(missing_cols, collapse=", ")))
}

if (!file.exists(SCAN_RES_FILE)) {
  stop(sprintf("isCORAL scan result file not found: %s", SCAN_RES_FILE))
}
scanning_results <- readRDS(SCAN_RES_FILE)

if (is.null(sample@misc$coral_scan_results)) sample@misc$coral_scan_results <- list()
if ("S_stab" %in% names(scanning_results)) {
  sample@misc$coral_scan_results[[BINSIZE_KEY]] <- scanning_results
} else {
  sample@misc$coral_scan_results <- c(sample@misc$coral_scan_results, scanning_results)
}

t_load_end <- Sys.time()

# ------------------------------------------------------------------------------
# 3. Stringent Preprocessing and Baseline DimReduc Construction
# ------------------------------------------------------------------------------
t_pre_start <- Sys.time()

sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)

var_genes <- VariableFeatures(sample)
mt_ribo_genes <- grep(pattern = "^mt-|^rps|^rpl", x = var_genes, value = TRUE, ignore.case = TRUE)
var_genes_filtered <- setdiff(var_genes, mt_ribo_genes)

VariableFeatures(sample) <- var_genes_filtered


features_to_scale <- VariableFeatures(sample)
#features_to_scale <- unique(c(VariableFeatures(sample), c("Hlf", "Mpo", "Pf4")))
# useless combination

sample <- ScaleData(sample, features = features_to_scale)

sample <- RunPCA(sample, features = VariableFeatures(sample), npcs = 50, reduction.name = "pca.unfilt", reduction.key = "PCAUNFILT_")
sample <- FindNeighbors(sample, dims = 1:50, reduction = "pca.unfilt", k.param = 20, graph.name = c("unfilt_nn", "unfilt_snn"))
sample <- RunUMAP(sample, dims = 1:50, reduction = "pca.unfilt", n.neighbors = 30, reduction.name = "umap.unfilt", reduction.key = "UMAPUNFILT_")

s.genes.mouse <- str_to_title(cc.genes$s.genes)
g2m.genes.mouse <- str_to_title(cc.genes$g2m.genes)
sample <- CellCycleScoring(sample, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse, set.ident = FALSE)

VariableFeatures(sample) <- var_genes_filtered
sample <- ScaleData(sample, features = features_to_scale, vars.to.regress = c("S.Score", "G2M.Score"))

sample <- RunPCA(sample, features = VariableFeatures(sample), npcs = 50, reduction.name="pca.raw")
sample <- FindNeighbors(sample, dims = 1:50, reduction="pca.raw", k.param = 20, graph.name = c("raw_nn", "raw_snn"))
sample <- RunUMAP(sample, dims = 1:50, reduction = "pca.raw", n.neighbors = 30, reduction.name = "umap.raw", reduction.key = "UMAPRAW_")

t_pre_end <- Sys.time()

# ------------------------------------------------------------------------------
# 4. Pipeline Execution via Pre-loaded isCORAL Functions
# ------------------------------------------------------------------------------
t_coral_start <- Sys.time()

sample <- CORALembedding(seurat_obj = sample, binsize_key = BINSIZE_KEY, n_dims = 30, workers_glmnet = 5)
sample <- SmoothEmbedding(seurat_obj = sample, lambda = 0.5)

t_coral_end <- Sys.time()

# ------------------------------------------------------------------------------
# 5. Scientific Visualization and PDF Archiving
# ------------------------------------------------------------------------------
t_plot_start <- Sys.time()

if(CELL_TYPE_COL %in% colnames(sample@meta.data)){
  sample$lineage_group <- as.character(sample[[CELL_TYPE_COL]][[1]])
  sample$lineage_group[!sample$lineage_group %in% c("HSC", "GMP", "Mk")] <- "Other"
  sample$lineage_group <- factor(sample$lineage_group, levels = c("HSC", "GMP", "Mk", "Other"))
}

# ✅ [Core update]: Add is_hollow to convert points to hollow circles (shape = 1).
create_base_plot_enhanced <- function(red_name, group_col, colors, title, type = "normal", is_hollow = FALSE) {
  p <- DimPlot(sample, reduction = red_name, group.by = group_col, pt.size = 0.4, label = (type=="normal"), repel = TRUE, label.size = 3.5) +
    scale_color_manual(values = colors) + theme_void() + ggtitle(title) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  
  if (is_hollow) {
    p$layers[[1]]$aes_params$shape <- 1  # Force hollow circles in the base scatter layer.
  }
  return(p)
}

# [Output 1] 4-panel manifold grid (solid points above, hollow points below)
p_grid_A <- create_base_plot_enhanced("umap.raw", "lineage_group", highlight_colors, "Raw UMAP", type="none", is_hollow=FALSE)
p_grid_B <- create_base_plot_enhanced("CORAL_cell_embed", "lineage_group", highlight_colors, "CORAL Cell Space", type="none", is_hollow=FALSE)
p_grid_C <- create_base_plot_enhanced("CORAL_smooth_lambda_05", "lineage_group", highlight_colors, "CORAL Smooth (λ=0.5)", type="none", is_hollow=TRUE)
p_grid_D <- create_base_plot_enhanced("CORAL_clone_embed", "lineage_group", highlight_colors, "CORAL Clone Space", type="none", is_hollow=TRUE)

p_eg_grid <- (p_grid_A | p_grid_B) / (p_grid_C | p_grid_D) +
  plot_annotation(title = "Native Lineage Communities Across Embeddings", theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
ggsave(filename = file.path(export_dir, "Manifold_4Panel_Grid.pdf"), plot = p_eg_grid, width = 12, height = 10, device = cairo_pdf)

# [Output 2] 1x3 Triptych Manifold Comparison
plot_A <- create_base_plot_enhanced("umap.raw", "lineage_group", highlight_colors, "A. Raw UMAP", is_hollow=FALSE)
plot_B <- create_base_plot_enhanced("CORAL_cell_embed",  "lineage_group", highlight_colors, "B. CORAL Cell Space", is_hollow=FALSE)
plot_C <- create_base_plot_enhanced("CORAL_clone_embed",  "lineage_group", highlight_colors, "C. CORAL Clone Space", is_hollow=TRUE) # Use hollow points here.

triptych_manifold <- plot_A + plot_B + plot_C + plot_layout(ncol = 3) + plot_annotation(title = "Adaptive Spatial Projections of HSC Lineages", theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
ggsave(filename = file.path(export_dir, "Manifold_Triptych_Comparison.pdf"), plot = triptych_manifold, width = 15, height = 5, device = cairo_pdf)

# [Output 3] Empirical Clonal Continuity Scatter Tracking
if(BARCODE_COL %in% colnames(sample@meta.data)){
  valid_barcodes <- sample[[BARCODE_COL]][[1]]
  valid_barcodes <- valid_barcodes[!is.na(valid_barcodes) & valid_barcodes != "NA"]
  
  if (length(valid_barcodes) > 0) {
    bc_counts     <- table(valid_barcodes)
    robust_clones <- names(bc_counts[bc_counts >= 4]) 
    
    if (length(robust_clones) > 0) {
      target_bc     <- base::sample(robust_clones, 1)
      target_cells  <- colnames(sample)[which(sample[[BARCODE_COL]][[1]] == target_bc)]
      target_count  <- length(target_cells)
      message(sprintf("  -> Clonal tracking locked Barcode: '%s' (n=%d)", target_bc, target_count))
      
      # ✅ [Highlight adaptation]: Use prominent hollow circles for highlighted cells in clone space.
      p_bc_base <- function(red, tit, is_hollow = FALSE) {
        p <- DimPlot(sample, reduction = red, cells.highlight = target_cells, 
                     cols.highlight = "#FF1E27", sizes.highlight = 2.2, cols = "#EDEDED", pt.size = 0.4) + 
          theme_void() + ggtitle(tit) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
        
        if (is_hollow) {
          p$layers[[1]]$aes_params$shape <- 1 # Render gray background cells as hollow circles.
          if (length(p$layers) >= 2) {
            p$layers[[2]]$aes_params$shape <- 1 # Render red highlighted cells as hollow circles.
            p$layers[[2]]$aes_params$stroke <- 1.2 # Thicken the highlighted outlines for visibility.
          }
        }
        return(p)
      }
      
      triptych_clonal <- p_bc_base("umap.raw", "A. Raw UMAP Clone Scatter", is_hollow=FALSE) +
        p_bc_base("CORAL_cell_embed",  "B. CORAL Cell Space Clone Scatter", is_hollow=FALSE) +
        p_bc_base("CORAL_clone_embed",  "C. CORAL Clone Space Clone Scatter", is_hollow=TRUE) + # Use hollow points here.
        plot_layout(ncol = 3) + 
        plot_annotation(
          title = sprintf("Empirical Clonal Continuity Benchmark (Barcode: %s | n=%d)", target_bc, target_count), 
          theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
        )
      
      ggsave(filename = file.path(export_dir, sprintf("Clonal_Triptych_%s.pdf", target_bc)), plot = triptych_clonal, width = 15, height = 5, device = cairo_pdf)
    }
  }
}

# [Output 4] Specific Hematopoietic Lineage Highlights
p_highlight_clean <- DimPlot(
  sample, reduction = "umap.raw", group.by = "lineage_group",
  cols = highlight_colors, order = c("HSC", "GMP", "Mk"), pt.size = 0.5
) +  
  ggtitle("Early HSC Expansion") +
  theme(
    axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15, face = "plain"),
    legend.text = element_text(size = 14), legend.title = element_blank(), legend.key.size = unit(1.5, "lines")   
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(filename = file.path(export_dir, "HSC_original_UMAP.pdf"), plot = p_highlight_clean, device = cairo_pdf)

# [Output 5] Core Marker Density Layer
feature_highlight <- scCustomize::Plot_Density_Custom(sample, c("Hlf", "Mpo", "Pf4"), reduction = "umap.raw")
ggsave(filename = file.path(export_dir, "HSC_original_featurePlot.pdf"), plot = feature_highlight, device = cairo_pdf)

t_plot_end <- Sys.time()

# ------------------------------------------------------------------------------
# 6. Pipeline Execution Summary
# ------------------------------------------------------------------------------
global_end_time <- Sys.time()
total_time_mins <- as.numeric(difftime(global_end_time, global_start_time, units = "mins"))

message("\n==============================================================================")
message("✅ Pipeline executed successfully!")
message(sprintf("🏁 Total Runtime: %.2f minutes", total_time_mins))
message("==============================================================================")
