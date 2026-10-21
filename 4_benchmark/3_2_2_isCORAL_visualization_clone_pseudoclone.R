# ==============================================================================
# --- [B] Matrix-level plotting with automatic cache adaptation, template-size detection, and PDF export ---
# ==============================================================================
message("\n--- Starting the B3 matrix visualization renderer with automatic cache and template discovery ---")

library(patchwork)
library(ggplot2)
library(dplyr)
library(Seurat)
library(ggdist)

library(extrafont)
# Load fonts for the PDF device
loadfonts(device = "pdf")

# -----------------------------------------------------------------
# 🟢 [0] Automatically load the latest B1 evaluation results in an isolated environment
# -----------------------------------------------------------------
# Set ISCORAL_PROJECT_ROOT and ISCORAL_SAMPLE_ID, or replace the sample placeholder below.
BASE_DIR      <- normalizePath(Sys.getenv("ISCORAL_PROJECT_ROOT", unset = "."), mustWork = TRUE)
SAMPLE_ID     <- Sys.getenv("ISCORAL_SAMPLE_ID", unset = "<SAMPLE_ID>")
if (identical(SAMPLE_ID, "<SAMPLE_ID>")) {
  stop("Set ISCORAL_SAMPLE_ID or replace <SAMPLE_ID> before running this template.")
}
RESULTS_DIR   <- file.path(BASE_DIR, "4_benchmark", "results")

PACKAGE_DIR <- Sys.getenv(
  "ISCORAL_PACKAGE_DIR",
  unset = file.path(BASE_DIR, "packages", "isCORAL_v6.1.0")
)
devtools::load_all(PACKAGE_DIR)

suppressPackageStartupMessages({
  library(ggplot2); library(scCustomize); library(ComplexHeatmap); library(circlize)
  library(ggridges); library(dplyr); library(tidyr); library(ggrepel)
  library(patchwork); library(svglite)
})



# 🎯 Discover and load the latest B1 cache automatically
all_b1_caches <- list.files(path = RESULTS_DIR, pattern = sprintf("^Eval_Cache_B1_%s_.*\\.rds$", SAMPLE_ID), recursive = TRUE, full.names = TRUE)
if (length(all_b1_caches) == 0) stop("❌ No B1 evaluation cache was found for the current sample. Run Script B1 first.")
latest_b1_cache_path <- rownames(file.info(all_b1_caches))[which.max(file.info(all_b1_caches)$mtime)]

message(sprintf("  🚀 Loading the latest evaluation data:\n     %s", latest_b1_cache_path))
b1_cache <- readRDS(latest_b1_cache_path)

# 🎯 Restore the plotting variables required by this script
plot_roc_df <- b1_cache$plot_roc_df
plot_pr_df  <- b1_cache$plot_pr_df
baseline_raw_matrices <- readRDS(b1_cache$run_parameters$CACHE_FILE)$baseline_raw_matrices

# Add the isCORAL matrix to the shared matrix list
baseline_raw_matrices[["isCORAL"]] <- b1_cache$matrix_raw_iscoral # Downstream extraction provides a fallback when B1 did not retain this matrix

# Load the Seurat object
seurat_obj <- readRDS(b1_cache$run_parameters$CACHE_FILE)$seurat_obj 

# Attach isCORAL results for UMAP visualization
SCAN_RES_FILE <- b1_cache$run_parameters$SCAN_RES_FILE
ISCORAL_BIN_SIZE <- b1_cache$run_parameters$ISCORAL_BIN_SIZE
ISCORAL_BINSIZE_NAME <- paste0("binsize_", ISCORAL_BIN_SIZE)

scanning_results <- readRDS(SCAN_RES_FILE)
seurat_obj@misc$CORAL_scan_results <- scanning_results
if (!is.null(scanning_results$S_stab) && is.null(scanning_results[[ISCORAL_BINSIZE_NAME]])) {
  tmp_list <- list(); tmp_list[[ISCORAL_BINSIZE_NAME]] <- scanning_results
  seurat_obj@misc$CORAL_scan_results <- tmp_list
}
S_stab_mat_fresh <- seurat_obj@misc$CORAL_scan_results[[ISCORAL_BINSIZE_NAME]]$S_stab

# Extract shared variables
GLOBAL_CLONE_COL <- b1_cache$run_parameters$GLOBAL_CLONE_COL
GLOBAL_INVALID_CLONES <- c("0","NA", "NaN", "", "None", "Unknown", "unknown", "Unassigned")
CHOSEN_CLONE_ID <- "11" 
UMAP_REDUCTION <- "bench_umap" 

export_dir <- dirname(latest_b1_cache_path)


# ==============================================================================
# 📏 Automatic template-size detection and PDF export using the raw PDF grep approach from B2
# ==============================================================================
TEMPLATE_DIR <- Sys.getenv(
  "ISCORAL_PLOT_TEMPLATE_DIR",
  unset = file.path(BASE_DIR, "4_benchmark", "plothelper", "Benchmark_graph_template")
)

get_template_dims <- function(template_dir, template_name) {
  base_name <- tools::file_path_sans_ext(template_name)
  pdf_path <- file.path(template_dir, paste0(base_name, ".pdf"))
  if (file.exists(pdf_path)) {
    tryCatch({
      cmd <- sprintf("grep -a '/MediaBox' %s | head -n 1", shQuote(pdf_path))
      box_line <- system(cmd, intern = TRUE)
      if (length(box_line) > 0) {
        nums <- regmatches(box_line, gregexpr("[0-9]+\\.?[0-9]*", box_line))[[1]]
        if (length(nums) >= 4) {
          w_in <- as.numeric(nums[3]) / 72
          h_in <- as.numeric(nums[4]) / 72
          message(sprintf("  📏 [Source probe succeeded] PDF template size: %.2f x %.2f in", w_in, h_in))
          return(list(width = w_in, height = h_in))
        }
      }
    }, error = function(e) {})
  }
  return(NULL)
}

smart_ggsave_pdf <- function(filename, plot_obj, dims, default_w, default_h, ...) {
  if (!is.null(dims) && dims$width > 4) {
    ggsave(filename = filename, plot = plot_obj, device = cairo_pdf, width = dims$width, height = dims$height, ...)
  } else {
    message(sprintf("  ⚠️ Template missing or unreadable; using safe default dimensions: %.1f x %.1f in", default_w, default_h))
    ggsave(filename = filename, plot = plot_obj, device = cairo_pdf, width = default_w, height = default_h, ...)
  }
}

# -----------------------------------------------------------------
# Theme and color palette definitions
# -----------------------------------------------------------------
theme_benchmark_vis <- function(base_size = 14) { 
  theme_bw(base_size = base_size, base_family = "Arial") + 
    theme(
      plot.title = element_text(size = rel(1.1), face = "bold", hjust = 0.5, family = "Arial"), 
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5, family = "Arial"), 
      legend.position = "bottom", 
      legend.box = "vertical", 
      strip.background = element_rect(fill = "grey90"), 
      strip.text = element_text(face = "bold", family = "Arial")
    ) 
}

method_colors <- c(
  "isCORAL" = "#C74B46", 
  "isCORAL (Probability Field)" = "#C74B46",
  "GEMLI" = "#3D5A80", 
  "GEMLI (Prediction Score)" = "#3D5A80",
  "CC-Regressed SNN" = "#2A9D8F", 
  "Standard SNN" = "#A8DADC",
  "CC-Regressed KNN" = "#E9C46A", 
  "Standard KNN" = "#F4E285"
)


# -----------------------------------------------------------------
# 📈 Render 5: Matrix-level ROC and PR curves without the theoretical optimum
# -----------------------------------------------------------------
message("\n  -> Render 5: Matrix-level ROC and PR curves without the theoretical optimum...")

# 🎯 0. Remove the theoretical optimum before plotting
plot_roc_df_clean <- plot_roc_df %>% 
  filter(!grepl("Theoretical Best", Method))

plot_pr_df_clean <- plot_pr_df %>% 
  filter(!grepl("Theoretical Best", Method))

# 🎯 1. Extract ordered breaks and labels containing metric values
roc_summary_mat <- plot_roc_df_clean %>% 
  group_by(Method) %>% summarise(auc = max(auc, na.rm=TRUE), .groups="drop") %>% arrange(desc(auc))
roc_breaks_mat <- roc_summary_mat$Method
roc_labels_mat <- setNames(paste0(roc_summary_mat$Method, " (AUC = ", round(roc_summary_mat$auc, 3), ")"), roc_summary_mat$Method)

pr_summary_mat <- plot_pr_df_clean %>% 
  group_by(Method) %>% summarise(auprc = max(auprc, na.rm=TRUE), .groups="drop") %>% arrange(desc(auprc))
pr_breaks_mat <- pr_summary_mat$Method
pr_labels_mat <- setNames(paste0(pr_summary_mat$Method, " (AUPRC = ", round(pr_summary_mat$auprc, 3), ")"), pr_summary_mat$Method)

# 🎯 2. Define the plotting theme
theme_curve_vis_mat <- function() {
  theme_bw(base_size = 14, base_family = "Arial") + 
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, family = "Arial", margin = margin(b = 15)), 
      panel.border = element_rect(color = "#333333", linewidth = 1.5), 
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.5, linetype = "dashed"), 
      panel.grid.minor = element_blank(), 
      axis.text = element_text(color = "black", size = 13, family = "Arial"), 
      axis.title = element_text(color = "black", size = 15,  family = "Arial"), 
      legend.text = element_text(size = 10, family = "Arial"), 
      legend.title = element_blank(), 
      legend.key.width = unit(0.8, "cm"),  
      legend.key.height = unit(0.5, "cm"), 
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA), 
      legend.box.background = element_rect(color = "#333333", linewidth = 0.5), 
      legend.margin = margin(t = 2, r = 4, b = 2, l = 4) 
    )
}

# 🎯 3. Draw the ROC curve
p_mat_roc <- ggplot(plot_roc_df_clean, aes(x = fpr, y = tpr, color = Method)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#A0A0A0", linewidth = 1) + 
  geom_line(aes(linewidth = ifelse(Method == "isCORAL", 1.5, 0.8)), alpha = 0.9) + 
  scale_linewidth_identity() + 
  scale_color_manual(values = method_colors, breaks = roc_breaks_mat, labels = roc_labels_mat) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Matrix Cell-Pair ROC Curve", x = "False Positive Rate", y = "True Positive Rate") + 
  theme_curve_vis_mat() + 
  theme(legend.justification = c(1, 0), legend.position = c(0.97, 0.03))

# 🎯 4. Draw the PR curve
p_mat_pr <- ggplot(plot_pr_df_clean, aes(x = recall, y = precision, color = Method)) + 
  geom_line(aes(linewidth = ifelse(Method == "isCORAL", 1.5, 0.8)), alpha = 0.9) + 
  scale_linewidth_identity() + 
  scale_color_manual(values = method_colors, breaks = pr_breaks_mat, labels = pr_labels_mat) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Matrix Cell-Pair Precision-Recall Curve", x = "Recall", y = "Precision") + 
  theme_curve_vis_mat() + 
  theme(legend.justification = c(1, 1), legend.position = c(0.97, 0.97))

# 🎯 5. Combine and export the panels
p_mat_curves <- p_mat_roc + p_mat_pr + plot_layout(nrow = 2)

print(p_mat_curves)

# 👉 Export using detected template dimensions
out_file_mat_curves <- file.path(export_dir, "Fig.Matrix_ROC_PR_Curves.pdf")
dims_5 <- get_template_dims(TEMPLATE_DIR, "Fig.Matrix_ROC_PR_Curves.pdf")
smart_ggsave_pdf(filename = out_file_mat_curves, plot = p_mat_curves, dims = dims_5, default_w = 12, default_h = 10)
message(sprintf("   ✅ Matrix-curve PDF exported to:\n      %s", out_file_mat_curves))


# ==============================================================================
# 🗺️ Render 8: UMAP probability fields with an absolute enrichment scale, a 500-fold range, and log mapping
# ==============================================================================
message("\n  -> Render 8: UMAP probability fields on a shared absolute fold-enrichment scale...")

# -----------------------------------------------------------------
# 🎯 0. Define the global absolute scale
# -----------------------------------------------------------------
# Extend the enrichment scale to 500-fold to accommodate highly confident biological trajectories
UNIVERSAL_FC_MAX <- 500 

# -----------------------------------------------------------------
# 1. Select and validate the seed cell
# -----------------------------------------------------------------
if (!exists("CHOSEN_SEED_CELL")) CHOSEN_SEED_CELL <- NULL 

if (all(grepl("^[0-9]+$", head(rownames(S_stab_mat_fresh), 50)))) {
  if (nrow(S_stab_mat_fresh) == ncol(seurat_obj)) {
    rownames(S_stab_mat_fresh) <- colnames(seurat_obj); colnames(S_stab_mat_fresh) <- colnames(seurat_obj)
  } else {
    stop("❌ The matrix cell count does not match the Seurat object.")
  }
}

valid_umap_cells <- intersect(rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data[[GLOBAL_CLONE_COL]]) & !(as.character(seurat_obj@meta.data[[GLOBAL_CLONE_COL]]) %in% GLOBAL_INVALID_CLONES)], rownames(S_stab_mat_fresh))

if (!is.null(CHOSEN_CLONE_ID) && CHOSEN_CLONE_ID %in% seurat_obj@meta.data[[GLOBAL_CLONE_COL]]) { 
  gt_label_str <- CHOSEN_CLONE_ID
  gt_cells <- rownames(seurat_obj@meta.data)[which(as.character(seurat_obj@meta.data[[GLOBAL_CLONE_COL]]) == gt_label_str)]
  valid_gt_cells <- intersect(gt_cells, valid_umap_cells)
  if (length(valid_gt_cells) == 0) stop(sprintf("❌ The selected clone %s contains no valid cells.", gt_label_str))
  
  if (!is.null(CHOSEN_SEED_CELL) && CHOSEN_SEED_CELL %in% valid_gt_cells) {
    seed_cell <- CHOSEN_SEED_CELL
  } else {
    seed_cell <- sample(valid_gt_cells, 1) 
  }
} else { 
  valid_labels <- as.character(seurat_obj@meta.data[valid_umap_cells, GLOBAL_CLONE_COL])
  clone_counts <- table(valid_labels)
  candidate_clones <- names(clone_counts[clone_counts > 3])
  if (length(candidate_clones) == 0) candidate_clones <- names(clone_counts)
  gt_label_str <- sample(candidate_clones, 1)
  gt_cells <- rownames(seurat_obj@meta.data)[which(as.character(seurat_obj@meta.data[[GLOBAL_CLONE_COL]]) == gt_label_str)] 
  valid_gt_cells <- intersect(gt_cells, valid_umap_cells)
  seed_cell <- sample(valid_gt_cells, 1)
}

message(sprintf("  🌟 Selected clone ID '%s'; seed cell: %s", gt_label_str, seed_cell))
seurat_obj@misc$last_umap_seed_info <- list(clone_id = gt_label_str, seed_cell = seed_cell)

# -----------------------------------------------------------------
# 2. Extract the seed coordinates and draw the ground-truth panel
# -----------------------------------------------------------------
umap_embed <- Embeddings(seurat_obj, UMAP_REDUCTION)
seed_df <- as.data.frame(umap_embed[seed_cell, , drop=FALSE])
colnames(seed_df) <- c("UMAP_1", "UMAP_2") 

seurat_obj$GT_Label <- factor(ifelse(rownames(seurat_obj[[]]) %in% gt_cells, "Lineage Member", "Others"), levels = c("Others", "Lineage Member"))
p_gt_umap <- DimPlot(seurat_obj, reduction = UMAP_REDUCTION, group.by = "GT_Label", cols = c("Others"="#E0E0E0", "Lineage Member"="#94C83D"), order = "Lineage Member") + 
  # Mark the seed cell with the same gold diamond used in the other panels
  geom_point(data = seed_df, aes(x = UMAP_1, y = UMAP_2), shape = 23, size = 3.5, fill = "#FFCC00", color = "black", stroke = 1) + 
  theme_void() + theme(plot.title = element_text(size=14, hjust=0.5, family="Arial", face="bold")) + labs(title = paste("Ground Truth (ID:", gt_label_str, ")"))

# -----------------------------------------------------------------
# 3. 🎯 Core renderer based on fold enrichment and log1p visual enhancement
# -----------------------------------------------------------------
plot_umap_enrichment_seurat <- function(mat, title, feature_name) { 
  if(is.null(mat)) return(ggplot() + theme_void() + labs(title=title))
  if (all(grepl("^[0-9]+$", head(rownames(mat), 50))) && nrow(mat) == ncol(seurat_obj)) { 
    rownames(mat) <- colnames(seurat_obj); colnames(mat) <- colnames(seurat_obj) 
  }
  if(!(seed_cell %in% rownames(mat))) return(ggplot() + theme_void() + labs(title=paste(title, "(Cell Excluded)")))
  
  # --- a. Strict Markov row normalization ---
  P_raw <- mat
  if (!inherits(P_raw, "CsparseMatrix")) P_raw <- as(P_raw, "CsparseMatrix")
  
  row_sums <- Matrix::rowSums(P_raw)
  row_sums[row_sums == 0] <- 1  
  
  P_norm <- Matrix::Diagonal(x = 1 / row_sums) %*% P_raw
  rownames(P_norm) <- rownames(P_raw); colnames(P_norm) <- colnames(P_raw)
  
  # --- b. Convert probabilities to a shared fold-enrichment scale ---
  N_cells <- ncol(P_norm)
  E_matrix <- P_norm * (N_cells - 1) 
  
  # --- c. Extract enrichment propagated from the seed cell ---
  enrich_vector <- E_matrix[seed_cell, ]
  enrich_vector[is.na(enrich_vector)] <- 0
  
  # Remove the dominant self-retention probability so it does not mask the propagated gradient
  if (seed_cell %in% names(enrich_vector)) enrich_vector[seed_cell] <- NA 
  
  seurat_obj[[feature_name]] <- enrich_vector[match(rownames(seurat_obj[[]]), colnames(E_matrix))]
  
  # --- d. Draw the field on the absolute scale with log1p correction for diffuse decay ---
  p <- FeaturePlot(seurat_obj, reduction = UMAP_REDUCTION, features = feature_name, order = TRUE) + 
    scale_color_gradientn(
      name = "Enrichment\n(Fold Change)", 
      colors = c("#E0E0E0", "#9ECAE1", "#4292C6", "#08519C", "#08306B"), 
      na.value = "#E0E0E0", 
      limits = c(0, UNIVERSAL_FC_MAX),  # 🔐 Fix the upper bound of the expanded absolute scale
      oob = scales::squish,             # Clamp connections above the upper bound to the highest color
      
      # 🎯 Use log mapping to retain diffuse signals
      trans = "log1p",
      
      # 🎯 Reposition log-scale breaks for the expanded 500-fold legend range
      breaks = c(0, 10, 50, 200, UNIVERSAL_FC_MAX)
    ) + 
    # Show the seed explicitly as a gold coordinate marker rather than encoding it as field data
    geom_point(data = seed_df, aes(x = UMAP_1, y = UMAP_2), shape = 23, size = 3.5, fill = "#FFCC00", color = "black", stroke = 1) + 
    theme_void() + 
    theme(
      plot.title = element_text(size=14, hjust=0.5, family="Arial", face="bold"), 
      legend.title = element_text(family="Arial", size=10, face="bold"), 
      legend.text = element_text(family="Arial", size=9)
    ) + 
    labs(title = title)
  
  return(p)
}

# -----------------------------------------------------------------
# 4. Generate, combine, and export the panels
# -----------------------------------------------------------------
# Load scales to support squish-based clamping
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

p_iscoral_umap <- plot_umap_enrichment_seurat(S_stab_mat_fresh, "isCORAL", "fc_iscoral")
p_gemli_umap   <- plot_umap_enrichment_seurat(baseline_raw_matrices[["GEMLI"]], "GEMLI", "fc_gemli")
p_snn_cc_umap  <- plot_umap_enrichment_seurat(baseline_raw_matrices[["CC-Regressed SNN"]], "CC-Regressed SNN", "fc_snn_cc")

p_all_umaps <- wrap_plots(p_gt_umap, p_iscoral_umap, p_gemli_umap, p_snn_cc_umap, ncol = 2) + 
  plot_annotation(
    title = "Markov Transition Enrichment Field (Log Scale)", 
    theme = theme(plot.title = element_text(size = 18, hjust = 0.5, family = "Arial", face = "bold"))
  )

print(p_all_umaps)

# 👉 Export using the shared template-aware scaling function
out_file_benchmark <- file.path(export_dir, sprintf("Fig.Benchmark_UMAP_%s.pdf", SAMPLE_ID))
dims_8 <- get_template_dims(TEMPLATE_DIR, "Fig.Benchmark_UMAP.pdf") 
smart_ggsave_pdf(filename = out_file_benchmark, plot = p_all_umaps, dims = dims_8, default_w = 12, default_h = 10, dpi = 300)
message(sprintf("   ✅ Benchmark UMAP panel with an absolute scale exported to:\n      %s", out_file_benchmark))

# ==============================================================================
# 📊 Render 9: Matrix-level core metrics with native split axes and absolute grid alignment
# ==============================================================================
message("\n  -> Render 9: Matrix-level core metrics with native split axes and cross-panel alignment...")

alg_levels_strict <- c(
  "isCORAL", "GEMLI", "Standard SNN", "CC-Regressed SNN", 
  "Standard KNN", "CC-Regressed KNN", "Theoretical Best", "Uniform (Worst)"
)

mat_group_mapping <- c(
  "Standard SNN" = "SNN baseline", 
  "CC-Regressed SNN" = "SNN baseline", 
  "Standard KNN" = "SNN baseline",
  "CC-Regressed KNN" = "SNN baseline",
  "GEMLI" = "competitors", 
  "isCORAL" = "isCORAL",
  "Theoretical Best" = "Boundaries", 
  "Uniform (Worst)" = "Boundaries"   
)
mat_group_ordered_levels <- c("isCORAL", "competitors", "SNN baseline", "Boundaries")

method_colors_extended <- c(
  method_colors, 
  "Theoretical Best" = "#81C784", 
  "Uniform (Worst)" = "#B0BEC5"
)

plot_auprc_df <- b1_cache$matrix_summary %>% 
  filter(!Algorithm %in% c("Theoretical Best (Bootstrap Floor)", "Uniform Allocation (Worst Case)")) %>%
  select(Algorithm, Score = AUPRC_Enrichment) %>% 
  mutate(Algorithm = dplyr::case_when(
    grepl("isCORAL", Algorithm) ~ "isCORAL", grepl("GEMLI", Algorithm) ~ "GEMLI", 
    grepl("Standard SNN", Algorithm) ~ "Standard SNN", grepl("CC-Regressed SNN", Algorithm) ~ "CC-Regressed SNN", 
    grepl("Standard KNN", Algorithm) ~ "Standard KNN", grepl("CC-Regressed KNN", Algorithm) ~ "CC-Regressed KNN", 
    TRUE ~ as.character(Algorithm)
  ))

current_base_rate <- mean(b1_cache$matrix_summary$Random_AUPRC, na.rm = TRUE)
boundaries_df <- data.frame(Algorithm = c("Theoretical Best", "Uniform (Worst)"), Score = c(1.0 / current_base_rate, 1.0))

plot_auprc_df <- bind_rows(plot_auprc_df, boundaries_df) %>% 
  mutate(
    Algorithm = factor(Algorithm, levels = rev(alg_levels_strict)),
    Algorithm_Group = factor(mat_group_mapping[as.character(Algorithm)], levels = mat_group_ordered_levels)
  )

clean_raw_dist_data <- function(raw_df, score_col) {
  raw_df %>%
    mutate(Algorithm = dplyr::case_when(
      grepl("isCORAL", Method) ~ "isCORAL", grepl("GEMLI", Method) ~ "GEMLI", 
      grepl("Standard SNN", Method) ~ "Standard SNN", grepl("CC-Regressed SNN", Method) ~ "CC-Regressed SNN", 
      grepl("Standard KNN", Method) ~ "Standard KNN", grepl("CC-Regressed KNN", Method) ~ "CC-Regressed KNN", 
      grepl("Theoretical Best", Method) ~ "Theoretical Best", grepl("Uniform Allocation", Method) ~ "Uniform (Worst)",
      TRUE ~ as.character(Method)
    )) %>%
    select(Algorithm, Score = !!sym(score_col)) %>%
    mutate(
      Algorithm = factor(Algorithm, levels = rev(alg_levels_strict)),
      Algorithm_Group = factor(mat_group_mapping[as.character(Algorithm)], levels = mat_group_ordered_levels)
    )
}

plot_ed_dist_df <- clean_raw_dist_data(b1_cache$df_all_ed_raw, "Energy_Distance")
plot_pb_dist_df <- clean_raw_dist_data(b1_cache$df_all_pb_raw, "PB_L1_Distance")

max_val <- max(plot_auprc_df$Score, na.rm = TRUE)
second_max <- max(plot_auprc_df$Score[plot_auprc_df$Score < max_val], na.rm = TRUE)

p_base_bar <- ggplot(plot_auprc_df, aes(y = Algorithm, x = Score, fill = Algorithm)) + 
  geom_col(aes(linewidth = ifelse(Algorithm == "isCORAL", 0.8, 0.4)), width = 0.75, color = "black") + 
  scale_linewidth_identity() +
  geom_text(aes(label = sprintf(ifelse(Score > 10, "%.1f", "%.2f"), Score)), 
            hjust = -0.15, size = 4.5, family = "Arial", color = "black") + 
  scale_fill_manual(values = method_colors_extended) + 
  scale_y_discrete(expand = expansion(add = c(0.7, 0.7))) + 
  facet_grid(Algorithm_Group ~ ., scales = "free_y", space = "free_y", switch = "y")

t_base <- theme_bw(base_size = 14, base_family = "Arial") + 
  theme(
    legend.position = "none", panel.border = element_rect(color = "#333333", linewidth = 0.8), 
    panel.grid = element_blank(), axis.text.x = element_text(color = "black", size = 11, family = "Arial"), 
    axis.ticks.y = element_blank(), strip.placement = "outside", 
    strip.background.y = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.8), 
    strip.text.y.left = element_text(size = 13, angle = 0, family = "Arial")
  )

p_enrich_left <- p_base_bar + 
  coord_cartesian(xlim = c(0, second_max * 1.15)) +  
  scale_x_continuous(breaks = c(0, 5, 10)) + 
  labs(title = "AUPRC Enrichment\n(Higher is Better)", x = NULL, y = NULL) + 
  t_base + 
  theme(
    plot.margin = margin(t = 10, r = 2, b = 10, l = 5), 
    plot.title = element_text(size = 15, hjust = 0.5, family = "Arial"),
    axis.text.y = element_text(color = "black", size = 13, family = "Arial")
  )

p_enrich_right <- p_base_bar + 
  coord_cartesian(xlim = c(max_val * 0.92, max_val * 1.15)) + 
  scale_x_continuous(breaks = c(round(max_val, 1))) +  
  labs(title = " \n ", x = NULL, y = NULL) + 
  t_base + 
  theme(
    plot.margin = margin(t = 10, r = 5, b = 10, l = 2), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_blank(), 
    strip.text.y.left = element_blank(), 
    strip.background.y = element_blank()
  )

make_raincloud_plot <- function(data, clean_title, show_y_text = FALSE) {
  title_suffix <- "\n(Lower is Better)"
  df_mean <- data %>% group_by(Algorithm, Algorithm_Group) %>% summarise(mean_val = mean(Score, na.rm=TRUE), .groups="drop")
  
  p <- ggplot(data, aes(y = Algorithm, x = Score, fill = Algorithm)) + 
    ggdist::stat_halfeye(adjust = 0.6, width = 0.55, .width = 0, justification = -0.15, point_colour = NA, color = "black", alpha = 0.85) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5, color = "black", linewidth = 0.6) +
    geom_jitter(aes(color = Algorithm), width = 0, height = 0.08, size = 0.6, alpha = 0.25) +
    geom_point(data = df_mean, aes(x = mean_val, y = Algorithm), color = "black", shape = 23, fill = "white", size = 2.5, stroke = 1, inherit.aes = FALSE) +
    geom_text(data = df_mean, aes(x = mean_val, y = Algorithm, label = sprintf("%.3f", mean_val)), 
              nudge_y = -0.35, size = 4.2, family = "Arial", color = "black", inherit.aes = FALSE) +
    scale_fill_manual(values = method_colors_extended) + scale_color_manual(values = method_colors_extended) + 
    scale_y_discrete(expand = expansion(add = c(0.7, 0.7))) + 
    facet_grid(Algorithm_Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = paste0(clean_title, title_suffix), x = NULL, y = NULL) + 
    theme_bw(base_size = 14, base_family = "Arial") + 
    theme(
      legend.position = "none", plot.title = element_text(size = 15, hjust = 0.5, family = "Arial"), 
      panel.border = element_rect(color = "#333333", linewidth = 0.8), panel.grid.major.x = element_line(color = "#E5E5E5", linetype = "dashed"), 
      panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(color = "black", size = 11, family = "Arial"), 
      axis.ticks.y = element_blank(), strip.placement = "outside", strip.background.y = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.8), 
      strip.text.y.left = element_text(size = 13, angle = 0, family = "Arial"), 
      plot.margin = margin(t = 10, r = 10, b = 10, l = 5)
    )
  
  if (show_y_text) p <- p + theme(axis.text.y = element_text(color = "black", size = 13, family = "Arial")) 
  else p <- p + theme(axis.text.y = element_blank(), strip.text.y.left = element_blank(), strip.background.y = element_blank()) 
  
  return(p)
}

p_mat_ed   <- make_raincloud_plot(plot_ed_dist_df, "Energy Distance", FALSE)
p_mat_pbl1 <- make_raincloud_plot(plot_pb_dist_df, "Pseudobulk L1", FALSE)

p_matrix_metrics_panel <- p_enrich_left + p_enrich_right + p_mat_ed + p_mat_pbl1 + 
  plot_layout(nrow = 1, widths = c(1.3, 0.4, 1.4, 1.4))

print(p_matrix_metrics_panel)

# 👉 Export using the sample-independent template name
out_file_mat_barplots_9 <- file.path(export_dir, sprintf("Fig.Matrix_Metrics_Rainclouds_%s.pdf", SAMPLE_ID))
dims_9 <- get_template_dims(TEMPLATE_DIR, "Fig.Matrix_Metrics_Rainclouds.pdf")
smart_ggsave_pdf(filename = out_file_mat_barplots_9, plot = p_matrix_metrics_panel, dims = dims_9, default_w = 18, default_h = 6, dpi = 300)
message(sprintf("   ✅ Matrix-level core metric panels exported to:\n      %s", out_file_mat_barplots_9))


# ==============================================================================
# 📊 Render 10: Simplified matrix-level core metrics for isCORAL, GEMLI, and CC-Regressed SNN
# ==============================================================================
message("\n  -> Render 10: Simplified matrix-level core metrics...")

alg_levels_strict <- c("isCORAL", "GEMLI", "CC-Regressed SNN")
mat_group_mapping <- c("CC-Regressed SNN" = "SNN baseline", "GEMLI" = "competitors", "isCORAL" = "isCORAL")
mat_group_ordered_levels <- c("isCORAL", "competitors", "SNN baseline")

plot_auprc_df <- b1_cache$matrix_summary %>% 
  filter(Algorithm %in% c("isCORAL", "GEMLI", "CC-Regressed SNN", "isCORAL (Probability Field)", "GEMLI (Prediction Score)")) %>%
  select(Algorithm, Score = AUPRC_Enrichment) %>% 
  mutate(Algorithm = dplyr::case_when(
    grepl("isCORAL", Algorithm) ~ "isCORAL", grepl("GEMLI", Algorithm) ~ "GEMLI", 
    grepl("CC-Regressed SNN", Algorithm) ~ "CC-Regressed SNN", TRUE ~ as.character(Algorithm)
  )) %>%
  mutate(
    Algorithm = factor(Algorithm, levels = rev(alg_levels_strict)),
    Algorithm_Group = factor(mat_group_mapping[as.character(Algorithm)], levels = mat_group_ordered_levels)
  )

clean_raw_dist_data <- function(raw_df, score_col) {
  raw_df %>%
    filter(grepl("isCORAL|GEMLI|CC-Regressed SNN", Method)) %>%
    mutate(Algorithm = dplyr::case_when(
      grepl("isCORAL", Method) ~ "isCORAL", grepl("GEMLI", Method) ~ "GEMLI", 
      grepl("CC-Regressed SNN", Method) ~ "CC-Regressed SNN", TRUE ~ as.character(Method)
    )) %>%
    select(Algorithm, Score = !!sym(score_col)) %>%
    mutate(
      Algorithm = factor(Algorithm, levels = rev(alg_levels_strict)),
      Algorithm_Group = factor(mat_group_mapping[as.character(Algorithm)], levels = mat_group_ordered_levels)
    )
}

plot_ed_dist_df <- clean_raw_dist_data(b1_cache$df_all_ed_raw, "Energy_Distance")
plot_pb_dist_df <- clean_raw_dist_data(b1_cache$df_all_pb_raw, "PB_L1_Distance")

p_mat_enrich <- ggplot(plot_auprc_df, aes(y = Algorithm, x = Score, fill = Algorithm)) + 
  geom_col(aes(linewidth = ifelse(Algorithm == "isCORAL", 0.8, 0.4)), width = 0.75, color = "black") + 
  scale_linewidth_identity() +
  geom_text(aes(label = sprintf(ifelse(Score > 10, "%.1f", "%.2f"), Score)), 
            hjust = -0.15, size = 4.5, family = "Arial", color = "black") + 
  scale_fill_manual(values = method_colors_extended) + 
  scale_y_discrete(expand = expansion(add = c(0.7, 0.7))) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) + 
  facet_grid(Algorithm_Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "AUPRC Enrichment\n(Higher is Better)", x = NULL, y = NULL) + 
  theme_bw(base_size = 14, base_family = "Arial") + 
  theme(
    legend.position = "none", plot.title = element_text(size = 15, hjust = 0.5, family = "Arial"), 
    panel.border = element_rect(color = "#333333", linewidth = 0.8), panel.grid = element_blank(), 
    axis.text.x = element_text(color = "black", size = 11, family = "Arial"), 
    axis.text.y = element_text(color = "black", size = 13, family = "Arial"),
    axis.ticks.y = element_blank(), strip.placement = "outside", 
    strip.background.y = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.8), 
    strip.text.y.left = element_text(size = 13, angle = 0, family = "Arial"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 5)
  )

p_mat_ed   <- make_raincloud_plot(plot_ed_dist_df, "Energy Distance", FALSE)
p_mat_pbl1 <- make_raincloud_plot(plot_pb_dist_df, "Pseudobulk L1", FALSE)

p_matrix_metrics_panel_simp <- p_mat_enrich + p_mat_ed + p_mat_pbl1 + 
  plot_layout(nrow = 1, widths = c(1.4, 1.4, 1.4))

print(p_matrix_metrics_panel_simp)

# 👉 Export using the sample-independent template name
out_file_mat_barplots_10 <- file.path(export_dir, sprintf("Fig.Matrix_Metrics_Rainclouds_Simplified_%s.pdf", SAMPLE_ID))
dims_10 <- get_template_dims(TEMPLATE_DIR, "Fig.Matrix_Metrics_Rainclouds_Simplified.pdf")
smart_ggsave_pdf(filename = out_file_mat_barplots_10, plot = p_matrix_metrics_panel_simp, dims = dims_10, default_w = 16, default_h = 6, dpi = 300)
message(sprintf("   ✅ Simplified matrix-level core metric panels exported to:\n      %s", out_file_mat_barplots_10))

message("\n🎉 B3 visualization rendering completed.")
