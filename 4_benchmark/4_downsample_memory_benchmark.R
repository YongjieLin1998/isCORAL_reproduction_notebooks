# ==============================================================================
# 🧮 Analysis Script 04: isCORAL Downsampling Benchmark Evaluation 
#   (Final Version: Absolute Scale + Defensive + Spearman + High-Contrast Plots)
# ==============================================================================

message("--- Loading Dependencies ---")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(Matrix)
  library(PRROC)
  library(pROC)
})

library(extrafont)
# Load fonts for the PDF device
loadfonts(device = "pdf")

# -----------------------------------------------------------------
# 🟢 [0] Unified Global Configuration Hub with strict input alignment
# -----------------------------------------------------------------
# Set ISCORAL_PROJECT_ROOT, ISCORAL_SAMPLE_ID, and ISCORAL_DOWNSAMPLE_RDS as needed.
BASE_DIR      <- normalizePath(Sys.getenv("ISCORAL_PROJECT_ROOT", unset = "."), mustWork = TRUE)
SAMPLE_ID     <- Sys.getenv("ISCORAL_SAMPLE_ID", unset = "<SAMPLE_ID>")
if (identical(SAMPLE_ID, "<SAMPLE_ID>")) {
  stop("Set ISCORAL_SAMPLE_ID or replace <SAMPLE_ID> before running this template.")
}

# Path to the downsampling result file containing all 18 runs
DOWNSAMPLE_RDS <- Sys.getenv(
  "ISCORAL_DOWNSAMPLE_RDS",
  unset = file.path(BASE_DIR, "3_results_new", "downsampling_analysis", "<DOWNSAMPLING_RESULTS_FILE>.rds")
)

# Core Metadata and System-Level Filtering
GLOBAL_CLONE_COL      <- "barcode" 
GLOBAL_INVALID_CLONES <- c("0","NA", "NaN", "", "None", "Unknown", "unknown", "Unassigned")

# Gene Engine and Display Parameters
TOP_N_MEMORY_GENES   <- 2000
FILTER_MEAN_GT       <- 0.1

ISCORAL_BIN_SIZE     <- 200  
ISCORAL_BINSIZE_NAME <- paste0("binsize_", ISCORAL_BIN_SIZE)

# -----------------------------------------------------------------
# --- STEP 1: Auto-Locate Cache & Load Results
# -----------------------------------------------------------------
message("--- Step 1: Loading Dependencies & Auto-Locating Latest Baseline Cache ---")

BASELINE_DIR <- file.path(BASE_DIR, "4_benchmark/baseline")
all_baseline_files <- list.files(
  path = BASELINE_DIR, 
  pattern = sprintf("^Baseline_Cache_%s_.*\\.rds$", SAMPLE_ID), 
  full.names = TRUE
)

if (length(all_baseline_files) == 0) {
  stop(sprintf("❌ No Baseline Cache found for SAMPLE_ID '%s'! Please check directory.", SAMPLE_ID))
}

file_info <- file.info(all_baseline_files)
CACHE_FILE <- rownames(file_info)[which.max(file_info$mtime)]

message(sprintf("  🚀 Auto-loading the latest baseline cache:\n     %s", CACHE_FILE))
cache <- readRDS(CACHE_FILE)
seurat_obj <- cache$seurat_obj

message("  🚀 Loading downsampling results grid...")
if (!file.exists(DOWNSAMPLE_RDS)) stop(sprintf("❌ Downsampling result file not found: %s", DOWNSAMPLE_RDS))
ds_master <- readRDS(DOWNSAMPLE_RDS)

# -----------------------------------------------------------------
# --- STEP 2: Extract Ground Truth & Calculate Baseline Specs
# -----------------------------------------------------------------
message("--- Step 2: Extracting Ground Truth & Baseline Metrics ---")

safe_get_matrix <- function(obj, assay_name, layer_name, force_dense = FALSE) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") mat <- LayerData(obj, assay = assay_name, layer = layer_name) 
  else mat <- GetAssayData(obj, assay = assay_name, slot = layer_name)
  if (force_dense) return(as.matrix(mat)) else return(mat)
}

exp_matrix_full <- safe_get_matrix(seurat_obj, DefaultAssay(seurat_obj), "data")
universal_gene_means <- data.frame(name = rownames(exp_matrix_full), mean_expr = Matrix::rowMeans(exp_matrix_full))

if (is.null(seurat_obj@misc$groundtruth$heritable_genes_df)) {
  stop("❌ CRITICAL ERROR: The seurat_obj in the baseline cache does not contain ground-truth data.")
}

gt_df <- seurat_obj@misc$groundtruth$heritable_genes_df %>% 
  dplyr::select(name, ground_truth_score = Omega_square) %>% 
  dplyr::filter(!is.na(ground_truth_score)) %>% 
  dplyr::inner_join(universal_gene_means, by = "name") %>% 
  dplyr::filter(mean_expr > FILTER_MEAN_GT) %>% 
  dplyr::arrange(desc(ground_truth_score)) %>% 
  dplyr::mutate(is_high_truth = dplyr::row_number() <= TOP_N_MEMORY_GENES)

# Calculate absolute baseline sizes
total_cells_raw <- ncol(seurat_obj)
total_umis_raw  <- sum(seurat_obj$nCount_RNA) # Assumes the RNA assay provides nCount_RNA
umis_per_cell_raw <- total_umis_raw / total_cells_raw

message(sprintf("  >>> Baseline size - total cells: %d | mean UMIs per cell: %.1f", total_cells_raw, umis_per_cell_raw))

# -----------------------------------------------------------------
# --- STEP 3: Downsampling Metric Calculation (Defensive Engine)
# -----------------------------------------------------------------
message("--- Step 3: Evaluating Downsampling Runs against Ground Truth ---")

calc_pr <- function(labels, scores) { P <- sum(labels, na.rm=T); N <- sum(!labels, na.rm=T); if(P==0 || N==0) return(list(df=data.frame(Recall=c(0,1), Precision=c(1,0), Threshold=c(Inf,-Inf)), auprc=0, max_f1=0, best_thresh=0)); df <- data.frame(score = scores, label = labels)[order(scores, decreasing = TRUE, na.last = NA), ]; TP <- 0; FP <- 0; pts <- data.frame(Recall=0, Precision=1, Threshold=Inf); for(i in 1:nrow(df)) { if(df$label[i]) TP <- TP+1 else FP <- FP+1; if(i < nrow(df) && df$score[i] == df$score[i+1]) next; pts <- rbind(pts, c(TP/P, ifelse(TP+FP==0,0,TP/(TP+FP)), df$score[i])) }; if(TP!=P || FP!=N) pts <- rbind(pts, c(1, P/(P+N), -Inf)); pr_auc_df <- pts[!duplicated(pts$Recall), ]; auprc <- sum(diff(pr_auc_df$Recall) * (pr_auc_df$Precision[-1] + pr_auc_df$Precision[-nrow(pr_auc_df)])/2, na.rm=T); return(list(auprc = auprc)) }

eval_list <- list()
fractions_to_test <- c(0.2, 0.5, 1.0)
modes_to_test <- c("Cell", "UMI")

for (mode in modes_to_test) {
  for (frac in fractions_to_test) {
    for (rep in 1:3) {
      
      run_id <- paste0(mode, "_", frac, "_rep", rep)
      test_res <- ds_master$results[[run_id]]
      
      if (is.null(test_res)) {
        message(sprintf("⚠️ Skipping [%s]: result data not found.", run_id))
        next
      }
      
      iscoral_raw <- test_res[[ISCORAL_BINSIZE_NAME]]$ss_filter_average
      
      if (is.null(iscoral_raw)) {
        message(sprintf("⚠️ Skipping [%s]: score table for %s not found.", run_id, ISCORAL_BINSIZE_NAME))
        next
      }
      
      score_col_name <- if("omega_square_infer" %in% colnames(iscoral_raw)) "omega_square_infer" else "score"
      
      gene_collated <- iscoral_raw %>% 
        dplyr::rename(score = !!sym(score_col_name)) %>%
        dplyr::select(name, score) %>% 
        dplyr::inner_join(gt_df, by="name") %>% 
        dplyr::filter(!is.na(is_high_truth), is.finite(ground_truth_score), is.finite(score))
      
      if(nrow(gene_collated) < 10) {
        message(sprintf("⚠️ Skipping [%s]: fewer than 10 shared genes were matched.", run_id))
        next
      }
      
      p_res <- calc_pr(gene_collated$is_high_truth, gene_collated$score)
      # 🌟 Use Spearman correlation
      spearman_r <- cor(gene_collated$ground_truth_score, gene_collated$score, method="spearman")
      
      eval_list[[run_id]] <- data.frame(
        Mode = mode, Fraction = as.factor(frac), Repeat = rep,
        AUPRC = p_res$auprc, Spearman_R = spearman_r # 🌟 Spearman correlation
      )
    }
  }
}

df_eval <- dplyr::bind_rows(eval_list)
if(nrow(df_eval) == 0) stop("❌ No valid data were collected because every run was skipped. Review the preceding errors.")

# -----------------------------------------------------------------
# --- STEP 4: Aggregation & Absolute X-Axis Mapping
# -----------------------------------------------------------------
message("--- Step 4: Data Aggregation & Coordinate Mapping ---")

df_summary <- df_eval %>%
  dplyr::group_by(Mode, Fraction) %>%
  dplyr::summarise(
    AUPRC_Mean = mean(AUPRC, na.rm=TRUE), AUPRC_SD = sd(AUPRC, na.rm=TRUE),
    Spearman_Mean = mean(Spearman_R, na.rm=TRUE), Spearman_SD = sd(Spearman_R, na.rm=TRUE),
    .groups = "drop"
  )

# 1. Map the summary table to absolute coordinates for trend lines and error bars
df_eval_plot <- df_summary %>%
  mutate(Fraction_num = as.numeric(as.character(Fraction))) %>%
  mutate(Absolute_X = case_when(
    Mode == "Cell" ~ Fraction_num * total_cells_raw,
    Mode == "UMI"  ~ Fraction_num * umis_per_cell_raw
  ))

# 2. Map the raw table to absolute coordinates for individual data points
df_raw_plot <- df_eval %>%
  mutate(Fraction_num = as.numeric(as.character(Fraction))) %>%
  mutate(Absolute_X = case_when(
    Mode == "Cell" ~ Fraction_num * total_cells_raw,
    Mode == "UMI"  ~ Fraction_num * umis_per_cell_raw
  ))

message("\n★★★ Downsampling Absolute Benchmark Scorecard ★★★")
print(knitr::kable(df_eval_plot, digits = 3))

export_dir = file.path(BASE_DIR, "4_benchmark/results/downsampling_reports")
if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

write.csv(df_summary, file.path(export_dir, sprintf("isCORAL_Downsampling_Absolute_Summary_%s.csv", SAMPLE_ID)), row.names = FALSE)
write.csv(df_eval, file.path(export_dir, sprintf("isCORAL_Downsampling_Raw_Metrics_%s.csv", SAMPLE_ID)), row.names = FALSE)


# -----------------------------------------------------------------
# --- STEP 5: Generate the 2x2 benchmark figure
# -----------------------------------------------------------------
message("--- Step 5: Generating Clean Academic Plots (Fixing Overlaps & Layers) ---")

# Publication theme with adjusted font sizes and spacing to prevent overlap
theme_clean_academic <- theme_bw() + 
  theme(
    text = element_text(size = 16),                                        
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)), # Add spacing below the title
    axis.text = element_text(size = 14, color = "black"),  
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.8, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),                                  
    legend.position = "none"
  )

create_panel_clean <- function(summary_data, raw_data, mode_filter, metric_col, sd_col, title_text) {
  
  p_summary <- summary_data %>% filter(Mode == mode_filter)
  p_raw     <- raw_data %>% filter(Mode == mode_filter)
  
  raw_metric_col <- if(startsWith(metric_col, "AUPRC")) "AUPRC" else "Spearman_R"
  main_color <- ifelse(mode_filter == "Cell", "#E64B35", "#4DBBD5")
  
  ggplot() +
    # Layer 1: Trend line connecting group means
    geom_line(data = p_summary, aes(x = Absolute_X, y = !!sym(metric_col)), 
              color = main_color, size = 1.2) +
    
    # Layer 2: Standard-deviation error bars
    geom_errorbar(data = p_summary, 
                  aes(x = Absolute_X, 
                      ymin = !!sym(metric_col) - !!sym(sd_col), 
                      ymax = !!sym(metric_col) + !!sym(sd_col)), 
                  width = (max(p_raw$Absolute_X) - min(p_raw$Absolute_X)) * 0.04, 
                  size = 0.9, color = "#444444", alpha = 0.7) +
    
    # Layer 3: Individual raw data points
    # Keep x jitter at zero to preserve exact continuous coordinates
    # Apply only minimal y jitter to make overlapping points distinguishable
    geom_jitter(data = p_raw, aes(x = Absolute_X, y = !!sym(raw_metric_col)),
                shape = 21, size = 3.5, fill = main_color, color = "white", stroke = 0.5,
                alpha = 0.6, width = 0, height = 0.005) +
    
    # Layer 4: Mean anchors shown as black crosses
    geom_point(data = p_summary, aes(x = Absolute_X, y = !!sym(metric_col)), 
               color = "#222222", size = 4, shape = 3, stroke = 1.5) + 
    # Layer 5: Compact mean-value labels
    # nudge_y = 0.04 places each label slightly above its point
    # vjust = 0 bottom-aligns labels for stable placement
    geom_label(data = p_summary, 
               aes(x = Absolute_X, y = !!sym(metric_col), 
                   label = sprintf("%.3f", !!sym(metric_col))),
               color = "#222222", size = 4.0, fontface = "bold", # Bold labels improve readability
               fill = "white", alpha = 0.85,                  # Semi-transparent background
               label.size = 0.1, label.padding = unit(0.12, "lines"),
               nudge_y = 0.04, vjust = 0) +
    
    # Retain a strict 0-1 metric scale while leaving room for the nudged labels
    coord_cartesian(ylim = c(0, 1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    coord_cartesian(ylim = c(0, 1.1)) + # Extend the y limit to prevent clipping of top labels
    scale_y_continuous(breaks = seq(0, 1, 0.5)) + 
    labs(title = title_text, x = NULL, y = NULL) +
    theme_clean_academic
}

# Draw the four panels
p1 <- create_panel_clean(df_eval_plot, df_raw_plot, "Cell", "AUPRC_Mean", "AUPRC_SD", "AUPRC vs Absolute Cell Number")
p2 <- create_panel_clean(df_eval_plot, df_raw_plot, "UMI", "AUPRC_Mean", "AUPRC_SD", "AUPRC vs UMIs per Cell")
p3 <- create_panel_clean(df_eval_plot, df_raw_plot, "Cell", "Spearman_Mean", "Spearman_SD", "Spearman R vs Absolute Cell Number")
p4 <- create_panel_clean(df_eval_plot, df_raw_plot, "UMI", "Spearman_Mean", "Spearman_SD", "Spearman R vs UMIs per Cell")

# Combine panels with additional vertical spacing and margins
final_grid <- (p1 | p2) / (p3 | p4) + 
  plot_annotation(
    title = sprintf("isCORAL Stress Test (%s)", SAMPLE_ID),
    subtitle = sprintf("Baseline: %d Cells | %.1f UMIs per Cell\n(Dark dots represent 3 independent replicates)", total_cells_raw, umis_per_cell_raw),
    theme = theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(t = 15, b = 10)),    
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30", lineheight = 1.2),
      plot.margin = margin(15, 15, 15, 15) # Add outer whitespace around the combined figure
    )
  )

print(final_grid)

# 💾 Export the figure as a 15 x 12 inch PDF
out_file_benchmark <- file.path(export_dir, sprintf("Fig.Benchmark_Downsampling_Fixed_%s.pdf", SAMPLE_ID))
ggsave(filename = out_file_benchmark, plot = final_grid, device = cairo_pdf)
message("   ✅ The PDF panel with individual points was exported successfully.")
