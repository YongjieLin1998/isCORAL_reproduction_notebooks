# ==============================================================================
# 🎨 Script B2: Plot-only module with automatic B1 cache discovery, PDF export, and source-level size detection
# ==============================================================================


library(extrafont)
# Load fonts for the PDF device
loadfonts(device = "pdf")

message("\n--- Starting the B2 visualization renderer with automatic cache discovery and PDF export ---")

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

# 🎯 Restore the required variables into the current environment
gene_collated           <- b1_cache$gene_collated
gene_summary_table      <- b1_cache$gene_summary_table
all_method_scatter_data <- b1_cache$all_method_scatter_data
gene_roc_df             <- b1_cache$gene_roc_df
gene_pr_df              <- b1_cache$gene_pr_df
df_all_ed_raw           <- b1_cache$df_all_ed_raw
matrix_roc_all          <- b1_cache$matrix_roc_all
matrix_pr_all           <- b1_cache$matrix_pr_all


# ==============================================================================
# 🛠️ [CORE FIX 1] Standardize baseline method names to prevent filtering losses caused by spelling or case differences
# ==============================================================================
message("  -> Aligning and standardizing evaluation method names...")

standardize_method_names <- function(df, col_name = "method") {
  if (!col_name %in% colnames(df)) return(df)
  df[[col_name]] <- dplyr::case_when(
    grepl("(?i)seurat.*mvp", df[[col_name]]) ~ "Seurat-MVP",      # Match MVP first
    grepl("(?i)seurat", df[[col_name]]) ~ "Seurat-vst",           # Then match VST and other Seurat labels
    grepl("(?i)scran", df[[col_name]]) ~ "scran",
    grepl("(?i)standard.*snn", df[[col_name]]) ~ "Standard SNN",
    grepl("(?i)cc.*snn|regressed.*snn", df[[col_name]]) ~ "CC-Regressed SNN",
    grepl("(?i)powerseek", df[[col_name]]) ~ "PowerSeek",
    grepl("(?i)gemli", df[[col_name]]) ~ "GEMLI",
    grepl("(?i)iscoral", df[[col_name]]) ~ "isCORAL",
    TRUE ~ as.character(df[[col_name]])
  )
  return(df)
}

# Apply the naming patch across all evaluation tables
gene_summary_table      <- standardize_method_names(gene_summary_table, "method")
gene_roc_df             <- standardize_method_names(gene_roc_df, "method")
gene_pr_df              <- standardize_method_names(gene_pr_df, "method")
all_method_scatter_data <- standardize_method_names(all_method_scatter_data, "method")
matrix_roc_all          <- standardize_method_names(matrix_roc_all, "Method")
matrix_pr_all           <- standardize_method_names(matrix_pr_all, "Method")


# ==============================================================================
# 🎯 [CORE FIX 2] Attach baseline labels and isCORAL results to the base Seurat object dynamically
# ==============================================================================
message("  -> Attaching the base Seurat object and all algorithm results...")
baseline_cache <- readRDS(b1_cache$run_parameters$CACHE_FILE)
seurat_obj <- baseline_cache$seurat_obj

# 🚀 Attach cached baseline labels directly to metadata for downstream violin plots
if (!is.null(baseline_cache$baseline_raw_labels)) {
  if (!is.null(baseline_cache$baseline_raw_labels[["GEMLI"]])) {
    seurat_obj$gemli_clones_inferred <- baseline_cache$baseline_raw_labels[["GEMLI"]]
  }
  if (!is.null(baseline_cache$baseline_raw_labels[["CC-Regressed SNN"]])) {
    seurat_obj$snn_cc_clones <- baseline_cache$baseline_raw_labels[["CC-Regressed SNN"]]
  }
  message("  ✅ GEMLI and CC-SNN clone labels were added to Seurat metadata.")
}

SCAN_RES_FILE <- b1_cache$run_parameters$SCAN_RES_FILE
ISCORAL_BIN_SIZE <- b1_cache$run_parameters$ISCORAL_BIN_SIZE
ISCORAL_BINSIZE_NAME <- paste0("binsize_", ISCORAL_BIN_SIZE)

scanning_results <- readRDS(SCAN_RES_FILE)
seurat_obj@misc$CORAL_scan_results <- scanning_results

# Support the legacy list structure
if (!is.null(scanning_results$S_stab) && is.null(scanning_results[[ISCORAL_BINSIZE_NAME]])) {
  tmp_list <- list()
  tmp_list[[ISCORAL_BINSIZE_NAME]] <- scanning_results
  seurat_obj@misc$CORAL_scan_results <- tmp_list
}

# 🎯 Recover plotting thresholds from the cached data frames
TOP_N_MEMORY_GENES <- b1_cache$run_parameters$TOP_N_MEMORY_GENES
gt_thresh <- min(gene_collated$ground_truth_score[gene_collated$is_high_truth], na.rm = TRUE)
iscoral_best_thresh <- all_method_scatter_data %>% filter(method == "isCORAL") %>% pull(best_thresh) %>% .[1]


# -----------------------------------------------------------------
# 🟢 [1] Dynamic visualization parameters
# -----------------------------------------------------------------
GENE_TO_PLOT       <- b1_cache$run_parameters$GENE_TO_PLOT       
CHOSEN_CLONE_ID    <- "73"         
MIN_CLONE_SIZE     <- 2            
N_CLONES_TO_SAMPLE <- 12           
UMAP_REDUCTION     <- "umap"        


# ==============================================================================
# 📏 Automatic template-size detection and PDF export using raw PDF grep
# ==============================================================================
TEMPLATE_DIR <- Sys.getenv(
  "ISCORAL_PLOT_TEMPLATE_DIR",
  unset = file.path(BASE_DIR, "4_benchmark", "plothelper", "Benchmark_graph_template")
)

get_template_dims <- function(template_dir, file_name) {
  base_name <- tools::file_path_sans_ext(file_name)
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

smart_ggsave_pdf <- function(filename, plot_obj, dims, default_w, default_h) {
  if (!is.null(dims) && dims$width > 4) {
    ggsave(filename = filename, plot = plot_obj, device = cairo_pdf, width = dims$width, height = dims$height)
  } else {
    message(sprintf("  ⚠️ Template missing or unreadable; using safe default dimensions: %.1f x %.1f in", default_w, default_h))
    ggsave(filename = filename, plot = plot_obj, device = cairo_pdf, width = default_w, height = default_h)
  }
}

export_dir <- dirname(latest_b1_cache_path)

# -----------------------------------------------------------------
# Shared plotting theme and color palette
# -----------------------------------------------------------------
theme_benchmark_vis <- function(base_size = 14) { 
  theme_bw(base_size = base_size) + theme(plot.title = element_text(size = rel(1.1), hjust = 0.5), plot.subtitle = element_text(size = rel(1.0), hjust = 0.5), legend.position = "bottom", legend.box = "vertical", strip.background = element_rect(fill = "grey90"), strip.text = element_text()) 
}

method_colors_clean <- c(
  "isCORAL" = "#C74B46", "GEMLI" = "#3D5A80", "PowerSeek" = "#98C1D9",
  "CC-Regressed SNN" = "#2A9D8F", "Standard SNN" = "#A8DADC",
  "scran" = "#8E8E8E", "Seurat-vst" = "#D9D9D9", "Seurat-MVP" = "#BDBDBD"
)

nature_colors <- c("TP" = "#D37648", "TN" = "#E0E0E0", "FP" = "#91B35C", "FN" = "#4C7AA5")

# 🎯 Data filtering
plot_roc_df <- matrix_roc_all %>% filter(!grepl("KNN", Method))
plot_pr_df  <- matrix_pr_all %>% filter(!grepl("KNN", Method))
plot_ed_df  <- df_all_ed_raw %>% filter(!grepl("KNN", Method))


# ==============================================================================
# 🎨 Render 10: isCORAL-specific publication-style scatter plot
# ==============================================================================
message("\n  -> Render 10: isCORAL-specific publication-style scatter plot...")

df_iscoral <- all_method_scatter_data %>% filter(method == "isCORAL")
n_genes <- nrow(df_iscoral)

TP <- sum(df_iscoral$class == "TP"); TN <- sum(df_iscoral$class == "TN")
FP <- sum(df_iscoral$class == "FP"); FN <- sum(df_iscoral$class == "FN")
Acc  <- (TP + TN) / n_genes
Prec <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
Rec  <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
F1   <- ifelse((Prec + Rec) == 0, 0, 2 * Prec * Rec / (Prec + Rec))
Pear <- cor(df_iscoral$ground_truth_score, df_iscoral$score, method = "pearson")

top_genes <- df_iscoral %>% filter(class == "TP") %>% arrange(desc(score)) %>% head(5)
max_x <- max(df_iscoral$ground_truth_score, na.rm = TRUE)
max_y <- max(df_iscoral$score, na.rm = TRUE)

legend_df <- data.frame(
  x = c(max_x * 0.75, max_x * 0.90, max_x * 0.75, max_x * 0.90),
  y = c(max_y * 0.25, max_y * 0.25, max_y * 0.10, max_y * 0.10),
  label = c("FN", "FP", "TN", "TP"),
  hex_color = unname(nature_colors[c("FN", "FP", "TN", "TP")])
)

p_main <- ggplot(df_iscoral, aes(x = ground_truth_score, y = score)) +
  geom_abline(slope = 1, intercept = 0, color = "#D0D0D0", linetype = "dashed", linewidth = 1.5) +
  geom_smooth(method = "lm", color = "#C1463D", linewidth = 2.5, se = FALSE) +
  geom_point(aes(color = class), size = 2.5, alpha = 0.9) +
  scale_color_manual(values = nature_colors) +
  geom_vline(xintercept = gt_thresh, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_hline(yintercept = iscoral_best_thresh, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_text_repel(data = top_genes, aes(label = name), color = "#D37648", size = 6, family = "Arial", segment.color = "black", segment.linetype = "dashed", segment.size = 1.2, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf) +
  annotate("text", x = max_x * 0.05, y = max_y * 0.95, label = paste(n_genes, "genes"), size = 7, hjust = 0, family = "Arial") +
  geom_text(data = legend_df, aes(x = x, y = y, label = label), color = legend_df$hex_color, size = 6, family = "Arial") +
  labs(title = "real v.s. inferred omega squared", x = "Groundtruth Score", y = "isCORAL inference") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b = 15)), axis.text = element_text(color = "black", size = 16), axis.title = element_text(color = "black", size = 18), panel.border = element_rect(color = "#333333", linewidth = 2.5), panel.grid.major = element_line(color = "#F0F0F0", linewidth = 1.5), panel.grid.minor = element_blank(), plot.margin = margin(t=10, r=5, b=10, l=10))

metrics_text <- paste0("isCORAL memory\ngenes detection\n\n", "Pearson r = ", sprintf("%.2f", Pear), "\n", "Accuracy = ", sprintf("%.3f", Acc), "\n", "Precision = ", sprintf("%.3f", Prec), "\n", "Recall = ", sprintf("%.2f", Rec), "\n", "F1-score = ", sprintf("%.3f", F1))
p_metrics <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = metrics_text, size = 6.5, hjust = 0.5, vjust = 0.5, family = "Arial") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_void() + theme(panel.border = element_rect(color = "#C1463D", linewidth = 4, fill = NA), plot.margin = margin(t=40, r=10, b=40, l=10))

p_final_showcase <- p_main + p_metrics + plot_layout(widths = c(3, 1))

out_10 <- file.path(export_dir, "Fig.isCORALD.pdf")
dims_10 <- get_template_dims(TEMPLATE_DIR, "Fig.isCORALD.pdf")
smart_ggsave_pdf(filename = out_10, plot_obj = p_final_showcase, dims = dims_10, default_w = 12, default_h = 6)
message(sprintf("✅ Fig.isCORALD.pdf exported to:\n   👉 %s", out_10))


# ==============================================================================
# 🎨 Render 11: Publication-style overview scatter plots for baseline methods with thin lines
# ==============================================================================
message("\n  -> Render 11: Publication-style overview scatter plots for baseline methods...")

df_baselines <- all_method_scatter_data %>% filter(method != "isCORAL")

p_baselines <- ggplot(df_baselines, aes(x = ground_truth_score, y = score)) +
  geom_abline(slope = 1, intercept = 0, color = "#D0D0D0", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "lm", color = "#C1463D", linewidth = 0.8, se = FALSE) +
  geom_point(aes(color = class), size = 0.8, alpha = 0.7) +
  scale_color_manual(values = nature_colors) +
  geom_vline(xintercept = gt_thresh, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(aes(yintercept = best_thresh), linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~method, scales = "free_y", ncol = 3) +
  labs(title = "Baseline Methods: real v.s. inferred omega squared", x = "Groundtruth Score", y = "Inferred Score") +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, margin = margin(b = 15), family = "Arial"), axis.text = element_text(color = "black", size = 10, family = "Arial"), axis.title = element_text(color = "black", size = 13, family = "Arial"), panel.border = element_rect(color = "#333333", linewidth = 0.5), panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.4), panel.grid.minor = element_blank(), strip.background = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.5), strip.text = element_text(size = 10, family = "Arial", color = "black", margin = margin(t=4, b=4)), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12, family = "Arial"), legend.key.size = unit(0.8, "cm"))

out_11 <- file.path(export_dir, "Fig.Baselines_Scatters.pdf")
dims_11 <- get_template_dims(TEMPLATE_DIR, "Fig.Baselines_Scatters.pdf")
smart_ggsave_pdf(filename = out_11, plot_obj = p_baselines, dims = dims_11, default_w = 12, default_h = 9)
message(sprintf("✅ Fig.Baselines_Scatters.pdf exported to:\n   👉 %s", out_11))


# ==============================================================================
# 📊 Render 2: Gene-level evaluation metrics with standardized method names
# ==============================================================================
message("\n  -> Render 2: Gene-level bar plots...")

group_ordered_levels <- c("isCORAL", "competitors", "SNN baseline", "HVG baseline")
method_ordered_levels <- c("Seurat-MVP", "Seurat-vst", "scran", "Standard SNN", "CC-Regressed SNN", "PowerSeek", "GEMLI", "isCORAL")
group_mapping <- c("Seurat-MVP"="HVG baseline", "Seurat-vst"="HVG baseline", "scran"="HVG baseline", "Standard SNN"="SNN baseline", "CC-Regressed SNN"="SNN baseline", "PowerSeek"="competitors", "GEMLI"="competitors", "isCORAL"="isCORAL")

target_metrics <- c("AUROC", "AUPRC", "Max_F1", "Spearman_R", "Pearson_R")
plot_data_long <- gene_summary_table %>% 
  tidyr::pivot_longer(cols = -method, names_to = "Metric", values_to = "Score") %>% 
  filter(Metric %in% target_metrics) %>% 
  filter(method %in% method_ordered_levels) %>% 
  mutate(Metric = factor(Metric, levels = target_metrics), method = factor(method, levels = method_ordered_levels), Algorithm_Group = factor(group_mapping[as.character(method)], levels = group_ordered_levels))

make_metric_plot_clean <- function(data, metric_name, show_y_text = FALSE) { 
  df_sub <- data %>% filter(Metric == metric_name)
  p <- ggplot(df_sub, aes(x = method, y = Score, fill = method)) + 
    geom_col(aes(linewidth = ifelse(method == "isCORAL", 0.8, 0.4)), width = 0.75, color = "black") + 
    scale_linewidth_identity() +
    geom_text(aes(label = sprintf("%.3f", Score), color = ifelse(method %in% c("Seurat-MVP", "Seurat-vst", "PowerSeek", "Standard SNN"), "black", "white")), hjust = 1.15, size = 4.5, family = "Arial") + 
    scale_color_identity() + scale_fill_manual(values = method_colors_clean) + 
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 1), labels = c("0", "1"), expand = c(0, 0)) +
    coord_flip(clip = "off") + 
    facet_grid(Algorithm_Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = metric_name, x = NULL, y = NULL) + 
    theme_bw(base_size = 14, base_family = "Arial") + 
    theme(legend.position = "none", plot.title = element_text(size = 16, hjust = 0.5, family = "Arial"), panel.border = element_rect(color = "#333333", linewidth = 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(color = "black", size = 13, family = "Arial"), axis.ticks.y = element_blank(), strip.placement = "outside", strip.background.y = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.8), strip.text.y.left = element_text(size = 14, angle = 0, family = "Arial"), panel.spacing = unit(0.1, "lines"), plot.margin = margin(t = 10, r = 10, b = 10, l = 5)) 
  if (show_y_text) { p <- p + theme(axis.text.y = element_text(color = "black", size = 14, family = "Arial")) } else { p <- p + theme(axis.text.y = element_blank(), strip.text.y.left = element_blank(), strip.background.y = element_blank()) }
  return(p) 
}

p_auroc    <- make_metric_plot_clean(plot_data_long, "AUROC", TRUE)
p_auprc    <- make_metric_plot_clean(plot_data_long, "AUPRC", FALSE)
p_max_f1   <- make_metric_plot_clean(plot_data_long, "Max_F1", FALSE)
p_spearman <- make_metric_plot_clean(plot_data_long, "Spearman_R", FALSE)
p_pearson  <- make_metric_plot_clean(plot_data_long, "Pearson_R", FALSE)

p_metrics_panel <- p_auroc + p_auprc + p_max_f1 + p_spearman + p_pearson + plot_layout(nrow = 1, widths = c(1.8, 1, 1, 1, 1))

out_2 <- file.path(export_dir, "Fig.Metrics_Grouped_Barplots.pdf")
dims_2 <- get_template_dims(TEMPLATE_DIR, "Fig.Metrics_Grouped_Barplots.pdf")
smart_ggsave_pdf(filename = out_2, plot_obj = p_metrics_panel, dims = dims_2, default_w = 15, default_h = 6)
message(sprintf("✅ Fig.Metrics_Grouped_Barplots.pdf exported to:\n   👉 %s", out_2))


# ==============================================================================
# 📈 Render 3: Gene-level ROC and PR curves
# ==============================================================================
message("\n  -> Render 3: Gene-level ROC and PR curves...")

gene_roc_clean <- gene_roc_df %>% filter(method %in% method_ordered_levels) %>% mutate(method = factor(method, levels = method_ordered_levels))
gene_pr_clean <- gene_pr_df %>% filter(method %in% method_ordered_levels) %>% mutate(method = factor(method, levels = method_ordered_levels))

roc_breaks <- (gene_roc_clean %>% group_by(method) %>% summarise(auc = max(auc), .groups = "drop") %>% arrange(desc(auc)))$method
pr_breaks <- (gene_pr_clean %>% group_by(method) %>% summarise(auprc = max(auprc), .groups = "drop") %>% arrange(desc(auprc)))$method
iscoral_f1_pt <- gene_pr_clean %>% filter(method == "isCORAL") %>% mutate(current_f1 = ifelse(Precision + Recall == 0, 0, 2 * (Precision * Recall) / (Precision + Recall))) %>% arrange(desc(current_f1)) %>% head(1)

theme_curve_vis <- function() {
  theme_bw(base_size = 14, base_family = "Arial") + theme(plot.title = element_text(size = 16, hjust = 0.5, family = "Arial", margin = margin(b = 15)), panel.border = element_rect(color = "#333333", linewidth = 1.5), panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.5, linetype = "dashed"), panel.grid.minor = element_blank(), axis.text = element_text(color = "black", size = 13, family = "Arial"), axis.title = element_text(color = "black", size = 15, family = "Arial"), legend.text = element_text(size = 13, family = "Arial"), legend.title = element_blank(), legend.key.width = unit(1.5, "cm"), legend.background = element_rect(fill = alpha("white", 0.8), color = NA), legend.box.background = element_rect(color = "#333333", linewidth = 0.5), legend.margin = margin(t = 4, r = 8, b = 4, l = 8))
}

p_gene_roc <- ggplot(gene_roc_clean, aes(x = FPR, y = TPR, color = method)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#A0A0A0", linewidth = 1) + geom_line(aes(linewidth = ifelse(method == "isCORAL", 1.5, 0.8)), alpha = 0.9) + scale_linewidth_identity() + scale_color_manual(values = method_colors_clean, breaks = roc_breaks) + scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) + labs(title = "Gene Level ROC Curve", x = "False Positive Rate", y = "True Positive Rate") + theme_curve_vis() + theme(legend.position = c(0.72, 0.25))

p_gene_pr <- ggplot(gene_pr_clean, aes(x = Recall, y = Precision, color = method)) + geom_line(aes(linewidth = ifelse(method == "isCORAL", 1.5, 0.8)), alpha = 0.9) + scale_linewidth_identity() + scale_color_manual(values = method_colors_clean, breaks = pr_breaks) + scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) + labs(title = "Gene Level Precision-Recall Curve", x = "Recall", y = "Precision") + theme_curve_vis() + theme(legend.position = c(0.28, 0.25))

if (nrow(iscoral_f1_pt) > 0) { p_gene_pr <- p_gene_pr + geom_point(data = iscoral_f1_pt, aes(x = Recall, y = Precision), color = method_colors_clean["isCORAL"], size = 5, shape = 18, inherit.aes = FALSE) + geom_text(data = iscoral_f1_pt, aes(x = Recall, y = Precision, label = "Max F1"), color = method_colors_clean["isCORAL"], size = 4.5, family = "Arial", nudge_y = 0.05, nudge_x = 0.03, inherit.aes = FALSE) }

p_curves_panel <- p_gene_roc + p_gene_pr + plot_layout(nrow = 1)

out_3 <- file.path(export_dir, "Fig.ROC_PR_Curves.pdf")
dims_3 <- get_template_dims(TEMPLATE_DIR, "Fig.ROC_PR_Curves.pdf")
smart_ggsave_pdf(filename = out_3, plot_obj = p_curves_panel, dims = dims_3, default_w = 12, default_h = 6)
message(sprintf("✅ Fig.ROC_PR_Curves.pdf exported to:\n   👉 %s", out_3))


# ==============================================================================
# 🎻 Render 4: Combined gene-level plots
# ==============================================================================
message("\n  -> Render 4: Combined gene-level plots...")

theme_violin_vis <- function() {
  theme_bw(base_size = 14, base_family = "Arial") + theme(legend.position = "none", plot.title = element_text(size = 16, hjust = 0.5, family = "Arial"), plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5, color = "#555555", family = "Arial", margin = margin(b = 8)), panel.border = element_rect(color = "#333333", linewidth = 1.2), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9, family = "Arial"), axis.text.y = element_text(color = "black", size = 11, family = "Arial"), axis.title = element_text(size = 13, family = "Arial"), plot.margin = margin(10, 10, 10, 10))
}

render_showcase_gene <- function(seurat_obj, gene_name) {
  if (!gene_name %in% rownames(seurat_obj)) { message(sprintf("  ⚠️ Gene '%s' was not found and will be skipped.", gene_name)); return(NULL) }
  message(sprintf("  -> Rendering showcase gene: %s", gene_name))
  
  gt_info <- gene_collated %>% distinct(name, ground_truth_score) %>% mutate(rank = min_rank(-ground_truth_score)) %>% filter(name == gene_name) %>% mutate(method = "Groundtruth", score = ground_truth_score) %>% select(method, rank, score)
  rank_df <- gene_collated %>% group_by(method) %>% mutate(rank = min_rank(-score)) %>% filter(name == gene_name) %>% select(method, rank, score)
  all_ranks_df <- bind_rows(gt_info, rank_df) %>% mutate(method = ifelse(method == "CC-Regressed SNN", "SNN-CC", method))
  get_anno <- function(m_name) { row <- all_ranks_df %>% filter(method == m_name) %>% head(1); if(nrow(row) > 0) return(sprintf("Score: %.3f  |  Rank: %d", row$score, row$rank)); return("") }
  
  p_feature <- scCustomize::FeaturePlot_scCustom(seurat_obj, features = gene_name, na_cutoff = NULL) + labs(title = paste(gene_name, "Expression"), subtitle = NULL) + theme_bw(base_size = 14, base_family = "Arial") + theme(plot.title = element_text(size = 18, hjust = 0.5, family = "Arial", margin = margin(b=15)), panel.border = element_rect(color = "#333333", linewidth = 1.5, fill = NA), panel.grid = element_blank(), axis.text = element_text(color = "black", size = 12, family = "Arial"), axis.title = element_text(size = 14, family = "Arial"), aspect.ratio = 1, legend.position = "right", legend.title = element_blank(), legend.key.height = unit(0.8, "cm"), legend.text = element_text(size = 11, family = "Arial"))
  
  ISCORAL_K_VALUE <- b1_cache$run_parameters$ISCORAL_K_VALUE
  ISCORAL_BIN_SIZE <- b1_cache$run_parameters$ISCORAL_BIN_SIZE
  iscoral_clone_col_name <- paste0("metaclone_k", ISCORAL_K_VALUE, "_bs", ISCORAL_BIN_SIZE)
  
  if (!iscoral_clone_col_name %in% colnames(seurat_obj@meta.data)) { 
    seurat_obj <- SimulateClone(seurat_obj, bin_size = ISCORAL_BIN_SIZE, k_values = ISCORAL_K_VALUE) 
  }
  
  global_exp <- Seurat::FetchData(seurat_obj, vars = gene_name)[[1]]
  y_min <- min(global_exp, na.rm = TRUE); y_max <- max(global_exp, na.rm = TRUE)
  
  make_violin <- function(col, title, hex_light, hex_dark, method_name_for_anno) { 
    c_ids <- as.character(seurat_obj@meta.data[[col]])
    GLOBAL_CLONE_COL <- b1_cache$run_parameters$GLOBAL_CLONE_COL
    GLOBAL_INVALID_CLONES <- c("0","NA", "NaN", "", "None", "Unknown", "unknown", "Unassigned")
    valid <- if (col == GLOBAL_CLONE_COL) names(which(table(c_ids[!(c_ids %in% GLOBAL_INVALID_CLONES) & !is.na(c_ids)]) >= MIN_CLONE_SIZE)) else names(which(table(c_ids[!is.na(c_ids)]) >= MIN_CLONE_SIZE))
    anno_text <- get_anno(method_name_for_anno)
    if (length(valid) == 0) { return(ggplot() + theme_void() + labs(title = title, subtitle = paste(anno_text, "(No clones)", sep="\n")) + theme(plot.title=element_text(family="Arial", hjust=0.5), plot.subtitle=element_text(family="Arial", face="italic", color="#555555", hjust=0.5))) }
    clone_counts <- sort(table(c_ids[c_ids %in% valid]), decreasing = TRUE)
    samp <- head(names(clone_counts), N_CLONES_TO_SAMPLE)
    df <- data.frame(id = c_ids, exp = global_exp) %>% filter(id %in% samp)
    ord <- df %>% group_by(id) %>% summarise(m = mean(exp)) %>% arrange(desc(m)) %>% pull(id)
    df$id <- factor(df$id, levels = ord) 
    cols <- colorRampPalette(c(hex_light, hex_dark))(length(ord))
    ggplot(df, aes(x = id, y = exp, fill = id)) + geom_violin(scale = "width", trim = TRUE, alpha = 0.9, color = "#333333", linewidth = 0.5) + scale_fill_manual(values = setNames(cols, ord)) + labs(title = title, subtitle = anno_text, x = "Top Clones", y = "Expression") + theme_violin_vis()
  }
  
  # Render using the dynamically attached standardized metadata columns
  p_gt      <- make_violin(b1_cache$run_parameters$GLOBAL_CLONE_COL, "Groundtruth", "#E0E0E0", "#404040", "Groundtruth")
  p_iscoral <- make_violin(iscoral_clone_col_name,  "isCORAL",        "#F5C6C5", "#C74B46", "isCORAL")
  p_gemli   <- make_violin("gemli_clones_inferred", "GEMLI",          "#D5DFE8", "#3D5A80", "GEMLI")
  p_snn     <- make_violin("snn_cc_clones",         "SNN-CC",         "#D0EAE5", "#2A9D8F", "SNN-CC")
  
  p_violins_matrix <- wrap_plots(p_gt, p_iscoral, p_gemli, p_snn, ncol = 2) & scale_y_continuous(limits = c(y_min, y_max))
  p_right_side <- wrap_elements(p_violins_matrix + plot_annotation(title = "Clonal Expression Distribution", theme = theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Arial", margin = margin(b = 15)))))
  p_final <- wrap_elements(p_feature) | p_right_side
  p_final <- p_final + plot_layout(widths = c(1, 1.4)) & theme(text = element_text(family = "Arial"))
  
  file_4 <- sprintf("Fig.Combined_Showcase_%s.pdf", gene_name)
  out_4 <- file.path(export_dir, file_4)
  dims_4 <- get_template_dims(TEMPLATE_DIR, file_4)
  smart_ggsave_pdf(filename = out_4, plot_obj = p_final, dims = dims_4, default_w = 16, default_h = 7.5)
  message(sprintf("    ✅ %s exported to:\n       👉 %s", file_4, out_4))
}

# 🚀 Render all selected genes
target_genes <- c("NNMT", "HIST1H1E")
for (g in target_genes) {
  render_showcase_gene(seurat_obj, g)
}
message("\n🎉 B2 visualization rendering completed.")
