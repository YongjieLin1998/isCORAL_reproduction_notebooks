# ==============================================================================
#   # 🧮 Script B1: isCORAL Rapid Evaluation - Core Computation & Data Prep
#   #              (Refactored: True PB L1 + Distributions + Bootstrap AUPRC)
# ==============================================================================

# -----------------------------------------------------------------
# 🟢 [0] Unified Global Configuration Hub
# -----------------------------------------------------------------
# Set ISCORAL_PROJECT_ROOT, ISCORAL_SAMPLE_ID, and ISCORAL_SCAN_RESULTS as needed.
BASE_DIR      <- normalizePath(Sys.getenv("ISCORAL_PROJECT_ROOT", unset = "."), mustWork = TRUE)
TIMESTAMP     <- format(Sys.time(), "%Y%m%d_%H%M%S")
SAMPLE_ID     <- Sys.getenv("ISCORAL_SAMPLE_ID", unset = "<SAMPLE_ID>")
if (identical(SAMPLE_ID, "<SAMPLE_ID>")) {
  stop("Set ISCORAL_SAMPLE_ID or replace <SAMPLE_ID> before running this template.")
}
SCAN_RES_FILE <- Sys.getenv(
  "ISCORAL_SCAN_RESULTS",
  unset = file.path(BASE_DIR, "3_results_new", "<SCAN_RESULTS_FILE>.rds")
)

# Core Metadata and System-Level Filtering
GLOBAL_CLONE_COL      <- "barcode"
GLOBAL_INVALID_CLONES <- c("0", "NA", "NaN", "", "None", "Unknown", "unknown", "Unassigned")

# =================================================================
# 📸 [T=0 Memory Snapshot Lock]: Automatically fetch path and lock code!
# =================================================================
get_running_script_path <- function() {
  tryCatch({
    if (requireNamespace("this.path", quietly = TRUE)) {
      return(normalizePath(this.path::this.path()))
    } else {
      return("")
    }
  }, error = function(e) {
    return("")
  })
}

current_script_path <- get_running_script_path()
SNAPSHOT_CODE <- "Source unavailable."

if (current_script_path != "" && file.exists(current_script_path)) {
  SNAPSHOT_CODE <- paste(readLines(current_script_path, warn = FALSE), collapse = "\n")
  message(sprintf("📸 T=0 Snapshot Successful: Code locked in from [%s]. You can now freely modify the original file!", basename(current_script_path)))
} else {
  warning("⚠️ Auto-path detection failed. Snapshot unavailable.")
}
# =================================================================

# Gene Engine and Display Parameters
TOP_N_MEMORY_GENES   <- 2000
FILTER_MEAN_GT       <- 0.1
ISCORAL_BIN_SIZE     <- 100
ISCORAL_BINSIZE_NAME <- paste0("binsize_", ISCORAL_BIN_SIZE)
ISCORAL_METHOD_NAME  <- "isCORAL"
MIN_CLONE_SIZE       <- 2
ISCORAL_K_VALUE      <- 50
N_CLONES_TO_SAMPLE   <- 12

# -----------------------------------------------------------------
# --- STEP 1: Auto-Locate Cache, Load Results & Verify Meta-parameters
# -----------------------------------------------------------------
message("--- Step 1: Loading Dependencies & Auto-Locating Latest Baseline Cache ---")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(pbapply)
  library(PRROC)
  library(pROC)
  library(scran)
  library(SingleCellExperiment) 
  library(irlba)
  library(matrixStats)
})

pboptions(type = "timer", style = 3, char = "=")

BASELINE_DIR <- file.path(BASE_DIR, "4_benchmark/baseline")
all_baseline_files <- list.files(
  path = BASELINE_DIR, 
  pattern = sprintf("^Baseline_Cache_%s_.*\\.rds$", SAMPLE_ID), 
  full.names = TRUE
)

if (length(all_baseline_files) == 0) {
  stop(sprintf("❌ No Baseline Cache found for SAMPLE_ID '%s'! Please run Script A first.", SAMPLE_ID))
}

file_info <- file.info(all_baseline_files)
CACHE_FILE <- rownames(file_info)[which.max(file_info$mtime)]

message(sprintf("  🚀 Auto-loading the latest baseline cache:\n     %s", CACHE_FILE))
cache <- readRDS(CACHE_FILE)

if (!is.null(cache$run_parameters)) {
  params_to_check <- list(SAMPLE_ID = SAMPLE_ID, GLOBAL_CLONE_COL = GLOBAL_CLONE_COL)
  mismatches <- c()
  for (p_name in names(params_to_check)) {
    val_b1 <- params_to_check[[p_name]]
    val_cache <- cache$run_parameters[[p_name]]
    if (!is.null(val_cache) && !identical(val_b1, val_cache)) {
      mismatches <- c(mismatches, sprintf("    - %s: Current B1 is set to '%s', but Cache is '%s'", 
                                          p_name, paste(val_b1, collapse=","), paste(val_cache, collapse=",")))
    }
  }
  if (length(mismatches) > 0) {
    warning(paste(c("\n⚠️ CRITICAL WARNING: Parameters in script B1 do not match Baseline Cache parameters!", mismatches), collapse = "\n"))
    if (SAMPLE_ID != cache$run_parameters$SAMPLE_ID) stop("❌ Sample ID completely mismatched! Execution forcefully stopped.")
  } else { 
    message("  [✔] Perfect! Current parameter configuration completely matches the Baseline Cache history.") 
  }
}

seurat_obj <- cache$seurat_obj
scanning_results <- readRDS(SCAN_RES_FILE)
seurat_obj@misc$CORAL_scan_results <- scanning_results

if (!is.null(scanning_results$S_stab) && is.null(scanning_results[[ISCORAL_BINSIZE_NAME]])) {
  seurat_obj@misc$CORAL_scan_results <- list()
  seurat_obj@misc$CORAL_scan_results[[ISCORAL_BINSIZE_NAME]] <- scanning_results
}
if (is.null(seurat_obj@misc$CORAL_scan_results[[ISCORAL_BINSIZE_NAME]])) {
  stop("❌ CRITICAL ERROR: isCORAL extracted data not found! Please check BIN SIZE naming.")
}

GLOBAL_HVG_LIST <- VariableFeatures(seurat_obj)
if(length(GLOBAL_HVG_LIST) == 0) GLOBAL_HVG_LIST <- rownames(seurat_obj)[1:2000] else GLOBAL_HVG_LIST <- head(GLOBAL_HVG_LIST, 2000)

# -----------------------------------------------------------------
# --- STEP 2: Evaluate Gene-Level Scores (Centralized)
# -----------------------------------------------------------------
message("--- Step 2: Evaluating Gene Scores (Centralized) ---")

safe_get_matrix <- function(obj, assay_name, layer_name, force_dense = FALSE) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") mat <- LayerData(obj, assay = assay_name, layer = layer_name) 
  else mat <- GetAssayData(obj, assay = assay_name, slot = layer_name)
  if (force_dense) return(as.matrix(mat)) else return(mat)
}

calc_omega_sq <- function(expr_mat, groups) {
  valid <- !is.na(groups); grp_counts <- table(groups[valid]); valid_groups <- names(grp_counts[grp_counts > 1])
  valid <- valid & (as.character(groups) %in% valid_groups); expr <- expr_mat[, valid, drop=FALSE]
  grp <- as.factor(as.character(groups[valid])); k <- length(levels(grp)); n_cells <- ncol(expr)
  if (k <= 1) return(rep(0, nrow(expr)))
  sst <- Matrix::rowSums((expr - Matrix::rowMeans(expr))^2)
  ssw <- Matrix::rowSums((expr - t(rowsum(t(as.matrix(expr)), grp) / as.vector(table(grp)))[, as.numeric(grp)])^2)
  omega <- ((sst - ssw) - (k - 1) * (ssw / max(1, n_cells - k))) / (sst + (ssw / max(1, n_cells - k)))
  omega[is.na(omega) | omega < 0] <- 0
  return(setNames(omega, rownames(expr)))
}

exp_matrix_full <- safe_get_matrix(seurat_obj, DefaultAssay(seurat_obj), "data")
universal_gene_means <- data.frame(name = rownames(exp_matrix_full), mean_expr = Matrix::rowMeans(exp_matrix_full))

gene_results_list <- list()
if (!is.null(cache$baseline_raw_labels)) {
  for (method_name in names(cache$baseline_raw_labels)) {
    gene_results_list[[method_name]] <- data.frame(name = rownames(exp_matrix_full), score = calc_omega_sq(exp_matrix_full, cache$baseline_raw_labels[[method_name]]), method = method_name)
  }
}

tryCatch({ 
  obj_vst <- FindVariableFeatures(seurat_obj, assay=DefaultAssay(seurat_obj), selection.method="vst", nfeatures=nrow(seurat_obj), verbose=FALSE)
  hvf_info <- HVFInfo(obj_vst, assay = DefaultAssay(seurat_obj)); scores <- hvf_info[[grep("variance.standardized", colnames(hvf_info), value = TRUE)[1]]]
  scores[is.na(scores) | !is.finite(scores)] <- min(scores, na.rm = TRUE) - 1e-5
  gene_results_list[["Seurat-VST"]] <- data.frame(name = rownames(hvf_info), score = scores, method = "Seurat-VST")
}, error = function(e) message("  ⚠️ Seurat-VST failed: ", e$message))

tryCatch({ 
  obj_mvp <- FindVariableFeatures(seurat_obj, assay=DefaultAssay(seurat_obj), selection.method="mean.var.plot", verbose=FALSE)
  hvf_info_mvp <- HVFInfo(obj_mvp, assay = DefaultAssay(seurat_obj)); score_col <- grep("dispersion.scaled", colnames(hvf_info_mvp), value = TRUE)[1]
  if (is.na(score_col)) score_col <- grep("dispersion", colnames(hvf_info_mvp), value = TRUE)[1]
  scores <- hvf_info_mvp[[score_col]]; scores[is.na(scores) | !is.finite(scores)] <- min(scores, na.rm = TRUE) - 1e-5
  gene_results_list[["Seurat-MVP"]] <- data.frame(name = rownames(hvf_info_mvp), score = scores, method = "Seurat-MVP")
}, error = function(e) message("  ⚠️ Seurat-MVP failed: ", e$message))

tryCatch({ 
  scran_df <- as.data.frame(scran::modelGeneVar(as.SingleCellExperiment(seurat_obj))) %>% tibble::rownames_to_column("name"); scores <- scran_df$bio
  scores[is.na(scores) | !is.finite(scores)] <- min(scores, na.rm = TRUE) - 1e-5
  gene_results_list[["Scran"]] <- data.frame(name = scran_df$name, score = scores, method="Scran") 
}, error = function(e) message("  ⚠️ Scran failed: ", e$message))

tryCatch({
  tmp_scaled <- safe_get_matrix(ScaleData(seurat_obj, features=rownames(seurat_obj), verbose=FALSE), DefaultAssay(seurat_obj), "scale.data", force_dense=FALSE)
  S_filt_mat <- as.matrix(tmp_scaled); rv <- matrixStats::rowVars(S_filt_mat); S_filt <- S_filt_mat[rv > 0, , drop=FALSE]
  k_calc <- min(30, min(nrow(S_filt), ncol(S_filt)) - 1); set.seed(42)
  svd_res <- irlba::irlba(S_filt, nv=k_calc, nu=k_calc)
  scores_all <- rowSums(sweep(svd_res$u ^ 2, 2, (svd_res$d ^ 2) / (ncol(S_filt) - 1), "*"))
  min_score <- min(scores_all, na.rm=TRUE); df_ps1 <- data.frame(name = rownames(S_filt), score = scores_all)
  diff_genes <- setdiff(rownames(seurat_obj), rownames(S_filt))
  if (length(diff_genes) > 0) {
    df_ps2 <- data.frame(name = diff_genes, score = min_score - abs(min_score) * 0.05)
    gene_results_list[["PowerSeek"]] <- bind_rows(df_ps1, df_ps2) %>% mutate(method="PowerSeek")
  } else { gene_results_list[["PowerSeek"]] <- df_ps1 %>% mutate(method="PowerSeek") }
}, error = function(e) message("  ⚠️ PowerSeek failed: ", e$message))

iscoral_raw <- seurat_obj@misc$CORAL_scan_results[[ISCORAL_BINSIZE_NAME]]$ss_filter_average
gene_results_list[[ISCORAL_METHOD_NAME]] <- iscoral_raw %>% dplyr::select(name, score = omega_square_infer) %>% mutate(method = ISCORAL_METHOD_NAME)

gt_df <- seurat_obj@misc$groundtruth$heritable_genes_df %>% dplyr::select(name, ground_truth_score = Omega_square) %>% 
  dplyr::filter(!is.na(ground_truth_score)) %>% dplyr::inner_join(universal_gene_means, by = "name") %>% 
  dplyr::filter(mean_expr > FILTER_MEAN_GT) %>% dplyr::arrange(desc(ground_truth_score)) %>% 
  dplyr::mutate(is_high_truth = dplyr::row_number() <= TOP_N_MEMORY_GENES)

gene_collated <- dplyr::bind_rows(gene_results_list) %>% dplyr::inner_join(gt_df, by="name") %>% 
  dplyr::filter(!is.na(is_high_truth), is.finite(ground_truth_score), is.finite(score))

GLOBAL_MEMORY_GENES <- head(gt_df$name[gt_df$is_high_truth], 2000)

calc_roc <- function(labels, scores) { 
  P <- sum(labels, na.rm=T); N <- sum(!labels, na.rm=T); if(P == 0 || N == 0) return(list(df=data.frame(FPR=c(0,1), TPR=c(0,1)), auc=0.5))
  thresh <- sort(unique(c(Inf, scores, -Inf)), decreasing = TRUE); FPRs <- numeric(length(thresh)); TPRs <- numeric(length(thresh))
  for(i in seq_along(thresh)) { preds <- scores >= thresh[i]; TPRs[i] <- sum(labels & preds, na.rm=T) / P; FPRs[i] <- sum(!labels & preds, na.rm=T) / N }
  ord <- order(FPRs); FPRs <- FPRs[ord]; TPRs <- TPRs[ord]; uniq <- !duplicated(FPRs, fromLast=TRUE); FPRs <- FPRs[uniq]; TPRs <- TPRs[uniq]
  auc <- sum(diff(FPRs) * (TPRs[-1] + TPRs[-length(TPRs)])/2, na.rm=T)
  return(list(df = data.frame(FPR=FPRs, TPR=TPRs), auc = auc)) 
}

calc_pr <- function(labels, scores) { 
  P <- sum(labels, na.rm=T); N <- sum(!labels, na.rm=T); if(P == 0 || N == 0) return(list(df=data.frame(Recall=c(0,1), Precision=c(1,0), Threshold=c(Inf,-Inf)), auprc=0, max_f1=0, best_thresh=0))
  df <- data.frame(score = scores, label = labels)[order(scores, decreasing = TRUE, na.last = NA), ]; TP <- 0; FP <- 0; pts <- data.frame(Recall=0, Precision=1, Threshold=Inf)
  for(i in 1:nrow(df)) { if(df$label[i]) TP <- TP + 1 else FP <- FP + 1; if(i < nrow(df) && df$score[i] == df$score[i+1]) next; pts <- rbind(pts, c(TP/P, ifelse(TP+FP==0, 0, TP/(TP+FP)), df$score[i])) }
  if(TP != P || FP != N) pts <- rbind(pts, c(1, P/(P+N), -Inf)); pr_auc_df <- pts[!duplicated(pts$Recall), ]
  auprc <- sum(diff(pr_auc_df$Recall) * (pr_auc_df$Precision[-1] + pr_auc_df$Precision[-nrow(pr_auc_df)])/2, na.rm=T)
  f1 <- ifelse(pts$Precision + pts$Recall == 0, 0, 2*(pts$Precision*pts$Recall)/(pts$Precision+pts$Recall)); f1[1] <- 0
  return(list(df = pts, auprc = auprc, max_f1 = max(f1, na.rm=T), best_thresh = pts$Threshold[which.max(f1)])) 
}

gene_roc_list <- list(); gene_pr_list <- list()
for (m_name in unique(gene_collated$method)) { 
  sub_d <- filter(gene_collated, method == m_name); r_res <- calc_roc(sub_d$is_high_truth, sub_d$score); p_res <- calc_pr(sub_d$is_high_truth, sub_d$score)
  gene_roc_list[[m_name]] <- r_res$df %>% mutate(method = m_name, auc = r_res$auc); gene_pr_list[[m_name]] <- p_res$df %>% mutate(method = m_name, auprc = p_res$auprc, max_f1 = p_res$max_f1, best_thresh = p_res$best_thresh) 
}

gene_roc_df <- bind_rows(gene_roc_list); gene_pr_df <- bind_rows(gene_pr_list)
omega_methods <- c("isCORAL", "GEMLI", "Standard SNN", "CC-Regressed SNN")

gene_summary_table <- gene_collated %>% group_by(method) %>% 
  summarise(Pearson_R = cor(ground_truth_score, score, method="pearson"), Spearman_R = cor(ground_truth_score, score, method="spearman"), RMSE = sqrt(mean((ground_truth_score - score)^2, na.rm = TRUE)), .groups="drop") %>% 
  mutate(RMSE = ifelse(method %in% omega_methods, RMSE, NA)) %>%
  left_join(gene_roc_df %>% group_by(method) %>% summarise(AUROC=max(auc)), by="method") %>% 
  left_join(gene_pr_df %>% group_by(method) %>% summarise(AUPRC=max(auprc), Max_F1=max(max_f1)), by="method") %>% select(method, AUROC, AUPRC, Max_F1, RMSE, Spearman_R, Pearson_R)

all_best_thresholds <- lapply(names(gene_pr_list), function(m_name) { 
  pr_data <- gene_pr_list[[m_name]]; if (is.null(pr_data) || nrow(pr_data) == 0) return(NULL)
  f1 <- ifelse(pr_data$Precision + pr_data$Recall == 0, 0, 2 * (pr_data$Precision * pr_data$Recall) / (pr_data$Precision + pr_data$Recall)); data.frame(method = m_name, best_thresh = pr_data$Threshold[which.max(f1)], stringsAsFactors = F) 
}) %>% dplyr::bind_rows()

all_method_scatter_data <- gene_collated %>% dplyr::left_join(all_best_thresholds, by = "method") %>% dplyr::filter(!is.na(best_thresh)) %>% 
  dplyr::mutate(class = case_when(is_high_truth & score >= best_thresh ~ "TP", !is_high_truth & score < best_thresh ~ "TN", !is_high_truth & score >= best_thresh ~ "FP", TRUE ~ "FN")) %>% 
  dplyr::mutate(method = factor(method, levels = c(ISCORAL_METHOD_NAME, setdiff(unique(method), ISCORAL_METHOD_NAME))))

# -----------------------------------------------------------------
# --- STEP 3: Evaluate Matrix-Level Scores 
# -----------------------------------------------------------------
message("--- Step 3: Evaluating Matrix Scores (Rigorous ED, AUPRC Enrichment & True PB L1 Distance) ---")

align_matrix_names <- function(mat, obj) { 
  if (is.null(mat)) return(NULL); cells <- colnames(obj)
  if (is.null(rownames(mat)) || length(intersect(rownames(mat), cells)) == 0) { 
    if (nrow(mat) == length(cells) && ncol(mat) == length(cells)) { rownames(mat) <- cells; colnames(mat) <- cells } else return(NULL) 
  }
  if (length(mat) > 0 && max(mat, na.rm=TRUE) > 0) diag(mat) <- max(mat, na.rm=TRUE); return(mat) 
}

S_stab_mat_fresh <- align_matrix_names(seurat_obj@misc$CORAL_scan_results[[ISCORAL_BINSIZE_NAME]]$S_stab, seurat_obj)
if(is.null(cache$baseline_raw_matrices)) cache$baseline_raw_matrices <- list()
cache$baseline_raw_matrices[["isCORAL"]] <- S_stab_mat_fresh

benchmark_lineage_engine_raw <- function(seurat_obj, prob_matrix, method_name, pca_dims=1:30) {
  if (is.null(prob_matrix)) return(NULL)
  metadata <- seurat_obj@meta.data; valid_cells <- rownames(metadata)[!is.na(metadata[[GLOBAL_CLONE_COL]]) & !as.character(metadata[[GLOBAL_CLONE_COL]]) %in% GLOBAL_INVALID_CLONES]
  is_uniform <- is.character(prob_matrix) && prob_matrix[1] == "UNIFORM"
  
  if (!is_uniform) { 
    valid_cells <- intersect(valid_cells, rownames(prob_matrix)); prob_matrix_sub <- as.matrix(prob_matrix[valid_cells, valid_cells, drop=FALSE]); prob_matrix_sub[is.na(prob_matrix_sub)] <- 0 
  }
  if (length(valid_cells) < 2) return(NULL)
  
  D_global <- as.matrix(stats::dist(Seurat::Embeddings(seurat_obj, "bench_pca")[valid_cells, pca_dims]))
  gt_labels_vec <- as.character(metadata[valid_cells, GLOBAL_CLONE_COL]); names(gt_labels_vec) <- valid_cells
  
  gt_stats <- lapply(unique(gt_labels_vec), function(gt_label) {
    bc_cells <- valid_cells[gt_labels_vec == gt_label]; n_c <- length(bc_cells); if (n_c <= 1) return(list(cells=bc_cells, eps_q=NA)) 
    corr2 <- n_c/(n_c-1); e_w2_in <- mean(D_global[bc_cells, bc_cells]) * corr2; cross_dist_vector <- rowMeans(D_global[, bc_cells, drop=FALSE])
    eps_q <- quantile(replicate(100, { samp_idx <- sample(bc_cells, n_c, replace=TRUE); max(0, 2*mean(cross_dist_vector[samp_idx]) - mean(D_global[samp_idx, samp_idx])*corr2 - e_w2_in) }), 0.9)
    list(cells=bc_cells, eps_q=eps_q, e_w2_in=e_w2_in, cross_dist_vector=cross_dist_vector)
  }); names(gt_stats) <- unique(gt_labels_vec)
  
  N_cells <- length(valid_cells)
  if (is_uniform) { D_row_sums <- rowSums(D_global); D_total_sum <- sum(D_row_sums) }
  
  res_list <- pbapply::pblapply(valid_cells, function(cell_i) {
    stats <- gt_stats[[gt_labels_vec[cell_i]]]
    if(is.null(stats$e_w2_in)) return(list(seed_cell=cell_i, Method=method_name, model_error_ed=NA, baseline_error_ed=stats$eps_q))
    
    if (is_uniform) { 
      cross_sum_excluding_i <- sum(stats$cross_dist_vector) - stats$cross_dist_vector[cell_i]
      E_cross <- cross_sum_excluding_i / (N_cells - 1); sum_D_minus_i <- D_total_sum - (2 * D_row_sums[cell_i]); E_internal <- sum_D_minus_i / ((N_cells - 1)^2)
      model_ed <- max(0, 2 * E_cross - E_internal - stats$e_w2_in) 
    } else {
      prob_vec <- prob_matrix_sub[cell_i, ]; if (cell_i %in% names(prob_vec)) { prob_vec[cell_i] <- 0 }; new_prob_sum <- sum(prob_vec)
      if(new_prob_sum == 0) return(list(seed_cell=cell_i, Method=method_name, model_error_ed=NA, baseline_error_ed=stats$eps_q))
      prob_vec <- prob_vec / new_prob_sum; nz_idx <- which(prob_vec > 1e-6); if (length(nz_idx) > 400) nz_idx <- nz_idx[order(prob_vec[nz_idx], decreasing = TRUE)[1:400]]
      if(length(nz_idx) == 0) return(list(seed_cell=cell_i, Method=method_name, model_error_ed=NA, baseline_error_ed=stats$eps_q))
      w_clean <- prob_vec[nz_idx] / sum(prob_vec[nz_idx]); model_ed <- max(0, 2 * sum(w_clean * stats$cross_dist_vector[nz_idx]) - sum(w_clean * (D_global[nz_idx, nz_idx, drop=FALSE] %*% w_clean)) - stats$e_w2_in)
    }
    list(seed_cell=cell_i, Method=method_name, model_error_ed=model_ed, baseline_error_ed=stats$eps_q)
  })
  if(!is_uniform) { rm(prob_matrix_sub); rm(D_global); gc() }
  return(dplyr::bind_rows(res_list))
}

benchmark_pairwise_dual <- function(prob_matrix, method_name) { 
  if (is.null(prob_matrix)) return(NULL)
  gt_labels <- as.character(seurat_obj@meta.data[[GLOBAL_CLONE_COL]]); valid_idx <- !is.na(gt_labels) & !(gt_labels %in% GLOBAL_INVALID_CLONES); valid_cells <- rownames(seurat_obj@meta.data)[valid_idx]
  clone_counts <- table(gt_labels[valid_idx]); valid_clones <- names(clone_counts[clone_counts > 1]); valid_cells <- intersect(valid_cells[gt_labels[valid_idx] %in% valid_clones], rownames(prob_matrix))
  if (length(valid_cells) < 2) return(NULL)
  
  gt_clean <- outer(gt_labels[match(valid_cells, rownames(seurat_obj@meta.data))], gt_labels[match(valid_cells, rownames(seurat_obj@meta.data))], "==")
  gt_clean <- gt_clean[upper.tri(matrix(0, length(valid_cells), length(valid_cells)))]
  score_clean <- as.matrix(prob_matrix[valid_cells, valid_cells])
  
  
  # === Diagnostic output for the current method ===
  message(sprintf("\n🔍 Inspecting method: %s", method_name))
  message(sprintf("   - Valid cells evaluated: %d", length(valid_cells)))
  message(sprintf("   - NA values in the subset matrix: %d", sum(is.na(score_clean))))
  # ==========================
  
  score_clean <- score_clean[upper.tri(matrix(0, length(valid_cells), length(valid_cells)))]
  
  v <- !is.na(gt_clean) & !is.na(score_clean); gt_clean <- gt_clean[v]; score_clean <- score_clean[v]
  if (length(gt_clean) > 2000000) { set.seed(42); idx <- sample(length(gt_clean), 2000000); gt_clean <- gt_clean[idx]; score_clean <- score_clean[idx] }
  random_auprc_base <- max(sum(gt_clean) / length(gt_clean), 1e-8)
  
  if (length(unique(gt_clean)) < 2) {
    return(list(roc = data.frame(fpr=c(0,1), tpr=c(0,1), Method=method_name, auc=NA), pr = data.frame(recall=c(0,1), precision=c(1,0), Method=method_name, auprc=NA, Random_AUPRC=random_auprc_base, AUPRC_Enrichment=NA)))
  }
  r <- pROC::roc(response=gt_clean, predictor=score_clean, levels=c(FALSE, TRUE), direction="<", quiet=TRUE)
  p <- PRROC::pr.curve(scores.class0=score_clean[gt_clean==TRUE], scores.class1=score_clean[gt_clean==FALSE], curve=TRUE)
  list(roc = data.frame(fpr=1-r$specificities, tpr=r$sensitivities, Method=method_name, auc=as.numeric(r$auc)), 
       pr = data.frame(recall=p$curve[,1], precision=p$curve[,2], Method=method_name, auprc=p$auc.integral, Random_AUPRC=random_auprc_base, AUPRC_Enrichment=p$auc.integral/random_auprc_base)) 
}

# 💡 [MODIFICATION] True Pseudobulk (Counts-then-Normalize) L1 Distance
benchmark_pseudobulk_l1 <- function(seurat_obj, prob_matrix, method_name) {
  if (is.null(prob_matrix)) return(data.frame(seed_cell=NA_character_, Method=method_name, PB_L1_Distance=NA, stringsAsFactors=FALSE))
  
  metadata <- seurat_obj@meta.data
  valid_idx <- !is.na(metadata[[GLOBAL_CLONE_COL]]) & !(metadata[[GLOBAL_CLONE_COL]] %in% GLOBAL_INVALID_CLONES)
  gt_labels <- as.character(metadata[[GLOBAL_CLONE_COL]])[valid_idx]
  valid_cells <- rownames(metadata)[valid_idx]
  
  clone_counts <- table(gt_labels); valid_clones <- names(clone_counts[clone_counts > 1])
  valid_cells <- valid_cells[gt_labels %in% valid_clones]; gt_labels <- gt_labels[gt_labels %in% valid_clones]
  is_uniform <- is.character(prob_matrix) && prob_matrix[1] == "UNIFORM"
  if (!is_uniform) valid_cells <- intersect(valid_cells, rownames(prob_matrix))
  
  if (length(valid_cells) < 2) return(data.frame(seed_cell=NA_character_, Method=method_name, PB_L1_Distance=NA, stringsAsFactors=FALSE))
  gt_labels <- as.character(metadata[valid_cells, GLOBAL_CLONE_COL])
  
  C_all <- safe_get_matrix(seurat_obj, DefaultAssay(seurat_obj), "counts")[, valid_cells, drop=FALSE]
  G <- outer(gt_labels, gt_labels, "==") * 1.0; diag(G) <- 0; G_rs <- rowSums(G); valid_G <- G_rs > 0; G[valid_G, ] <- G[valid_G, ] / G_rs[valid_G]
  
  if (is_uniform) {
    P <- matrix(1/(length(valid_cells)-1), nrow=length(valid_cells), ncol=length(valid_cells), dimnames=list(valid_cells, valid_cells)); diag(P) <- 0
  } else {
    P <- as.matrix(prob_matrix[valid_cells, valid_cells, drop=FALSE]); P[is.na(P)] <- 0; diag(P) <- 0; P_rs <- rowSums(P); valid_P <- P_rs > 1e-6; P[valid_P, ] <- P[valid_P, ] / P_rs[valid_P]
  }
  
  valid_eval <- if(is_uniform) valid_G else (valid_G & valid_P)
  if(sum(valid_eval) < 2) return(data.frame(seed_cell=NA_character_, Method=method_name, PB_L1_Distance=NA, stringsAsFactors=FALSE))
  
  C_true_bulk <- C_all %*% t(G[valid_eval, , drop=FALSE]); C_pred_bulk <- C_all %*% t(P[valid_eval, , drop=FALSE])
  normalize_pseudobulk <- function(count_mat) { col_sums <- colSums(count_mat); col_sums[col_sums == 0] <- 1; return(log1p(sweep(count_mat, 2, col_sums, "/") * 10000)) }
  Y_true_norm <- normalize_pseudobulk(C_true_bulk); Y_pred_norm <- normalize_pseudobulk(C_pred_bulk)
  l1_all <- sapply(1:ncol(Y_pred_norm), function(i) mean(abs(Y_pred_norm[,i] - Y_true_norm[,i]), na.rm=TRUE))
  
  return(data.frame(seed_cell = colnames(Y_pred_norm), Method = method_name, PB_L1_Distance = l1_all, stringsAsFactors = FALSE))
}

all_ed_results <- list(); all_roc_results <- list(); all_pr_results <- list(); all_pb_results <- list()
valid_gt_cells <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data[[GLOBAL_CLONE_COL]]) & !as.character(seurat_obj@meta.data[[GLOBAL_CLONE_COL]]) %in% GLOBAL_INVALID_CLONES]
gt_labels_clean <- as.character(seurat_obj@meta.data[[GLOBAL_CLONE_COL]][match(valid_gt_cells, rownames(seurat_obj@meta.data))])

gt_matrix_perfect <- outer(gt_labels_clean, gt_labels_clean, "==") * 1.0; rownames(gt_matrix_perfect) <- valid_gt_cells; colnames(gt_matrix_perfect) <- valid_gt_cells

set.seed(42)
gt_matrix_bootstrap <- matrix(0, nrow=length(valid_gt_cells), ncol=length(valid_gt_cells), dimnames=list(valid_gt_cells, valid_gt_cells))
for (clone in unique(gt_labels_clean)) {
  idx <- which(gt_labels_clean == clone)
  if (length(idx) > 1) {
    for (i in idx) {
      peers <- setdiff(idx, i)
      if (length(peers) == 1) samp <- peers else samp <- sample(peers, size=length(peers), replace=TRUE)
      counts <- table(samp); gt_matrix_bootstrap[i, as.integer(names(counts))] <- as.numeric(counts)
    }
  }
}

# 🔴 [CRITICAL FIX] Symmetrize the bootstrap matrix so upper-triangle AUPRC evaluation retains one-way sampling probabilities
gt_matrix_bootstrap <- (gt_matrix_bootstrap + t(gt_matrix_bootstrap)) / 2

all_ed_results[["Theoretical Best (Bootstrap Floor)"]] <- benchmark_lineage_engine_raw(seurat_obj, gt_matrix_bootstrap, "Theoretical Best (Bootstrap Floor)")
all_pb_results[["Theoretical Best (Bootstrap Floor)"]] <- benchmark_pseudobulk_l1(seurat_obj, gt_matrix_bootstrap, "Theoretical Best (Bootstrap Floor)")

# 🔴 [MODIFIED] Explicitly calculate and retain ROC and PR scores for the theoretical optimum
dual_res_best <- benchmark_pairwise_dual(gt_matrix_bootstrap, "Theoretical Best (Bootstrap Floor)")
if (!is.null(dual_res_best)) {
  all_roc_results[["Theoretical Best (Bootstrap Floor)"]] <- dual_res_best$roc
  all_pr_results[["Theoretical Best (Bootstrap Floor)"]] <- dual_res_best$pr
}

all_ed_results[["Uniform Allocation (Worst Case)"]] <- benchmark_lineage_engine_raw(seurat_obj, "UNIFORM", "Uniform Allocation (Worst Case)")
all_pb_results[["Uniform Allocation (Worst Case)"]] <- benchmark_pseudobulk_l1(seurat_obj, "UNIFORM", "Uniform Allocation (Worst Case)")

for (method_name in names(cache$baseline_raw_matrices)) {
  mat <- cache$baseline_raw_matrices[[method_name]]
  if (!is.null(mat)) {
    all_ed_results[[method_name]] <- benchmark_lineage_engine_raw(seurat_obj, mat, method_name)
    all_pb_results[[method_name]] <- benchmark_pseudobulk_l1(seurat_obj, mat, method_name)
    dual_res <- benchmark_pairwise_dual(mat, method_name)
    if (!is.null(dual_res)) { all_roc_results[[method_name]] <- dual_res$roc; all_pr_results[[method_name]] <- dual_res$pr }
  }
}

df_all_ed_raw <- dplyr::bind_rows(all_ed_results) %>% dplyr::rename(Energy_Distance = model_error_ed)
matrix_roc_all <- dplyr::bind_rows(all_roc_results)
matrix_pr_all <- dplyr::bind_rows(all_pr_results)
df_all_pb_raw <- dplyr::bind_rows(all_pb_results)

# -----------------------------------------------------------------
# --- STEP 4: Data Aggregation & Export 
# -----------------------------------------------------------------
message("--- Step 4: Data Aggregation & Exporting Excel ---")

pb_summary_table <- df_all_pb_raw %>% group_by(Method) %>% summarise(Mean_PB_L1 = mean(PB_L1_Distance, na.rm=TRUE), Median_PB_L1 = median(PB_L1_Distance, na.rm=TRUE), SD_PB_L1 = sd(PB_L1_Distance, na.rm=TRUE), .groups="drop") %>% arrange(Mean_PB_L1)
ed_summary_table <- df_all_ed_raw %>% group_by(Method) %>% summarise(Mean_ED = mean(Energy_Distance, na.rm=T), Median_ED = median(Energy_Distance, na.rm=T), SD_ED = sd(Energy_Distance, na.rm=T), Mean_Baseline_ED_Floor = mean(baseline_error_ed, na.rm=T), .groups="drop") %>% arrange(Mean_ED)

if(nrow(matrix_roc_all) > 0) matrix_roc_agg <- matrix_roc_all %>% dplyr::group_by(Method) %>% dplyr::summarise(AUROC = max(auc, na.rm=TRUE), .groups="drop") else matrix_roc_agg <- data.frame(Method = character(), AUROC = numeric()) 
if(nrow(matrix_pr_all) > 0) matrix_pr_agg <- matrix_pr_all %>% dplyr::group_by(Method) %>% dplyr::summarise(AUPRC = max(auprc, na.rm=TRUE), Random_AUPRC = mean(Random_AUPRC, na.rm=TRUE), AUPRC_Enrichment = max(AUPRC_Enrichment, na.rm=TRUE), .groups="drop") else matrix_pr_agg <- data.frame(Method = character(), AUPRC = numeric(), Random_AUPRC=numeric(), AUPRC_Enrichment=numeric()) 

matrix_summary <- df_all_ed_raw %>% 
  dplyr::filter(!Method %in% c("Theoretical Best (Bootstrap Floor)", "Uniform Allocation (Worst Case)")) %>% 
  dplyr::group_by(Method) %>% dplyr::summarise(Mean_Energy_Distance = mean(Energy_Distance, na.rm=T), .groups="drop") %>% 
  dplyr::mutate(Method_Short = dplyr::case_when(grepl("isCORAL", Method) ~ "isCORAL", grepl("GEMLI", Method) ~ "GEMLI", grepl("Standard SNN", Method) ~ "Standard SNN", grepl("CC-Regressed SNN", Method) ~ "CC-Regressed SNN", grepl("Standard KNN", Method) ~ "Standard KNN", grepl("CC-Regressed KNN", Method) ~ "CC-Regressed KNN", TRUE ~ as.character(Method))) %>% 
  dplyr::left_join(matrix_roc_agg %>% dplyr::rename(Method_Short = Method), by="Method_Short") %>% 
  dplyr::left_join(matrix_pr_agg %>% dplyr::rename(Method_Short = Method), by="Method_Short") %>% 
  dplyr::left_join(pb_summary_table %>% dplyr::rename(PB_L1_Distance = Mean_PB_L1) %>% dplyr::mutate(Method_Short = dplyr::case_when(grepl("isCORAL", Method) ~ "isCORAL", grepl("GEMLI", Method) ~ "GEMLI", grepl("Standard SNN", Method) ~ "Standard SNN", grepl("CC-Regressed SNN", Method) ~ "CC-Regressed SNN", grepl("Standard KNN", Method) ~ "Standard KNN", grepl("CC-Regressed KNN", Method) ~ "CC-Regressed KNN", TRUE ~ as.character(Method))) %>% dplyr::select(-Method), by="Method_Short") %>%
  dplyr::select(Algorithm = Method_Short, Mean_Energy_Distance, AUROC, AUPRC, Random_AUPRC, AUPRC_Enrichment, PB_L1_Distance) %>% dplyr::arrange(desc(AUPRC_Enrichment))

message("\n★★★ [Table 1] Gene Engine Scorecard ★★★")
print(knitr::kable(gene_summary_table, digits = 3, align = "lcccccc"))

message("\n★★★ [Table 2] Matrix Engine Comprehensive Scorecard (With True Pseudobulk L1 Distance) ★★★")
print(knitr::kable(matrix_summary, digits = 3, align = "lcccccc"))

run_uid <- paste0(SAMPLE_ID, "_bin", ISCORAL_BIN_SIZE, "_", TIMESTAMP)
export_dir <- file.path(BASE_DIR, "4_benchmark/results", paste0("Benchmark_Results_", run_uid))
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

excel_sheets_list <- list()
run_params <- seurat_obj@misc[["CORAL_scan_results"]][["run_parameters"]]
if (!is.null(run_params)) { 
  params_clean <- lapply(run_params, function(x) { if(is.null(x)) return("NULL"); if(length(x) > 1) return(paste(x, collapse = " | ")); return(as.character(x)) })
  params_clean$benchmark_completion_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  excel_sheets_list[["Run_Parameters"]] <- data.frame(Parameter = names(params_clean), Value = unlist(params_clean), stringsAsFactors = FALSE) 
}

excel_sheets_list[["Gene_Level_Metrics"]]    <- gene_summary_table
excel_sheets_list[["Matrix_Level_Metrics"]]  <- matrix_summary
excel_sheets_list[["Energy_Distance_Stats"]] <- ed_summary_table
excel_sheets_list[["PB_L1_Distance_Stats"]]  <- pb_summary_table

openxlsx::write.xlsx(excel_sheets_list, file = file.path(export_dir, paste0("Comprehensive_Benchmark_", run_uid, ".xlsx")), rowNames = FALSE)

#plot_roc_df <- if(nrow(matrix_roc_all) > 0) matrix_roc_all %>% filter(!grepl("KNN", Method)) else matrix_roc_all
#plot_pr_df  <- if(nrow(matrix_pr_all) > 0) matrix_pr_all %>% filter(!grepl("KNN", Method)) else matrix_pr_all
#plot_ed_df  <- if(nrow(df_all_ed_raw) > 0) df_all_ed_raw %>% filter(!grepl("KNN", Method)) else df_all_ed_raw


plot_roc_df <- matrix_roc_all
plot_pr_df  <- matrix_pr_all
plot_ed_df  <- df_all_ed_raw

# ==============================================================================
# 🗂️ STEP 5: Comprehensive Archiving
# ==============================================================================
message("--- Step 5: Archiving B1 Results & Script ---")

backup_script_name <- sprintf("Eval_Script_B1_%s.R", run_uid)
target_script_path <- file.path(export_dir, backup_script_name)
writeLines(SNAPSHOT_CODE, target_script_path)
message(sprintf("  ✅ Script source file (T=0 Snapshot) successfully backed up to: %s", target_script_path))

session_info_str <- paste(capture.output(print(sessionInfo())), collapse = "\n")

run_parameters_b1 <- list(BASE_DIR = BASE_DIR, SAMPLE_ID = SAMPLE_ID, TIMESTAMP = TIMESTAMP, SCAN_RES_FILE = SCAN_RES_FILE, CACHE_FILE = CACHE_FILE, GLOBAL_CLONE_COL = GLOBAL_CLONE_COL, TOP_N_MEMORY_GENES = TOP_N_MEMORY_GENES, FILTER_MEAN_GT = FILTER_MEAN_GT, ISCORAL_BIN_SIZE = ISCORAL_BIN_SIZE, ISCORAL_K_VALUE = ISCORAL_K_VALUE, K_SCAN_RANGE = tryCatch(K_SCAN_RANGE, error=function(e) NULL), SCRIPT_PATH = current_script_path, SOURCE_CODE = SNAPSHOT_CODE, SESSION_INFO = session_info_str)

b1_cache <- list(run_parameters = run_parameters_b1, gene_collated = gene_collated, gene_roc_df = gene_roc_df, gene_pr_df = gene_pr_df, gene_summary_table = gene_summary_table, all_method_scatter_data = all_method_scatter_data, df_all_ed_raw = df_all_ed_raw, ed_summary_table = ed_summary_table, df_all_pb_raw = df_all_pb_raw, pb_summary_table = pb_summary_table, matrix_roc_all = matrix_roc_all, matrix_pr_all = matrix_pr_all, matrix_summary = matrix_summary, plot_roc_df = plot_roc_df, plot_pr_df = plot_pr_df, plot_ed_df = plot_ed_df)

b1_cache_file <- file.path(export_dir, paste0("Eval_Cache_B1_", run_uid, ".rds"))
saveRDS(b1_cache, file = b1_cache_file)
message(sprintf("  ✅ B1 evaluation results and plotting data fully archived to: %s", b1_cache_file))
