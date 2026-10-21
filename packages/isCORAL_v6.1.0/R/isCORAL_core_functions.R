# ==============================================================================
# --- Section 0: Helper Functions (Injected into Background Workers) ---
# ==============================================================================

#' @export
CalculateGeneOmegaSquared <- function(expression_matrix, 
                                      groups, 
                                      chunk_size = 2000, 
                                      verbose = FALSE) {
  valid_idx <- which(!is.na(groups))
  n_cells_valid <- length(valid_idx)
  
  if (n_cells_valid < 2) {
    if(verbose) message("Warning: Too few valid cells (n<2). Returning NAs.")
    return(setNames(rep(NA, nrow(expression_matrix)), rownames(expression_matrix)))
  }
  
  if (n_cells_valid == ncol(expression_matrix)) {
      expr_valid <- expression_matrix
      groups_valid <- droplevels(as.factor(groups))
  } else {
      expr_valid <- expression_matrix[, valid_idx, drop = FALSE]
      groups_valid <- droplevels(as.factor(groups[valid_idx]))
  }
  
  n_genes <- nrow(expr_valid)
  n_cells <- ncol(expr_valid)
  n_groups <- length(levels(groups_valid))
  
  df_between <- n_groups - 1
  df_within <- n_cells - n_groups
  
  if (n_groups <= 1 || df_within <= 0) {
    if(verbose) message("Warning: Not enough groups or degrees of freedom. Returning 0s.")
    return(setNames(rep(0, n_genes), rownames(expr_valid)))
  }
  
  design_mat <- Matrix::sparse.model.matrix(~ 0 + groups_valid)
  n_k <- Matrix::colSums(design_mat) 
  
  global_sum_all <- Matrix::rowSums(expr_valid)
  global_mean_all <- global_sum_all / n_cells
  
  num_chunks <- ceiling(n_genes / chunk_size)
  omega_sq_list <- vector("list", num_chunks)
  
  for (i in 1:num_chunks) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, n_genes)
    
    expr_chunk <- expr_valid[start_row:end_row, , drop = FALSE]
    global_mean_chunk <- global_mean_all[start_row:end_row]
    
    group_sums <- expr_chunk %*% design_mat 
    group_sums_dense <- as.matrix(group_sums)
    
    group_means <- group_sums_dense / rep(n_k, each = nrow(group_sums_dense))
    
    diff_means <- group_means - global_mean_chunk
    weighted_sq_diff <- (diff_means^2) * rep(n_k, each = nrow(diff_means))
    ssb_chunk <- rowSums(weighted_sq_diff)
    
    chunk_sq <- expr_chunk
    chunk_sq@x <- chunk_sq@x^2
    sum_x_sq <- Matrix::rowSums(chunk_sq)
    
    sst_chunk <- sum_x_sq - (n_cells * global_mean_chunk^2)
    sst_chunk[sst_chunk < 0] <- 0 
    
    ssw_chunk <- sst_chunk - ssb_chunk
    ssw_chunk[ssw_chunk < 0] <- 0
    
    msw_chunk <- ssw_chunk / df_within
    numerator <- ssb_chunk - (df_between * msw_chunk)
    denominator <- sst_chunk + msw_chunk
    
    omega_sq_chunk <- numerator / denominator
    omega_sq_chunk[denominator == 0 | !is.finite(omega_sq_chunk)] <- 0
    
    omega_sq_list[[i]] <- omega_sq_chunk
    
    rm(expr_chunk, group_sums, group_sums_dense, group_means, diff_means, weighted_sq_diff, chunk_sq, sum_x_sq, sst_chunk, ssw_chunk, numerator, denominator, omega_sq_chunk)
  }
  
  final_omega_sq <- do.call(c, omega_sq_list)
  names(final_omega_sq) <- rownames(expr_valid)
  
  rm(expr_valid, design_mat, global_sum_all, global_mean_all)
  return(final_omega_sq)
}

# ------------------------------------------------------------------------------
# RSpectra Clustering Wrapper
# ------------------------------------------------------------------------------
get_cell_bins_IN_MEMORY_ROBUST <- function(expression_cache, feature_genes, k, gene_names_vec, kernel = "laplacedot", cell_indices_to_use = NULL, fallback_sigma = 0.1) {
  
  if (!is.null(cell_indices_to_use)) {
    expression_cache <- expression_cache[, cell_indices_to_use, drop = FALSE]
  }
  num_cells_in_run <- ncol(expression_cache)
  
  local_k <- k
  if (local_k >= num_cells_in_run) local_k <- max(2, num_cells_in_run - 1) 
  if (local_k <= 1 && num_cells_in_run > 1) local_k <- 2
  
  result <- tryCatch({
    if (num_cells_in_run <= 1) stop("Cannot cluster 1 or fewer cells.")
    feature_genes_intersect <- intersect(feature_genes, gene_names_vec)
    if (length(feature_genes_intersect) < 2) stop("Error: Fewer than 2 feature genes found.")
    
    data_subset_local <- expression_cache[feature_genes_intersect, , drop = FALSE]
    data_for_specc <- as.matrix(t(data_subset_local)) 
    data_for_specc[!is.finite(data_for_specc)] <- 0
    
    estimated_sigma_raw <- tryCatch({
      kernlab::sigest(data_for_specc, scaled = FALSE)["0.5"]
    }, error = function(e) fallback_sigma)
    sigma_to_use <- as.numeric(estimated_sigma_raw)
    if (is.null(sigma_to_use) || !is.finite(sigma_to_use) || sigma_to_use <= 1e-9) sigma_to_use <- fallback_sigma
    
    kernel_function_obj <- NULL
    if (is.character(kernel)) {
      if (kernel == "vanilladot") { kernel_function_obj <- kernlab::vanilladot()
      } else if (kernel == "rbfdot") { kernel_function_obj <- kernlab::rbfdot(sigma = sigma_to_use)
      } else if (kernel == "polydot") { kernel_function_obj <- kernlab::polydot(degree = 2, scale = 1, offset = 1)
      } else if (kernel == "laplacedot") { kernel_function_obj <- kernlab::laplacedot(sigma = sigma_to_use)
      } else if (kernel == "tanhdot") { kernel_function_obj <- kernlab::tanhdot(scale = 1, offset = 1)
      } else {
        kernel_function_obj <- kernlab::rbfdot(sigma = sigma_to_use)
      }
    } else {
      kernel_function_obj <- kernel
    }
    
    W <- kernlab::kernelMatrix(kernel_function_obj, x = data_for_specc)
    if (inherits(W, "kernelMatrix")) W <- W@.Data 
    
    D <- rowSums(W)
    D_inv_sqrt <- 1 / sqrt(D + 1e-9)
    
    L_sym <- t(t(W * D_inv_sqrt) * D_inv_sqrt)
    rm(W) 
    
    U_norm <- NULL
    
    k_calc <- min(nrow(L_sym), max(30, ceiling(local_k * 2.5)))
    if (k_calc >= nrow(L_sym)) k_calc <- nrow(L_sym) - 1
    res <- RSpectra::eigs_sym(L_sym, k = k_calc, which = "LM", tol = 1e-13, maxitr = 5000)
    U <- res$vectors[, 1:local_k, drop = FALSE]
    
    rm(L_sym, data_subset_local, data_for_specc) 
    
    row_norms <- sqrt(rowSums(U^2))
    U_norm <- U / (row_norms + 1e-9)
    
    km_res <- stats::kmeans(U_norm, centers = local_k, nstart = 100, iter.max = 100)
    
    return(as.factor(km_res$cluster))
  }, error = function(e) return(rep(NA, num_cells_in_run)))
  
  return(result)
}

# ------------------------------------------------------------------------------
# RcppML Wrapper
# ------------------------------------------------------------------------------
run_gemli_nmf_RcppML_ROBUST <- function(input_matrix_subset, feature_n, k_val, nrun = 1, threads = 1) { 
  best_result <- NULL; best_rss <- Inf 
  
  if (inherits(input_matrix_subset, "dgCMatrix")) {
      mat_for_nmf <- input_matrix_subset
      mat_for_nmf@x[is.na(mat_for_nmf@x) | is.infinite(mat_for_nmf@x)] <- 0
      mat_for_nmf@x[mat_for_nmf@x < 0] <- 0
      max_val <- max(mat_for_nmf@x)
  } else {
      mat_for_nmf <- as.matrix(input_matrix_subset)
      mat_for_nmf[is.na(mat_for_nmf) | is.infinite(mat_for_nmf)] <- 0
      mat_for_nmf[mat_for_nmf < 0] <- 0 
      max_val <- max(mat_for_nmf)
  }
  
  if (max_val > 50) {
      if (inherits(mat_for_nmf, "dgCMatrix")) mat_for_nmf@x <- log1p(mat_for_nmf@x)
      else mat_for_nmf <- log1p(mat_for_nmf)
  }

  if (inherits(mat_for_nmf, "dgCMatrix")) {
      col_norms <- sqrt(Matrix::colSums(mat_for_nmf^2))
      col_norms[col_norms == 0] <- 1
      mat_for_nmf <- mat_for_nmf %*% Matrix::Diagonal(x = 1/col_norms)
  } else {
      col_norms <- sqrt(colSums(mat_for_nmf^2))
      col_norms[col_norms == 0] <- 1
      mat_for_nmf <- mat_for_nmf / rep(col_norms, each = nrow(mat_for_nmf))
  }

  for (i in 1:nrun) {
    nmf_result <- RcppML::nmf(mat_for_nmf, k = k_val, verbose = FALSE, L1 = 0, L2 = 0, threads = threads) 
    
    W_mat <- nmf_result$w
    H_mat <- nmf_result$h
    current_rss <- 0
    for (c in 1:ncol(mat_for_nmf)) {
        diff_col <- mat_for_nmf[, c] - (W_mat %*% H_mat[, c, drop = FALSE])
        current_rss <- current_rss + sum(diff_col^2)
    }
    
    if (current_rss < best_rss) { 
        best_rss <- current_rss
        best_result <- nmf_result 
    }
    rm(W_mat, H_mat, diff_col) 
  }
  
  loadings_matrix <- best_result$w
  rownames(loadings_matrix) <- rownames(input_matrix_subset)
  
  feature_genes <- apply(loadings_matrix, 2, function(x) {
    rownames(loadings_matrix)[order(x, decreasing = TRUE)[1:feature_n]]
  })
  
  return(feature_genes)
}

# ------------------------------------------------------------------------------
# Consensus Matrix Computation
# ------------------------------------------------------------------------------
run_simple_consensus <- function(intermediate_results, num_cells, cell_names) {
  valid_intermediate_results <- Filter(function(res) {
    !is.null(res$final_bins) && length(res$final_bins) == num_cells &&
      !any(is.na(res$final_bins)) && length(unique(res$final_bins[!is.na(res$final_bins)])) > 1
  }, intermediate_results)
  
  if (length(valid_intermediate_results) < 2) stop("Fewer than 2 fully valid intermediate runs found.", call. = FALSE)
  
  indicator_list <- lapply(valid_intermediate_results, function(res) {
    fb <- res$final_bins
    names(fb) <- cell_names
    Matrix::sparse.model.matrix(~ 0 + fb)
  })
  
  combined_indicators <- do.call(cbind, indicator_list)
  rm(indicator_list)
  
  co_assoc_matrix <- Matrix::tcrossprod(combined_indicators)
  rm(combined_indicators)
  
  S_stab <- co_assoc_matrix / length(valid_intermediate_results)
  diag(S_stab) <- 1
  rm(co_assoc_matrix)

  return(list(S_stab = S_stab))
}

# ==============================================================================
# --- Section 3: Main isCORAL Function ---
# ==============================================================================
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom RSpectra eigs_sym
#' @importFrom matrixStats colRanks 
#' @importFrom kernlab kernelMatrix rbfdot sigest specc
#' @importFrom stats kmeans
#' @import RcppML
#' @import callr
#' @importFrom Matrix Matrix sparse.model.matrix t tcrossprod
isCORAL_v6_1_0 <- function(seurat_obj,
                           assay = "RNA",
                           results_dir = "iscoral_results_raw",
                           species = "mouse",
                           random_seed = NULL,
                           downsampleUMI = NULL,
                           bin_size = 50,
                           n_repeats = 100,
                           parallel_cores = 30,
                           max_cells = NULL,
                           gemli_mean_q1 = 0.98,
                           gemli_mean_q2 = 0.90,
                           gemli_fano_q2 = 0.40,
                           gemli_mean_q3 = 0.60,
                           gemli_fano_q3 = 0.90,
                           genes_to_remove = NULL,
                           gene_min_mean_for_ranking = 0.05,
                           nmf_feature_selection_ratio = 0.2,
                           nmf_nrun = 1, 
                           nmf_feature_number = 5,
                           nmf_k_value = 50,
                           downsample_gemli_ratio = 0.75,
                           specc_kernel = "laplacedot",
                           nmf_eval_downsample_cells = 2000,
                           custom_log_file = NULL) {
  
  start_time_total <- Sys.time()
  global_run_uid <- paste0("iscoral_opt_", format(start_time_total, "%Y%m%d%H%M%S"), "_",
                           paste(sample(c(letters, 8:9), 8), collapse = ""))
  
  # ----------------------------------------------------------------------------
  # Initialize Logging System
  # ----------------------------------------------------------------------------
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  log_dir <- file.path(results_dir, "logs") 
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  if (!is.null(custom_log_file)) {
      progress_log_file <- custom_log_file 
      if(!dir.exists(dirname(progress_log_file))) dir.create(dirname(progress_log_file), recursive = TRUE)
  } else {
      progress_log_file <- normalizePath(file.path(log_dir, paste0(global_run_uid, "_progress.log")), mustWork = FALSE)
  }
  
  cat(sprintf("=== isCORAL HPC Run Started at %s ===\n", start_time_total), file = progress_log_file)
  cat(sprintf("Target: %d repeats, %d cores. Wait for workers to report...\n", n_repeats, parallel_cores), file = progress_log_file, append = TRUE)
  
  message(paste("--- [isCORAL v6.1.0 CALLR HPC EDITION] Run UID:", global_run_uid, "---"))
  message(">>> REAL-TIME PROGRESS LOG: ", progress_log_file, " <<<")
  message(">>> TIP: In your Linux terminal, run: tail -f ", progress_log_file, " <<<")
  # ----------------------------------------------------------------------------

  run_parameters <- as.list(environment())
  run_parameters$seurat_obj <- NULL
  
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
  seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)
  
  if (!is.null(downsampleUMI) && is.numeric(downsampleUMI)) {
    message(paste("Performing UMI downsampling to", downsampleUMI, "..."))
    original_counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
    if (is.null(original_counts) || prod(dim(original_counts)) == 0) stop("No valid 'counts'.")
    downsampled_counts <- Seurat::SampleUMI(original_counts, max.umi = downsampleUMI, upsample = FALSE, verbose = FALSE)
    seurat_obj[[assay]]$counts <- downsampled_counts
  }
  
  if (!is.null(max_cells) && ncol(seurat_obj) > max_cells) {
    message(paste("Downsampling to", max_cells, "cells."))
    seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), size = max_cells, replace = FALSE))
  }
  num_cells_for_analysis <- ncol(seurat_obj)
  cell_names_for_analysis <- colnames(seurat_obj)
  
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  sc_expression_full <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data") 
  if (!inherits(sc_expression_full, "dgCMatrix")) {
    sc_expression_full <- methods::as(sc_expression_full, "dgCMatrix")
  }
  
  message("Calculating gene statistics...")
  the_means <- Matrix::rowMeans(sc_expression_full)
  
  temp_sq <- sc_expression_full
  temp_sq@x <- temp_sq@x^2
  the_sst <- Matrix::rowSums(temp_sq) - (ncol(sc_expression_full) * the_means^2)
  the_sst[the_sst < 0] <- 0
  rm(temp_sq) 
  
  expression_stats_base <- data.frame(mean = the_means, sst = the_sst, name = rownames(sc_expression_full))
  expression_stats_base$fano <- (expression_stats_base$sst / ncol(sc_expression_full)) / expression_stats_base$mean
  expression_stats_base[is.na(expression_stats_base)] <- 0
  
  if (tolower(species) == "human") { mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
  } else if (tolower(species) == "mouse") { mt_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE)
  } else { mt_genes <- c() }
  
  s_genes_human <- Seurat::cc.genes$s.genes; g2m_genes_human <- Seurat::cc.genes$g2m.genes
  if (tolower(species) == "mouse") {
    to_mouse_case <- function(genes) { return(paste0(substr(genes, 1, 1), tolower(substr(genes, 2, nchar(genes))))) }
    s_genes <- to_mouse_case(s_genes_human); g2m_genes <- to_mouse_case(g2m_genes_human)
  } else { s_genes <- s_genes_human; g2m_genes <- g2m_genes_human }
  
  combined_genes_to_remove <- unique(c(mt_genes, s_genes, g2m_genes, "MALAT1", genes_to_remove))
  gene_names_vec <- rownames(sc_expression_full)
  
  gemli_mean_cond1 <- expression_stats_base$mean >= quantile(expression_stats_base$mean, gemli_mean_q1)
  gemli_mean_cond2 <- expression_stats_base$mean >= quantile(expression_stats_base$mean, gemli_mean_q2) & 
                      expression_stats_base$fano >= quantile(expression_stats_base$fano, gemli_fano_q2)
  gemli_mean_cond3 <- expression_stats_base$mean >= quantile(expression_stats_base$mean, gemli_mean_q3) & 
                      expression_stats_base$fano >= quantile(expression_stats_base$fano, gemli_fano_q3)
                      
  gemli_genes_unfiltered <- rownames(expression_stats_base[gemli_mean_cond1 | gemli_mean_cond2 | gemli_mean_cond3, ])
  
  gemli_genes <- setdiff(gemli_genes_unfiltered, combined_genes_to_remove)
  
  ss_filter <- subset(expression_stats_base, mean > gene_min_mean_for_ranking)
  
  bin_size_values <- bin_size
  all_scan_results <- list(run_parameters = run_parameters)
  
  message("Packaging payload to avoid excessive memory usage...")
  
  helper_funcs <- list(
    CalculateGeneOmegaSquared = CalculateGeneOmegaSquared,
    get_cell_bins_IN_MEMORY_ROBUST = get_cell_bins_IN_MEMORY_ROBUST,
    run_gemli_nmf_RcppML_ROBUST = run_gemli_nmf_RcppML_ROBUST
  )
  
  payload <- list(
    expr_matrix = sc_expression_full,
    gemli_genes = gemli_genes,
    gene_names_vec = gene_names_vec,
    ss_filter = ss_filter,
    funcs = helper_funcs,
    random_seed = random_seed,
    progress_log = progress_log_file, 
    params = list(
      n_tasks = n_repeats,
      parallel_cores = parallel_cores,
      bin_size_values = bin_size_values,
      num_cells = num_cells_for_analysis,
      specc_kernel = specc_kernel,
      nmf_feature_number = nmf_feature_number,
      nmf_k_value = nmf_k_value,
      nmf_nrun = nmf_nrun,
      downsample_gemli_ratio = downsample_gemli_ratio,
      nmf_eval_downsample_cells = nmf_eval_downsample_cells,
      nmf_feature_selection_ratio = nmf_feature_selection_ratio
    )
  )

  message(paste("Launching isolated callr backend with", payload$params$parallel_cores, "cores. Main process standing by..."))
  
  # ============================================================================
  # Launch callr Background Process
  # ============================================================================
  background_results <- callr::r(function(data) {
    
    list2env(data$funcs, envir = globalenv())
    
    suppressPackageStartupMessages({
      library(Matrix)
      library(doMC)
      library(foreach)
      library(RSpectra)
      library(kernlab)
      library(matrixStats)
      library(RcppML)
      library(doRNG)
    })
    
    # ----------------------------------------------------------------------
    # Suppress underlying implicit multi-threading 
    # to prevent thread collisions with doMC
    # ----------------------------------------------------------------------
    Sys.setenv(OMP_NUM_THREADS = "1")
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    Sys.setenv(MKL_NUM_THREADS = "1")
    Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
    Sys.setenv(NUMEXPR_NUM_THREADS = "1")
    Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
    
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
    # ----------------------------------------------------------------------
    
    registerDoMC(cores = data$params$parallel_cores)
    if (!is.null(data$random_seed)) doRNG::registerDoRNG(data$random_seed)
    
    results_internal <- list()
    
    for (current_bin_size in data$params$bin_size_values) {
      num_bins <- max(2, floor(data$params$num_cells / current_bin_size))
      
      iteration_results <- foreach::foreach(
        i = 1:data$params$n_tasks,
        .errorhandling = "pass"
      ) %dopar% {
        
        gemli_subset <- sample(data$gemli_genes, floor(length(data$gemli_genes) * data$params$downsample_gemli_ratio))
        nmf_gene_indices <- which(data$gene_names_vec %in% gemli_subset)
        nmf_matrix <- data$expr_matrix[nmf_gene_indices, , drop = FALSE]
        
        NMF_features <- run_gemli_nmf_RcppML_ROBUST(
          nmf_matrix, feature_n = data$params$nmf_feature_number, 
          k_val = data$params$nmf_k_value, nrun = data$params$nmf_nrun, threads = 1 
        )
        
        eval_cell_indices <- NULL
        k_to_use <- num_bins
        if (!is.null(data$params$nmf_eval_downsample_cells) && is.numeric(data$params$nmf_eval_downsample_cells)) {
           n_to_sample <- min(data$params$nmf_eval_downsample_cells, ncol(data$expr_matrix))
           eval_cell_indices <- sample(1:ncol(data$expr_matrix), size = n_to_sample, replace = FALSE)
           k_to_use <- max(2, min(num_bins, floor(n_to_sample / 5)))
        }
        
        omega_eval_matrix <- matrix(NA, nrow = length(data$gene_names_vec), ncol = ncol(NMF_features))
        rownames(omega_eval_matrix) <- data$gene_names_vec
        
        for (j in 1:ncol(NMF_features)) {
             bins_temp_raw <- get_cell_bins_IN_MEMORY_ROBUST(
               data$expr_matrix, NMF_features[, j], k_to_use, data$gene_names_vec, data$params$specc_kernel, 
               cell_indices_to_use = eval_cell_indices
             )
             
             full_bins_vector <- rep(NA, ncol(data$expr_matrix))
             if (!is.null(eval_cell_indices)) {
                 if (length(bins_temp_raw) == length(eval_cell_indices)) full_bins_vector[eval_cell_indices] <- bins_temp_raw
             } else {
                 full_bins_vector <- bins_temp_raw
             }
             
             omega_eval_matrix[, j] <- CalculateGeneOmegaSquared(data$expr_matrix, as.factor(full_bins_vector))
             rm(bins_temp_raw, full_bins_vector)
        }
        
        valid_genes_rank <- intersect(rownames(omega_eval_matrix), data$ss_filter$name)
        omega_subset <- omega_eval_matrix[valid_genes_rank, , drop=FALSE]
        omega_subset[is.na(omega_subset)] <- 0
        
        ranks_matrix <- matrixStats::colRanks(-omega_subset, preserveShape = TRUE)
        rownames(ranks_matrix) <- valid_genes_rank
        
        eval_genes <- intersect(valid_genes_rank, gemli_subset)
        final_feature_genes <- unique(as.vector(NMF_features))
        
        if(length(eval_genes) > 0) {
             mod_scores <- colMeans(ranks_matrix[eval_genes, , drop = FALSE], na.rm = TRUE)
             n_keep <- max(2, floor(data$params$nmf_feature_selection_ratio * length(mod_scores)))
             best_indices <- order(mod_scores)[1:n_keep]
             final_feature_genes <- unique(as.vector(NMF_features[, best_indices]))
        }

        # Perform final spectral clustering
        final_bins <- get_cell_bins_IN_MEMORY_ROBUST(
           data$expr_matrix, final_feature_genes, num_bins, data$gene_names_vec, data$params$specc_kernel,
           cell_indices_to_use = NULL 
        )
        final_omega <- CalculateGeneOmegaSquared(data$expr_matrix, final_bins)
        
        res_out <- list(final_bins = final_bins, omega_squared_run = final_omega)
        
        # ----------------------------------------------------------------------
        # Timestamp and Progress Logging
        # ----------------------------------------------------------------------
        cat(sprintf("[%s] [Worker PID: %d] SUCCESS - Finished Repeat %d / %d (bin_size: %d)\n", 
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                    Sys.getpid(), 
                    i, 
                    data$params$n_tasks, 
                    current_bin_size), 
            file = data$progress_log, append = TRUE)
        # ----------------------------------------------------------------------

        # Memory cleanup
        rm(list = setdiff(ls(), "res_out"))
        
        return(res_out)
      }
      results_internal[[paste0("binsize_", current_bin_size)]] <- iteration_results
    }
    
    return(results_internal)
    
  }, args = list(data = payload), show = TRUE) 
  # ============================================================================
  
  message("\nBackground processing complete. Integrating results into Seurat...")
  
  rm(payload, helper_funcs)
  gc(verbose = FALSE, reset = TRUE)

  for (current_bin_size in bin_size_values) {
    iteration_results <- background_results[[paste0("binsize_", current_bin_size)]]
    
    valid_results <- Filter(function(x) !is.null(x) && !inherits(x, "error") && !is.null(x$final_bins) && !all(is.na(x$final_bins)), iteration_results)
    
    if (length(valid_results) > 1) {
      consensus_results <- run_simple_consensus(
        intermediate_results = valid_results,
        num_cells = num_cells_for_analysis,
        cell_names = cell_names_for_analysis
      )
      
      sum_omega <- rep(0, length(gene_names_vec))
      count_omega <- rep(0, length(gene_names_vec))
      
      for(res in valid_results) {
         if(!is.null(res$omega_squared_run)) {
            clean_res <- res$omega_squared_run
            clean_res[is.na(clean_res)] <- 0
            sum_omega <- sum_omega + clean_res
            count_omega <- count_omega + 1
         }
      }
      avg_omega <- sum_omega / pmax(1, count_omega)
      names(avg_omega) <- gene_names_vec
      avg_omega[avg_omega < 0] <- 0
      scores_df <- cbind(expression_stats_base, omega_square_infer = avg_omega[expression_stats_base$name])
      
      all_scan_results[[paste0("binsize_", current_bin_size)]] <- list(
        ss_filter_average = scores_df,
        S_stab = consensus_results$S_stab
      )
    } else {
      warning(paste("Not enough valid runs for bin_size", current_bin_size))
    }
  }
  
  seurat_obj@misc$CORAL_scan_results <- all_scan_results
  total_time_secs <- as.numeric(difftime(Sys.time(), start_time_total, units = "secs"))
  message(paste("\nFull analysis scan complete. Total time:", round(total_time_secs, 2), "seconds")) 
  
  suppressWarnings({
    if (exists("sc_expression_full")) rm(sc_expression_full)
    if (exists("background_results")) rm(background_results)
    if (exists("all_scan_results")) rm(all_scan_results)
    if (exists("the_sst")) rm(the_sst)
    if (exists("the_means")) rm(the_means)
    if (exists("expression_stats_base")) rm(expression_stats_base)
  })
  gc(verbose = FALSE, reset = TRUE) 
  
  return(seurat_obj)
}

# ==============================================================================
# --- Section 4: Downstream Utilities ---
# ==============================================================================
#' Simulate Meta-Clones from isCORAL consensus matrix
#'
#' @param seurat_obj A Seurat object containing isCORAL scan results.
#' @param bin_size The bin size used in the isCORAL scan.
#' @param k_values A numeric vector of k values for cutting the hierarchical tree.
#' @param hclust_method The agglomeration method to be used (default is "ward.D2").
#' @return A modified Seurat object with meta-clone labels added to metadata.
#' @export
SimulateClone <- function(seurat_obj, bin_size, k_values, hclust_method = "ward.D2") {
  
  # Check if isCORAL results exist in the object
  if (!"CORAL_scan_results" %in% names(seurat_obj@misc)) {
    stop("Error: isCORAL results not found in seurat_obj@misc$CORAL_scan_results.")
  }
  
  result_name <- paste0("binsize_", bin_size)
  if (!result_name %in% names(seurat_obj@misc$CORAL_scan_results)) {
    stop(paste("Error: No isCORAL results found for bin_size =", bin_size))
  }
  
  coral_result <- seurat_obj@misc$CORAL_scan_results[[result_name]]
  
  if (!"S_stab" %in% names(coral_result)) {
    stop(paste("Error: 'S_stab' consensus matrix not found for bin_size =", bin_size))
  }
  
  message("Computing hierarchical clustering from consensus matrix (S_stab)...")
  S_mat <- as.matrix(coral_result$S_stab)
  
  # Automatically fill in cell names if the matrix lacks them but dimensions match the Seurat object
  if (is.null(rownames(S_mat)) || length(intersect(rownames(S_mat), colnames(seurat_obj))) == 0) {
    if (nrow(S_mat) == ncol(seurat_obj) && ncol(S_mat) == ncol(seurat_obj)) {
      rownames(S_mat) <- colnames(seurat_obj)
      colnames(S_mat) <- colnames(seurat_obj)
    } else {
      stop("Error: S_stab dimensions do not match the Seurat object, and no matching cell names were found.")
    }
  }

  dist_matrix <- stats::as.dist(1 - S_mat)
  hclust_obj <- stats::hclust(dist_matrix, method = hclust_method)
  
  message(paste("Simulating Meta-Clones from bin_size =", bin_size, "at k =", paste(k_values, collapse = ", ")))
  
  for (k in k_values) {
    clusters <- as.factor(stats::cutree(hclust_obj, k = as.integer(k)))
    
    # Forcibly assign cell names to clusters to pass Seurat safety checks
    names(clusters) <- rownames(S_mat)
    
    col_name <- paste0("metaclone_k", k, "_bs", bin_size)
    seurat_obj[[col_name]] <- clusters
  }
  
  return(seurat_obj)
}