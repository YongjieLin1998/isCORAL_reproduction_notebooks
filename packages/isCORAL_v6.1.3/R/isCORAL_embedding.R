#' @title Unified Integrated isCORAL Dual-Layer Manifold Embedding Engine (CORALembedding)
#' @description 
#' Directly fits the latent space of pseudoclones using raw transcriptomic features, 
#' then applies full-weight graph diffusion at the output space.
#' Permanently mounts "CORAL_cell_embed" (pure gene cell space) and "CORAL_clone_embed" (fully diffused clone space),
#' while automatically caching the generated normalized graph probability, MDS scaffold, and initial projections 
#' into @misc$isCORAL_embedding for persistent reuse.
#' 
#' @param seurat_obj A preprocessed (Normalize/Scale) Seurat object compatible with V4 and V5.
#' @param binsize_key Character key for locating the graph transition matrix inside @misc. Default is "binsize_200".
#' @param n_dims Target dimensionality of the latent manifold space (typically 20 or 30).
#' @param top_n_genes Number of core driving genes filtered by inferred Omega-Square ranking.
#' @param max_cells_lasso Maximum number of cells subsampled for training the regression mapper. 
#'        Default is `NULL` (subsampling disabled, uses all cells). Set to an integer (e.g., 3000) to enable subsampling for very large datasets to prevent memory exhaustion.
#' @param alpha_elasticnet Elastic net penalty coefficient. Default is 0.05.
#' @param workers_glmnet Number of parallel threads used during cross-validation training, hard-capped at 5 internally.
#' 
#' @return A Seurat object containing two standard dimensional reductions and the persistent @misc topological cache.
#' @export
CORALembedding <- function(
    seurat_obj, 
    binsize_key = "binsize_200", 
    n_dims = 30, 
    top_n_genes = 3000, 
    max_cells_lasso = NULL,          
    alpha_elasticnet = 0.05, 
    workers_glmnet = 5
) {
  old_omp <- Sys.getenv("OMP_NUM_THREADS")
  Sys.setenv(OMP_NUM_THREADS = 1) 
  on.exit({ Sys.setenv(OMP_NUM_THREADS = old_omp) }, add = TRUE)
  
  message("🚀 [isCORAL] Initiating unified manifold embedding pipeline...")
  
  coral_misc_key <- ifelse("CORAL_scan_results" %in% names(seurat_obj@misc), "CORAL_scan_results", "coral_scan_results")
  coral_res <- seurat_obj@misc[[coral_misc_key]][[binsize_key]]
  if (is.null(coral_res)) stop(sprintf("🚨 Fatal Error: Unable to locate network scan results for key '%s' inside @misc!", binsize_key))
  
  prob_matrix <- coral_res$S_stab
  ss_stats <- coral_res$ss_filter_average
  
  E_mat <- tryCatch({ Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data") }, 
                    error = function(e) { Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data") })
  
  seurat_cells <- colnames(E_mat)
  N_cells_raw <- length(seurat_cells)
  
  if (is.null(rownames(prob_matrix)) || length(intersect(rownames(prob_matrix), seurat_cells)) == 0) {
    if (nrow(prob_matrix) == N_cells_raw && ncol(prob_matrix) == N_cells_raw) {
      rownames(prob_matrix) <- colnames(prob_matrix) <- seurat_cells
    } else {
      stop("🚨 Dimension mismatch between expression matrix and probabilistic matrix! Please align cell barcodes beforehand.")
    }
  }
  
  valid_cells <- intersect(seurat_cells, rownames(prob_matrix))
  N_cells <- length(valid_cells)
  if (N_cells == 0) stop("🚨 Fatal Error: No valid overlapping cells found between expression matrix and prob_matrix.")
  
  ss_stats_sorted <- ss_stats[order(ss_stats$omega_square_infer, decreasing = TRUE), , drop = FALSE]
  top_genes <- intersect(rownames(ss_stats_sorted)[1:min(top_n_genes, nrow(ss_stats_sorted))], rownames(E_mat))
  
  E_mat <- E_mat[top_genes, valid_cells, drop = FALSE]
  
  P <- prob_matrix[valid_cells, valid_cells, drop = FALSE]
  if (!inherits(P, "CsparseMatrix")) {
    P <- as(P, "CsparseMatrix")
  }
  
  row_sums <- Matrix::rowSums(P)
  row_sums[row_sums == 0] <- 1  
  P_norm <- Matrix::Diagonal(x = 1 / row_sums) %*% P
  
  v_norms <- Matrix::colSums(E_mat^2)
  chunk_size <- ifelse(N_cells < 15000, 5000, 2000) 
  starts <- seq(1, N_cells, by = chunk_size)
  E_cross <- matrix(0, nrow = N_cells, ncol = N_cells)
  
  for (s in starts) {
    idx <- s:min(s + chunk_size - 1, N_cells)
    cross_chunk <- as.matrix(Matrix::crossprod(E_mat, E_mat[, idx, drop = FALSE])) 
    D_sq_chunk <- -2 * cross_chunk + v_norms + rep(v_norms[idx], each = N_cells)
    D_sq_chunk[D_sq_chunk < 0] <- 0 
    E_cross <- E_cross + as.matrix(P_norm %*% sqrt(D_sq_chunk)) %*% as.matrix(t(P_norm[, idx, drop = FALSE]))
  }
  
  E_in <- diag(E_cross)
  Dist_E_MMD2 <- 2 * E_cross
  Dist_E_MMD2 <- sweep(Dist_E_MMD2, 1, E_in, "-")
  Dist_E_MMD2 <- sweep(Dist_E_MMD2, 2, E_in, "-")
  Dist_E_MMD2[Dist_E_MMD2 < 0] <- 0 
  rm(E_cross)
  
  row_mean <- rowMeans(Dist_E_MMD2)
  mean_all <- mean(row_mean)
  B_mat <- -0.5 * (sweep(sweep(Dist_E_MMD2, 1, row_mean, "-"), 2, row_mean, "-") + mean_all)
  
  safe_n_dims <- min(n_dims, N_cells - 5)
  if (safe_n_dims < 2) stop("🚨 Not enough cells to perform dimensional reduction. N_cells too small.")
  
  svd_res <- irlba::irlba(B_mat, nv = safe_n_dims, fastpath = TRUE, tol = 1e-4)
  
  original_MDS_coords <- sweep(svd_res$u, 2, sqrt(pmax(svd_res$d[1:safe_n_dims], 0)), "*")
  rownames(original_MDS_coords) <- valid_cells
  
  rm(B_mat, Dist_E_MMD2, svd_res); gc(verbose = FALSE)
  
  X_raw <- t(as.matrix(E_mat))
  train_cell_norms <- sqrt(rowSums(X_raw^2) + 1e-8)
  X_l2 <- sweep(X_raw, 1, train_cell_norms, "/")
  ref_gene_means <- matrixStats::colMeans2(X_l2)
  ref_gene_sds <- matrixStats::colSds(X_l2)
  ref_gene_sds[ref_gene_sds == 0] <- 1e-8
  X_scaled <- as.matrix(scale(X_l2, center = ref_gene_means, scale = ref_gene_sds))
  X_scaled[is.na(X_scaled)] <- 0
  
  do_subsample <- TRUE
  if (is.null(max_cells_lasso) || is.infinite(max_cells_lasso)) {
    do_subsample <- FALSE
  } else if (nrow(X_scaled) <= max_cells_lasso) {
    do_subsample <- FALSE
  }
  
  if (do_subsample) {
    train_idx <- base::sample(1:nrow(X_scaled), max_cells_lasso)
    message(sprintf("📉 [isCORAL] Subsampling enabled: Randomly selected %d cells for ElasticNet training.", max_cells_lasso))
  } else {
    train_idx <- 1:nrow(X_scaled)
    message(sprintf("📈 [isCORAL] Subsampling bypassed: Using all %d cells for ElasticNet training.", length(train_idx)))
  }
  
  actual_workers <- min(workers_glmnet, 5, parallel::detectCores() - 1) 
  cl <- parallel::makeCluster(actual_workers)
  doParallel::registerDoParallel(cl)
  
  cv_fit <- tryCatch({
    glmnet::cv.glmnet(x = X_scaled[train_idx, ], y = original_MDS_coords[train_idx, ], 
                      family = "mgaussian", alpha = alpha_elasticnet, nfolds = 5, parallel = TRUE)
  }, finally = {
    parallel::stopCluster(cl)
  })
  
  Y_ref_raw <- predict(cv_fit, newx = X_scaled, s = "lambda.min")[, , 1]
  Y_ref_raw <- sweep(Y_ref_raw, 2, matrixStats::colSds(original_MDS_coords) / pmax(matrixStats::colSds(Y_ref_raw), 1e-8), "*")
  rownames(Y_ref_raw) <- valid_cells
  colnames(Y_ref_raw) <- paste0("CORALcell_", 1:safe_n_dims)
  
  seurat_obj[["CORAL_cell_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = Y_ref_raw, 
    key = "CORALcell_", 
    assay = Seurat::DefaultAssay(seurat_obj)
  )
  
  Y_smoothed <- as.matrix(P_norm %*% Y_ref_raw)
  rownames(Y_smoothed) <- valid_cells
  colnames(Y_smoothed) <- paste0("CORALclone_", 1:safe_n_dims)
  
  seurat_obj[["CORAL_clone_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = Y_smoothed, 
    key = "CORALclone_", 
    assay = Seurat::DefaultAssay(seurat_obj)
  )
  
  seurat_obj@misc$isCORAL_embedding <- list(
    P_norm = P_norm,                           
    original_MDS_coords = original_MDS_coords,   
    Y_ref_raw = Y_ref_raw,                     
    valid_cells = valid_cells,                 
    top_genes = top_genes,                     
    source_binsize = binsize_key,              
    lasso_model = cv_fit                       
  )
  
  return(seurat_obj)
}

#' @title Dynamic Manifold Spectral Smoothing
#' @export
SmoothEmbedding <- function(seurat_obj, lambda = 0.5) {
  
  if (is.null(seurat_obj@misc$isCORAL_embedding)) {
    stop("🚨 Error: 'isCORAL_embedding' persistent cache not detected in the object! Please run CORALembedding() first.")
  }
  
  if (lambda < 0 || lambda > 1) {
    stop("🚨 Parameter Out of Bounds: lambda smoothing coefficient must be strictly restricted within the [0, 1] interval!")
  }
  
  cache       <- seurat_obj@misc$isCORAL_embedding
  P_norm      <- cache$P_norm
  Y_ref_raw   <- cache$Y_ref_raw
  valid_cells <- cache$valid_cells
  n_dims      <- ncol(Y_ref_raw)
  
  message(sprintf("⚙️ [isCORAL] Invoking persistent basis, dynamically compiling manifold spectrum for lambda = %s...", lambda))
  
  Y_smooth_max <- as.matrix(P_norm %*% Y_ref_raw)
  Y_final      <- (1 - lambda) * Y_ref_raw + lambda * Y_smooth_max
  
  lambda_tag     <- gsub("\\.", "", as.character(lambda))
  reduction_name <- paste0("CORAL_smooth_lambda_", lambda_tag)
  reduction_key  <- paste0("CORALsmooth", lambda_tag, "_")
  
  rownames(Y_final) <- valid_cells
  colnames(Y_final) <- paste0(reduction_key, 1:n_dims)
  
  seurat_obj[[reduction_name]] <- Seurat::CreateDimReducObject(
    embeddings = Y_final,
    key = reduction_key,
    assay = Seurat::DefaultAssay(seurat_obj)
  )
  return(seurat_obj)
}