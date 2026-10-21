# ==============================================================================
# 🧱 Script A: Single-Cell Lineage Inference - Heavy Baseline Generator Factory
#              (Refactored: Pure Topology Generation, ZERO Evaluation)
# ==============================================================================

# -----------------------------------------------------------------
# 🟢 [0] Unified Global Configuration Console
# -----------------------------------------------------------------
# Set ISCORAL_PROJECT_ROOT and ISCORAL_SAMPLE_ID, or replace the placeholders below.
BASE_DIR      <- normalizePath(Sys.getenv("ISCORAL_PROJECT_ROOT", unset = "."), mustWork = TRUE)
TIMESTAMP     <- format(Sys.time(), "%Y%m%d_%H%M%S")

# =================================================================
# 📸 [T=0 Memory Snapshot Lock]: Automatically fetch path and lock code!
# =================================================================
# Reverted to your original dynamic path fetching function
get_running_script_path <- function() {
  tryCatch({
    # Force load this.path namespace and get path, normalizePath ensures absolute path
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

PACKAGE_DIR <- Sys.getenv(
  "ISCORAL_PACKAGE_DIR",
  unset = file.path(BASE_DIR, "packages", "isCORAL_v6.1.0")
)
devtools::load_all(PACKAGE_DIR)

SAMPLE_ID <- Sys.getenv("ISCORAL_SAMPLE_ID", unset = "<SAMPLE_ID>")
if (identical(SAMPLE_ID, "<SAMPLE_ID>")) {
  stop("Set ISCORAL_SAMPLE_ID or replace <SAMPLE_ID> before running this template.")
}
RAW_DATA_FILE <- Sys.getenv(
  "ISCORAL_RAW_DATA_FILE",
  unset = file.path(BASE_DIR, "1_data", "raw", paste0(SAMPLE_ID, ".RData"))
)

# Core Metadata
GLOBAL_CLONE_COL      <- "barcode_day0" 
GLOBAL_INVALID_CLONES <- c("0","NA", "NaN", "", "None", "Unknown", "unknown", "Unassigned")

# Parameter Console
GLOBAL_HVG         <- 2000         
GLOBAL_PC          <- 30  
GLOBAL_RESOLUTION  <- 1.5    
GEMLI_CUTOFF       <- 50         

# -----------------------------------------------------------------
# --- STEP 0: Loading Libraries & Initializing Data
# -----------------------------------------------------------------
message(sprintf("--- Step 0: Initializing Data for Run [%s] ---", TIMESTAMP))

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(irlba); library(pbapply); library(matrixStats); library(Matrix)
})
pboptions(type = "timer", style = 3, char = "=")

load(RAW_DATA_FILE)

# =================================================================
# 🐛 [FIX]: Convert the data to a sparse matrix before calculate_heritable_genes accesses the underlying @x slot
# =================================================================
message("   -> [Fix] Ensuring data layer is a sparse matrix before core calculation...")
assay_name <- DefaultAssay(sample)
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # Seurat v5-compatible access
  LayerData(sample, assay = assay_name, layer = "data") <- as(LayerData(sample, assay = assay_name, layer = "data"), "CsparseMatrix")
} else {
  # Seurat v4-compatible access
  sample <- SetAssayData(sample, assay = assay_name, slot = "data", new.data = as(GetAssayData(sample, assay = assay_name, slot = "data"), "CsparseMatrix"))
}
# =================================================================

sample <- calculate_heritable_genes(sample, group_by = GLOBAL_CLONE_COL, invalid_barcodes = GLOBAL_INVALID_CLONES, permutation_repeats = 100, n_cores = 20)
seurat_obj <- sample
assay_to_use <- DefaultAssay(seurat_obj) 

safe_get_matrix <- function(obj, assay_name, layer_name, force_dense = FALSE) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") { mat <- LayerData(obj, assay = assay_name, layer = layer_name) } 
  else { mat <- GetAssayData(obj, assay = assay_name, slot = layer_name) }
  if (force_dense) return(as.matrix(mat)) else return(mat)
}

# -----------------------------------------------------------------
# --- STEP 1: Preprocessing & Running Heavy Algorithms
# -----------------------------------------------------------------
message("--- Step 1: Preprocessing & Executing Baseline Engines ---")
seurat_obj <- FindVariableFeatures(seurat_obj, assay=assay_to_use, nfeatures=GLOBAL_HVG, verbose=FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs=GLOBAL_PC, verbose=FALSE, reduction.name = "bench_pca", reduction.key = "benchPC_")
seurat_obj <- RunUMAP(seurat_obj, dims=1:GLOBAL_PC, reduction = "bench_pca", reduction.name = "bench_umap", reduction.key = "benchUMAP_", verbose=FALSE)

# [A] GEMLI
message("  -> Calculating GEMLI Lineages...")
exp_matrix_sparse <- safe_get_matrix(seurat_obj, assay_to_use, "data")

gemli_items <- tryCatch({ GEMLI::predict_lineages(list(gene_expression = exp_matrix_sparse)) }, 
                        error = function(e) { GEMLI::predict_lineages(list(gene_expression = as.matrix(exp_matrix_sparse))) })
gemli_inferred_labels <- if(!is.null(gemli_items)) GEMLI::prediction_to_lineage_information(gemli_items, cutoff = GEMLI_CUTOFF)$predicted_lineages else rep(NA, ncol(seurat_obj))
gemli_mat_fresh <- if(!is.null(gemli_items)) pmax(gemli_items$prediction, t(gemli_items$prediction), na.rm = TRUE) else NULL
if(!is.null(gemli_mat_fresh)) { diag(gemli_mat_fresh) <- 1 } 

# [B] Standard SNN & KNN
message("  -> Calculating Standard SNN & KNN...")
obj_snn <- FindNeighbors(seurat_obj, reduction = "bench_pca", dims=1:GLOBAL_PC, verbose=FALSE)
snn_labels_fresh <- FindClusters(obj_snn, resolution=GLOBAL_RESOLUTION, verbose=FALSE)$seurat_clusters
snn_mat_fresh <- obj_snn@graphs[[grep("_snn", names(obj_snn@graphs), value=TRUE)[1]]]
knn_mat_fresh <- obj_snn@graphs[[grep("_nn$", names(obj_snn@graphs), value=TRUE)[1]]] 

# [C] CC-Regressed SNN & KNN
message("  -> Calculating CC-Regressed SNN & KNN...")
human_s_genes <- cc.genes.updated.2019$s.genes; human_g2m_genes <- cc.genes.updated.2019$g2m.genes
mouse_s_genes <- paste0(toupper(substr(human_s_genes, 1, 1)), tolower(substr(human_s_genes, 2, nchar(human_s_genes))))
mouse_g2m_genes <- paste0(toupper(substr(human_g2m_genes, 1, 1)), tolower(substr(human_g2m_genes, 2, nchar(human_g2m_genes))))

s_genes_use <- if(sum(human_s_genes %in% rownames(seurat_obj)) > sum(mouse_s_genes %in% rownames(seurat_obj))) human_s_genes else mouse_s_genes
g2m_genes_use <- if(sum(human_g2m_genes %in% rownames(seurat_obj)) > sum(mouse_g2m_genes %in% rownames(seurat_obj))) human_g2m_genes else mouse_g2m_genes 

obj_cc <- CellCycleScoring(seurat_obj, s.features=s_genes_use, g2m.features=g2m_genes_use, set.ident=FALSE)
obj_cc <- RunPCA(ScaleData(obj_cc, vars.to.regress=c("S.Score", "G2M.Score"), verbose=FALSE), npcs=GLOBAL_PC, verbose=FALSE, reduction.name = "bench_pca_cc")
obj_cc <- FindNeighbors(obj_cc, reduction = "bench_pca_cc", dims=1:GLOBAL_PC, verbose=FALSE)
snn_cc_labels_fresh <- FindClusters(obj_cc, resolution=GLOBAL_RESOLUTION, verbose=FALSE)$seurat_clusters
snn_cc_mat_fresh <- obj_cc@graphs[[grep("_snn", names(obj_cc@graphs), value=TRUE)[1]]]
knn_cc_mat_fresh <- obj_cc@graphs[[grep("_nn$", names(obj_cc@graphs), value=TRUE)[1]]]

rm(obj_snn, obj_cc); gc()

# -----------------------------------------------------------------
# --- STEP 2: Matrix Alignment & Packaging
# -----------------------------------------------------------------
message("--- Step 2: Packaging Raw Predictions ---")

align_matrix_names <- function(mat, obj) { 
  if (is.null(mat)) return(NULL); cells <- colnames(obj)
  if (is.null(rownames(mat)) || length(intersect(rownames(mat), cells)) == 0) {
    if (nrow(mat) == length(cells) && ncol(mat) == length(cells)) { rownames(mat)<-cells; colnames(mat)<-cells } else return(NULL)
  }
  if (length(mat) > 0 && max(mat, na.rm=TRUE) > 0) diag(mat) <- max(mat, na.rm=TRUE); return(mat) 
}

baseline_raw_labels <- list(
  "Standard SNN"     = snn_labels_fresh,
  "CC-Regressed SNN" = snn_cc_labels_fresh,
  "GEMLI"            = gemli_inferred_labels
)

baseline_raw_matrices <- list(
  "GEMLI"            = align_matrix_names(gemli_mat_fresh, seurat_obj),
  "Standard SNN"     = align_matrix_names(snn_mat_fresh, seurat_obj),
  "CC-Regressed SNN" = align_matrix_names(snn_cc_mat_fresh, seurat_obj),
  "Standard KNN"     = align_matrix_names(knn_mat_fresh, seurat_obj),     
  "CC-Regressed KNN" = align_matrix_names(knn_cc_mat_fresh, seurat_obj)   
)

# -----------------------------------------------------------------
# --- STEP 3: Save Unique Cache & Archiving (T=0 Memory Snapshot)
# -----------------------------------------------------------------
message("--- Step 3: Saving Raw Baseline Cache ---")
cache_dir <- file.path(BASE_DIR, "4_benchmark", "baseline")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

# 1. Physical Backup: Save the code stored in memory at T=0 as an independent .R archive file
backup_name <- sprintf("Baseline_Script_%s_%s.R", SAMPLE_ID, TIMESTAMP)
backup_path <- file.path(cache_dir, backup_name)
writeLines(SNAPSHOT_CODE, backup_path)
message(sprintf("✅ Script source file (T=0 Snapshot) successfully backed up to: %s", backup_path))

# 2. Record SessionInfo
session_info_str <- paste(capture.output(print(sessionInfo())), collapse = "\n")

# 3. Package RDS Cache
cache_filename <- sprintf("Baseline_Cache_%s_%s.rds", SAMPLE_ID, TIMESTAMP)
cache_file_path <- file.path(cache_dir, cache_filename)

baseline_cache <- list(
  run_parameters        = list(
    SAMPLE_ID    = SAMPLE_ID, 
    TIMESTAMP    = TIMESTAMP, 
    SCRIPT_PATH  = current_script_path, # Records the dynamically fetched path
    SOURCE_CODE  = SNAPSHOT_CODE,       # Inject T=0 memory snapshot
    SESSION_INFO = session_info_str
  ), 
  seurat_obj            = seurat_obj,
  baseline_raw_labels   = baseline_raw_labels,    
  baseline_raw_matrices = baseline_raw_matrices   
)

saveRDS(baseline_cache, file = cache_file_path)
message(sprintf("✅ Success! Raw Baseline Data cached to: %s", cache_file_path))
