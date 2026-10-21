################################################################################
#
#       isCORAL v6.1.0 Analysis Script 03: Downsampling Analysis
#       [Updated: RAM Anchoring, OS-Aware Parallelism & Grid Logger Bridge]
#       [Storage-Optimized: Consolidates 18 runs into a SINGLE .rds file]
#
################################################################################


# ===================================================================

# --- Step 1: Environment Preparation ---
message("--- Loading required packages... ---")
library(devtools)
library(yaml)
library(here)
library(foreach)

# [CRITICAL FIX]: Force doMC for Copy-on-Write memory sharing on Linux
if (.Platform$OS.type != "windows") {
  if (!requireNamespace("doMC", quietly = TRUE)) {
    stop("CRITICAL ERROR: 'doMC' package is missing! Running without doMC on Linux will cause massive memory duplication (OOM). Please install doMC.")
  }
  library(doMC)
  message("Memory Sharing Backend (doMC) successfully loaded.")
} else {
  library(doSNOW)
  warning("Running on Windows. Memory sharing via Copy-on-Write is disabled.")
}

library(Seurat)
library(Matrix)
library(stringr)
message("--- All packages loaded. ---")

# --- Step 2: Load Config & Setup Logger Bridge ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("CRITICAL ERROR: Please provide a config file path!")
config_path <- args[1]

# [Logger Bridge]: Receive Archive directory passed from Bash
archive_dir <- NULL
if (length(args) >= 2) {
  archive_dir <- args[2]
  message(">>> Connected to Bash Archive Directory: ", archive_dir)
}

if (!file.exists(config_path) && file.exists(here::here(config_path))) config_path <- here::here(config_path)

config <- read_yaml(config_path)
package_path <- here::here("packages", config$run_settings$package_version)
load_all(package_path)

total_run_start_time <- Sys.time()
message("========================================================")
message("Starting DOWNSAMPLING ANALYSIS with isCORAL v6.1.0...")

# --- Step 3: Load Original Data ---
input_data_path <- config$paths$input_seurat_object
load_env <- new.env()
load(input_data_path, envir = load_env)
loaded_objects <- ls(envir = load_env)

if (length(loaded_objects) == 1) {
    seurat_obj_name <- loaded_objects[1]
} else {
    seurat_obj_name <- config$run_settings$input_obj_name
}
seurat_obj_original <- get(seurat_obj_name, envir = load_env)
rm(load_env); gc()
message(paste("Original data loaded. Total cells:", ncol(seurat_obj_original)))

# -------------------------------------------------------------------
# [Physical RAM Anchoring Patch]
# Prevents the kernel from reclaiming Page Cache during the heavy 
# 18-iteration loop, which causes lazy-load corruption errors.
# -------------------------------------------------------------------
if (requireNamespace("doRNG", quietly = TRUE)) {
  message("--- [PATCH] Anchoring parallel libraries into physical RAM ---")
  library(doRNG)
  library(rngtools)
  invisible(as.list(getNamespace("doRNG")))
  invisible(as.list(getNamespace("rngtools")))
  message("--- [PATCH] Libraries anchored. IO-related corruptions blocked. ---")
}
# -------------------------------------------------------------------

# --- Step 4: Define Downsampling Grid & Initialize Master List ---
test_modes <- c("Cell", "UMI")
fractions <- c(0.2, 0.5, 1)
repeats <- 1:3

temp_dir_ssd <- file.path(config$paths$ssd_temp_root, "downsampling_analysis_run")
if(!dir.exists(temp_dir_ssd)) dir.create(temp_dir_ssd, recursive = TRUE)

params_from_config <- config$iscoral_params
params_from_config$seurat_obj <- NULL
params_from_config$results_dir <- NULL
params_from_config$random_seed <- NULL
assay_to_use <- if(!is.null(params_from_config$assay)) params_from_config$assay else "RNA"

benchmark_master_list <- list()
benchmark_master_list$metadata <- list(
  source_dataset = basename(input_data_path),
  run_time = Sys.time(),
  base_parameters = params_from_config,
  grid_setup = list(modes = test_modes, fractions = fractions, repeats = repeats)
)
benchmark_master_list$results <- list()

# --- Step 5: Execute Downsampling Loop ---
message("\n########################################################")
message("Starting Grid Iteration: 18 Independent Runs")

for (mode in test_modes) {
  for (frac in fractions) {
    for (rep in repeats) {
      
      run_id <- paste0(mode, "_", frac, "_rep", rep)
      message(paste("\n>>> Now Running:", run_id))
      
      current_seed <- config$run_settings$master_seed + rep * 100 
      set.seed(current_seed)
      
      test_obj <- seurat_obj_original
      
      if (mode == "Cell") {
        target_cells <- max(10, round(ncol(test_obj) * frac))
        sampled_cells <- sample(colnames(test_obj), size = target_cells, replace = FALSE)
        test_obj <- subset(test_obj, cells = sampled_cells)
        message(paste("    - Subsampled to", target_cells, "cells."))
        
      } else if (mode == "UMI") {
        counts_mat <- Seurat::GetAssayData(test_obj, assay = assay_to_use, layer = "counts")
        new_counts <- counts_mat
        new_counts@x <- as.numeric(rbinom(n = length(counts_mat@x), size = counts_mat@x, prob = frac))
        new_counts <- Matrix::drop0(new_counts)
        test_obj[[assay_to_use]]$counts <- new_counts
        message(paste("    - UMI counts downsampled to ~", frac*100, "%."))
      }
      
      iscoral_args <- c(
        list(seurat_obj = test_obj, results_dir = temp_dir_ssd, random_seed = current_seed),
        params_from_config
      )
      
      # [Grid Logger Logic]: Generate independent log files for each iteration to prevent overwrites
      if (!is.null(archive_dir)) {
          grid_log_file <- file.path(archive_dir, paste0("progress_", run_id, ".log"))
          iscoral_args$custom_log_file <- grid_log_file
          message("    - Real-time progress routed to: ", grid_log_file)
      }

      # Execute isCORAL
      seurat_obj_processed <- do.call(isCORAL_v6_1_0, iscoral_args)
      
      if (!"CORAL_scan_results" %in% names(seurat_obj_processed@misc)) {
        stop("    - [ERROR] isCORAL results not generated for this run.")
      }
      
      benchmark_master_list$results[[run_id]] <- seurat_obj_processed@misc$CORAL_scan_results
      message(paste("    - Run", run_id, "completed and saved to memory list."))
      
      # Aggressively clear memory for the current iteration
      rm(test_obj, seurat_obj_processed, iscoral_args)
      gc(verbose = FALSE, reset = TRUE)
    }
  }
}

# --- Step 6: Save the Consolidated Master List ---
message("\n########################################################")
message("Grid iteration finished. Consolidating and saving final object.")

base_filename <- config$paths$final_rdata_filename
if (is.null(base_filename) || tolower(base_filename) == "auto") {
  input_basename <- tools::file_path_sans_ext(basename(config$paths$input_seurat_object))
  current_bin <- config$iscoral_params$bin_size
  base_filename_core <- paste0("coral_res_", input_basename, "_bin", current_bin)
} else {
  base_filename_core <- sub("\\.RData$|\\.rda$|\\.rds$", "", base_filename, ignore.case = TRUE)
}

final_save_dir <- file.path(config$paths$hdd_final_root, "downsampling_analysis")
if (!dir.exists(final_save_dir)) dir.create(final_save_dir, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
final_filename <- paste0(base_filename_core, "_DownsamplingGrid_", timestamp, ".rds")
final_save_path <- file.path(final_save_dir, final_filename)

saveRDS(benchmark_master_list, file = final_save_path)

total_run_end_time <- Sys.time()
total_duration <- difftime(total_run_end_time, total_run_start_time, units = "hours")

message("\n========================================================")
message("🎉 DOWNSAMPLING ANALYSIS COMPLETELY FINISHED!")
message(paste(">> Total Analysis Time:", round(as.numeric(total_duration), 2), "hours"))
message(paste(">> Single consolidated .rds file saved at:\n  ", final_save_path))
message("========================================================")


# ===================================================================
# --- [NEW] Step 7: Environment Snapshot Archiving                ---
# ===================================================================
if (!is.null(archive_dir) && dir.exists(archive_dir)) {
  message("--- Saving strict environment dependencies snapshot... ---")
  
  # Save as a machine-readable .rds object (core for future automated reproduction)
  saveRDS(sessionInfo(), file.path(archive_dir, "session_info_snapshot.rds"))
  
  # Save as a human-readable .txt file
  writeLines(capture.output(print(sessionInfo())), 
             file.path(archive_dir, "session_info_snapshot.txt"))
  
  # -----------------------------------------------------------------
  # [NEW]: Snapshot the actual isCORAL package source code
  # -----------------------------------------------------------------
  if (exists("package_path") && dir.exists(package_path)) {
    message("--- Archiving isCORAL package source code... ---")
    
    # Define the archive filename
    tar_filename <- file.path(archive_dir, "isCORAL_source_snapshot.tar.gz")
    
    # Compress and archive the package_path directory
    tar(tarfile = tar_filename, files = package_path, compression = "gzip")
    
    message("Package source successfully archived to: ", tar_filename)
  } else {
    warning("Package path not found in memory. Cannot archive source code.")
  }
  # -----------------------------------------------------------------
  
  message("Environment snapshot successfully anchored to archive.")
} else {
  warning("Archive directory not found or not provided. Skipping environment snapshot.")
}
message("========================================================")