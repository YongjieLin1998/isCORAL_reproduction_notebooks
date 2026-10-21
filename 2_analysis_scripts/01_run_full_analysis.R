################################################################################
#
#       isCORAL v6.1.0 Full Dataset Analysis Script (Production)
#       [Updated: Storage-Optimized, RAM Anchoring & Real-Time Logger]
#
################################################################################



# --- Step 1: Environment Preparation ---
message("--- Loading all required packages... ---")
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

library(data.table)
library(Seurat)
library(dplyr)
library(stringr)
library(purrr)
library(matrixStats)
library(kernlab)
library(Matrix)
library(RcppML)
message("--- All packages loaded. ---")

# --- Step 2: Load Config File and Setup Logger Bridge ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) { 
  stop("CRITICAL ERROR: Please provide a config file path!") 
}
config_path <- args[1]

# [Logger Bridge]: Receive Archive directory from Bash (second argument)
archive_dir <- NULL
if (length(args) >= 2) {
  archive_dir <- args[2]
  message(">>> Connected to Bash Archive Directory: ", archive_dir)
}

if (!file.exists(config_path) && file.exists(here::here(config_path))) {
   config_path <- here::here(config_path)
}

config <- read_yaml(config_path)
message("Using config file: ", config_path)

package_path <- here::here("packages", config$run_settings$package_version)
message("Loading algorithm from: ", package_path)
load_all(package_path)

# --- Step 3: Script Startup and Timing ---
total_run_start_time <- Sys.time()
message("========================================================")
message("Starting FULL DATASET ANALYSIS with isCORAL v6.1.0...")

# --- Step 4: Load Data (Robust Dynamic Loading) ---
message("\n--- Step 4: Loading Seurat Input Data ---")
input_data_path <- config$paths$input_seurat_object
load_env <- new.env()
load(input_data_path, envir = load_env)
loaded_objects <- ls(envir = load_env)

if (length(loaded_objects) == 1) {
    seurat_obj_name <- loaded_objects[1]
} else {
    if (!is.null(config$run_settings$input_obj_name) && config$run_settings$input_obj_name %in% loaded_objects) {
        seurat_obj_name <- config$run_settings$input_obj_name
    } else {
        stop("Multiple objects found in input file and config is missing.")
    }
}
seurat_obj <- get(seurat_obj_name, envir = load_env)
rm(load_env); gc()
message("Data object loaded successfully.")

# --- Step 5: Execute Core Analysis ---
message("\n########################################################")
message("Step 5: Executing Core isCORAL v6.1.0 Analysis...")

# -------------------------------------------------------------------
# [Physical RAM Anchoring Patch]
# Prevents the Linux kernel from dropping Page Cache under heavy memory 
# pressure (e.g., 50GB), which causes lazy-load corruption errors.
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

temp_dir_ssd <- file.path(config$paths$ssd_temp_root, "full_data_run")

params_from_config <- config$iscoral_params
params_from_config$seurat_obj <- NULL
params_from_config$results_dir <- NULL
params_from_config$random_seed <- NULL

iscoral_args <- c(
  list(
    seurat_obj = seurat_obj,
    results_dir = temp_dir_ssd,
    random_seed = config$run_settings$master_seed
  ),
  params_from_config
)

# [Logger Bridge]: Route real-time progress to Archive directory if it exists
if (!is.null(archive_dir)) {
    iscoral_args$custom_log_file <- file.path(archive_dir, "realtime_progress.log")
    message(">>> Real-time progress will be routed to: ", iscoral_args$custom_log_file)
}

# Execute the analysis
seurat_obj_processed <- do.call(isCORAL_v6_1_0, iscoral_args)

message("Core isCORAL v6.1.0 analysis complete!")
message("########################################################\n")

# --- Step 6: Save Final Results (STORAGE OPTIMIZED & SMART NAMING) ---
message("\n========================================================")
message("Step 6: Analysis finished. Extracting and saving isCORAL results ONLY.")

if ("CORAL_scan_results" %in% names(seurat_obj_processed@misc)) {
  coral_results_only <- seurat_obj_processed@misc$CORAL_scan_results
} else {
  stop("CRITICAL ERROR: isCORAL results not found in the processed object.")
}

# [Smart Naming Logic]
base_filename <- config$paths$final_rdata_filename

if (is.null(base_filename) || tolower(base_filename) == "auto") {
  input_basename <- tools::file_path_sans_ext(basename(config$paths$input_seurat_object))
  current_bin <- config$iscoral_params$bin_size
  base_filename <- paste0("coral_res_", input_basename, "_bin", current_bin, ".rds")
} else {
  base_filename <- sub("\\.RData$|\\.rda$", ".rds", base_filename, ignore.case = TRUE)
  if (!grepl("\\.rds$", base_filename, ignore.case = TRUE)) base_filename <- paste0(base_filename, ".rds")
}

# [Append Timestamp]
if (!is.null(config$run_settings$add_timestamp_to_output) && config$run_settings$add_timestamp_to_output) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  final_filename <- sub("(\\.[^.]+)$", paste0("_", timestamp, "\\1"), base_filename)
} else {
  final_filename <- base_filename
}

final_save_dir <- config$paths$hdd_final_root
if (!dir.exists(final_save_dir)) dir.create(final_save_dir, recursive = TRUE)
final_save_path <- file.path(final_save_dir, final_filename)

message(paste("Saving standalone isCORAL results to:", final_filename))
saveRDS(coral_results_only, file = final_save_path)

total_run_end_time <- Sys.time()
total_duration <- difftime(total_run_end_time, total_run_start_time, units = "mins")
message(paste("Successfully saved light-weight results to:", final_save_path))
message(paste("Total running time:", round(as.numeric(total_duration), 2), "minutes"))
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
  
  message("Environment snapshot successfully anchored to archive.")
} else {
  warning("Archive directory not found or not provided. Skipping environment snapshot.")
}
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