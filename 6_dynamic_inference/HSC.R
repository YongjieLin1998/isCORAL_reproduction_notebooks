# Minimal HSC dynamic-inference reproduction.
# Run from the reproduction project root, or set ISCORAL_PROJECT_ROOT.
# Supply a processed Seurat RDS, or reuse `sample` from the embedding workflow.

BASE_DIR <- normalizePath(Sys.getenv("ISCORAL_PROJECT_ROOT", unset = "."), mustWork = TRUE)
WORKFLOW_DIR <- Sys.getenv("ISCORAL_DYNAMIC_DIR", unset = file.path(BASE_DIR, "6_dynamic_inference"))
SAMPLE_ID <- Sys.getenv("ISCORAL_SAMPLE_ID", unset = "<SAMPLE_ID>")

CONFIG <- list(
  profile = "HSC",
  sample_id = SAMPLE_ID,
  input_rds = Sys.getenv("ISCORAL_DYNAMIC_INPUT", unset = ""),
  package_dir = Sys.getenv("ISCORAL_PACKAGE_DIR", unset = file.path(BASE_DIR, "packages", "isCORAL_v6.1.4")),
  output_root = Sys.getenv("ISCORAL_DYNAMIC_OUTPUT_DIR", unset = file.path(WORKFLOW_DIR, "results")),
  run_id = Sys.getenv("ISCORAL_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S")),
  run_cytotrace2 = TRUE,              # FALSE reuses the existing score column.
  cytotrace2_species = "mouse",
  cytotrace2_slot = "counts",
  cytotrace2_seed = 42L,
  set_seed_before_cytotrace2 = FALSE, # Matches the original HSC entry script.
  potential_col = "CytoTRACE2_Score",
  metaclone_col = "metaclone",
  cell_type_col = "cluster_integrated",
  epsilon = 0.015,
  edge_quantile = 0.70,
  min_forward_flux = 0,
  tolerance = 1e-12,
  # These biological labels apply only to the original HSC metaclone numbering.
  # Use NULL or replace this map when analyzing a newly clustered object.
  annotations = c(
    "1" = "MK-biased HSC/MPP", "2" = "GMP-primed late MPP",
    "3" = "Early GMP", "4" = "MK-primed late MPP2",
    "5" = "MK-biased HSC", "6" = "Balanced early MPP",
    "7" = "Balanced MPP", "8" = "MK Progenitor",
    "9" = "GMP-biased HSC", "10" = "MK-biased MPP",
    "11" = "MK-biased MPP3", "12" = "MK-primed late MPP1",
    "13" = "GMP-primed early MPP", "14" = "MK-biased MPP2",
    "15" = "Balanced late MPP", "16" = "Late Mono-biased GMP"
  ),
  metaclone_colors = NULL,            # Use isCORAL::GetMetacloneColors().
  flow_roots = c("5", "9"),          # Historical HSC roots; edit for new labels.
  flow_sinks = NULL,                 # All terminal states of the denoised model.
  flow_method = "most_probable",
  auto_root_count = 2L,              # Used only when flow_roots is empty.
  sister_fate_seeds = c("5", "9", "15"),
  auto_sister_seed_count = 3L,
  timeline_reductions = c("umap.raw", "CORAL_cell_embed"),
  dag_reduction = "CORAL_clone_embed",
  pcoral_umap = NULL,                # HSC uses the stored embedding axes.
  font_family = "Arial"
)

# Run these lines in order. Each plotting block previews and then saves its plot.
source(file.path(WORKFLOW_DIR, "01_inference.R"), encoding = "UTF-8")
source(file.path(WORKFLOW_DIR, "02_plots.R"), encoding = "UTF-8")
