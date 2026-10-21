# Minimal human-embryo dynamic-inference reproduction.
# Input must already contain the intended cells, pCORAL embeddings and metaclones.
# Spatial reconstruction and upstream clustering are separate analyses.

BASE_DIR <- normalizePath(Sys.getenv("ISCORAL_PROJECT_ROOT", unset = "."), mustWork = TRUE)
WORKFLOW_DIR <- Sys.getenv("ISCORAL_DYNAMIC_DIR", unset = file.path(BASE_DIR, "6_dynamic_inference"))
SAMPLE_ID <- Sys.getenv("ISCORAL_SAMPLE_ID", unset = "<SAMPLE_ID>")

CONFIG <- list(
  profile = "embryo",
  sample_id = SAMPLE_ID,
  input_rds = Sys.getenv("ISCORAL_DYNAMIC_INPUT", unset = ""),
  package_dir = Sys.getenv("ISCORAL_PACKAGE_DIR", unset = file.path(BASE_DIR, "packages", "isCORAL_v6.1.4")),
  output_root = Sys.getenv("ISCORAL_DYNAMIC_OUTPUT_DIR", unset = file.path(WORKFLOW_DIR, "results")),
  run_id = Sys.getenv("ISCORAL_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S")),
  run_cytotrace2 = TRUE,              # FALSE reuses the existing score column.
  cytotrace2_species = "human",
  cytotrace2_slot = "counts",
  cytotrace2_seed = 42L,
  set_seed_before_cytotrace2 = TRUE,  # Matches the original embryo entry script.
  potential_col = "CytoTRACE2_Score",
  metaclone_col = "metaclone",
  cell_type_col = "celltype",
  epsilon = 0.015,
  edge_quantile = 0.70,
  min_forward_flux = 0,
  tolerance = 1e-12,
  annotations = NULL,
  # Fixed embryo palette from the original dataset configuration.
  metaclone_colors = c(
    "1" = "#D55E00", "2" = "#009E73", "3" = "#6F4E9C", "4" = "#E69F00",
    "5" = "#0072B2", "6" = "#B79F00", "7" = "#CC79A7", "8" = "#4D4D4D",
    "9" = "#56B4E9", "10" = "#8C564B", "11" = "#882255", "12" = "#66A61E",
    "13" = "#003F5C", "14" = "#B2182B", "15" = "#00A6A6", "16" = "#B39DDB",
    "17" = "#F2C14E", "18" = "#1B7837", "19" = "#E7298A", "20" = "#607D8B"
  ),
  # Empty roots reproduce the selection rule in the original 6.2b script:
  # the two transient metaclones with the highest mean CytoTRACE2 potential.
  # Supply explicit IDs here when reproducing a manually selected historical run.
  flow_roots = character(),
  flow_sinks = NULL,
  flow_method = "most_probable",
  auto_root_count = 2L,
  sister_fate_seeds = character(),
  auto_sister_seed_count = 3L,
  timeline_reductions = c("umap", "umap.pcoral.cell"),
  dag_reduction = "umap.pcoral.clone",
  # These projections affect visualization only. Defaults match original 6.0b.
  pcoral_umap = list(
    dims = 1:30, n_neighbors = 30L, min_dist = 0.3,
    metric = "cosine", seed = 42L, overwrite = TRUE,
    cell_source = "CORAL_cell_embed", cell_target = "umap.pcoral.cell",
    clone_source = "CORAL_clone_embed", clone_target = "umap.pcoral.clone"
  ),
  font_family = "Arial"
)

# Run these lines in order. Each plotting block previews and then saves its plot.
source(file.path(WORKFLOW_DIR, "01_inference.R"), encoding = "UTF-8")
source(file.path(WORKFLOW_DIR, "02_plots.R"), encoding = "UTF-8")
