# Shared numerical workflow. First run the configuration block in HSC.R or embryo.R.
# The numerical model and path functions are retained from the original workflow.

if (!exists("CONFIG") || !exists("WORKFLOW_DIR")) {
  stop("Run the configuration block in HSC.R or embryo.R first.", call. = FALSE)
}
if (!nzchar(CONFIG$sample_id) || grepl("[<>/\\\\]", CONFIG$sample_id)) {
  stop("Set ISCORAL_SAMPLE_ID or replace <SAMPLE_ID> with a sample identifier.", call. = FALSE)
}
WORKFLOW_DIR <- normalizePath(WORKFLOW_DIR, winslash = "/", mustWork = TRUE)
source(file.path(WORKFLOW_DIR, "R", "metaclone_dynamics_internal.R"), encoding = "UTF-8")
source(file.path(WORKFLOW_DIR, "R", "metaclone_dynamics_api.R"), encoding = "UTF-8")
.iscoral_require_packages(c("Seurat", "Matrix", "igraph"))

# Match the historical package version; no packages are installed by this workflow.
if (!"isCORAL" %in% loadedNamespaces()) {
  if (dir.exists(CONFIG$package_dir)) {
    .iscoral_require_packages("devtools")
    devtools::load_all(CONFIG$package_dir, quiet = TRUE)
  } else {
    .iscoral_require_packages("isCORAL")
    loadNamespace("isCORAL")
  }
}
if (as.character(utils::packageVersion("isCORAL")) != "6.1.4") {
  stop("This reproduction uses isCORAL 6.1.4. Restart R and load that version.", call. = FALSE)
}

# A supplied RDS takes precedence over an in-session object.
if (nzchar(CONFIG$input_rds)) {
  CONFIG$input_rds <- normalizePath(CONFIG$input_rds, winslash = "/", mustWork = TRUE)
  sample <- readRDS(CONFIG$input_rds)
  input_source <- CONFIG$input_rds
} else {
  if (!exists("sample", inherits = FALSE) || !inherits(sample, "Seurat")) {
    stop("Set ISCORAL_DYNAMIC_INPUT to a processed Seurat RDS, or supply `sample`.", call. = FALSE)
  }
  input_source <- "in-session Seurat object: sample"
}
.iscoral_assert_seurat(sample)
required_metadata <- c(CONFIG$metaclone_col, CONFIG$cell_type_col)
if (!all(required_metadata %in% colnames(sample@meta.data))) {
  stop(sprintf("Missing metadata: %s", paste(setdiff(required_metadata, colnames(sample@meta.data)), collapse = ", ")), call. = FALSE)
}
if (is.null(sample@misc$isCORAL_embedding$P_norm) ||
    is.null(sample@misc$isCORAL_embedding$valid_cells)) {
  stop("Input must contain misc$isCORAL_embedding$P_norm and $valid_cells.", call. = FALSE)
}
if (!is.numeric(CONFIG$auto_root_count) || length(CONFIG$auto_root_count) != 1L ||
    !is.finite(CONFIG$auto_root_count) || CONFIG$auto_root_count < 1 ||
    !is.numeric(CONFIG$auto_sister_seed_count) || length(CONFIG$auto_sister_seed_count) != 1L ||
    !is.finite(CONFIG$auto_sister_seed_count) || CONFIG$auto_sister_seed_count < 1) {
  stop("Automatic root and seed counts must be positive numbers.", call. = FALSE)
}

OUTPUT_DIR <- file.path(CONFIG$output_root, CONFIG$profile, CONFIG$sample_id, CONFIG$run_id)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUTPUT_DIR)) stop("Cannot create OUTPUT_DIR.", call. = FALSE)
message(sprintf("Dataset: %s / %s; input: %s", CONFIG$profile, CONFIG$sample_id, input_source))
message(sprintf("Output: %s", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE)))

# Optional embryo UMAP projections. Preserve the original Seurat UMAP and the
# high-dimensional pCORAL embeddings used to infer the upstream metaclones.
umap_settings <- list()
if (!is.null(CONFIG$pcoral_umap)) {
  u <- CONFIG$pcoral_umap
  projection_sources <- c(cell = u$cell_source, clone = u$clone_source)
  projection_targets <- c(cell = u$cell_target, clone = u$clone_target)
  if (anyDuplicated(projection_targets) ||
      any(projection_targets %in% c(projection_sources, CONFIG$timeline_reductions[[1L]]))) {
    stop("pCORAL UMAP targets must have distinct names, separate from source embeddings and the original UMAP.", call. = FALSE)
  }
  for (layer in names(projection_sources)) {
    source_reduction <- projection_sources[[layer]]
    target_reduction <- projection_targets[[layer]]
    if (!source_reduction %in% names(sample@reductions)) {
      stop(sprintf("Missing source reduction: %s", source_reduction), call. = FALSE)
    }
    coordinates <- Seurat::Embeddings(sample, source_reduction)
    dims <- u$dims[u$dims <= ncol(coordinates)]
    neighbors <- min(u$n_neighbors, nrow(coordinates) - 1L)
    if (length(dims) < 2L || neighbors < 2L) {
      stop("Insufficient embedding dimensions or cells for pCORAL UMAP.", call. = FALSE)
    }
    recompute <- isTRUE(u$overwrite) || !target_reduction %in% names(sample@reductions)
    if (recompute) {
      if (target_reduction %in% names(sample@reductions)) sample[[target_reduction]] <- NULL
      set.seed(u$seed)
      sample <- Seurat::RunUMAP(
        sample, reduction = source_reduction, dims = dims,
        n.neighbors = neighbors, min.dist = u$min_dist, metric = u$metric,
        seed.use = u$seed, reduction.name = target_reduction,
        reduction.key = if (layer == "cell") "PCORALCELLUMAP_" else "PCORALCLONEUMAP_",
        verbose = TRUE
      )
    }
    umap_settings[[layer]] <- list(
      source = source_reduction, target = target_reduction,
      dims = dims, n_neighbors = neighbors, min_dist = u$min_dist,
      metric = u$metric, seed = u$seed, recomputed = recompute
    )
  }
}

# CytoTRACE2 is required only when recomputing developmental potential.
if (isTRUE(CONFIG$run_cytotrace2)) {
  .iscoral_require_packages("CytoTRACE2")
  if (isTRUE(CONFIG$set_seed_before_cytotrace2)) set.seed(CONFIG$cytotrace2_seed)
  sample <- CytoTRACE2::cytotrace2(
    sample, is_seurat = TRUE, slot_type = CONFIG$cytotrace2_slot,
    species = CONFIG$cytotrace2_species, seed = CONFIG$cytotrace2_seed
  )
  cytotrace2_settings <- list(
    species = CONFIG$cytotrace2_species, slot_type = CONFIG$cytotrace2_slot,
    seed = CONFIG$cytotrace2_seed,
    version = as.character(utils::packageVersion("CytoTRACE2")),
    recomputed = TRUE
  )
} else {
  cytotrace2_settings <- list(recomputed = FALSE, score_column = CONFIG$potential_col)
  message(sprintf("Reusing stored scores from %s.", CONFIG$potential_col))
}

# This function also writes cell- and metaclone-level potential/pseudotime to metadata.
sample <- RunMetacloneDynamics(
  sample, cluster_col = CONFIG$metaclone_col, cell_type_col = CONFIG$cell_type_col,
  potential_col = CONFIG$potential_col, epsilon = CONFIG$epsilon,
  edge_quantile = CONFIG$edge_quantile, min_forward_flux = CONFIG$min_forward_flux,
  tolerance = CONFIG$tolerance, annotations = CONFIG$annotations,
  cytotrace2_settings = cytotrace2_settings
)
model <- sample@misc$isCORAL_dynamics$model
terminal <- sample@misc$isCORAL_dynamics$fate$terminal
print(model$state_summary)

# Resolve embryo defaults without editing the configuration source file.
transient_order <- terminal$transient_nodes[order(
  sample@misc$isCORAL_dynamics$potential$metaclone_potential[terminal$transient_nodes],
  decreasing = TRUE
)]
flow_roots <- as.character(CONFIG$flow_roots)
if (length(flow_roots) == 0L) {
  flow_roots <- head(transient_order, CONFIG$auto_root_count)
}
sister_seeds <- as.character(CONFIG$sister_fate_seeds)
if (length(sister_seeds) == 0L) {
  sister_seeds <- head(transient_order, CONFIG$auto_sister_seed_count)
}
if (!all(sister_seeds %in% rownames(model$P_dynamic))) {
  stop("One or more sister_fate_seeds are absent from this model. Update CONFIG.", call. = FALSE)
}

# Paths use the jump chain, whereas cumulative edge flow retains residency in P_dynamic.
if (length(flow_roots) > 0L) {
  sample <- RunMetacloneFlow(
    sample, roots = flow_roots, sinks = CONFIG$flow_sinks,
    method = CONFIG$flow_method, tolerance = CONFIG$tolerance
  )
} else {
  message("All metaclones are terminal; no transient root is available for path inference.")
}

# Expected edge traversals from a unit initial mass at each root (original step 9.2).
expected_flow <- list()
for (root in flow_roots) {
  flux <- matrix(0, nrow(model$P_dynamic), ncol(model$P_dynamic), dimnames = dimnames(model$P_dynamic))
  if (root %in% terminal$transient_nodes) {
    visits <- terminal$fundamental_matrix[match(root, terminal$transient_nodes), ]
    flux[terminal$transient_nodes, ] <- sweep(
      model$P_dynamic[terminal$transient_nodes, , drop = FALSE], 1, visits, "*"
    )
  }
  diag(flux) <- 0
  expected_flow[[root]] <- flux
}
sample@misc$isCORAL_dynamics$expected_flow <- list(
  source_model_id = model$model_id,
  definition = "off-diagonal expected edge traversals: N[root,i] * P_dynamic[i,j]",
  roots = flow_roots, matrices = expected_flow
)

# Record the requested configuration and the values resolved for this input.
source_files <- c(
  file.path(WORKFLOW_DIR, paste0(CONFIG$profile, ".R")),
  file.path(WORKFLOW_DIR, c("01_inference.R", "02_plots.R")),
  list.files(file.path(WORKFLOW_DIR, "R"), pattern = "[.]R$", full.names = TRUE)
)
source_files <- source_files[file.exists(source_files)]
sample@misc$isCORAL_dynamics$reproduction <- list(
  config = CONFIG, input_source = input_source,
  n_cells = ncol(sample), n_features = nrow(sample),
  default_assay = Seurat::DefaultAssay(sample),
  resolved_flow_roots = flow_roots, resolved_sister_seeds = sister_seeds,
  pcoral_umap = umap_settings,
  source_md5 = tools::md5sum(source_files),
  session_info = utils::sessionInfo()
)

# Export one full object plus a compact cache and model tables.
cache <- sample@misc$isCORAL_dynamics
saveRDS(sample, file.path(OUTPUT_DIR, "sample_with_dynamics.rds"))
saveRDS(cache, file.path(OUTPUT_DIR, "dynamic_cache.rds"))
dput(CONFIG, file = file.path(OUTPUT_DIR, "config_used.R"))
writeLines(capture.output(utils::sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
utils::write.csv(cache$model$state_summary, file.path(OUTPUT_DIR, "metaclone_state_summary.csv"), row.names = FALSE)
utils::write.csv(cache$potential$seed_core_identity, file.path(OUTPUT_DIR, "seed_core_identity.csv"), row.names = FALSE)
utils::write.csv(cache$potential$metaclone_potential_statistics, file.path(OUTPUT_DIR, "seed_core_potential_summary.csv"), row.names = FALSE)
matrix_outputs <- list(
  block_flux = model$block_flux,
  forward_flux_before_sparsification = model$forward_flux_formula,
  forward_flux_after_sparsification = model$forward_flux_model,
  P_dynamic_formula = model$P_dynamic_formula,
  P_dynamic = model$P_dynamic,
  sister_fate_forward_flux = cache$fate$sister$forward_mass,
  sister_fate_conditional = cache$fate$sister$conditional_probability,
  terminal_fate = terminal$fate_matrix
)
for (name in names(matrix_outputs)) {
  utils::write.csv(matrix_outputs[[name]], file.path(OUTPUT_DIR, paste0(name, ".csv")))
}
if (!is.null(cache$flow)) {
  flow_table <- do.call(rbind, lapply(names(cache$flow$results), function(name) {
    result <- cache$flow$results[[name]]
    data.frame(
      Comparison = name, Method = result$method, Root = result$root, Sink = result$sink,
      Reachable = isTRUE(result$reachable), Nodes = paste(result$nodes %||% character(), collapse = " -> "),
      PathProbability = result$path_probability %||% NA_real_,
      BottleneckProbability = result$bottleneck_probability %||% NA_real_,
      MaxFlowValue = result$max_flow_value %||% NA_real_, stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(flow_table, file.path(OUTPUT_DIR, "path_summary.csv"), row.names = FALSE)
  utils::write.csv(cache$flow$jump_chain, file.path(OUTPUT_DIR, "jump_chain.csv"))
}
if (length(expected_flow) > 0L) {
  expected_edges <- do.call(rbind, lapply(names(expected_flow), function(root) {
    edges <- .iscoral_matrix_edges(expected_flow[[root]], "ExpectedTraversals", CONFIG$tolerance)
    data.frame(Root = rep(root, nrow(edges)), edges, stringsAsFactors = FALSE)
  }))
  utils::write.csv(expected_edges, file.path(OUTPUT_DIR, "expected_edge_traversals.csv"), row.names = FALSE)
}
message(sprintf("Numerical results saved to %s", OUTPUT_DIR))
