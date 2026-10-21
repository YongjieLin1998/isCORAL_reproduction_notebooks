# Workflow-local numerical functions for isCORAL dynamic inference.
# Source metaclone_dynamics_internal.R before this file.

#' Infer a denoised metaclone dynamic model
#'
#' Builds seed-core potential summaries, directionally constrained and
#' sparsified metaclone transition matrices, sister-fate flux, and terminal
#' absorption probabilities. All terminal fates are computed from the denoised
#' P_dynamic matrix.
#'
#' @return The modified Seurat object.
RunMetacloneDynamics <- function(
    seurat_obj,
    cluster_col = "metaclone",
    cell_type_col = "cluster_integrated",
    potential_col = "CytoTRACE2_Score",
    epsilon = 0.015,
    edge_quantile = 0.70,
    min_forward_flux = 0,
    tolerance = 1e-12,
    annotations = NULL,
    cytotrace2_settings = NULL
) {
  .iscoral_assert_seurat(seurat_obj)
  .iscoral_require_packages(c("Matrix", "Seurat", "isCORAL"))
  loaded_version <- as.character(utils::packageVersion("isCORAL"))
  if (!identical(loaded_version, .iscoral_required_version)) {
    stop(
      sprintf("This workflow requires isCORAL %s; loaded version is %s.",
              .iscoral_required_version, loaded_version),
      call. = FALSE
    )
  }
  required_metadata <- c(cluster_col, potential_col)
  missing_metadata <- setdiff(required_metadata, colnames(seurat_obj@meta.data))
  if (length(missing_metadata) > 0L) {
    stop(sprintf("Missing metadata: %s", paste(missing_metadata, collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0) {
    stop("epsilon must be one finite non-negative number.", call. = FALSE)
  }
  if (is.null(edge_quantile) || !is.numeric(edge_quantile) ||
      length(edge_quantile) != 1L || !is.finite(edge_quantile) ||
      edge_quantile < 0 || edge_quantile > 1) {
    stop("edge_quantile must be one number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(min_forward_flux) || length(min_forward_flux) != 1L ||
      !is.finite(min_forward_flux) || min_forward_flux < 0) {
    stop("min_forward_flux must be one finite non-negative number.", call. = FALSE)
  }

  aggregate <- .iscoral_aggregate_metaclone_flux(seurat_obj, cluster_col)
  metaclone_ids <- rownames(aggregate$block_flux)
  metaclone_annotations <- .iscoral_resolve_annotations(metaclone_ids, annotations)

  cell_potential <- as.numeric(seurat_obj@meta.data[[potential_col]])
  names(cell_potential) <- rownames(seurat_obj@meta.data)
  if (anyNA(cell_potential) || any(!is.finite(cell_potential)) ||
      any(cell_potential < 0 | cell_potential > 1)) {
    stop(sprintf("Metadata '%s' must contain finite values in [0, 1].", potential_col), call. = FALSE)
  }

  core_identity <- do.call(rbind, lapply(metaclone_ids, function(metaclone_id) {
    cells <- aggregate$core_cells[[metaclone_id]]
    data.frame(
      Cell = cells,
      Metaclone = metaclone_id,
      CellType = if (cell_type_col %in% colnames(seurat_obj@meta.data)) {
        as.character(seurat_obj@meta.data[cells, cell_type_col])
      } else {
        NA_character_
      },
      DevelopmentalPotential = unname(cell_potential[cells]),
      Pseudotime = 1 - unname(cell_potential[cells]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(core_identity) <- core_identity$Cell

  metaclone_potential <- stats::setNames(
    vapply(aggregate$core_cells, function(cells) mean(cell_potential[cells]), numeric(1)),
    metaclone_ids
  )
  metaclone_potential_statistics <- do.call(rbind, lapply(
    metaclone_ids,
    function(metaclone_id) {
      values <- cell_potential[aggregate$core_cells[[metaclone_id]]]
      data.frame(
        Metaclone = metaclone_id,
        Mean = mean(values),
        Median = stats::median(values),
        SD = if (length(values) > 1L) stats::sd(values) else NA_real_,
        Q25 = as.numeric(stats::quantile(values, 0.25, names = FALSE, type = 8)),
        Q75 = as.numeric(stats::quantile(values, 0.75, names = FALSE, type = 8)),
        NSeedCells = length(values),
        stringsAsFactors = FALSE,
        row.names = metaclone_id
      )
    }
  ))
  metaclone_pseudotime <- 1 - metaclone_potential
  cell_type_counts <- NULL
  if (cell_type_col %in% colnames(seurat_obj@meta.data)) {
    cell_type_counts <- as.matrix(table(
      Metaclone = factor(core_identity$Metaclone, levels = metaclone_ids),
      CellType = core_identity$CellType,
      useNA = "ifany"
    ))
  }

  block_flux <- aggregate$block_flux
  residency_flux <- diag(block_flux)
  names(residency_flux) <- metaclone_ids
  potential_drop <- outer(metaclone_potential, metaclone_potential, "-")
  dimnames(potential_drop) <- list(metaclone_ids, metaclone_ids)
  forward_mask <- potential_drop > epsilon
  diag(forward_mask) <- FALSE

  forward_flux_formula <- block_flux
  diag(forward_flux_formula) <- 0
  forward_flux_formula[!forward_mask] <- 0
  positive_flux <- forward_flux_formula[forward_flux_formula > 0]
  quantile_threshold <- if (length(positive_flux) > 0L) {
    as.numeric(stats::quantile(
      positive_flux,
      probs = edge_quantile,
      names = FALSE,
      type = 8
    ))
  } else {
    0
  }
  effective_threshold <- max(quantile_threshold, min_forward_flux)
  forward_flux_model <- forward_flux_formula
  forward_flux_model[forward_flux_model < effective_threshold] <- 0

  formula_dynamic <- .iscoral_make_dynamic_matrix(
    residency_flux,
    forward_flux_formula,
    tolerance
  )
  model_dynamic <- .iscoral_make_dynamic_matrix(
    residency_flux,
    forward_flux_model,
    tolerance
  )
  P_dynamic_formula <- formula_dynamic$matrix
  P_dynamic <- model_dynamic$matrix
  absorption <- .iscoral_compute_absorption(P_dynamic)
  sister_fate <- .iscoral_make_sister_fate(forward_flux_model, residency_flux)

  off_diagonal_formula <- P_dynamic_formula
  off_diagonal_model <- P_dynamic
  diag(off_diagonal_formula) <- 0
  diag(off_diagonal_model) <- 0
  sinks_formula <- metaclone_ids[rowSums(off_diagonal_formula) <= tolerance]
  sinks_model <- metaclone_ids[rowSums(off_diagonal_model) <= tolerance]
  discarded_by_direction <- rowSums(block_flux) - residency_flux -
    rowSums(forward_flux_formula)
  removed_by_sparsification <- rowSums(forward_flux_formula - forward_flux_model)

  state_summary <- data.frame(
    Metaclone = metaclone_ids,
    Annotation = unname(metaclone_annotations[metaclone_ids]),
    SeedCoreSize = unname(aggregate$core_size[metaclone_ids]),
    DevelopmentalPotential = unname(metaclone_potential[metaclone_ids]),
    MedianDevelopmentalPotential = metaclone_potential_statistics[metaclone_ids, "Median"],
    SDDevelopmentalPotential = metaclone_potential_statistics[metaclone_ids, "SD"],
    Q25DevelopmentalPotential = metaclone_potential_statistics[metaclone_ids, "Q25"],
    Q75DevelopmentalPotential = metaclone_potential_statistics[metaclone_ids, "Q75"],
    Pseudotime = unname(metaclone_pseudotime[metaclone_ids]),
    ResidencyFlux = unname(residency_flux[metaclone_ids]),
    ForwardFluxBeforeSparsification = rowSums(forward_flux_formula),
    ForwardFluxAfterSparsification = rowSums(forward_flux_model),
    DirectionallyDiscardedFlux = unname(discarded_by_direction),
    SparsificationRemovedFlux = unname(removed_by_sparsification),
    IsFormulaSink = metaclone_ids %in% sinks_formula,
    IsModelSink = metaclone_ids %in% sinks_model,
    stringsAsFactors = FALSE,
    row.names = metaclone_ids
  )

  settings <- list(
    cluster_col = cluster_col,
    cell_type_col = cell_type_col,
    potential_col = potential_col,
    epsilon = epsilon,
    edge_quantile = edge_quantile,
    min_forward_flux = min_forward_flux,
    quantile_threshold = quantile_threshold,
    effective_forward_flux_threshold = effective_threshold,
    tolerance = tolerance
  )
  model_id <- .iscoral_model_id(list(
    P_dynamic = P_dynamic,
    settings = settings,
    core_cells = aggregate$core_cells
  ))
  created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

  sister_fate$source_model_id <- model_id
  terminal_fate <- c(
    list(
      definition = "sparsified_model_only",
      source_matrix = "model$P_dynamic",
      source_model_id = model_id
    ),
    absorption
  )

  cell_ids <- colnames(seurat_obj)
  if (!all(cell_ids %in% rownames(seurat_obj@meta.data))) {
    stop("Seurat metadata cannot be aligned to all cells by barcode.", call. = FALSE)
  }
  metaclone_by_cell <- as.character(
    seurat_obj@meta.data[cell_ids, cluster_col]
  )
  metadata_update <- data.frame(
    cell_developmental_potential = unname(cell_potential[cell_ids]),
    cell_pseudotime = 1 - unname(cell_potential[cell_ids]),
    metaclone_developmental_potential = unname(
      metaclone_potential[metaclone_by_cell]
    ),
    metaclone_pseudotime = unname(
      metaclone_pseudotime[metaclone_by_cell]
    ),
    stringsAsFactors = FALSE,
    row.names = cell_ids
  )
  if (anyNA(metadata_update)) {
    stop(
      "Developmental metadata contains NA after cell/metaclone alignment.",
      call. = FALSE
    )
  }
  seurat_obj <- Seurat::AddMetaData(
    object = seurat_obj,
    metadata = metadata_update
  )

  seurat_obj@misc$isCORAL_dynamics <- list(
    schema_version = "2.0",
    iscoral_version = loaded_version,
    created_at = created_at,
    potential = list(
      method = "CytoTRACE2 mean over P_norm-aligned seed cores",
      cytotrace2_settings = cytotrace2_settings,
      score_column = potential_col,
      seed_core_cells = aggregate$core_cells,
      seed_core_identity = core_identity,
      seed_core_cell_type_counts = cell_type_counts,
      metaclone_potential = metaclone_potential,
      metaclone_potential_statistics = metaclone_potential_statistics,
      metaclone_pseudotime = metaclone_pseudotime
    ),
    model = list(
      model_id = model_id,
      created_at = created_at,
      formula_orientation = list(
        rows = "seed-core metaclone",
        columns = "relevant sister-cell metaclone",
        potential = "larger CytoTRACE2 score is more upstream"
      ),
      settings = settings,
      annotations = metaclone_annotations,
      core_cells = aggregate$core_cells,
      core_size = aggregate$core_size,
      cell_to_metaclone = aggregate$cell_to_metaclone,
      block_sum = aggregate$block_sum,
      block_flux = block_flux,
      residency_flux = residency_flux,
      potential_drop = potential_drop,
      forward_mask = forward_mask,
      forward_flux_formula = forward_flux_formula,
      forward_flux_model = forward_flux_model,
      P_dynamic_formula = P_dynamic_formula,
      P_dynamic = P_dynamic,
      sinks_formula = sinks_formula,
      sinks_model = sinks_model,
      state_summary = state_summary,
      p_norm_row_mass = aggregate$p_norm_row_mass,
      zero_mass_nodes_formula = formula_dynamic$zero_mass_nodes,
      zero_mass_nodes_model = model_dynamic$zero_mass_nodes
    ),
    fate = list(
      sister = sister_fate,
      terminal = terminal_fate
    ),
    annotation = NULL,
    flow = NULL
  )

  message(sprintf(
    "RunMetacloneDynamics completed: %d metaclones, %d denoised sinks, model_id=%s.",
    length(metaclone_ids), length(sinks_model), model_id
  ))
  seurat_obj
}

#' Infer metaclone paths or graph-theoretic flow
#'
#' @return The modified Seurat object.
RunMetacloneFlow <- function(
    seurat_obj,
    roots,
    sinks = NULL,
    method = c("most_probable", "widest", "max_flow"),
    tolerance = 1e-12
) {
  .iscoral_assert_seurat(seurat_obj)
  .iscoral_require_packages("igraph")
  method <- match.arg(method)
  cache <- .iscoral_get_cache(seurat_obj, c("potential", "model", "fate"))
  model <- cache$model
  terminal <- cache$fate$terminal
  .iscoral_validate_source_model(terminal, model$model_id, "terminal fate")

  metaclone_ids <- rownames(model$P_dynamic)
  roots <- unique(as.character(roots))
  if (length(roots) == 0L) stop("At least one root metaclone is required.", call. = FALSE)
  invalid_roots <- setdiff(roots, metaclone_ids)
  if (length(invalid_roots) > 0L) {
    stop(sprintf("Unknown root metaclones: %s", paste(invalid_roots, collapse = ", ")), call. = FALSE)
  }
  sinks <- unique(as.character(sinks %||% terminal$sink_nodes))
  invalid_sinks <- setdiff(sinks, terminal$sink_nodes)
  if (length(invalid_sinks) > 0L) {
    stop(
      sprintf("Requested sinks are not denoised terminal states: %s",
              paste(invalid_sinks, collapse = ", ")),
      call. = FALSE
    )
  }

  jump_chain <- .iscoral_make_jump_chain(model$P_dynamic, tolerance)
  root_summary <- model$state_summary[roots, , drop = FALSE]
  root_identity_counts <- cache$potential$seed_core_cell_type_counts
  if (!is.null(root_identity_counts)) {
    root_identity_counts <- root_identity_counts[roots, , drop = FALSE]
  }
  results <- list()
  for (root in roots) {
    for (sink in sinks) {
      key <- sprintf("root_%s__sink_%s", root, sink)
      result <- switch(
        method,
        most_probable = .iscoral_most_probable_path(jump_chain, root, sink, tolerance),
        widest = .iscoral_widest_path(jump_chain, root, sink, tolerance),
        max_flow = .iscoral_true_max_flow(model$forward_flux_model, root, sink, tolerance)
      )
      result$method <- method
      result$root <- root
      result$sink <- sink
      result$root_annotation <- unname(model$annotations[root])
      result$sink_annotation <- unname(model$annotations[sink])
      if (!is.null(result$nodes) && length(result$nodes) > 0L) {
        result$state_summary <- model$state_summary[result$nodes, , drop = FALSE]
      }
      results[[key]] <- result
      if (!isTRUE(result$reachable)) {
        warning(sprintf("No %s connection from metaclone %s to %s.", method, root, sink), call. = FALSE)
      }
    }
  }

  seurat_obj@misc$isCORAL_dynamics$flow <- list(
    source_model_id = model$model_id,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    definition = switch(
      method,
      most_probable = "most probable path in the denoised jump chain",
      widest = "path maximizing the minimum denoised jump-chain probability",
      max_flow = "graph-theoretic maximum flow using denoised forward flux as capacity"
    ),
    settings = list(
      method = method,
      roots = roots,
      sinks = sinks,
      tolerance = tolerance
    ),
    root_summary = root_summary,
    root_identity_counts = root_identity_counts,
    jump_chain = jump_chain,
    results = results
  )
  message(sprintf(
    "RunMetacloneFlow completed: method=%s, %d root-sink comparisons.",
    method, length(results)
  ))
  seurat_obj
}
