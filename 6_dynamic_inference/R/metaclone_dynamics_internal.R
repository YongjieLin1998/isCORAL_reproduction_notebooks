# Internal utilities for the portable isCORAL dynamic-inference workflow.
# Functions in this file are implementation details, not user-facing APIs.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.iscoral_required_version <- "6.1.4"

.iscoral_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      sprintf("Required package(s) not installed: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.iscoral_assert_seurat <- function(seurat_obj) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  invisible(TRUE)
}

.iscoral_sort_ids <- function(ids) {
  ids <- unique(as.character(ids))
  numeric_ids <- suppressWarnings(as.numeric(ids))
  if (length(ids) > 0L && all(is.finite(numeric_ids))) {
    ids[order(numeric_ids, ids)]
  } else {
    sort(ids)
  }
}

.iscoral_validate_square_matrix <- function(x, label) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(sprintf("%s must be a matrix.", label), call. = FALSE)
  }
  if (nrow(x) != ncol(x)) {
    stop(sprintf("%s must be square.", label), call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      !identical(rownames(x), colnames(x))) {
    stop(sprintf("%s must have identical row and column names.", label), call. = FALSE)
  }
  invisible(TRUE)
}

.iscoral_validate_named_numeric <- function(x, required_names, label) {
  if (!is.numeric(x) || is.null(names(x))) {
    stop(sprintf("%s must be a named numeric vector.", label), call. = FALSE)
  }
  missing <- setdiff(required_names, names(x))
  if (length(missing) > 0L) {
    stop(sprintf("%s is missing: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
  }
  x <- x[required_names]
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("%s contains missing or non-finite values.", label), call. = FALSE)
  }
  x
}

.iscoral_model_id <- function(payload) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(payload, path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

.iscoral_get_cache <- function(seurat_obj, require = character()) {
  .iscoral_assert_seurat(seurat_obj)
  cache <- seurat_obj@misc$isCORAL_dynamics
  if (is.null(cache)) {
    stop("isCORAL dynamic cache was not found. Run RunMetacloneDynamics() first.", call. = FALSE)
  }
  missing <- require[!vapply(require, function(x) !is.null(cache[[x]]), logical(1))]
  if (length(missing) > 0L) {
    stop(sprintf("Dynamic cache is missing: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  cache
}

.iscoral_validate_source_model <- function(component, model_id, label) {
  if (is.null(component)) return(invisible(FALSE))
  source_id <- component$source_model_id
  if (is.null(source_id) || !identical(source_id, model_id)) {
    stop(sprintf("Cached %s does not belong to the current dynamic model.", label), call. = FALSE)
  }
  invisible(TRUE)
}

.iscoral_align_p_norm <- function(seurat_obj) {
  .iscoral_require_packages("Matrix")
  embedding_cache <- seurat_obj@misc$isCORAL_embedding
  if (is.null(embedding_cache$P_norm) || is.null(embedding_cache$valid_cells)) {
    stop(
      "P_norm or valid_cells was not found in misc$isCORAL_embedding. Run CORALembedding() first.",
      call. = FALSE
    )
  }

  P_norm <- embedding_cache$P_norm
  valid_cells <- as.character(embedding_cache$valid_cells)
  if (length(valid_cells) == 0L || anyNA(valid_cells) || any(!nzchar(valid_cells)) ||
      anyDuplicated(valid_cells)) {
    stop("misc$isCORAL_embedding$valid_cells is invalid.", call. = FALSE)
  }
  if (nrow(P_norm) != length(valid_cells) || ncol(P_norm) != length(valid_cells)) {
    stop("P_norm dimensions do not match valid_cells.", call. = FALSE)
  }

  if (is.null(rownames(P_norm)) && is.null(colnames(P_norm))) {
    rownames(P_norm) <- colnames(P_norm) <- valid_cells
  } else if (is.null(rownames(P_norm))) {
    rownames(P_norm) <- colnames(P_norm)
  } else if (is.null(colnames(P_norm))) {
    colnames(P_norm) <- rownames(P_norm)
  }
  .iscoral_validate_square_matrix(P_norm, "P_norm")

  if (!setequal(valid_cells, rownames(P_norm))) {
    stop("P_norm and valid_cells do not contain the same barcodes.", call. = FALSE)
  }
  if (!all(valid_cells %in% colnames(seurat_obj))) {
    stop("Some P_norm cells are absent from the Seurat object.", call. = FALSE)
  }

  P_norm <- P_norm[valid_cells, valid_cells, drop = FALSE]
  values <- if (inherits(P_norm, "sparseMatrix")) P_norm@x else as.numeric(P_norm)
  if (anyNA(values) || any(!is.finite(values)) || any(values < -1e-12)) {
    stop("P_norm contains invalid probability values.", call. = FALSE)
  }

  row_mass <- as.numeric(Matrix::rowSums(P_norm))
  names(row_mass) <- valid_cells
  positive <- row_mass > 1e-12
  if (any(positive) && max(abs(row_mass[positive] - 1)) > 1e-6) {
    warning("Non-zero P_norm rows are not row-stochastic within tolerance.", call. = FALSE)
  }

  list(P_norm = P_norm, valid_cells = valid_cells, row_mass = row_mass)
}

.iscoral_aggregate_metaclone_flux <- function(seurat_obj, cluster_col) {
  aligned <- .iscoral_align_p_norm(seurat_obj)
  valid_cells <- aligned$valid_cells
  metaclone <- as.character(seurat_obj@meta.data[valid_cells, cluster_col])
  if (anyNA(metaclone) || any(!nzchar(metaclone))) {
    stop("One or more P_norm-aligned cells lack a metaclone label.", call. = FALSE)
  }

  metaclone_ids <- .iscoral_sort_ids(metaclone)
  metaclone_factor <- factor(metaclone, levels = metaclone_ids)
  membership <- Matrix::sparseMatrix(
    i = seq_along(valid_cells),
    j = as.integer(metaclone_factor),
    x = 1,
    dims = c(length(valid_cells), length(metaclone_ids)),
    dimnames = list(valid_cells, metaclone_ids)
  )
  core_size <- as.numeric(Matrix::colSums(membership))
  names(core_size) <- metaclone_ids
  if (any(core_size <= 0)) {
    stop("At least one metaclone has an empty seed core.", call. = FALSE)
  }

  block_sum <- as.matrix(Matrix::t(membership) %*% aligned$P_norm %*% membership)
  dimnames(block_sum) <- list(metaclone_ids, metaclone_ids)
  block_flux <- sweep(block_sum, 1, core_size, "/")
  dimnames(block_flux) <- list(metaclone_ids, metaclone_ids)
  core_cells <- split(valid_cells, metaclone_factor, drop = TRUE)[metaclone_ids]

  list(
    valid_cells = valid_cells,
    cell_to_metaclone = stats::setNames(metaclone, valid_cells),
    core_cells = core_cells,
    core_size = core_size,
    block_sum = block_sum,
    block_flux = block_flux,
    p_norm_row_mass = aligned$row_mass
  )
}

.iscoral_make_dynamic_matrix <- function(residency, forward_flux, tolerance) {
  metaclone_ids <- rownames(forward_flux)
  denominator <- residency + rowSums(forward_flux)
  P_dynamic <- matrix(
    0,
    nrow = length(metaclone_ids),
    ncol = length(metaclone_ids),
    dimnames = list(metaclone_ids, metaclone_ids)
  )
  zero_mass_nodes <- character()

  for (metaclone_id in metaclone_ids) {
    if (denominator[[metaclone_id]] > tolerance) {
      P_dynamic[metaclone_id, ] <-
        forward_flux[metaclone_id, ] / denominator[[metaclone_id]]
      P_dynamic[metaclone_id, metaclone_id] <-
        residency[[metaclone_id]] / denominator[[metaclone_id]]
    } else {
      P_dynamic[metaclone_id, metaclone_id] <- 1
      zero_mass_nodes <- c(zero_mass_nodes, metaclone_id)
    }
  }

  if (any(abs(rowSums(P_dynamic) - 1) > 1e-8)) {
    stop("The dynamic transition matrix is not row-stochastic.", call. = FALSE)
  }
  list(matrix = P_dynamic, zero_mass_nodes = zero_mass_nodes)
}

.iscoral_compute_absorption <- function(P_dynamic, tolerance = 1e-10) {
  .iscoral_validate_square_matrix(P_dynamic, "P_dynamic")
  metaclone_ids <- rownames(P_dynamic)
  if (anyNA(P_dynamic) || any(!is.finite(P_dynamic)) || any(P_dynamic < -tolerance)) {
    stop("P_dynamic contains invalid probabilities.", call. = FALSE)
  }
  if (any(abs(rowSums(P_dynamic) - 1) > 1e-8)) {
    stop("P_dynamic must be row-stochastic.", call. = FALSE)
  }

  off_diagonal <- P_dynamic
  diag(off_diagonal) <- 0
  sink_nodes <- metaclone_ids[rowSums(off_diagonal) <= tolerance]
  transient_nodes <- setdiff(metaclone_ids, sink_nodes)
  if (length(sink_nodes) == 0L) {
    stop("No absorbing metaclone was found in denoised P_dynamic.", call. = FALSE)
  }

  fate_matrix <- matrix(
    0,
    nrow = length(metaclone_ids),
    ncol = length(sink_nodes),
    dimnames = list(metaclone_ids, sink_nodes)
  )
  for (sink_node in sink_nodes) {
    fate_matrix[sink_node, sink_node] <- 1
  }
  Q <- matrix(numeric(), nrow = 0, ncol = 0)
  A <- matrix(numeric(), nrow = 0, ncol = length(sink_nodes))
  N <- matrix(numeric(), nrow = 0, ncol = 0)
  expected_steps <- stats::setNames(rep(0, length(metaclone_ids)), metaclone_ids)

  if (length(transient_nodes) > 0L) {
    Q <- P_dynamic[transient_nodes, transient_nodes, drop = FALSE]
    A <- P_dynamic[transient_nodes, sink_nodes, drop = FALSE]
    system_matrix <- diag(length(transient_nodes)) - Q
    reciprocal_condition <- rcond(system_matrix)
    if (!is.finite(reciprocal_condition) || reciprocal_condition < 1e-12) {
      stop("The transient AMC system is singular or nearly singular.", call. = FALSE)
    }
    N <- solve(system_matrix)
    absorption <- N %*% A
    if (any(absorption < -1e-8) || any(absorption > 1 + 1e-8)) {
      stop("Absorption probabilities fall outside [0, 1].", call. = FALSE)
    }
    absorption <- pmin(pmax(absorption, 0), 1)
    absorption_mass <- rowSums(absorption)
    if (any(abs(absorption_mass - 1) > 1e-6)) {
      stop("One or more metaclones do not fully absorb into denoised sinks.", call. = FALSE)
    }
    fate_matrix[transient_nodes, ] <- sweep(absorption, 1, absorption_mass, "/")
    expected_steps[transient_nodes] <- rowSums(N)
  }

  list(
    sink_nodes = sink_nodes,
    transient_nodes = transient_nodes,
    fate_matrix = fate_matrix,
    transient_matrix = Q,
    absorbing_matrix = A,
    fundamental_matrix = N,
    expected_steps_to_absorption = expected_steps
  )
}

.iscoral_make_sister_fate <- function(forward_flux_model, residency_flux) {
  forward_mass <- forward_flux_model
  diag(forward_mass) <- 0
  total_forward_mass <- rowSums(forward_mass)
  conditional_probability <- matrix(
    0,
    nrow = nrow(forward_mass),
    ncol = ncol(forward_mass),
    dimnames = dimnames(forward_mass)
  )
  positive <- total_forward_mass > 0
  if (any(positive)) {
    conditional_probability[positive, ] <-
      sweep(forward_mass[positive, , drop = FALSE], 1, total_forward_mass[positive], "/")
  }
  list(
    definition = "denoised_directionally_admissible_outward_flux",
    forward_mass = forward_mass,
    conditional_probability = conditional_probability,
    total_forward_mass = total_forward_mass,
    residency_flux = residency_flux
  )
}

.iscoral_resolve_annotations <- function(metaclone_ids, annotations = NULL) {
  if (is.null(annotations)) {
    return(stats::setNames(paste("Metaclone", metaclone_ids), metaclone_ids))
  }
  annotations <- as.character(annotations)
  if (is.null(names(annotations))) {
    if (length(annotations) != length(metaclone_ids)) {
      stop("Unnamed annotations must match the number of metaclones.", call. = FALSE)
    }
    names(annotations) <- metaclone_ids
  }
  missing <- setdiff(metaclone_ids, names(annotations))
  if (length(missing) > 0L) {
    stop(sprintf("Annotations are missing metaclones: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  annotations[metaclone_ids]
}

.iscoral_get_metaclone_colors <- function(metaclone_ids) {
  configured_colors <- getOption("iscoral.dynamic.metaclone_colors", NULL)
  if (!is.null(configured_colors)) {
    if (is.null(names(configured_colors)) ||
        !all(metaclone_ids %in% names(configured_colors))) {
      stop("The configured metaclone palette does not cover all observed IDs.", call. = FALSE)
    }
    return(configured_colors[metaclone_ids])
  }
  if (!requireNamespace("isCORAL", quietly = TRUE)) {
    stop("isCORAL v6.1.4 must be loaded before plotting.", call. = FALSE)
  }
  colors <- isCORAL::GetMetacloneColors(
    n = length(metaclone_ids),
    cluster_names = metaclone_ids
  )
  if (is.null(names(colors))) names(colors) <- metaclone_ids
  if (!all(metaclone_ids %in% names(colors))) {
    stop("GetMetacloneColors() returned an incomplete palette.", call. = FALSE)
  }
  colors[metaclone_ids]
}

.iscoral_reduction_data <- function(seurat_obj, reduction, extra = NULL) {
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' was not found.", reduction), call. = FALSE)
  }
  embedding <- Seurat::Embeddings(seurat_obj, reduction)
  metadata <- seurat_obj@meta.data
  if (is.null(rownames(embedding)) || anyDuplicated(rownames(embedding)) ||
      !all(rownames(embedding) %in% rownames(metadata))) {
    stop(
      sprintf("Reduction '%s' cannot be aligned to metadata by cell barcode.", reduction),
      call. = FALSE
    )
  }
  output <- data.frame(
    Cell = rownames(embedding),
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    stringsAsFactors = FALSE,
    row.names = rownames(embedding)
  )
  if (!is.null(extra)) output <- cbind(output, extra[rownames(output), , drop = FALSE])
  output
}

.iscoral_matrix_edges <- function(x, value_name = "probability", tolerance = 0,
                                  include_diagonal = FALSE) {
  .iscoral_validate_square_matrix(x, deparse(substitute(x)))
  mask <- x > tolerance
  if (!include_diagonal) diag(mask) <- FALSE
  index <- which(mask, arr.ind = TRUE)
  if (nrow(index) == 0L) {
    output <- data.frame(from = character(), to = character(), value = numeric())
  } else {
    output <- data.frame(
      from = rownames(x)[index[, "row"]],
      to = colnames(x)[index[, "col"]],
      value = x[index],
      stringsAsFactors = FALSE
    )
  }
  names(output)[3] <- value_name
  output
}

.iscoral_make_jump_chain <- function(P_dynamic, tolerance = 1e-12) {
  jump_chain <- P_dynamic
  diag(jump_chain) <- 0
  outgoing <- rowSums(jump_chain)
  positive <- outgoing > tolerance
  if (any(positive)) {
    jump_chain[positive, ] <-
      sweep(jump_chain[positive, , drop = FALSE], 1, outgoing[positive], "/")
  }
  jump_chain[!positive, ] <- 0
  jump_chain
}

.iscoral_path_edge_probabilities <- function(matrix, nodes) {
  if (length(nodes) < 2L) return(numeric())
  row_index <- match(nodes[-length(nodes)], rownames(matrix))
  col_index <- match(nodes[-1L], colnames(matrix))
  matrix[cbind(row_index, col_index)]
}

.iscoral_most_probable_path <- function(jump_chain, from, to, tolerance = 1e-12) {
  .iscoral_require_packages("igraph")
  ids <- rownames(jump_chain)
  if (identical(from, to)) {
    return(list(
      reachable = TRUE, nodes = from, edge_probability = numeric(),
      path_probability = 1, bottleneck_probability = 1
    ))
  }
  edges <- .iscoral_matrix_edges(jump_chain, "probability", tolerance)
  graph <- igraph::graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = ids, stringsAsFactors = FALSE)
  )
  igraph::E(graph)$cost <- -log(pmax(igraph::E(graph)$probability, .Machine$double.xmin))
  path <- igraph::shortest_paths(
    graph,
    from = from,
    to = to,
    mode = "out",
    weights = igraph::E(graph)$cost,
    output = "vpath"
  )$vpath[[1]]
  if (length(path) == 0L) {
    return(list(reachable = FALSE, nodes = character(), edge_probability = numeric()))
  }
  nodes <- as.character(igraph::as_ids(path))
  probability <- .iscoral_path_edge_probabilities(jump_chain, nodes)
  list(
    reachable = TRUE,
    nodes = nodes,
    edge_probability = probability,
    path_probability = prod(probability),
    bottleneck_probability = if (length(probability)) min(probability) else 1
  )
}

.iscoral_widest_path <- function(jump_chain, from, to, tolerance = 1e-12) {
  ids <- rownames(jump_chain)
  if (identical(from, to)) {
    return(list(
      reachable = TRUE, nodes = from, edge_probability = numeric(),
      path_probability = 1, bottleneck_probability = 1
    ))
  }
  width <- stats::setNames(rep(0, length(ids)), ids)
  previous <- stats::setNames(rep(NA_character_, length(ids)), ids)
  visited <- stats::setNames(rep(FALSE, length(ids)), ids)
  width[[from]] <- Inf

  repeat {
    candidates <- ids[!visited]
    if (length(candidates) == 0L) break
    current <- candidates[which.max(width[candidates])]
    if (width[[current]] <= tolerance || identical(current, to)) break
    visited[[current]] <- TRUE
    neighbors <- ids[jump_chain[current, ] > tolerance]
    for (neighbor in neighbors) {
      candidate_width <- min(width[[current]], jump_chain[current, neighbor])
      if (candidate_width > width[[neighbor]]) {
        width[[neighbor]] <- candidate_width
        previous[[neighbor]] <- current
      }
    }
  }

  if (width[[to]] <= tolerance) {
    return(list(reachable = FALSE, nodes = character(), edge_probability = numeric()))
  }
  nodes <- to
  while (!identical(nodes[[1]], from)) {
    parent <- previous[[nodes[[1]]]]
    if (is.na(parent)) {
      return(list(reachable = FALSE, nodes = character(), edge_probability = numeric()))
    }
    nodes <- c(parent, nodes)
  }
  probability <- .iscoral_path_edge_probabilities(jump_chain, nodes)
  list(
    reachable = TRUE,
    nodes = nodes,
    edge_probability = probability,
    path_probability = prod(probability),
    bottleneck_probability = min(probability)
  )
}

.iscoral_true_max_flow <- function(forward_flux_model, from, to, tolerance = 1e-12) {
  .iscoral_require_packages("igraph")
  ids <- rownames(forward_flux_model)
  if (identical(from, to)) {
    return(list(
      reachable = TRUE,
      max_flow_value = 0,
      flow_edges = data.frame(
        from = character(), to = character(), capacity = numeric(), flow = numeric(),
        stringsAsFactors = FALSE
      )
    ))
  }
  edges <- .iscoral_matrix_edges(forward_flux_model, "capacity", tolerance)
  graph <- igraph::graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = ids, stringsAsFactors = FALSE)
  )
  result <- igraph::max_flow(
    graph,
    source = from,
    target = to,
    capacity = igraph::E(graph)$capacity
  )
  edges$flow <- as.numeric(result$flow)
  flow_edges <- edges[edges$flow > tolerance, , drop = FALSE]
  list(
    reachable = result$value > tolerance,
    max_flow_value = as.numeric(result$value),
    flow_edges = flow_edges,
    cut = result$cut,
    partition1 = as.character(igraph::as_ids(result$partition1)),
    partition2 = as.character(igraph::as_ids(result$partition2))
  )
}
