# Workflow-local plotting functions for isCORAL dynamic inference.
# Source internal and API files before this file.

.iscoral_base_family <- function() {
  getOption("iscoral.base_family", "Arial")
}

# Rasterize a ggplot layer with ggrastr when available (keeps high-density
# cell-point PDFs small and responsive). Falls back to the vector layer with a
# one-time warning; never installs packages. DPI defaults to
# getOption("iscoral.ggrastr_dpi", 300L).
.iscoral_rasterize_points <- function(
    layer,
    dpi = getOption("iscoral.ggrastr_dpi", 300L)
) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    return(ggrastr::rasterise(layer, dpi = dpi))
  }
  if (is.null(getOption("iscoral.ggrastr_warned"))) {
    options(iscoral.ggrastr_warned = TRUE)
    warning(
      "ggrastr is not installed; point layers stay vector ",
      "(no automatic installation).",
      call. = FALSE
    )
  }
  layer
}

.iscoral_theme_void <- function() {
  ggplot2::theme_void(base_family = .iscoral_base_family(), base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "plain", size = 14),
      legend.text = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.7)
    )
}

.iscoral_theme_grid <- function() {
  ggplot2::theme_bw(base_family = .iscoral_base_family(), base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.text = ggplot2::element_text(face = "plain", size = 14),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "black")
    )
}

.iscoral_file_token <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

.iscoral_embedding_value_plot <- function(
    seurat_obj, reduction, values, title, legend_title,
    low = "grey90", high = "#D33F6A", midpoint = NULL,
    limits = NULL
) {
  data <- .iscoral_reduction_data(seurat_obj, reduction)
  data$Value <- unname(values[data$Cell])
  data <- data[order(data$Value, na.last = TRUE), , drop = FALSE]
  plot <- ggplot2::ggplot(data, ggplot2::aes(Dim1, Dim2, color = Value)) +
    .iscoral_rasterize_points(
      ggplot2::geom_point(size = 0.8, alpha = 0.85, shape = 16)
    ) +
    .iscoral_theme_void() +
    ggplot2::labs(title = title, color = legend_title)
  if (is.null(midpoint)) {
    plot <- plot + ggplot2::scale_color_gradient(low = low, high = high, limits = limits)
  } else {
    plot <- plot + ggplot2::scale_color_gradient2(
      low = "#3B4CC0",
      mid = "grey95",
      high = "#B40426",
      midpoint = midpoint,
      limits = limits
    )
  }
  plot
}

.iscoral_centroid_data <- function(seurat_obj, reduction, cluster_col) {
  data <- .iscoral_reduction_data(seurat_obj, reduction)
  data$Metaclone <- as.character(seurat_obj@meta.data[data$Cell, cluster_col])
  stats::aggregate(cbind(Dim1, Dim2) ~ Metaclone, data = data, FUN = mean)
}

.iscoral_add_edge_coordinates <- function(edges, centroids) {
  if (nrow(edges) == 0L) {
    edges$x <- numeric()
    edges$y <- numeric()
    edges$xend <- numeric()
    edges$yend <- numeric()
    return(edges)
  }
  from_index <- match(edges$from, centroids$Metaclone)
  to_index <- match(edges$to, centroids$Metaclone)
  edges$x <- centroids$Dim1[from_index]
  edges$y <- centroids$Dim2[from_index]
  edges$xend <- centroids$Dim1[to_index]
  edges$yend <- centroids$Dim2[to_index]
  edges
}

.iscoral_terminal_plot <- function(
    seurat_obj, reduction, cluster_col, sink_nodes, annotations, metaclone_colors
) {
  data <- .iscoral_reduction_data(seurat_obj, reduction)
  data$Metaclone <- as.character(seurat_obj@meta.data[data$Cell, cluster_col])
  data$Status <- "Background Cells"
  data$PlotOrder <- 1L
  sink_nodes <- .iscoral_sort_ids(sink_nodes)
  for (node in sink_nodes) {
    label <- sprintf("%s\n(MC%s)", annotations[[node]], node)
    selected <- data$Metaclone == node
    data$Status[selected] <- label
    data$PlotOrder[selected] <- 2L
  }
  data <- data[order(data$PlotOrder), , drop = FALSE]
  levels <- c(
    "Background Cells",
    vapply(sink_nodes, function(node) sprintf("%s\n(MC%s)", annotations[[node]], node), character(1))
  )
  data$Status <- factor(data$Status, levels = levels)
  colors <- c("Background Cells" = "#E0E0E0")
  sizes <- c("Background Cells" = 0.3)
  alphas <- c("Background Cells" = 0.5)
  for (node in sink_nodes) {
    label <- sprintf("%s\n(MC%s)", annotations[[node]], node)
    colors[[label]] <- metaclone_colors[[node]]
    sizes[[label]] <- 1.2
    alphas[[label]] <- 1
  }
  ggplot2::ggplot(data, ggplot2::aes(Dim1, Dim2, color = Status)) +
    .iscoral_rasterize_points(
      ggplot2::geom_point(ggplot2::aes(size = Status, alpha = Status), shape = 16)
    ) +
    ggplot2::scale_color_manual(values = colors, name = "Terminal States") +
    ggplot2::scale_size_manual(values = sizes, guide = "none") +
    ggplot2::scale_alpha_manual(values = alphas, guide = "none") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1))) +
    .iscoral_theme_void() +
    ggplot2::labs(title = sprintf("Terminal Metaclones: %s", reduction))
}

.iscoral_probability_panels <- function(
    seurat_obj, reduction, cluster_col, value_matrix, columns,
    annotations, title_prefix, legend_title, metaclone_colors
) {
  data <- .iscoral_reduction_data(seurat_obj, reduction)
  data$Metaclone <- as.character(seurat_obj@meta.data[data$Cell, cluster_col])
  unknown <- setdiff(unique(data$Metaclone), rownames(value_matrix))
  if (length(unknown) > 0L) {
    stop(sprintf("Plot contains metaclones absent from the value matrix: %s",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }
  sink_colors <- metaclone_colors[columns]
  plots <- lapply(columns, function(column) {
    panel <- data
    panel$Probability <- value_matrix[panel$Metaclone, column]
    panel <- panel[order(panel$Probability), , drop = FALSE]
    ggplot2::ggplot(panel, ggplot2::aes(Dim1, Dim2, color = Probability)) +
      .iscoral_rasterize_points(
        ggplot2::geom_point(size = 0.8, alpha = 0.8, shape = 16)
      ) +
      # Fate-priming ramp: the full terminal-metaclone color at probability 1,
      # fading to white at probability 0.
      ggplot2::scale_color_gradient(
        low = "#FFFFFF",
        high = sink_colors[[column]],
        limits = c(0, 1),
        name = legend_title
      ) +
      .iscoral_theme_void() +
      ggplot2::labs(title = sprintf("%s: %s (MC%s)", title_prefix, annotations[[column]], column))
  })
  names(plots) <- columns
  plots
}

# Core plotting views retained from the embryo/HSC workflow.
PlotMetacloneDynamics <- function(
    seurat_obj,
    view = c("potential", "transition", "dag", "terminal", "terminal_fate"),
    reductions = c("umap", "CORAL_cell_embed"),
    seed_metaclones = NULL,
    sister_value = c("forward_mass", "conditional_probability"),
    dag_reduction = "CORAL_clone_embed"
) {
  .iscoral_assert_seurat(seurat_obj)
  .iscoral_require_packages(c("ggplot2", "patchwork", "Seurat", "isCORAL"))
  allowed <- c(
    "potential", "potential_distribution", "seed_identity", "transition",
    "model_diagnostics", "dag", "terminal", "terminal_fate", "sister_fate"
  )
  invalid <- setdiff(view, allowed)
  if (length(invalid) > 0L) {
    stop(sprintf("Unknown dynamic plot view(s): %s", paste(invalid, collapse = ", ")), call. = FALSE)
  }
  cache <- .iscoral_get_cache(seurat_obj, c("potential", "model", "fate"))
  model <- cache$model
  sister <- cache$fate$sister
  terminal <- cache$fate$terminal
  sister_value <- match.arg(sister_value)
  .iscoral_validate_source_model(sister, model$model_id, "sister fate")
  .iscoral_validate_source_model(terminal, model$model_id, "terminal fate")
  cluster_col <- model$settings$cluster_col
  metaclone_ids <- rownames(model$P_dynamic)
  annotations <- model$annotations
  metaclone_colors <- .iscoral_get_metaclone_colors(metaclone_ids)
  plots <- list()

  if ("potential" %in% view) {
    cell_metaclone <- as.character(seurat_obj@meta.data[[cluster_col]])
    metaclone_value <- model$state_summary$Pseudotime
    names(metaclone_value) <- model$state_summary$Metaclone
    metaclone_cell_value <- stats::setNames(
      metaclone_value[cell_metaclone],
      rownames(seurat_obj@meta.data)
    )
    cell_value <- stats::setNames(
      as.numeric(seurat_obj@meta.data$cell_pseudotime),
      rownames(seurat_obj@meta.data)
    )
    for (reduction in reductions) {
      key <- sprintf("potential_%s", .iscoral_file_token(reduction))
      cell_plot <- .iscoral_embedding_value_plot(
        seurat_obj,
        reduction,
        cell_value,
        "Cell-level CytoTRACE2 pseudotime",
        "Pseudotime",
        low = "#440154",
        high = "#FDE725",
        limits = c(0, 1)
      )
      metaclone_plot <- .iscoral_embedding_value_plot(
        seurat_obj,
        reduction,
        metaclone_cell_value,
        "Metaclone-mean pseudotime",
        "Pseudotime",
        low = "#440154",
        high = "#FDE725",
        limits = c(0, 1)
      )
      plots[[key]] <- patchwork::wrap_plots(
        cell_plot,
        metaclone_plot,
        ncol = 2,
        guides = "collect"
      ) + patchwork::plot_annotation(
        title = sprintf("Developmental timeline: %s", reduction)
      )
    }
  }

  if ("potential_distribution" %in% view) {
    distribution_data <- cache$potential$seed_core_identity
    distribution_data$Metaclone <- as.character(distribution_data$Metaclone)
    upstream_order <- model$state_summary$Metaclone[order(
      model$state_summary$DevelopmentalPotential,
      decreasing = TRUE
    )]
    metaclone_labels <- stats::setNames(
      sprintf("MC%s | %s", metaclone_ids, annotations[metaclone_ids]),
      metaclone_ids
    )
    distribution_data$MetacloneLabel <- factor(
      unname(metaclone_labels[distribution_data$Metaclone]),
      levels = rev(unname(metaclone_labels[upstream_order]))
    )
    plots$potential_distribution <- ggplot2::ggplot(
      distribution_data,
      ggplot2::aes(
        x = MetacloneLabel,
        y = DevelopmentalPotential,
        fill = Metaclone
      )
    ) +
      ggplot2::geom_violin(
        scale = "width", trim = FALSE, width = 0.85,
        color = "#333333", linewidth = 0.35
      ) +
      ggplot2::geom_boxplot(
        width = 0.16, outlier.shape = NA,
        fill = "white", color = "#333333", linewidth = 0.4
      ) +
      ggplot2::stat_summary(
        fun = mean, geom = "point", shape = 23,
        size = 2.4, fill = "black", color = "white", stroke = 0.4
      ) +
      ggplot2::scale_fill_manual(values = metaclone_colors, guide = "none") +
      ggplot2::scale_y_continuous(
        limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      ggplot2::coord_flip() +
      .iscoral_theme_grid() +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 12)
      ) +
      ggplot2::labs(
        title = "Seed-core CytoTRACE2 Developmental Potential",
        x = "Metaclone (ordered from downstream to upstream on the axis)",
        y = "CytoTRACE2 developmental potential"
      )
  }

  if ("seed_identity" %in% view) {
    identity_data <- cache$potential$seed_core_identity
    identity_data$Metaclone <- as.character(identity_data$Metaclone)
    identity_data$CellTypePlot <- as.character(identity_data$CellType)
    missing_identity <- is.na(identity_data$CellTypePlot) |
      !nzchar(identity_data$CellTypePlot)
    identity_data$CellTypePlot[missing_identity] <- "Missing metadata"
    identity_counts <- as.data.frame(table(
      Metaclone = identity_data$Metaclone,
      CellType = identity_data$CellTypePlot
    ), stringsAsFactors = FALSE)
    identity_counts <- identity_counts[identity_counts$Freq > 0, , drop = FALSE]
    identity_counts$Proportion <- ave(
      identity_counts$Freq,
      identity_counts$Metaclone,
      FUN = function(x) x / sum(x)
    )
    upstream_order <- model$state_summary$Metaclone[order(
      model$state_summary$DevelopmentalPotential,
      decreasing = TRUE
    )]
    identity_counts$Metaclone <- factor(
      identity_counts$Metaclone,
      levels = upstream_order
    )
    cell_type_levels <- sort(unique(identity_counts$CellType))
    cell_type_colors <- stats::setNames(
      grDevices::hcl.colors(length(cell_type_levels), palette = "Dark 3"),
      cell_type_levels
    )
    plots$seed_identity <- ggplot2::ggplot(
      identity_counts,
      ggplot2::aes(Metaclone, Proportion, fill = CellType)
    ) +
      ggplot2::geom_col(width = 0.78, color = "black", linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = cell_type_colors, name = "Seed-core identity") +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        labels = function(x) paste0(round(100 * x), "%"),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::theme_classic(base_family = .iscoral_base_family(), base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = ggplot2::element_text(color = "black"),
        legend.position = "right"
      ) +
      ggplot2::labs(
        title = "Seed-core Cell Identity Composition",
        x = "Seed-core metaclone",
        y = "Cell-type proportion"
      )
  }

  if ("transition" %in% view) {
    formula_data <- as.data.frame(as.table(model$P_dynamic_formula), stringsAsFactors = FALSE)
    names(formula_data) <- c("From", "To", "Probability")
    formula_data$Matrix <- "Direction-constrained formula"
    model_data <- as.data.frame(as.table(model$P_dynamic), stringsAsFactors = FALSE)
    names(model_data) <- c("From", "To", "Probability")
    model_data$Matrix <- "Sparsified model"
    transition_data <- rbind(formula_data, model_data)
    transition_data$From <- factor(transition_data$From, levels = rev(metaclone_ids))
    transition_data$To <- factor(transition_data$To, levels = metaclone_ids)
    plots$transition <- ggplot2::ggplot(
      transition_data,
      ggplot2::aes(To, From, fill = Probability)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.15) +
      ggplot2::facet_wrap(~Matrix, ncol = 2) +
      ggplot2::scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1)) +
      ggplot2::coord_equal() +
      .iscoral_theme_grid() +
      ggplot2::labs(
        title = "Metaclone Dynamic Transition Matrices",
        x = "Target metaclone", y = "Seed-core metaclone", fill = "Probability"
      )
  }

  if ("model_diagnostics" %in% view) {
    state_data <- model$state_summary
    metaclone_order <- state_data$Metaclone[order(state_data$Pseudotime)]
    state_data$Metaclone <- factor(state_data$Metaclone, levels = metaclone_order)
    flux_data <- rbind(
      data.frame(
        Metaclone = state_data$Metaclone,
        Component = "Residency",
        Flux = state_data$ResidencyFlux
      ),
      data.frame(
        Metaclone = state_data$Metaclone,
        Component = "Retained forward flux",
        Flux = state_data$ForwardFluxAfterSparsification
      ),
      data.frame(
        Metaclone = state_data$Metaclone,
        Component = "Removed by sparsification",
        Flux = state_data$SparsificationRemovedFlux
      ),
      data.frame(
        Metaclone = state_data$Metaclone,
        Component = "Discarded by direction",
        Flux = state_data$DirectionallyDiscardedFlux
      )
    )
    flux_data$Component <- factor(
      flux_data$Component,
      levels = c(
        "Residency", "Retained forward flux",
        "Removed by sparsification", "Discarded by direction"
      )
    )
    flux_colors <- c(
      "Residency" = "#4C78A8",
      "Retained forward flux" = "#59A14F",
      "Removed by sparsification" = "#F28E2B",
      "Discarded by direction" = "#B9B9B9"
    )
    flux_plot <- ggplot2::ggplot(
      flux_data,
      ggplot2::aes(Metaclone, Flux, fill = Component)
    ) +
      ggplot2::geom_col(width = 0.78, color = "black", linewidth = 0.25) +
      ggplot2::scale_fill_manual(values = flux_colors, name = "Flux component") +
      .iscoral_theme_grid() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      ggplot2::labs(
        title = "Metaclone Flux Budget",
        x = "Metaclone (upstream to downstream)",
        y = "Average P_norm probability mass"
      )

    diagnostic_data <- rbind(
      data.frame(
        Metaclone = metaclone_order,
        Diagnostic = "Formula sink",
        Flag = metaclone_order %in% model$sinks_formula
      ),
      data.frame(
        Metaclone = metaclone_order,
        Diagnostic = "Denoised-model sink",
        Flag = metaclone_order %in% model$sinks_model
      ),
      data.frame(
        Metaclone = metaclone_order,
        Diagnostic = "Zero retained denominator",
        Flag = metaclone_order %in% model$zero_mass_nodes_model
      )
    )
    diagnostic_data$Metaclone <- factor(
      diagnostic_data$Metaclone,
      levels = metaclone_order
    )
    diagnostic_data$Diagnostic <- factor(
      diagnostic_data$Diagnostic,
      levels = rev(c(
        "Formula sink", "Denoised-model sink", "Zero retained denominator"
      ))
    )
    sink_plot <- ggplot2::ggplot(
      diagnostic_data,
      ggplot2::aes(Metaclone, Diagnostic, fill = Flag)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.6) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "#E0E4E7", "TRUE" = "#C74B46"),
        name = "Flag"
      ) +
      .iscoral_theme_grid() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        legend.position = "right"
      ) +
      ggplot2::labs(x = "Metaclone", y = NULL)

    plots$model_diagnostics <- flux_plot / sink_plot +
      patchwork::plot_layout(heights = c(3, 1))
  }

  if ("dag" %in% view) {
    if (length(dag_reduction) != 1L || is.na(dag_reduction) ||
        !nzchar(dag_reduction)) {
      stop("dag_reduction must be one non-missing reduction name.",
           call. = FALSE)
    }
    # Dag view style constants (adjustable here).
    dag_point_size <- 9
    dag_label_size <- 5.5
    dag_curvature <- 0.08
    dag_edge_color <- "#888888"
    dag_arrow_length_mm <- 2.5
    dag_linewidth_range <- c(0.8, 3.5)
    dag_alpha_range <- c(0.25, 0.75)

    centroids <- .iscoral_centroid_data(
      seurat_obj,
      dag_reduction,
      cluster_col
    )
    edges <- .iscoral_matrix_edges(model$P_dynamic, "Probability", model$settings$tolerance)
    edges <- .iscoral_add_edge_coordinates(edges, centroids)
    centroids$Pseudotime <- model$state_summary[centroids$Metaclone, "Pseudotime"]
    plots$dag <- ggplot2::ggplot() +
      ggplot2::geom_curve(
        data = edges,
        ggplot2::aes(x, y, xend = xend, yend = yend, linewidth = Probability, alpha = Probability),
        curvature = dag_curvature,
        color = dag_edge_color,
        arrow = grid::arrow(length = grid::unit(dag_arrow_length_mm, "mm"), type = "closed")
      ) +
      ggplot2::geom_point(
        data = centroids,
        ggplot2::aes(Dim1, Dim2, fill = Metaclone),
        shape = 21, size = dag_point_size, color = "black", stroke = 0.7
      ) +
      ggplot2::geom_text(
        data = centroids,
        ggplot2::aes(Dim1, Dim2, label = Metaclone),
        size = dag_label_size, family = .iscoral_base_family()
      ) +
      ggplot2::scale_fill_manual(values = metaclone_colors, guide = "none") +
      ggplot2::scale_linewidth_continuous(range = dag_linewidth_range, guide = "none") +
      ggplot2::scale_alpha_continuous(range = dag_alpha_range, guide = "none") +
      .iscoral_theme_void() +
      ggplot2::labs(
        title = "Denoised Metaclone Dynamic Graph"
      )
  }

  if ("terminal" %in% view) {
    terminal_plots <- lapply(reductions, function(reduction) {
      .iscoral_terminal_plot(
        seurat_obj, reduction, cluster_col, terminal$sink_nodes,
        annotations, metaclone_colors
      )
    })
    names(terminal_plots) <- reductions
    plots$terminal <- patchwork::wrap_plots(terminal_plots, ncol = length(terminal_plots)) +
      patchwork::plot_annotation(title = "Denoised Terminal Metaclone Localization")
  }

  if ("terminal_fate" %in% view) {
    for (reduction in reductions) {
      panels <- .iscoral_probability_panels(
        seurat_obj, reduction, cluster_col, terminal$fate_matrix,
        terminal$sink_nodes, annotations, "Terminal fate", "Probability",
        metaclone_colors
      )
      key <- sprintf("terminal_fate_%s", .iscoral_file_token(reduction))
      plots[[key]] <- patchwork::wrap_plots(panels, ncol = min(3L, length(panels))) +
        patchwork::plot_annotation(
          title = sprintf("Metaclone-level Terminal Fate Landscapes: %s", reduction)
        )
    }
  }

  if ("sister_fate" %in% view) {
    if (is.null(seed_metaclones) || length(seed_metaclones) == 0L) {
      if (!is.null(cache$flow$settings$roots) &&
          length(cache$flow$settings$roots) > 0L) {
        seed_metaclones <- cache$flow$settings$roots
      } else {
        message(
          "[PlotMetacloneDynamics] sister_fate view skipped: no seed metaclones ",
          "configured (set CONFIG$sister_fate_seeds before inference)."
        )
        seed_metaclones <- NULL
      }
    }
  }
  if ("sister_fate" %in% view && !is.null(seed_metaclones) &&
      length(seed_metaclones) > 0L) {
    seed_metaclones <- unique(as.character(seed_metaclones))
    sister_matrix <- sister[[sister_value]]
    invalid_seed <- setdiff(seed_metaclones, rownames(sister_matrix))
    if (length(invalid_seed) > 0L) {
      stop(sprintf("Unknown seed metaclones: %s", paste(invalid_seed, collapse = ", ")), call. = FALSE)
    }
    for (reduction in reductions) {
      data <- .iscoral_reduction_data(seurat_obj, reduction)
      data$Metaclone <- as.character(seurat_obj@meta.data[data$Cell, cluster_col])
      panels <- lapply(seed_metaclones, function(seed) {
        panel <- data
        panel$Probability <- sister_matrix[seed, panel$Metaclone]
        panel <- panel[order(panel$Probability), , drop = FALSE]
        ggplot2::ggplot(panel, ggplot2::aes(Dim1, Dim2, color = Probability)) +
          .iscoral_rasterize_points(
            ggplot2::geom_point(size = 0.8, alpha = 0.8, shape = 16)
          ) +
          ggplot2::scale_color_viridis_c(
            option = "magma", limits = c(0, 1),
            name = if (identical(sister_value, "forward_mass")) "Outward\nflux" else "Conditional\nprobability"
          ) +
          .iscoral_theme_void() +
          ggplot2::labs(title = sprintf("Seed %s (MC%s)", annotations[[seed]], seed))
      })
      key <- sprintf("sister_fate_%s", .iscoral_file_token(reduction))
      plots[[key]] <- patchwork::wrap_plots(panels, ncol = min(3L, length(panels))) +
        patchwork::plot_annotation(
          title = sprintf("Direct Relevant/Sister Fate: %s", reduction)
        )
    }
  }

  plots
}

.iscoral_flow_nodes <- function(result) {
  if (identical(result$method, "max_flow")) {
    unique(c(result$root, result$flow_edges$from, result$flow_edges$to, result$sink))
  } else {
    result$nodes %||% character()
  }
}

.iscoral_flow_network_plot <- function(
    result, model, centroids, metaclone_colors
) {
  background_edges <- .iscoral_matrix_edges(model$P_dynamic, "value", model$settings$tolerance)
  background_edges <- .iscoral_add_edge_coordinates(background_edges, centroids)
  if (identical(result$method, "max_flow")) {
    selected_edges <- result$flow_edges
    if (nrow(selected_edges) > 0L) names(selected_edges)[names(selected_edges) == "flow"] <- "value"
  } else if (length(result$nodes) >= 2L) {
    selected_edges <- data.frame(
      from = result$nodes[-length(result$nodes)],
      to = result$nodes[-1L],
      value = result$edge_probability,
      stringsAsFactors = FALSE
    )
  } else {
    selected_edges <- data.frame(from = character(), to = character(), value = numeric())
  }
  if (!"value" %in% names(selected_edges)) selected_edges$value <- numeric(nrow(selected_edges))
  selected_edges <- .iscoral_add_edge_coordinates(selected_edges, centroids)
  selected_nodes <- .iscoral_flow_nodes(result)
  centroids$Selected <- centroids$Metaclone %in% selected_nodes

  ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = background_edges,
      ggplot2::aes(x, y, xend = xend, yend = yend),
      curvature = 0.10, color = "grey82", linewidth = 0.35,
      arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
    ) +
    ggplot2::geom_curve(
      data = selected_edges,
      ggplot2::aes(x, y, xend = xend, yend = yend, linewidth = value),
      curvature = 0.10, color = "#B2182B", alpha = 0.9,
      arrow = grid::arrow(length = grid::unit(3.0, "mm"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = centroids,
      ggplot2::aes(Dim1, Dim2, fill = Metaclone, size = Selected),
      shape = 21, color = "black", stroke = 0.6
    ) +
    ggplot2::geom_text(
      data = centroids,
      ggplot2::aes(Dim1, Dim2, label = Metaclone),
      size = 5, family = .iscoral_base_family()
    ) +
    ggplot2::scale_fill_manual(values = metaclone_colors, guide = "none") +
    ggplot2::scale_size_manual(values = c("FALSE" = 6, "TRUE" = 9), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.7, 3), guide = "none") +
    .iscoral_theme_void() +
    ggplot2::labs(
      title = sprintf("%s: MC%s to MC%s", result$method, result$root, result$sink)
    )
}
