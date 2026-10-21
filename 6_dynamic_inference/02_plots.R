# Core figures from the completed dynamic model.
# Run the setup below once, then run each print/save block manually in RStudio.
# Resize the Plot pane between print() and the following save block if needed.
# Sourcing the entire script previews and exports every available figure.

if (!exists("sample", inherits = FALSE) || !inherits(sample, "Seurat") ||
    is.null(sample@misc$isCORAL_dynamics$model) ||
    !exists("CONFIG") || !exists("OUTPUT_DIR") || !exists("WORKFLOW_DIR")) {
  stop("Run 01_inference.R first, or restore sample, CONFIG, OUTPUT_DIR and WORKFLOW_DIR.", call. = FALSE)
}
source(file.path(WORKFLOW_DIR, "R", "metaclone_dynamics_internal.R"), encoding = "UTF-8")
.iscoral_require_packages(c("ggplot2", "patchwork", "Seurat", "isCORAL"))
if (as.character(utils::packageVersion("isCORAL")) != "6.1.4") {
  stop("Load isCORAL 6.1.4 before reproducing these figures.", call. = FALSE)
}
suppressPackageStartupMessages(library(patchwork))
source(file.path(WORKFLOW_DIR, "R", "metaclone_dynamics_plots.R"), encoding = "UTF-8")
options(
  iscoral.base_family = CONFIG$font_family,
  iscoral.dynamic.metaclone_colors = CONFIG$metaclone_colors
)
missing_reductions <- setdiff(
  c(CONFIG$timeline_reductions, CONFIG$dag_reduction), names(sample@reductions)
)
if (length(missing_reductions) > 0L) {
  stop(sprintf("Missing plotting reductions: %s", paste(missing_reductions, collapse = ", ")), call. = FALSE)
}
if (length(CONFIG$timeline_reductions) != 2L) {
  stop("This minimal figure script expects two timeline_reductions.", call. = FALSE)
}

FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
cache <- sample@misc$isCORAL_dynamics
sister_seeds <- cache$reproduction$resolved_sister_seeds
dynamic_plots <- PlotMetacloneDynamics(
  sample,
  view = c("potential", "potential_distribution", "seed_identity", "transition",
           "model_diagnostics", "dag", "terminal", "terminal_fate", "sister_fate"),
  reductions = CONFIG$timeline_reductions,
  seed_metaclones = sister_seeds, sister_value = "forward_mass",
  dag_reduction = CONFIG$dag_reduction
)

# 1. Cell-level and metaclone-level pseudotime on the original stored UMAP.
p_timeline_original <- dynamic_plots[[paste0("potential_", .iscoral_file_token(CONFIG$timeline_reductions[[1L]]))]]
print(p_timeline_original)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "01_timeline_original_umap.pdf"), p_timeline_original,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 2. The same pseudotime values in the configured pCORAL visualization space.
p_timeline_coral <- dynamic_plots[[paste0("potential_", .iscoral_file_token(CONFIG$timeline_reductions[[2L]]))]]
print(p_timeline_coral)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "02_timeline_pcoral.pdf"), p_timeline_coral,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 3. Transition matrices before and after forward-flux sparsification.
p_transition <- dynamic_plots$transition
print(p_transition)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "03_transition_matrices.pdf"), p_transition,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 4. Denoised transition graph anchored by the configured metaclone centroids.
p_dag <- dynamic_plots$dag
print(p_dag)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "04_dynamic_graph.pdf"), p_dag,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 5. Seed-core developmental-potential distributions.
p_potential_distribution <- dynamic_plots$potential_distribution
print(p_potential_distribution)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "05_seed_core_potential.pdf"), p_potential_distribution,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 6. Observed cell-type composition of the P_norm-aligned seed cores.
p_identity <- dynamic_plots$seed_identity
print(p_identity)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "06_seed_core_identity.pdf"), p_identity,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 7. Retained, direction-discarded and sparsification-removed flux per metaclone.
p_diagnostics <- dynamic_plots$model_diagnostics
print(p_diagnostics)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "07_flux_and_sink_diagnostics.pdf"), p_diagnostics,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 8. Terminal metaclones in both visualization spaces.
p_terminal <- dynamic_plots$terminal
print(p_terminal)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "08_terminal_states.pdf"), p_terminal,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 9. Terminal absorption probabilities, mapped to cells by their metaclone labels.
p_terminal_fate <- dynamic_plots[[paste0("terminal_fate_", .iscoral_file_token(CONFIG$timeline_reductions[[1L]]))]]
print(p_terminal_fate)
plot_size <- grDevices::dev.size("in")
ggplot2::ggsave(file.path(FIGURE_DIR, "09_terminal_fate_original_umap.pdf"), p_terminal_fate,
                device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")

# 10. One-step denoised outward flux. This is distinct from terminal absorption.
if (length(sister_seeds) > 0L) {
  p_sister_fate <- dynamic_plots[[paste0("sister_fate_", .iscoral_file_token(CONFIG$timeline_reductions[[1L]]))]]
  print(p_sister_fate)
  plot_size <- grDevices::dev.size("in")
  ggplot2::ggsave(file.path(FIGURE_DIR, "10_direct_sister_fate_original_umap.pdf"), p_sister_fate,
                  device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")
}

# 11. Reachable root-to-terminal paths, or maximum-flow networks if configured.
if (!is.null(cache$flow)) {
  .iscoral_validate_source_model(cache$flow, cache$model$model_id, "flow inference")
  reachable_results <- cache$flow$results[vapply(cache$flow$results, function(x) isTRUE(x$reachable), logical(1))]
  if (length(reachable_results) > 0L) {
    centroids <- .iscoral_centroid_data(sample, CONFIG$dag_reduction, CONFIG$metaclone_col)
    metaclone_colors <- .iscoral_get_metaclone_colors(rownames(cache$model$P_dynamic))
    path_panels <- lapply(reachable_results, .iscoral_flow_network_plot,
                         model = cache$model, centroids = centroids, metaclone_colors = metaclone_colors)
    p_paths <- patchwork::wrap_plots(path_panels, ncol = min(2L, length(path_panels)))
    print(p_paths)
    plot_size <- grDevices::dev.size("in")
    ggplot2::ggsave(file.path(FIGURE_DIR, "11_root_to_terminal_paths.pdf"), p_paths,
                    device = pdf_device, width = plot_size[1], height = plot_size[2], units = "in")
  }
}
message(sprintf("Core figures exported to %s", FIGURE_DIR))
