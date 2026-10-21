# ==============================================================================
# 🚀 Metaclone API execution pipeline
# ==============================================================================
if (!exists("sample", inherits = FALSE) || !inherits(sample, "Seurat")) {
  stop("Run 5_1.R first, or provide a Seurat object named `sample`.")
}
if (!exists("export_dir", inherits = FALSE)) {
  stop("Run 5_1.R first, or define `export_dir` before running this script.")
}

suppressWarnings(loadfonts(device = "pdf", quiet = TRUE))

# Prerequisites: define the lineage palette and required metadata columns.
highlight_colors <- c(
  "HSC"   = "#FF52B7",
  "GMP"   = "#00D68B",
  "Mk"    = "#FFAA00",
  "Other" = "#707070"
)
BARCODE_COL <- if (exists("BARCODE_COL", inherits = FALSE)) {
  BARCODE_COL
} else {
  Sys.getenv("ISCORAL_BARCODE_COL", unset = "larry_sapphire_merged")
}
CLUSTER_COL <- Sys.getenv("ISCORAL_METACLONE_COL", unset = "metaclone")

message("\n📊 Running the core pipeline and computing visualization matrices...")

# 1. Run topology-based community detection and attach metaclone labels.
sample <- FindMetaclone(sample, resolution = 0.9, reduction = "CORAL_clone_embed", cluster_col_name = CLUSTER_COL)

# 2. Build multimodal validation maps with full-matrix community rendering.
plot_landscape <- PlotMetacloneLandscape(sample, reduction = "CORAL_clone_embed", cluster_col = CLUSTER_COL)

plot_landscape <- PlotMetacloneLandscape(sample, reduction = "CORAL_cell_embed", cluster_col = CLUSTER_COL)


# ⚠️ [Core adaptation 1]: Pass the required barcode_col and retain each returned plot object.
suite_cell <- PlotValidationSuite(sample, reduction = "CORAL_cell_embed", space_title = "CORAL Cell Space", cluster_col = CLUSTER_COL, barcode_col = BARCODE_COL)
suite_umap <- PlotValidationSuite(sample, reduction = "umap.raw", space_title = "UMAP Space", cluster_col = CLUSTER_COL, barcode_col = BARCODE_COL)


# 5. Export all plots to the configured embedding output directory.
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

message("💾 Writing high-resolution PDF files...")

# Save the community landscape.
ggsave(file.path(export_dir, "01_metaclone_landscape.pdf"), plot = plot_landscape, device = cairo_pdf)

# ⚠️ [Core adaptation 3]: Save suite_cell and suite_umap directly; do not access $facets.
ggsave(file.path(export_dir, "02_validation_facets_cell.pdf"), plot = suite_cell, device = cairo_pdf)
ggsave(file.path(export_dir, "03_validation_facets_umap.pdf"), plot = suite_umap, device = cairo_pdf)

