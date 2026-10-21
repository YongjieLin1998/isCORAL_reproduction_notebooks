# ==============================================================================
# 🌌 Step 6: Visualizing 'Probability Cloud Collapse' (Diptych Version)
# ==============================================================================
library(ggplot2)
library(patchwork)
library(stringr)
library(ggrepel)
library(extrafont) 

# This continuation script expects 5_1.R to have created `sample` and `export_dir`.
if (!exists("sample", inherits = FALSE) || !inherits(sample, "Seurat")) {
  stop("Run 5_1.R first, or provide a Seurat object named `sample`.")
}
if (!exists("export_dir", inherits = FALSE)) {
  stop("Run 5_1.R first, or define `export_dir` before running this script.")
}

# ------------------------------------------------------------------------------
# 0. Global Arial font and background-color configuration
# ------------------------------------------------------------------------------
suppressWarnings(loadfonts(device = "pdf", quiet = TRUE))

update_geom_defaults("text", list(family = "Arial"))
update_geom_defaults("label", list(family = "Arial"))
update_geom_defaults("text_repel", list(family = "Arial")) 

# High-saturation, high-contrast hematopoietic lineage palette
highlight_colors <- c(
  "HSC"   = "#FF52B7",  # Fluorescent magenta
  "GMP"   = "#00D68B",  # Bright emerald green
  "Mk"    = "#FFAA00",  # Saturated golden yellow
  "Other" = "#707070"   # Light gray
)

# Shared background color: a slightly darker cool gray
bg_gray <- "#E0E4E7"

# ------------------------------------------------------------------------------
# 1. Core diptych renderer for a minimal comparison
# ------------------------------------------------------------------------------
plot_probability_collapse_diptych <- function(seurat_obj, seed_cells, clone_col = "larry_sapphire_merged", point_size = 0.5) {
  
  cache <- seurat_obj@misc$isCORAL_embedding
  if (is.null(cache)) stop("🚨 Error: isCORAL_embedding cache not found. Run CORALembedding() first.")
  
  if (!clone_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("🚨 Error: ground-truth clone/lineage column '%s' was not found.", clone_col))
  }
  
  P_norm <- cache$P_norm
  valid_cells <- cache$valid_cells
  rownames(P_norm) <- valid_cells
  
  mat_umap <- Embeddings(seurat_obj, "umap.raw")[valid_cells, 1:2]
  all_clones <- as.character(seurat_obj@meta.data[valid_cells, clone_col])
  names(all_clones) <- valid_cells
  
  theme_collapse <- theme_void(base_family = "Arial") + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "plain"), 
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40", face = "plain"),
      legend.position = "none", 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5) 
    )
  
  plot_list <- list()
  
  for (i in seq_along(seed_cells)) {
    seed_cell <- seed_cells[i]
    if (!seed_cell %in% valid_cells) next
    
    seed_clone <- all_clones[seed_cell]
    if (is.na(seed_clone)) seed_clone <- "Unknown/NA"
    
    cell_probs <- as.numeric(P_norm[seed_cell, ])
    cell_probs[is.na(cell_probs)] <- 0
    if (max(cell_probs) > 0) {
      cell_probs <- cell_probs / max(cell_probs) 
    }
    
    raw_clusters <- as.character(seurat_obj@meta.data[valid_cells, "cluster_integrated"])
    mapped_clusters <- ifelse(raw_clusters %in% c("HSC", "GMP", "Mk"), raw_clusters, "Other")
    
    df <- data.frame(
      U1 = mat_umap[, 1], U2 = mat_umap[, 2],
      Prob = cell_probs,
      Cluster = mapped_clusters,
      IsSeed = (valid_cells == seed_cell),
      IsSibling = (all_clones == seed_clone) & (valid_cells != seed_cell) & (seed_clone != "Unknown/NA"),
      stringsAsFactors = FALSE
    )
    
    df_seed <- df[df$IsSeed, ]
    df_sibling <- df[df$IsSibling, ]
    df_cloud <- df[order(df$Prob), ] 
    
    # --- Left: isCORAL inference (probability cloud) ---
    p_cloud <- ggplot() +
      geom_point(data = subset(df_cloud, Prob == 0), aes(x = U1, y = U2), color = bg_gray, size = point_size) +
      geom_point(data = subset(df_cloud, Prob > 0), aes(x = U1, y = U2, color = Cluster, alpha = Prob), size = point_size * 1.6) +
      geom_point(data = df_seed, aes(x = U1, y = U2, fill = Cluster), shape = 21, size = 5.0, color = "black", stroke = 1.2) +
      scale_color_manual(values = highlight_colors) +
      scale_fill_manual(values = highlight_colors) +
      # Preserve the linear probability-to-alpha mapping.
      scale_alpha_continuous(range = c(0, 1), guide = "none") + 
      theme_collapse +
      labs(title = "isCORAL inference", subtitle = sprintf("seed cell type: %s", df_seed$Cluster))
    
    # --- Right: LARRY ground truth ---
    clean_clone_id <- sub("^[Ss]apphire_?", "", seed_clone)
    
    p_gt <- ggplot() +
      geom_point(data = subset(df, !IsSibling & !IsSeed), aes(x = U1, y = U2), color = bg_gray, size = point_size, alpha = 0.8) +
      geom_point(data = df_sibling, aes(x = U1, y = U2, fill = Cluster), shape = 21, size = 2.8, color = "black", stroke = 0.7) +
      geom_point(data = df_seed, aes(x = U1, y = U2, fill = Cluster), shape = 21, size = 5.0, color = "black", stroke = 1.2) +
      scale_fill_manual(values = highlight_colors) +
      theme_collapse +
      labs(title = "LARRY groundtruth", subtitle = sprintf("Clone ID: %s", clean_clone_id))
    
    plot_list[[i]] <- (p_cloud | p_gt)
  }
  
  final_composite <- wrap_plots(plot_list, ncol = 1)
  return(final_composite)
}

# ------------------------------------------------------------------------------
# 2. Specify representative seed cells for reproducibility
# ------------------------------------------------------------------------------
clone_col_name <- if (exists("BARCODE_COL", inherits = FALSE)) {
  BARCODE_COL
} else {
  Sys.getenv("ISCORAL_BARCODE_COL", unset = "larry_sapphire_merged")
}

message("\n📊 [Locking representative seed cells]")

demo_cells <- c(
  Sys.getenv("ISCORAL_HSC_SEED_CELL", unset = "<HSC_SEED_CELL>"),
  Sys.getenv("ISCORAL_GMP_SEED_CELL", unset = "<GMP_SEED_CELL>"),
  Sys.getenv("ISCORAL_MK_SEED_CELL", unset = "<MK_SEED_CELL>")
)
if (any(grepl("^<[^>]+>$", demo_cells))) {
  stop(
    "Set ISCORAL_HSC_SEED_CELL, ISCORAL_GMP_SEED_CELL, and ISCORAL_MK_SEED_CELL, ",
    "or replace the seed-cell placeholders before running this template."
  )
}

cache_global <- sample@misc$isCORAL_embedding
valid_cells_global <- cache_global$valid_cells

missing_cells <- setdiff(demo_cells, valid_cells_global)
if (length(missing_cells) > 0) {
  stop(sprintf("🚨 Fatal error: the following seed-cell barcodes are absent from the embedding: %s", paste(missing_cells, collapse=", ")))
}

message(sprintf("   ✅ HSC seed locked: %s", demo_cells[1]))
message(sprintf("   ✅ GMP seed locked: %s", demo_cells[2]))
message(sprintf("   ✅ Mk seed locked: %s", demo_cells[3]))

diptych_plot <- plot_probability_collapse_diptych(
  seurat_obj = sample, 
  seed_cells = demo_cells, 
  clone_col = clone_col_name
)


# ==============================================================================
# 🌌 Step 7: Master Anchor Map (Combined Clone Embedding with Fate Blending)
# ==============================================================================
message("\n🌟 [Step 7] Computing global fate-blended colors and preparing the combined manifold anchor map...")

hex_to_rgb <- function(hex) { as.numeric(col2rgb(hex)) / 255 }
rgb_HSC   <- hex_to_rgb(highlight_colors["HSC"])
rgb_GMP   <- hex_to_rgb(highlight_colors["GMP"])
rgb_Mk    <- hex_to_rgb(highlight_colors["Mk"])
rgb_Other <- hex_to_rgb(highlight_colors["Other"]) 

mat_clone_global <- Embeddings(sample, "CORAL_clone_embed")[valid_cells_global, 1:2]
raw_clusters_global <- as.character(sample@meta.data[valid_cells_global, "cluster_integrated"])
mapped_clusters_global <- ifelse(raw_clusters_global %in% c("HSC", "GMP", "Mk"), raw_clusters_global, "Other")

P_norm_global <- cache_global$P_norm
is_HSC   <- mapped_clusters_global == "HSC"
is_GMP   <- mapped_clusters_global == "GMP"
is_Mk    <- mapped_clusters_global == "Mk"
is_Other <- mapped_clusters_global == "Other"

prob_HSC_all   <- as.numeric(Matrix::rowSums(P_norm_global[, is_HSC, drop = FALSE]))
prob_GMP_all   <- as.numeric(Matrix::rowSums(P_norm_global[, is_GMP, drop = FALSE]))
prob_Mk_all    <- as.numeric(Matrix::rowSums(P_norm_global[, is_Mk, drop = FALSE]))
prob_Other_all <- as.numeric(Matrix::rowSums(P_norm_global[, is_Other, drop = FALSE]))

total_core_prob <- prob_HSC_all + prob_GMP_all + prob_Mk_all + prob_Other_all
total_core_prob[total_core_prob == 0] <- 1e-10 

w_HSC   <- prob_HSC_all   / total_core_prob
w_GMP   <- prob_GMP_all   / total_core_prob
w_Mk    <- prob_Mk_all    / total_core_prob
w_Other <- prob_Other_all / total_core_prob

df_global <- data.frame(
  P1 = mat_clone_global[, 1],
  P2 = mat_clone_global[, 2],
  Cluster = mapped_clusters_global,
  stringsAsFactors = FALSE
)
rownames(df_global) <- valid_cells_global

df_global$Blended_Color <- rgb(
  red   = pmin(pmax(w_HSC * rgb_HSC[1] + w_GMP * rgb_GMP[1] + w_Mk * rgb_Mk[1] + w_Other * rgb_Other[1], 0), 1),
  green = pmin(pmax(w_HSC * rgb_HSC[2] + w_GMP * rgb_GMP[2] + w_Mk * rgb_Mk[2] + w_Other * rgb_Other[2], 0), 1),
  blue  = pmin(pmax(w_HSC * rgb_HSC[3] + w_GMP * rgb_GMP[3] + w_Mk * rgb_Mk[3] + w_Other * rgb_Other[3], 0), 1)
)

df_seeds_combined <- df_global[demo_cells, ]
df_seeds_combined$Seed_ID <- paste0("Seed ", 1:nrow(df_seeds_combined))

p_master_anchor <- ggplot() +
  # 1. Background cells: fate-blended colors with thin hollow circles (shape = 1)
  geom_point(data = df_global, aes(x = P1, y = P2), 
             shape = 1, color = df_global$Blended_Color, size = 1.2, stroke = 0.5, alpha = 0.85) +
  
  # 2. Seed cells: thick outlined circles with white fill to mask the background
  geom_point(data = df_seeds_combined, aes(x = P1, y = P2, color = Cluster), 
             shape = 21, fill = "white", size = 6.5, stroke = 2.5) +
  
  geom_text_repel(data = df_seeds_combined, aes(x = P1, y = P2, label = Seed_ID),
                  size = 5, fontface = "plain", box.padding = 0.8, 
                  point.padding = 0.5, segment.color = "grey50") +
  
  scale_color_manual(values = highlight_colors) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"), 
    legend.position = "right", 
    legend.title = element_text(size = 12, face = "plain"), 
    legend.text = element_text(size = 11, face = "plain"),
    panel.border = element_blank(),
    plot.margin = ggplot2::margin(15, 15, 15, 15)
  ) +
  labs(
    title = "CORAL clone embedding",
    subtitle = "Fate blended pseudoclone space", 
    color = "seed cell lineage" 
  )


# ==============================================================================
# 🌟 Step 8 (enhanced): In-situ overlay of ground truth on the cell embedding
# ==============================================================================
message("\n🌟 [Step 8 enhanced] Plotting inferred probability clouds and LARRY ground truth in CORAL_cell_embed space...")

# Extract cell-level embedding coordinates.
mat_cell <- Embeddings(sample, "CORAL_cell_embed")[valid_cells_global, 1:2]
seed_lineages <- c("HSC", "GMP", "Mk")

# Extract clone labels globally to identify sibling cells.
all_clones_global <- as.character(sample@meta.data[valid_cells_global, clone_col_name])
names(all_clones_global) <- valid_cells_global

df_facet_combined <- data.frame()
rownames(cache_global$P_norm) <- valid_cells_global

for (i in seq_along(demo_cells)) {
  seed_cell <- demo_cells[i]
  lineage <- seed_lineages[i]
  
  # 1. Extract and normalize inferred probabilities for the current seed cell.
  probs <- as.numeric(cache_global$P_norm[seed_cell, ])
  probs[is.na(probs)] <- 0
  if (max(probs) > 0) {
    probs <- probs / max(probs) 
  }
  
  # 2. Identify true LARRY clone siblings, excluding the seed and unknown clones.
  seed_clone <- all_clones_global[seed_cell]
  if (is.na(seed_clone)) seed_clone <- "Unknown/NA"
  
  is_sibling <- (all_clones_global == seed_clone) & 
    (valid_cells_global != seed_cell) & 
    (seed_clone != "Unknown/NA")
  
  # 3. Assemble the data frame for the current facet.
  temp_df <- data.frame(
    C1 = mat_cell[, 1],
    C2 = mat_cell[, 2],
    Prob = probs,
    Seed_Lineage = factor(lineage, levels = c("HSC", "GMP", "Mk")), 
    Target_Cluster = mapped_clusters_global,                        
    IsSeed = (valid_cells_global == seed_cell),
    IsSibling = is_sibling,
    stringsAsFactors = FALSE
  )
  
  # Draw high-probability cells last so the zero-probability background does not obscure them.
  temp_df <- temp_df[order(temp_df$Prob), ]
  df_facet_combined <- rbind(df_facet_combined, temp_df)
}

# Extract highlight subsets for precise top-layer placement in ggplot2.
df_siblings_only <- subset(df_facet_combined, IsSibling == TRUE)
df_seeds_only    <- subset(df_facet_combined, IsSeed == TRUE)

# 4. Build the multilayer composite.
p_gradient_overlay <- ggplot(df_facet_combined, aes(x = C1, y = C2)) +
  # Layer 1: gray background cells with zero inferred probability
  geom_point(data = subset(df_facet_combined, Prob == 0), 
             color = bg_gray, size = 0.5) +
  
  # Layer 2: probability cloud; color denotes lineage tendency and alpha denotes probability
  geom_point(data = subset(df_facet_combined, Prob > 0), 
             aes(color = Target_Cluster, alpha = Prob), size = 0.8) +
  
  # Layer 3: LARRY ground-truth siblings as black-bordered circles
  geom_point(data = df_siblings_only, 
             aes(fill = Target_Cluster), shape = 21, color = "black", size = 2.5, stroke = 0.6) +
  
  # Layer 4: seed cells as large top-layer markers
  geom_point(data = df_seeds_only, 
             aes(fill = Target_Cluster), shape = 21, color = "black", size = 5.5, stroke = 1.2) +
  
  # Keep palette and mappings aligned.
  scale_color_manual(values = highlight_colors, name = "Cell lineage") +
  scale_fill_manual(values = highlight_colors, guide = "none") +
  scale_alpha_continuous(range = c(0.1, 1), name = "isCORAL probability") +
  
  # Arrange the three lineages in parallel facets.
  facet_wrap(~ Seed_Lineage, ncol = 3) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40", face = "plain", margin = margin(b = 15)),
    strip.text = element_text(size = 14, face = "plain", margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "plain", vjust = 1),
    legend.text = element_text(size = 12, face = "plain"),
    legend.key.width = unit(2.5, "line"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(1.5, "lines"),
    plot.margin = ggplot2::margin(15, 15, 15, 15)
  ) +
  labs(
    title = "Pseudoclone probability distributions on cell embedding",
    subtitle = "In-situ Overlay: Colored clouds = isCORAL inference | Black-bordered circles = LARRY groundtruth"
  )


# ==============================================================================
# 🎨 Step 9: Final rendering and PDF export
# ==============================================================================
message("\n🖥️ [Step 9] Rendering all visualization layers and exporting PDF files...")

print(diptych_plot)
print(p_master_anchor)
print(p_gradient_overlay) # Preview the enhanced overlay.

if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

filename_diptych <- file.path(export_dir, "Probability_Cloud_Collapse_Diptych.pdf")
ggsave(filename = filename_diptych, plot = diptych_plot, device = cairo_pdf)
message(sprintf("   ✅ Diptych comparison saved to: %s", filename_diptych))

filename_master <- file.path(export_dir, "Master_Anchor_Map_Combined.pdf")
ggsave(filename = filename_master, plot = p_master_anchor, device = cairo_pdf)
message(sprintf("   ✅ Master anchor map saved to: %s", filename_master))

# Save the enhanced overlay as a PDF.
filename_gradient <- file.path(export_dir, "Pseudoclone_Probability_Gradients_CellEmbed_Overlay.pdf")
ggsave(filename = filename_gradient, plot = p_gradient_overlay, device = cairo_pdf)
message(sprintf("   ✅ Probability-cloud and ground-truth overlay saved to: %s", filename_gradient))

message("\n==============================================================================")
message("🎉 All manifold visualizations have been generated.")
message("==============================================================================")
