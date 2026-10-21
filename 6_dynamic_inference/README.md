# Minimal isCORAL dynamic inference: HSC and embryo

This directory contains ordinary R scripts for the shared dynamic-inference
workflow. It starts from a processed Seurat object with pCORAL embeddings and
metaclone assignments, and produces developmental-potential summaries,
direction-filtered and sparsified transitions, terminal fate probabilities,
root-to-terminal paths, expected edge traversals, and core figures.

There is no archive to unpack and no dependency on the old workflow directories.
The isCORAL package and the input dataset are supplied separately.

## Files

```text
6_dynamic_inference/
  HSC.R                              # HSC configuration and execution order
  embryo.R                           # Human-embryo configuration and execution order
  01_inference.R                     # Load input, infer dynamics and export results
  02_plots.R                         # Explicit figure-by-figure preview and export
  R/
    metaclone_dynamics_internal.R     # Numerical helpers and palette selection
    metaclone_dynamics_api.R          # Original dynamics and path-inference functions
    metaclone_dynamics_plots.R        # Core plotting functions
  README.md
```

## Input and dependencies

Use the Seurat object from the completed embedding/metaclone analysis. It must
contain:

- `misc$isCORAL_embedding$P_norm` and its matching `valid_cells`;
- the configured metaclone and cell-type metadata columns;
- the reductions listed in the selected entry script;
- a suitable assay with the configured expression slot/layer for CytoTRACE2,
  or existing finite `CytoTRACE2_Score` values in `[0, 1]` when scores are reused.

The original CytoTRACE2 settings use `counts`. Check the object's default assay
and, for Seurat v5, prepare a single usable counts layer before this stage. This
workflow does not join layers, normalize expression or reconstruct upstream
metaclones automatically.

Required packages are **isCORAL 6.1.4**, Seurat, Matrix, igraph, ggplot2 and
patchwork. CytoTRACE2 is required when `run_cytotrace2 = TRUE`. Loading an isCORAL
source directory also requires devtools. Optional ggrastr rasterizes dense
point layers. No script installs software or fonts.

The historical workflows explicitly require isCORAL 6.1.4; this requirement is
retained. Point `ISCORAL_PACKAGE_DIR` to that source package, or make that exact
version available in the R library. The beta package is not substituted.

## Run in RStudio

Set the working directory to the reproduction project root, which contains
`1_data`, `5_embeddings` and this `6_dynamic_inference` directory. Replace the
sample placeholder and supply a processed Seurat RDS:

```r
Sys.setenv(
  ISCORAL_SAMPLE_ID = "<SAMPLE_ID>",
  ISCORAL_DYNAMIC_INPUT = file.path("5_embeddings", "<SAMPLE_ID>", "sample_with_embeddings.rds")
)

# Choose ONE entry script for the current dataset.
source("6_dynamic_inference/HSC.R", encoding = "UTF-8")
# source("6_dynamic_inference/embryo.R", encoding = "UTF-8")
```

An alternative is to continue in the session used for the embedding analysis:

```r
# `sample` must already be the processed Seurat object for this dataset.
Sys.setenv(ISCORAL_SAMPLE_ID = "<SAMPLE_ID>")
Sys.unsetenv("ISCORAL_DYNAMIC_INPUT")
source("6_dynamic_inference/HSC.R", encoding = "UTF-8")
```

If the embedding scripts have not saved a checkpoint, save the prepared object
once with `saveRDS(sample, "sample_with_embeddings.rds")` and use that file as
`ISCORAL_DYNAMIC_INPUT`. A supplied RDS always takes precedence over an existing
in-session `sample`.

For manual figure adjustment, open the entry script, run its configuration block
and the `source(..., "01_inference.R")` line, then open `02_plots.R`. Run the plot
setup, followed by each explicit `print()` and save block. Resize the Plot pane
between preview and save. Each PDF uses `grDevices::dev.size("in")` from the
active preview device. Sourcing the whole entry script also exports the figures
without pausing; figure exports are intended for an active RStudio Plot device.

To reuse previously computed CytoTRACE2 scores, set
`CONFIG$run_cytotrace2 = FALSE` in the entry script before running the inference
step. The reuse is recorded in the output cache.

## Dataset settings

| Setting | HSC | Embryo |
| --- | --- | --- |
| CytoTRACE2 species | `mouse` | `human` |
| CytoTRACE2 input and seed | `counts`, `42` | `counts`, `42` |
| Cell-type metadata | `cluster_integrated` | `celltype` |
| Metaclone metadata | `metaclone` | `metaclone` |
| Direction threshold | `0.015` | `0.015` |
| Positive forward-flux quantile | `0.70` | `0.70` |
| Minimum flux / numerical tolerance | `0` / `1e-12` | `0` / `1e-12` |
| Flow roots | MC5, MC9 | Top two transient states by mean potential |
| Sister-fate seeds | MC5, MC9, MC15 | Top three transient states by mean potential |
| Default path method | `most_probable` | `most_probable` |
| Original UMAP | `umap.raw` | `umap` |
| Second timeline view | `CORAL_cell_embed` | `umap.pcoral.cell` |
| DAG layout | `CORAL_clone_embed` | `umap.pcoral.clone` |

The HSC annotations and root/seed IDs belong to the original metaclone numbering.
For a different clustering, edit these fields and use `annotations = NULL` until
biological labels have been assigned. They are not inferred from the sample ID.

The embryo root-selection rule comes from the original `6.2b` helper. If a
historical run used manually chosen roots, put those IDs in `flow_roots` to
reproduce that selection. The resolved roots and seeds are saved in the cache;
the configuration source file is never rewritten.

Embryo pCORAL UMAPs use dimensions `1:30`, 30 neighbors, `min.dist = 0.3`, cosine
distance and seed 42. As in the original script, unavailable dimensions are
excluded and neighbors are capped at the number of cells minus one; the effective
values are recorded. `overwrite = TRUE` recomputes only the two configured
pCORAL UMAPs. Set it to `FALSE` to preserve stored coordinates. These projections
do not change `P_norm`, metaclone membership or the dynamic model. The original
Seurat UMAP is read as stored.

## Paths and outputs

All paths can be edited directly in the entry script. Optional environment
variables are:

| Variable | Default or meaning |
| --- | --- |
| `ISCORAL_PROJECT_ROOT` | Current working directory |
| `ISCORAL_DYNAMIC_DIR` | `<project_root>/6_dynamic_inference` |
| `ISCORAL_SAMPLE_ID` | Required sample identifier |
| `ISCORAL_DYNAMIC_INPUT` | Processed Seurat RDS; empty reuses `sample` |
| `ISCORAL_PACKAGE_DIR` | `<project_root>/packages/isCORAL_v6.1.4` |
| `ISCORAL_DYNAMIC_OUTPUT_DIR` | `<dynamic_dir>/results` |
| `ISCORAL_RUN_ID` | Current timestamp |

Results are written to
`results/<HSC_or_embryo>/<SAMPLE_ID>/<RUN_ID>/`. Use a fresh run ID for a new
analysis; reusing the same ID writes the same output filenames.

- `sample_with_dynamics.rds`: complete Seurat object, including model cache.
- `dynamic_cache.rds`: dynamic model, fate, paths, expected flow and provenance.
- `config_used.R`, `sessionInfo.txt`: requested parameters and software environment.
- `metaclone_state_summary.csv`, `seed_core_identity.csv`,
  `seed_core_potential_summary.csv`: seed-core membership and potential summaries.
- `block_flux.csv`, `forward_flux_*.csv`, `P_dynamic*.csv`: transition stages.
- `sister_fate_*.csv`, `terminal_fate.csv`: direct outward fate and terminal absorption.
- `path_summary.csv`, `jump_chain.csv`: root-to-terminal inference, when roots exist.
- `expected_edge_traversals.csv`: cumulative off-diagonal edge flow for each root.
- `figures/`: numbered core PDFs.

The cache also records resolved roots, actual UMAP settings, input identity,
default assay, source-file MD5 values and session information. Detailed
maximum-flow edge results remain available in `cache$flow$results` when
`flow_method = "max_flow"` is selected.

To redraw a completed run without repeating inference, restore these objects and
then source the plotting script:

```r
WORKFLOW_DIR <- "6_dynamic_inference"
OUTPUT_DIR <- file.path(WORKFLOW_DIR, "results", "HSC", "<SAMPLE_ID>", "<RUN_ID>")
sample <- readRDS(file.path(OUTPUT_DIR, "sample_with_dynamics.rds"))
CONFIG <- sample@misc$isCORAL_dynamics$reproduction$config
# Load isCORAL 6.1.4 in this R session before plotting.
source(file.path(WORKFLOW_DIR, "02_plots.R"), encoding = "UTF-8")
```

## Scope and model interpretation

The model averages `P_norm` over source seed cores and sums over target seed
cores. It retains an edge only when the source-minus-target mean potential is
strictly greater than `epsilon`, then applies a global quantile threshold to
positive admissible forward fluxes (quantile type 8). Residency contributes to
the final row normalization. Terminal fates come from absorption in this final
denoised `P_dynamic`.

The default path is the most probable path in the jump chain after removal of
self-transitions and row renormalization. `expected_edge_traversals.csv` instead
uses the absorbing-chain fundamental matrix and the original `P_dynamic`:
`N[root, i] * P_dynamic[i, j]`, with diagonal entries omitted. It is not a
maximum-flow score or a separately normalized edge probability. A root that is
already terminal has zero off-diagonal expected traversals. All-terminal models
retain their numerical outputs; automatic transient-root selection yields no paths.

The minimum workflow covers the numerical content of original steps 6.1, 6.2,
7.1 and the expected-flow calculation in 9.2, with selected 6.3/6.4/7.2 views.
The two numerical API bodies are retained from `dynamic_inference_RStudio` and
also match their embryo counterparts. Plot helpers derive from the shared HSC
and embryo code, with export handling separated, larger text and explanatory
subtitles removed.

This directory does not include upstream CORAL scans or metaclone clustering,
embryo spatial reconstruction, DEG/GO and TF analyses, route-gene comparisons,
terminal cell-type pies, or all exploratory figures. It starts from the exact
prepared cell set; it does not reconstruct a historical spatial subset.

Matching historical numerical results requires the same input cells, `P_norm`,
metaclone labels, CytoTRACE2 scores/settings, selected roots and software versions.
The provided source files have not been executed or rendered as part of this
reorganization; runtime or figure equivalence has not been established here.
