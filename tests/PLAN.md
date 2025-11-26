## Testing plan for `scoriapy`

This document summarizes the main ideas for testing the `scoriapy`
package with `pytest`.

### 1. Structure

- Use a `tests/` package alongside `src/`.
- Organize tests roughly by module:
  - `tests/test_utils_extract.py`
  - `tests/test_utils_markers_logging.py`
  - `tests/test_pp_hvg.py`
  - `tests/test_pl_embedding.py`
  - `tests/test_pl_heatmap.py`
  - `tests/test_tl_de_aggregate_stats.py`

### 2. Fixtures

Defined in `tests/conftest.py`:

- `adata_base` – synthetic AnnData with:
  - 50 cells × 10 genes
  - `obs` columns: `group`, `age`, `phase`, `sample`, `sample_name`
  - `.raw` set to a copy of the object
- `adata_with_umap` – `adata_base` with:
  - simple preprocessing (normalize, log1p, PCA, neighbors)
  - UMAP embedding in `.obsm["X_umap"]`
  - clustering labels in `obs["leiden"]` and corresponding colors
- `adata_for_heatmap` – `adata_base` with:
  - `.raw` set
  - a `var["gene_cat"]` annotation for heatmap tests

These fixtures keep tests fast and deterministic without requiring any
external data.

### 3. What to test

- **`scoriapy.utils.markers`**
  - `subset_by_markers`:
    - Different strategies (`any`, `count`, `threshold`) on small, hand‑crafted
      arrays.
    - Return shape and that a *new* AnnData is returned (original untouched).
    - Optional `add_obs_key` behavior.

- **`scoriapy.utils._extract`**
  - `create_cell_mask`:
    - Single value, list, and callable conditions; validate masks manually.
  - `subset_adata_by_mask`:
    - Shapes match expectations and `use_raw=True` raises if `.raw` is missing.
  - `get_mat_df` / `get_full_df`:
    - Column names and shapes; z‑score properties when `scale=True`.
  - `make_annotation_df` / `make_mat_dfs`:
    - Simple `adata_dict` and expected DataFrame shapes/content.
  - `scale_df`:
    - Output has ~zero mean and unit variance per column.

- **`scoriapy.pp.run_hvg`**
  - Uses `adata_base`:
    - Returned AnnData has `n_vars == n_top`.
    - Original AnnData still has all genes.

- **`scoriapy.pl._embedding`**
  - Plotting functions (`plot_umap_by_categories`, `plot_embeddings_with_labels`,
    `plot_umap_by_group`, `plot_pca_for_samples`,
    `plot_umap_grid_by_sample`, `plot_umap_param_grid`):
    - Smoke tests – run without exceptions and, where applicable, return a
      `Figure`.
    - Do not mutate input `adata` in unexpected ways (shape and `obs` preserved).
  - `calculate_median_position`:
    - On a tiny dataset, medians match expected coordinates.

- **`scoriapy.pl._heatmap`**
  - `subset_adata`:
    - Filters by group and genes, sets `var["gene_cat"]`.
  - `sort_cells_by_metadata`:
    - Returns a *new* AnnData with sorted cells, leaving the original unchanged.
  - `get_heatmap_dfs`:
    - Output keys and shapes as expected.
  - `plot_distribution` / `plot_cluster_heatmap`:
    - Smoke tests – run without exceptions.

- **`scoriapy.tl._de`, `_aggregate`, `_stats`**
  - `get_pairwise_deg`, `run_rank_genes_groups`, `get_de_genes`,
    `get_deg_df`, `detect_gene_flips`:
    - Use simple two‑group scenarios; check keys, basic columns, and
      that calls complete without errors.
  - `aggregate_single_cell_data` / `create_group_dicts`:
    - Expected shapes and group definitions.
  - `summarize_category_distribution`:
    - Counts align with manual crosstabs; normalized percentages sum to ~100.
  - `compare_cluster_phase`, `compare_gene_expression`,
    `compare_gene_expression_two`:
    - Smoke tests – return values of the right type and do not crash on
      small synthetic inputs.

### 4. In‑place vs new‑object behavior

For functions where non‑mutating behaviour is important, tests should:

- Capture a copy of `adata.obs` (and/or `adata.shape`) before the call.
- Assert that:
  - `id(adata)` is unchanged.
  - `obs` content is unchanged where we expect pure functions.
- For functions that are explicitly documented to *return* a sorted or
  transformed AnnData, verify that the returned object has the expected
  ordering/values and the original remains unchanged.

