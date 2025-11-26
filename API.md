# ScoriaPy API Reference

This document summarizes the public API of the `scoriapy` package.
It follows the Scanpy-style namespace:

- `scoriapy.pp` – preprocessing / IO utilities
- `scoriapy.tl` – analysis tools / statistics
- `scoriapy.pl` – plotting utilities
- `scoriapy.utils` – generic helpers

For implementation details, see the corresponding modules in `src/scoriapy`.

---

## Top-level package

```python
import scoriapy

scoriapy.pp      # preprocessing / IO
scoriapy.tl      # tools / analysis
scoriapy.pl      # plotting
scoriapy.utils   # generic helpers
```

---

## `scoriapy.pp` – preprocessing / IO

Module: `src/scoriapy/pp`

### `run_hvg(adata, n_top=3000, flavor="seurat", log=None)`

Compute and subset highly variable genes using `scanpy.pp.highly_variable_genes`.
Returns a copy of `adata` restricted to HVGs.

### `read_10x_h5_concat(samples, sample_info)`

Wrapper around `scanpy.read_10x_h5` to read and concatenate multiple 10x HDF5
datasets.

- `samples`: iterable of `(sample_name, Path)` pointing to filtered H5 files.
- `sample_info`: `pandas.DataFrame` indexed by `sample_name`, must contain at
  least a `"treatment"` column.
- Returns an `AnnData` with `sample_name` and `treatment` in `.obs`.

### `read_10x_mtx_concat(main_folder_path, batch_key="batch")`

Read 10x MTX folders inside a parent directory and concatenate them.

- Each direct subdirectory is treated as a separate batch.
- Returns an `AnnData` with a batch label column in `.obs` named `batch_key`.

---

## `scoriapy.tl` – analysis tools

Module: `src/scoriapy/tl`

### Differential expression (`_de.py`)

- `get_pairwise_deg(adata, comparisons, groupby_col="group", method="wilcoxon")`  
  Convenience wrapper around `scanpy.tl.rank_genes_groups` to compute pairwise
  DE between groups; returns a dict `{ "A_vs_B": DataFrame, ... }`.

- `run_rank_genes_groups(adata, groupby_attribute, rank_key=None)`  
  Run `scanpy.tl.rank_genes_groups` for a grouping attribute and display the
  standard Scanpy DE plots.

- `get_de_genes(adata, group_col="age", key="rank_genes_groups", query=None)`  
  Extract sets of DE genes per group from an existing `rank_genes_groups`
  result in `adata.uns[key]`. Returns `{group: set(genes)}`.

- `get_deg_df(adata, cluster_key, groupby_attribute, rank_key, query)`  
  Run DE, filter by an expression/statistics query, and return:
  - `test_df`: expression matrix for the top genes
  - `filter_df`: statistics table of the selected genes.

- `detect_gene_flips(one_down, one_up, two_down, two_up, return_df=True)`  
  Identify genes that switch from down→up or up→down between two conditions.

### Aggregation and categorical summaries (`_aggregate.py`)

- `combine_group_values(group_values)`  
  Join group values into a single name, e.g. `"sample_celltype"`.

- `check_group_columns_exist(adata, group_cols)`  
  Validate that `group_cols` exist in `adata.obs`.

- `aggregate_group_expression(adata, group_mask, group_name)`  
  Aggregate expression for a subset of cells into a single pseudo-bulk column.

- `aggregate_single_cell_data(adata, group_cols)`  
  Group by one or more `obs` columns and compute pseudo-bulk expression
  matrices. Returns a DataFrame whose columns are groups.

- `create_group_dicts(adata, group_cols)`  
  Return a list of dictionaries describing groups, useful for downstream
  tooling.

- `summarize_category_distribution(adata, category_col, group_col, selected_groups=None, normalize=False, sort_order=None)`  
  Crosstab counts of a categorical column across group(s), with optional
  prefix filtering, sorting, and normalization. Returns a styled DataFrame
  (for notebooks) or raw counts.

### Statistical tests (`_stats.py`)

- `compare_cluster_phase(cross_table, cluster_id, level=0, group_level="age")`  
  Apply Kruskal–Wallis tests over phases (G1, G2M, S) within a given cluster
  index in a cross-tab.

- `compare_gene_expression(adata, clustering_key, cluster_ids, gene, groupby_col)`  
  Kruskal–Wallis test of a single gene across multiple groups, with summary
  stats and distribution plots.

- `compare_gene_expression_two(adata, clustering_key, cluster_ids, gene, groupby_col)`  
  Mann–Whitney U test of a single gene between exactly two groups, plus
  violin/strip plots.

---

## `scoriapy.pl` – plotting utilities

Module: `src/scoriapy/pl`

### Embedding / UMAP plotting (`_embedding.py`)

- `plot_umap_by_categories(adata, obs_columns, basis="X_umap", color="leiden", ...)`  
  Plot UMAPs for each unique combination of categories in `obs_columns`,
  preserving Scanpy color mappings and optionally adding on-data labels and
  legends.

- `calculate_median_position(adata_subset, basis, overlay_attribute)`  
  Compute median embedding coordinates per label for later annotation.

- `add_text_labels(ax, median_positions, text_props_kw=None)`  
  Draw labels at median positions on a matplotlib Axes.

- `plot_embeddings_with_labels(adata, basis, color_attribute, overlay_attribute, group_col=None, ...)`  
  Plot embeddings with median-position labels, optionally split by a grouping
  column.

- `plot_umap_by_group(adata, group_column, **kwargs)`  
  Multiple UMAPs, one subplot per unique value in `group_column`.

- `plot_pca_for_samples(adata, sample_list)`  
  PCA plots for a list of sample IDs, checking presence in `adata.obs["sample"]`.

- `plot_umap_grid_by_sample(adata, col_per_row=2, figsize_per_plot=(4, 4))`  
  Grid of UMAPs, one per `sample_name` in `adata.obs["sample_name"]`.

- `plot_umap_param_grid(adata, min_dists=(0.1, 1, 2), spreads=(0.5, 1, 5))`  
  Explore UMAP layout by plotting all combinations of `min_dist` and `spread`.

### Heatmaps (`_heatmap.py`)

- `subset_adata(adata, group_col, keys, genes_map)`  
  Subset cells by `group_col` and genes by a mapping, and annotate gene
  categories in `.var["gene_cat"]`.

- `sort_cells_by_metadata(adata, primary_sort_col, primary_categories, second_sort_col, second_categories)`  
  Reorder categorical `obs` columns and sort cells accordingly.

- `get_heatmap_dfs(adata, obs_cols, var_cols, scale=False)`  
  Prepare:
  - `df_mat`: expression matrix (genes x cells), optionally z-scored
  - `col_df`: selected columns from `adata.obs`
  - `row_df`: selected columns from `adata.var`

- `plot_distribution(df, n=1000, x_intercept=-1.5)`  
  KDE plot for matrix values, with an optional threshold line.

- `plot_cluster_heatmap(df_mat, col_df, row_df, figsize=(15, 15), **kwargs)`  
  Clustered heatmap with annotations using PyComplexHeatmap, splitting columns
  by `age` and annotating genes by `gene_cat`.

- `create_heatmap_annotation(col_df)`  
  Build a `HeatmapAnnotation` from a column DataFrame, auto-selecting color
  maps for known fields (e.g. `phase`).

### Gene-level plots (`_gene.py`)

- `plot_group_embeddings(adata, group_col, color, basis="X_pca_umap", **kwargs)`  
  Embedding plots, one per group in `group_col`.

- `plot_violin_genes(adata1, adata2, common_genes, group_col_adata1, group_col_adata2, subplot_width=2)`  
  Grid of violin plots for a list of genes across groups from two datasets.

---

## `scoriapy.utils` – generic helpers

Module: `src/scoriapy/utils`

### Logging (`logging.py`)

- `setup_logger(log_file=None, level=logging.INFO)`  
  Configure a `logging.Logger` named `"scoriapy"` with stream/file handlers,
  avoiding duplicate handlers.

### Marker-based filtering (`markers.py`)

- `subset_by_markers(adata, markers, strategy="any", threshold=0.0, require_n_markers=1, min_count=1, use_raw=False, add_obs_key=None)`  
  Flexible marker-based cell filtering with different strategies:
  `"any"`, `"all"`, `"count"`, `"threshold"`. Optionally stores the mask in
  `adata.obs[add_obs_key]`.

### Extraction helpers (`_extract.py`)

- `create_cell_mask(adata, conditions)`  
  Build a boolean mask based on flexible conditions over `adata.obs`.

- `subset_adata_by_mask(adata, cell_mask, gene_list, use_raw=True)`  
  Subset cells and genes, optionally using `adata.raw`.

- `extract_cell_cluster(adata, col_obs, keys)`  
  Convenience wrapper to subset by discrete values in an `obs` column.

- `get_mat_df(adata, keys=None, use_raw=True)`  
  Return a `DataFrame` of expression for selected keys (rows = cells).

- `get_full_df(adata, keys=None, scale=True)`  
  Expression from `adata.raw`, optionally z-scored.

- `make_annotation_df(adata_dict)`  
  Build a combined annotation DataFrame from an `{"age": ..., "rapa": ...}`
  dict of AnnData objects.

- `make_mat_dfs(adata_dict, common_genes)`  
  Extract scaled expression matrices for a list of genes from multiple AnnData
  objects, returning `{name: DataFrame (genes x cells)}`.

- `scale_df(df)`  
  Column-wise z-score scaling of a DataFrame.

### Color utilities (`_colors.py`)

- `map_colors(reference_dict, target_dict, mapping_dict)`  
  Map colors from a reference cluster map to a target cluster map using a
  mapping dictionary; defaults unmapped entries to gray.

- `display_colors(color_list)`  
  Show color blocks inline in a notebook.

- `display_colors_with_keys(color_map)`  
  Show labeled color blocks inline in a notebook.

- `confirm_and_set_colors(adata, color_attribute, color_map)`  
  Interactively confirm and apply a color map to `adata.uns[f"{color_attribute}_colors"]`.

---

## Notes

- Many functions are designed to be notebook-friendly (using styled DataFrames
  and inline HTML/plots). For non-interactive pipelines, you can call the
  underlying computation helpers (e.g. `aggregate_single_cell_data`,
  `get_de_genes`, `get_mat_df`) and ignore the plotting/styling.
- For detailed provenance of each helper relative to the original notebook
  scripts in `src/inbox`, see `src/inbox/history.md`.

