"""
Heatmap plotting utilities based on PyComplexHeatmap.
"""

from __future__ import annotations

from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import zscore

import PyComplexHeatmap as pch


def subset_adata(
    adata: sc.AnnData,
    group_col: str,
    keys: List[str],
    genes_map: Dict[str, str],
) -> sc.AnnData:
    """
    Subset AnnData by selected groups and genes, assigning gene categories.

    Parameters
    ----------
    adata
        AnnData object with a populated ``raw`` attribute.
    group_col
        Column in ``adata.obs`` used to select cells.
    keys
        Values in ``group_col`` that define the subset of cells.
    genes_map
        Mapping from gene name to category label. Only genes present in
        ``adata.raw.var_names`` are used.

    Returns
    -------
    AnnData
        Subsetted AnnData with selected genes and cells; gene categories are
        stored in ``.var['gene_cat']``.
    """
    genes = [gene for gene in genes_map if gene in adata.raw.var_names]
    genes_cat = [genes_map[gene] for gene in genes]

    adata_sub = adata.raw[:, genes].to_adata()
    adata_sub.var["gene_cat"] = genes_cat

    cell_mask = adata.obs[group_col].isin(keys)
    return adata_sub[cell_mask, :]


def sort_cells_by_metadata(
    adata: sc.AnnData,
    primary_sort_col: str,
    primary_categories: List[str],
    second_sort_col: str,
    second_categories: List[str],
) -> sc.AnnData:
    """
    Sort cells by two categorical obs columns.

    Parameters
    ----------
    adata
        AnnData object to sort.
    primary_sort_col
        Name of the primary categorical column (e.g. ``\"age\"``).
    primary_categories
        Desired ordering of categories in ``primary_sort_col``.
    second_sort_col
        Name of the secondary categorical column.
    second_categories
        Desired ordering of categories in ``second_sort_col``.

    Returns
    -------
    AnnData
        AnnData object reordered according to the sorted index.
    """
    adata_sorted = adata.copy()
    adata_sorted.obs[primary_sort_col] = adata_sorted.obs[
        primary_sort_col
    ].cat.reorder_categories(primary_categories, ordered=True)
    adata_sorted.obs[second_sort_col] = adata_sorted.obs[
        second_sort_col
    ].cat.reorder_categories(second_categories, ordered=True)

    sorted_index = adata_sorted.obs.sort_values(
        by=[primary_sort_col, second_sort_col]
    ).index
    return adata_sorted[sorted_index]


def get_heatmap_dfs(
    adata: sc.AnnData,
    obs_cols: List[str],
    var_cols: List[str],
    scale: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Prepare matrices and annotation DataFrames for heatmap plotting.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    obs_cols
        List of observation columns to include as column annotations.
    var_cols
        List of variable columns to include as row annotations.
    scale
        If ``True``, z-score rows of the expression matrix.

    Returns
    -------
    dict
        Dictionary with keys ``\"df_mat\"``, ``\"col_df\"``, and ``\"row_df\"``.
    """
    # Expression matrix: genes x cells, stored as float for downstream scaling
    df_mat = sc.get.obs_df(adata, list(adata.var_names)).T.astype(float)
    if scale:
        # Apply z-score row-wise while preserving DataFrame shape
        numeric = df_mat.select_dtypes(include=[np.number])
        scaled = numeric.apply(
            lambda x: zscore(x, nan_policy="omit"), axis=1, result_type="expand"
        )
        scaled.columns = numeric.columns
        df_mat.loc[:, numeric.columns] = scaled

    col_df = sc.get.obs_df(adata, obs_cols)
    row_df = adata.var.loc[:, var_cols]
    return {"df_mat": df_mat, "col_df": col_df, "row_df": row_df}


def plot_distribution(df: pd.DataFrame, n: int = 1000, x_intercept: float = -1.5):
    """
    Plot a KDE of values in ``df`` with an optional vertical cutoff line.

    Parameters
    ----------
    df
        DataFrame whose numeric values to visualize.
    n
        Maximum number of columns to sample for plotting.
    x_intercept
        Location of the vertical reference line.
    """
    # Handle both 1D (Series) and 2D (DataFrame) inputs
    if hasattr(df, "ndim") and df.ndim == 1:
        data = df.to_numpy().ravel()
    else:
        if df.shape[1] > n:
            df = df.sample(n=n, axis=1)
        data = df.to_numpy().ravel()

    fig, ax = plt.subplots(figsize=(6, 3))
    sns.kdeplot(data, ax=ax, bw_adjust=0.5)
    ax.axvline(x=x_intercept, color="r", linestyle="--")
    plt.show()


def plot_cluster_heatmap(
    df_mat: pd.DataFrame,
    col_df: pd.DataFrame,
    row_df: pd.DataFrame,
    figsize: tuple[int, int] = (15, 15),
    **kwargs,
):
    """
    Plot a clustered heatmap with annotations for genes and cells.

    Parameters
    ----------
    df_mat
        Expression matrix (rows = genes, columns = cells).
    col_df
        Column annotation DataFrame (index should align with columns of
        ``df_mat``).
    row_df
        Row annotation DataFrame (index should align with rows of
        ``df_mat`` and include a ``\"gene_cat\"`` column).
    figsize
        Size of the matplotlib figure.
    **kwargs
        Additional keyword arguments forwarded to
        :class:`PyComplexHeatmap.ClusterMapPlotter`.
    """
    row_ha = pch.HeatmapAnnotation(
        group=pch.anno_simple(row_df.gene_cat, legend=True, cmap="Accent"),
        axis=0,
        plot=False,
        label_side="bottom",
    )

    col_ha = pch.HeatmapAnnotation(
        Age=pch.anno_simple(col_df.age, legend=True, cmap="Dark2"),
        Phase=pch.anno_simple(col_df.phase, add_text=False, legend=True, cmap="Set1"),
        legend=True,
        legend_gap=5,
        hgap=0.5,
        plot=False,
    )

    plt.figure(figsize=figsize)
    pch.ClusterMapPlotter(
        data=df_mat,
        top_annotation=col_ha,
        right_annotation=row_ha,
        col_split=col_df.age,
        z_score=0,
        vmin=-1.2,
        vmax=2.5,
        show_rownames=True,
        show_colnames=False,
        row_cluster=False,
        col_cluster=True,
        row_dendrogram=False,
        col_dendrogram=True,
        cmap="RdYlBu_r",
        **kwargs,
    )


def create_heatmap_annotation(col_df: pd.DataFrame):
    """
    Create a :class:`PyComplexHeatmap.HeatmapAnnotation` from a column DataFrame.

    Parameters
    ----------
    col_df
        Column annotation DataFrame whose columns will be turned into
        annotation tracks.

    Returns
    -------
    PyComplexHeatmap.HeatmapAnnotation
        Configured annotation object.
    """
    annotations = {}
    for col in col_df.columns:
        display_name = col.capitalize()
        cmap = "Set1" if col.lower() == "phase" else "Dark2"
        annotations[display_name] = pch.anno_simple(col_df[col], legend=True, cmap=cmap)

    return pch.HeatmapAnnotation(
        **annotations, legend=True, legend_gap=5, hgap=0.5, plot=False
    )
