"""
Embedding and UMAP plotting utilities.
"""

from __future__ import annotations

import math
from itertools import product
from typing import Iterable, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def plot_umap_by_categories_v1(
    adata,
    obs_columns: Sequence[str],
    basis: str = "X_umap",
    color: str = "leiden",
    ncol: int = 3,
    figsize: tuple[int, int] = (15, 5),
    sort_order: Sequence[str] | None = None,
    cluster_to_cell_mapping: Mapping[str, str] | None = None,
    legend_location: str = "on_data",
    point_size: float | None = None,
):
    """
    Plot UMAPs for each combination of categories in ``obs_columns``.

    Parameters
    ----------
    adata
        AnnData object containing an embedding in ``.obsm[basis]``.
    obs_columns
        Sequence of observation columns whose unique combinations define
        subplots (e.g. ``[\"age\", \"sample\"]``).
    basis
        Key in ``adata.obsm`` containing the embedding coordinates.
    color
        Observation column used to color cells (e.g. ``\"leiden\"``).
    ncol
        Number of subplots per row.
    figsize
        Base figure size; actual height scales with the number of rows.
    sort_order
        Optional ordering for the first component of the combined group name.
    cluster_to_cell_mapping
        Optional mapping from cluster IDs to cell type labels for the legend.
    legend_location
        Where to place legends (``\"on_data\"``, ``\"right\"``, or ``\"none\"``).
    point_size
        Optional marker size. If ``None``, a size proportional to the total
        cell count is used.
    """
    color_key = f"{color}_colors"
    if color_key not in adata.uns:
        raise ValueError(
            f"Color mapping not found in adata.uns['{color_key}']. "
            "Check if colors were set."
        )

    original_categories = adata.obs[color].astype("category").cat.categories
    original_colors = adata.uns[color_key]

    if len(original_colors) < len(original_categories):
        raise ValueError(
            f"Original colors for '{color}' are missing or incomplete "
            f"in adata.uns['{color_key}']."
        )

    if point_size is None:
        point_size = 120000 / adata.n_obs

    unique_combinations = adata.obs[obs_columns].drop_duplicates()
    unique_combinations["combined"] = unique_combinations.apply(
        lambda row: "_".join(map(str, row)), axis=1
    )

    if sort_order is not None:
        unique_combinations = unique_combinations.sort_values(
            by="combined",
            key=lambda x: [
                (
                    sort_order.index(i.split("_")[0])
                    if i.split("_")[0] in sort_order
                    else float("inf")
                )
                for i in x
            ],
        )

    sorted_combinations = unique_combinations.drop(columns=["combined"])

    num_plots = len(sorted_combinations)
    nrow = (num_plots + ncol - 1) // ncol

    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(figsize[0], nrow * figsize[1] // ncol),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    axes = axes.flat

    for ax in axes[num_plots:]:
        ax.set_visible(False)

    for idx, row in enumerate(sorted_combinations.itertuples(index=False)):
        subset_query = " & ".join(
            [f"{col} == '{getattr(row, col)}'" for col in obs_columns]
        )
        adata_sub = adata[adata.obs.eval(subset_query), :].copy()

        title = " | ".join([f"{col}: {getattr(row, col)}" for col in obs_columns])

        adata_sub.obs[color] = adata_sub.obs[color].astype("category")
        adata_sub.obs[color] = adata_sub.obs[color].cat.set_categories(
            original_categories
        )
        adata_sub.uns[color_key] = original_colors

        sc.pl.embedding(
            adata_sub,
            basis=basis,
            color=color,
            title=title,
            return_fig=False,
            show=False,
            ax=axes[idx],
            legend_loc="none",
            size=point_size,
        )

        if legend_location == "on_data":
            umap_coords = adata_sub.obsm[basis]
            clusters = adata_sub.obs[color].astype(str)
            cluster_medians = {
                c: np.median(umap_coords[clusters == c], axis=0)
                for c in clusters.unique()
            }
            for cluster, (x, y) in cluster_medians.items():
                axes[idx].text(
                    x,
                    y,
                    str(cluster),
                    fontsize=10,
                    weight="bold",
                    ha="center",
                    va="center",
                    color="black",
                )

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=original_colors[i],
            markersize=10,
            label=(
                f"Cluster {cat} - {cluster_to_cell_mapping.get(cat, 'Unknown')}"
                if cluster_to_cell_mapping
                else f"Cluster {cat}"
            ),
        )
        for i, cat in enumerate(original_categories)
    ]

    if legend_location != "right":
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            fontsize=8,
            ncol=min(ncol, len(handles)),
        )

    plt.tight_layout()
    plt.show()


def calculate_median_position(adata_subset, basis: str, overlay_attribute: str):
    """
    Calculate median positions of clusters in an embedding.

    Parameters
    ----------
    adata_subset
        AnnData object with an embedding in ``.obsm[basis]``.
    basis
        Key of the embedding in ``adata_subset.obsm``.
    overlay_attribute
        Observation column defining the groups whose medians are computed.

    Returns
    -------
    dict
        Mapping from label to a 2D median coordinate.
    """
    unique_labels = adata_subset.obs[overlay_attribute].unique()
    median_positions = {}
    for label in unique_labels:
        label_subset = adata_subset[adata_subset.obs[overlay_attribute] == label]
        median_positions[label] = np.median(label_subset.obsm[basis], axis=0)
    return median_positions


def add_text_labels(ax, median_positions, text_props_kw: dict | None = None):
    """
    Add text labels at median positions on a matplotlib Axes.

    Parameters
    ----------
    ax
        Target matplotlib Axes.
    median_positions
        Mapping from label to ``(x, y)`` coordinate.
    text_props_kw
        Optional dict of keyword arguments forwarded to ``ax.text``.
    """
    if text_props_kw is None:
        text_props_kw = {}
    for label, position in median_positions.items():
        ax.text(
            position[0], position[1], label, ha="center", va="center", **text_props_kw
        )


def _plot_one_embedding_with_label(
    adata,
    basis: str,
    color_attribute: str,
    overlay_attribute: str,
    ax,
    text_props_kw: dict | None = None,
):
    """
    Plot a single embedding and overlay text labels at median positions.

    Parameters
    ----------
    adata
        AnnData object containing the embedding.
    basis
        Key in ``adata.obsm`` with embedding coordinates.
    color_attribute
        Observation column used to color cells.
    overlay_attribute
        Observation column used to define label groups.
    ax
        Target matplotlib Axes.
    text_props_kw
        Optional dict of keyword arguments forwarded to ``ax.text``.
    """
    sc.pl.embedding(adata, basis=basis, color=color_attribute, ax=ax, show=False)
    median_positions = calculate_median_position(adata, basis, overlay_attribute)
    add_text_labels(ax, median_positions, text_props_kw=text_props_kw)


def plot_embeddings_with_labels(
    adata,
    basis: str,
    color_attribute: str,
    overlay_attribute: str,
    group_col: str | None = None,
    subplot_width: float = 5,
    text_props_kw: dict | None = None,
):
    """
    Plot embeddings with median-position labels, optionally split by a grouping column.

    Parameters
    ----------
    adata
        AnnData object with an embedding in ``.obsm[basis]``.
    basis
        Key in ``adata.obsm`` for the embedding (e.g. ``\"X_umap\"``).
    color_attribute
        Observation column used to color cells.
    overlay_attribute
        Observation column used to define label groups.
    group_col
        Optional column in ``adata.obs`` used to split the data into
        separate subplots.
    subplot_width
        Base width of each subplot (in inches).
    text_props_kw
        Optional dict of keyword arguments forwarded to ``ax.text``.
    """
    if group_col and group_col in adata.obs:
        groups = adata.obs[group_col].unique()
        fig, axes = plt.subplots(
            1,
            len(groups),
            figsize=(np.ceil(subplot_width * (len(groups) + 0.5)), subplot_width),
            squeeze=False,
        )
        axes = axes.flatten()
        for ax, group in zip(axes, groups):
            subset = adata[adata.obs[group_col] == group]
            _plot_one_embedding_with_label(
                subset, basis, color_attribute, overlay_attribute, ax, text_props_kw
            )
            ax.set_title(f"{group_col} = {group}")
    else:
        fig, ax = plt.subplots(figsize=(subplot_width + 1, subplot_width))
        _plot_one_embedding_with_label(
            adata, basis, color_attribute, overlay_attribute, ax, text_props_kw
        )
        ax.set_title("All data")

    plt.tight_layout()
    plt.show()


def plot_umap_by_group(adata, group_column: str, **kwargs):
    """
    Plot UMAPs for each unique value in ``group_column``.

    Parameters
    ----------
    adata
        AnnData object with a computed UMAP embedding.
    group_column
        Column in ``adata.obs`` whose unique values define subplots.
    **kwargs
        Additional keyword arguments forwarded to :func:`scanpy.pl.umap`.
    """
    if group_column not in adata.obs.columns:
        raise ValueError(f"'{group_column}' not found in adata.obs.")

    unique_groups = adata.obs[group_column].unique()
    num_groups = len(unique_groups)

    fig, axs = plt.subplots(
        1, num_groups, figsize=(4 * num_groups, 4), sharex=True, sharey=True
    )
    if num_groups == 1:
        axs = [axs]

    for i, group in enumerate(unique_groups):
        sc.pl.umap(
            adata[adata.obs[group_column] == group],
            color=group_column,
            ax=axs[i],
            show=False,
            **kwargs,
        )
        axs[i].set_title(f"{group_column} = {group}")

    fig.tight_layout(pad=2.0)
    plt.show()


def plot_pca_for_samples(adata, sample_list):
    """
    Plot PCA embeddings for a list of sample IDs.

    Parameters
    ----------
    adata
        AnnData object containing a ``\"sample\"`` column in ``obs``.
    sample_list
        Iterable of sample IDs to visualize.
    """
    available_samples_set = set(adata.obs["sample"].cat.categories)

    for sample_id in sample_list:
        if sample_id not in available_samples_set:
            print(f"Sample {sample_id} is not present in the provided adata.")
            continue

        sample_subset = adata[adata.obs["sample"] == sample_id]
        if sample_subset.n_obs > 0:
            sc.pl.pca(sample_subset, color="sample", title=sample_id)
        else:
            print(f"No data available for sample: {sample_id}")


def plot_umap_grid_by_sample(
    adata,
    col_per_row: int = 2,
    figsize_per_plot: tuple[int, int] = (4, 4),
):
    """
    Plot a grid of UMAPs, one per sample in ``adata.obs['sample_name']``.

    Parameters
    ----------
    adata
        AnnData object with a UMAP embedding.
    col_per_row
        Number of UMAP subplots per row.
    figsize_per_plot
        Size of each individual UMAP subplot.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the UMAP grid.
    """
    sample_ids = adata.obs["sample_name"].cat.categories.tolist()
    num_samples = len(sample_ids)

    rows = math.ceil(num_samples / col_per_row)
    cols = min(col_per_row, num_samples)

    fig_width = cols * figsize_per_plot[0]
    fig_height = rows * figsize_per_plot[1]

    fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height))

    if rows == 1 and cols == 1:
        axes = [[axes]]
    elif rows == 1:
        axes = [axes]
    elif cols == 1:
        axes = [[ax] for ax in axes]

    for i, sid in enumerate(sample_ids):
        row, col = divmod(i, cols)
        sc.pl.umap(
            adata[adata.obs["sample_name"] == sid],
            color="leiden",
            title=sid,
            ax=axes[row][col],
            legend_loc="on data",
            show=False,
        )

    fig.tight_layout()
    plt.close(fig)
    return fig


def _validate_umap_params(min_dists: Sequence[float], spreads: Sequence[float]):
    """
    Validate UMAP parameter grids and construct index/value combinations.

    Parameters
    ----------
    min_dists
        Sequence of ``min_dist`` values.
    spreads
        Sequence of ``spread`` values.

    Returns
    -------
    tuple
        ``(min_dists, spreads, param_combinations)`` where
        ``param_combinations`` is a list of ``((i, min_dist), (j, spread))``
        index/value pairs.
    """
    if not all(isinstance(x, (int, float)) for x in list(min_dists) + list(spreads)):
        raise ValueError(
            "All entries in min_dists and spreads must be integers or floats."
        )
    if not min_dists or not spreads:
        raise ValueError("min_dists and spreads must be non-empty sequences.")

    param_combinations = list(product(enumerate(min_dists), enumerate(spreads)))
    return min_dists, spreads, param_combinations


def _compute_and_plot_umap(adata, min_dist: float, spread: float, ax):
    """
    Compute a UMAP embedding with given parameters and plot it on ``ax``.

    Notes
    -----
    This helper mutates ``adata`` in-place by adding the UMAP result.
    """
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)

    # Prefer a categorical obs column for coloring if available
    if "leiden" in adata.obs:
        color_key = "leiden"
    elif "group" in adata.obs:
        color_key = "group"
    else:
        color_key = None

    sc.pl.umap(
        adata,
        color=color_key,
        title=f"min_dist = {min_dist}, spread = {spread}",
        s=40,
        ax=ax,
        show=False,
    )


def _setup_umap_param_grid(min_dists: Sequence[float], spreads: Sequence[float]):
    """
    Create a matplotlib figure and axes grid for a UMAP parameter sweep.

    The grid size is ``len(min_dists) x len(spreads)``.
    """
    fig, axes = plt.subplots(
        len(min_dists),
        len(spreads),
        figsize=(len(spreads) * 3 + 2, len(min_dists) * 3),
    )
    # Ensure a 2D array of axes for consistent indexing
    axes = np.atleast_2d(axes)
    return fig, axes


def plot_umap_param_grid(
    adata,
    min_dists: Sequence[float] = (0.1, 1, 2),
    spreads: Sequence[float] = (0.5, 1, 5),
):
    """
    Plot a grid of UMAPs for all combinations of ``min_dist`` and ``spread``.

    Parameters
    ----------
    adata
        AnnData object with precomputed neighbors/PCA, ready for UMAP.
    min_dists
        Sequence of ``min_dist`` values to evaluate.
    spreads
        Sequence of ``spread`` values to evaluate.
    """
    min_dists, spreads, param_combinations = _validate_umap_params(min_dists, spreads)
    adata_temp = adata.copy()
    fig, axes = _setup_umap_param_grid(min_dists, spreads)

    for (i, min_dist), (j, spread) in param_combinations:
        ax = axes[i, j]
        _compute_and_plot_umap(adata_temp, min_dist, spread, ax)

    plt.tight_layout()
    plt.show()
    plt.close(fig)
    del adata_temp
