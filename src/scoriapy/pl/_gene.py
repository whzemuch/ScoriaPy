"""
Plotting utilities focused on gene-level summaries.
"""

from __future__ import annotations

import math
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def plot_group_embeddings(
    adata,
    group_col: str,
    color: str,
    basis: str = "X_pca_umap",
    **kwargs,
):
    """
    Create embedding plots for each group in ``group_col``.

    Parameters
    ----------
    adata
        AnnData object with an embedding in ``.obsm[basis]``.
    group_col
        Column in ``adata.obs`` used to define groups.
    color
        Observation column used to color cells.
    basis
        Key in ``adata.obsm`` containing the embedding coordinates.
    **kwargs
        Additional keyword arguments forwarded to :func:`scanpy.pl.embedding`.
    """
    uniq_groups = adata.obs[group_col].unique()
    num_groups = uniq_groups.size

    fig, axes = plt.subplots(1, num_groups, figsize=(5 * num_groups, 4), squeeze=False)
    axes = axes.flatten()

    for ax, value in zip(axes, uniq_groups):
        adata_subset = adata[adata.obs[group_col] == value]
        sc.pl.embedding(
            adata_subset,
            basis=basis,
            color=color,
            ax=ax,
            show=False,
            **kwargs,
        )
        ax.set_title(f"{group_col} = {value}")

    fig.tight_layout(pad=2.0)
    plt.show()


def plot_violin_genes(
    adata1,
    adata2,
    common_genes: List[str],
    group_col_adata1: str,
    group_col_adata2: str,
    subplot_width: float = 2,
):
    """
    Plot violin plots for genes across specified groups from two AnnData objects.

    Parameters
    ----------
    adata1
        First AnnData object.
    adata2
        Second AnnData object.
    common_genes
        List of genes to visualize.
    group_col_adata1
        Grouping column in ``adata1.obs``.
    group_col_adata2
        Grouping column in ``adata2.obs``.
    subplot_width
        Width of each subplot in inches.
    """
    groups_adata1 = adata1.obs[group_col_adata1].unique()
    groups_adata2 = adata2.obs[group_col_adata2].unique()
    comb_groups = np.concatenate([groups_adata1, groups_adata2])

    num_groups = len(comb_groups)
    num_genes = len(common_genes)

    fig, axes = plt.subplots(
        num_genes,
        num_groups,
        figsize=(math.ceil(subplot_width * num_groups), subplot_width * num_genes),
        squeeze=False,
    )

    for i, gene in enumerate(common_genes):
        for j, group in enumerate(comb_groups):
            ax = axes[i, j]
            subset = adata1 if group in groups_adata1 else adata2
            subset = subset[
                subset.obs[
                    group_col_adata1 if group in groups_adata1 else group_col_adata2
                ]
                == group
            ]

            sc.pl.violin(subset, [gene], groupby="pca-leiden", ax=ax, show=False)
            ax.set_title(f"{gene} - {group}" if j == 0 else group)

    plt.tight_layout()
    plt.show()
