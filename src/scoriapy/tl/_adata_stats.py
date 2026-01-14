"""
Statistical utilities for comparing clusters, phases, and gene expression.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
import scikit_posthocs as sp
from IPython.display import display
from scipy.stats import kruskal, mannwhitneyu


def compare_cluster_phase(
    cross_table: pd.DataFrame,
    cluster_id: str,
    level: int = 0,
    group_level: str = "age",
):
    """
    Apply a Kruskal-Wallis test to each phase for a given cluster.

    Parameters
    ----------
    cross_table
        Cross-tabulated table (e.g. from :func:`pandas.crosstab`) whose index
        includes the cluster identifier.
    cluster_id
        Specific cluster label to analyze.
    level
        Level of the index in ``cross_table`` corresponding to ``cluster_id``.
    group_level
        Index level used to define groups for the test (e.g. ``\"age\"``).

    Returns
    -------
    dict
        Mapping from phase (e.g. ``\"G1\"``, ``\"G2M\"``, ``\"S\"``) to
        Kruskal-Wallis test results.
    """
    cluster_data = cross_table.xs(cluster_id, level=level)
    results = {
        phase: kruskal(
            *[data for _, data in cluster_data.groupby(level=group_level)[phase]]
        )
        for phase in ["G1", "G2M", "S"]
    }
    return results


def compare_gene_expression(
    adata: sc.AnnData,
    clustering_key: str,
    cluster_ids,
    gene: str,
    groupby_col: str,
):
    """
    Compare expression of a single gene across multiple groups with Kruskal-Wallis.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    clustering_key
        Column in ``adata.obs`` used to subset to specific clusters.
    cluster_ids
        One or more cluster identifiers to include.
    gene
        Gene name to test.
    groupby_col
        Column in ``adata.obs`` defining the groups to compare.

    Returns
    -------
    pandas.DataFrame or None
        Dunn's test results as a DataFrame if the global test is significant,
        otherwise ``None``.
    """
    adata_sub = adata[adata.obs[clustering_key].isin(cluster_ids)]

    dff = sc.get.obs_df(adata_sub, keys=[groupby_col, gene], use_raw=True).reset_index(
        drop=True
    )

    groups = dff.groupby(groupby_col)[gene]
    summary_stats = groups.agg(["mean", "median", "std"]).reset_index()
    display(summary_stats)
    groups_list = groups.apply(list)

    h_statistic, p_value = kruskal(*groups_list)
    print(f"Kruskal-Wallis H-Statistic for {gene}: {h_statistic}")
    print(f"Kruskal-Wallis P-Value for {gene}: {p_value}")

    fig, axes = plt.subplots(1, 2, figsize=(20, 6))

    sns.histplot(dff[gene], kde=False, ax=axes[0])
    axes[0].set_title(f"Expression of {gene}")
    axes[0].set_xlabel(gene)
    axes[0].set_ylabel("Density")

    sns.violinplot(x=groupby_col, y=gene, data=dff, inner=None, ax=axes[1])
    sns.stripplot(x=groupby_col, y=gene, data=dff, color="black", alpha=0.5, ax=axes[1])
    axes[1].set_title(f"Expression of {gene} by {groupby_col}")
    axes[1].set_xlabel(groupby_col)
    axes[1].set_ylabel("Expression")

    plt.tight_layout()
    plt.show()

    if p_value < 0.05:
        dunn_results = sp.posthoc_dunn(
            dff, val_col=gene, group_col=groupby_col, p_adjust="bonferroni"
        )
        display(dunn_results)
        return dunn_results

    print("No significant differences found.")
    return None


def compare_gene_expression_two(
    adata: sc.AnnData,
    clustering_key: str,
    cluster_ids,
    gene: str,
    groupby_col: str,
):
    """
    Compare expression of a single gene between exactly two groups using Mann-Whitney U.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    clustering_key
        Column in ``adata.obs`` used to subset to specific clusters.
    cluster_ids
        One or more cluster identifiers to include.
    gene
        Gene name to test.
    groupby_col
        Column in ``adata.obs`` defining the two groups being compared.

    Returns
    -------
    dict
        Dictionary with keys ``\"gene\"``, ``\"statistic\"``, and ``\"p_value\"``.
    """
    adata_sub = adata[adata.obs[clustering_key].isin(cluster_ids)]

    dff = sc.get.obs_df(adata_sub, keys=[groupby_col, gene], use_raw=True).reset_index(
        drop=True
    )

    groups = dff.groupby(groupby_col)[gene]
    display(groups.agg(["mean", "median", "std"]))

    unique_groups = dff[groupby_col].unique()
    if len(unique_groups) != 2:
        raise ValueError(
            f"Expected exactly 2 groups for comparison, but found {len(unique_groups)}."
        )

    group1, group2 = unique_groups
    group1_data = dff[dff[groupby_col] == group1][gene]
    group2_data = dff[dff[groupby_col] == group2][gene]

    stat, p_value = mannwhitneyu(group1_data, group2_data)
    print(f"Mann-Whitney U Statistic for {gene}: {stat}")
    print(f"Mann-Whitney U P-Value for {gene}: {p_value}")

    plt.figure(figsize=(10, 6))
    sns.violinplot(x=groupby_col, y=gene, data=dff)
    sns.stripplot(x=groupby_col, y=gene, data=dff, color="black", alpha=0.5)
    plt.title(f"Expression of {gene} by {groupby_col}")
    plt.xlabel(groupby_col)
    plt.ylabel("Expression")
    plt.show()

    return {"gene": gene, "statistic": stat, "p_value": p_value}
