"""
Differential expression utilities built on top of scanpy.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from tqdm.auto import tqdm


def run_and_store_pairwise_rg(
    adata,
    groupby,
    comparisons,
    method="wilcoxon",
    n_genes=20,
    logger=None
):
    """
    Run multiple pairwise differential expression (DE) comparisons using
    Scanpy's rank_genes_groups, automatically generate plots, and store
    results tables inside ``adata.uns``.

    This is a convenience wrapper that iterates over a dictionary of
    comparisons and performs three steps for each:

    1. Compute DE with ``sc.tl.rank_genes_groups``.
    2. Save the standard Scanpy DE plot as PDF.
    3. Extract results into a pandas DataFrame and store it at:
       ``adata.uns[key]["df"]``.

    Logging messages are emitted to track progress and outputs.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must contain the grouping column in
        ``adata.obs[groupby]``.

    groupby : str
        Column name in ``adata.obs`` used for grouping cells
        (e.g., "age", "condition", "cluster", "leiden").

    comparisons : dict[str, tuple[str, str]]
        Dictionary mapping result keys → (group, reference).

        Example:
        {
            "OvsY": ("Old", "Young"),
            "MvsY": ("Middle", "Young")
        }

        Each key:
        - becomes the Scanpy result key_added
        - names the saved plot file: ``_{key}.pdf``
        - stores the DEG table in ``adata.uns[key]["df"]``

    method : str, default="wilcoxon"
        Statistical test for ``rank_genes_groups``.
        Common options:
        {"wilcoxon", "t-test", "logreg"}

    n_genes : int, default=20
        Number of top genes shown in the generated plot.

    logger : logging.Logger or None, default=None
        Logger for progress reporting. If None, uses
        ``logging.getLogger(__name__)``.

    Returns
    -------
    None
        Results are stored in-place inside ``adata``:

        - DEG statistics: ``adata.uns[key]`` (Scanpy native format)
        - DataFrame table: ``adata.uns[key]["df"]``
        - Plots saved to disk: ``rank_genes_groups_{key}.pdf``

    Side Effects
    ------------
    Writes PDF files to the current Scanpy figure directory
    (controlled by ``sc.settings.figdir``).

    Examples
    --------
    Basic usage:

    >>> import logging
    >>> logging.basicConfig(level=logging.INFO)

    >>> comparisons = {
    ...     "OvsY": ("O", "Y"),
    ...     "MvsY": ("M", "Y")
    ... }

    >>> run_and_store_pairwise_rg(
    ...     adata,
    ...     groupby="age",
    ...     comparisons=comparisons,
    ...     method="wilcoxon",
    ...     n_genes=20
    ... )

    Access DEG table:

    >>> df = adata.uns["OvsY"]["df"]
    >>> df.head()

    Notes
    -----
    This function is particularly useful for:
    - automated cell-type specific DE
    - running many contrasts reproducibly
    - keeping plots + tables synchronized
    - logging long Scanpy workflows
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Starting pairwise rank_genes_groups")
    logger.info("Groupby: %s | Method: %s", groupby, method)
    logger.info("Comparisons: %s", list(comparisons.keys()))

    for key, (group, ref) in comparisons.items():

        logger.info("Running %s (%s vs %s)", key, group, ref)

        # 1. Run DE
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            groups=[group],
            reference=ref,
            method=method,
            key_added=key
        )

        logger.info(
            "Finished rank_genes_groups for %s", key
        )

        # 2. Plot
        sc.pl.rank_genes_groups(
            adata,
            key=key,
            groups=[group],
            n_genes=n_genes,
            show=False,
            save=f"_{key}.pdf"
        )

        logger.info(
            "Saved rank_genes_groups plot for %s", key
        )

        # 3. Extract and store DataFrame
        df = sc.get.rank_genes_groups_df(
            adata, key=key, group=group
        )

        adata.uns[key]["df"] = df

        logger.info(
            "Stored DEG table for %s (%d genes)",
            key, df.shape[0]
        )

    logger.info("Completed all pairwise comparisons")






def get_pairwise_deg(
    adata: sc.AnnData,
    comparisons: Iterable[Tuple[str, str]],
    groupby_col: str = "group",
    method: str = "wilcoxon",
) -> Dict[str, pd.DataFrame]:
    """
    Perform pairwise differential expression analysis for specified comparisons.

    This is a convenience wrapper around :func:`scanpy.tl.rank_genes_groups`.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    comparisons
        Iterable of pairs ``(group1, group2)`` to compare.
    groupby_col
        Column in ``adata.obs`` used for grouping.
    method
        Differential expression method passed to Scanpy (e.g. ``\"wilcoxon\"``).

    Returns
    -------
    dict
        Mapping from comparison label (e.g. ``\"A_vs_B\"``) to a
        :class:`pandas.DataFrame` of DE statistics.
    """
    if "log1p" not in adata.uns or "base" not in adata.uns["log1p"]:
        print("Warning: log1p base not found in adata.uns. Setting default base e.")
        adata.uns["log1p"] = {"base": np.e}

    results_dict: Dict[str, pd.DataFrame] = {}
    pbar = tqdm(total=len(list(comparisons)), desc="Processing comparisons")

    for group1, group2 in comparisons:
        key = f"{group1}_vs_{group2}"
        pbar.set_description(f"Processing {key}")

        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby_col,
            groups=[group1],
            reference=group2,
            method=method,
            key_added=key,
            pts=True,
        )

        df = sc.get.rank_genes_groups_df(adata, group=group1, key=key)
        df["comparison"] = key
        results_dict[key] = df
        pbar.update(1)

    pbar.close()
    return results_dict


def run_rank_genes_groups(
    adata: sc.AnnData,
    groupby_attribute: str,
    rank_key: str | None = None,
) -> None:
    """
    Run :func:`scanpy.tl.rank_genes_groups` and plot the result.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    groupby_attribute
        Column in ``adata.obs`` used to define groups.
    rank_key
        Optional key under which to store results in ``adata.uns``. If
        ``None``, defaults to ``f\"{groupby_attribute}_rank\"``.
    """
    # Use Scanpy's default key name for easier downstream use
    if rank_key is None:
        rank_key = "rank_genes_groups"

    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby_attribute,
        method="wilcoxon",
        pts=True,
        key_added=rank_key,
    )
    sc.pl.rank_genes_groups(adata, sharey=False, key=rank_key)


def get_de_genes(
    adata: sc.AnnData,
    group_col: str = "age",
    key: str = "rank_genes_groups",
    query: str | None = None,
) -> Dict[str, set]:
    """
    Extract sets of differentially expressed genes for each group.

    Assumes that differential expression results are already present in
    ``adata.uns[key]`` as created by :func:`scanpy.tl.rank_genes_groups`.

    Parameters
    ----------
    adata
        AnnData object with DE results.
    group_col
        Column in ``adata.obs`` that defines the groups.
    key
        Key under which DE results are stored in ``adata.uns``.
    query
        Optional pandas query string to filter the DE table for each group.

    Returns
    -------
    dict
        Mapping from group label to a set of gene names.
    """
    if key not in adata.uns:
        raise ValueError(
            f"The key '{key}' does not exist in the provided AnnData object. "
            "Please run the differential expression analysis first."
        )

    if query is None:
        query = "(logfoldchanges > 1) | (logfoldchanges < -1) & (pvals_adj < 0.01)"

    gene_dict: Dict[str, set] = {}
    for group in adata.obs[group_col].unique():
        df_rank = sc.get.rank_genes_groups_df(adata, group=group, key=key)
        filtered_genes = df_rank.query(query)
        gene_dict[group] = set(filtered_genes["names"].values)

    return gene_dict


def get_deg_df(
    adata: sc.AnnData,
    cluster_key: str,
    groupby_attribute: str,
    rank_key: str,
    query: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute a DEG table for a given query and return expression for selected genes.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    cluster_key
        Column in ``adata.obs`` indicating clusters (currently not used in
        the implementation but kept for API compatibility).
    groupby_attribute
        Column in ``adata.obs`` used for grouping in DE analysis.
    rank_key
        Key under which to store DE results in ``adata.uns``.
    query
        Pandas query string used to filter the DE table.

    Returns
    -------
    test_df
        Expression matrix (genes x cells) for the selected genes.
    filter_df
        DataFrame with statistics for the selected genes.
    """
    from ..utils._extract import get_full_df

    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby_attribute,
        method="wilcoxon",
        pts=True,
        key_added=rank_key,
    )
    sc.pl.rank_genes_groups(adata, sharey=False, key=rank_key)

    # Choose the first available group from the grouping column
    # to make this helper work on generic synthetic data.
    groups = list(adata.obs[groupby_attribute].unique())
    if not groups:
        raise ValueError(
            f"No groups found in adata.obs['{groupby_attribute}'] for DEG computation."
        )
    group_for_df = groups[0]

    deg_df = sc.get.rank_genes_groups_df(adata, group=group_for_df, key=rank_key)
    filter_df = deg_df.query(query).sort_values("pvals_adj")
    deg_genes: List[str] = list(filter_df["names"].iloc[:150].values)

    full_df = get_full_df(adata, scale=True)
    test_df = full_df.loc[:, deg_genes].T
    return test_df, filter_df


def detect_gene_flips(
    one_down: List[str],
    one_up: List[str],
    two_down: List[str],
    two_up: List[str],
    return_df: bool = True,
):
    """
    Identify genes that flip direction between two conditions.

    Parameters
    ----------
    one_down, one_up, two_down, two_up
        Gene lists for the two conditions.
    return_df
        If ``True``, return a DataFrame, otherwise a dict.
    """
    down_2_up = set(one_down) & set(two_up)
    up_2_down = set(one_up) & set(two_down)

    results = {"down_to_up": down_2_up, "up_to_down": up_2_down}

    if return_df:
        df = pd.DataFrame(
            [
                (gene, condition)
                for condition, genes in results.items()
                for gene in genes
            ],
            columns=["Gene", "Condition"],
        )
        return df

    return results
