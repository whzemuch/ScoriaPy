"""
Aggregation utilities: pseudo-bulk expression and categorical summaries.
"""

from __future__ import annotations

from typing import List

import pandas as pd


def combine_group_values(group_values) -> str:
    """
    Combine group column values into a single string, separated by underscores.

    Parameters
    ----------
    group_values
        Iterable of values describing a group.

    Returns
    -------
    str
        Combined group name.
    """
    return "_".join(str(col) for col in group_values)


def check_group_columns_exist(adata, group_cols: List[str]) -> None:
    """
    Ensure that all requested grouping columns exist in ``adata.obs``.

    Parameters
    ----------
    adata
        AnnData object whose ``obs`` will be checked.
    group_cols
        List of column names that must be present in ``adata.obs``.

    Raises
    ------
    ValueError
        If any of the requested columns are missing.
    """
    missing_cols = [col for col in group_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Columns {missing_cols} not found in adata.obs")


def aggregate_group_expression(adata, group_mask, group_name: str) -> pd.DataFrame:
    """
    Aggregate gene expression for a specific group.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    group_mask
        Index or boolean mask selecting cells belonging to the group.
    group_name
        Name to use for the aggregated pseudo-bulk sample.

    Returns
    -------
    pandas.DataFrame
        Single-column DataFrame indexed by gene names.
    """
    adata_filtered = adata[group_mask, :]
    return adata_filtered.to_df().T.mean(axis=1).to_frame(name=group_name)


def aggregate_single_cell_data(adata, group_cols: List[str]) -> pd.DataFrame:
    """
    Aggregate single-cell expression data by the given group columns.

    Parameters
    ----------
    adata
        AnnData object with expression data.
    group_cols
        List of ``obs`` columns to group by (e.g. ``[\"sample\", \"cell_type\"]``).

    Returns
    -------
    pandas.DataFrame
        Pseudo-bulk expression matrix with groups as columns.
    """
    check_group_columns_exist(adata, group_cols)

    grouped = adata.obs.groupby(group_cols)
    df_list = []

    for group_values, group in grouped:
        group_name = combine_group_values(group_values)
        group_mask = group.index
        df = aggregate_group_expression(adata, group_mask, group_name)
        df_list.append(df)

    return pd.concat(df_list, axis=1)


def create_group_dicts(adata, group_cols: List[str]) -> list[dict]:
    """
    Create a list of dictionaries defining groups based on specified columns.

    Parameters
    ----------
    adata
        AnnData object whose ``obs`` is used to define groups.
    group_cols
        List of column names to group by.

    Returns
    -------
    list of dict
        List of dictionaries describing each group combination.
    """
    check_group_columns_exist(adata, group_cols)

    group_dict_list = (
        adata.obs.groupby(group_cols, observed=False)
        .size()
        .reset_index(name="size")
        .drop(columns="size")
        .to_dict("records")
    )
    return group_dict_list


def summarize_category_distribution(
    adata,
    category_col: str,
    group_col,
    selected_groups: list[str] | None = None,
    normalize: bool = False,
    sort_order: list[str] | None = None,
):
    """
    Summarize the distribution of a categorical variable across groups.

    Parameters
    ----------
    adata
        AnnData object containing the relevant categorical data in ``obs``.
    category_col
        Name of the column in ``obs`` whose distribution is summarized.
    group_col
        Column name or list of column names in ``obs`` defining groups.
    selected_groups
        Optional list of column-prefixes to keep (e.g. ``[\"Y\", \"M\", \"O\"]``).
    normalize
        If ``True``, express counts as percentages per column and add a
        ``\"Total\"`` row.
    sort_order
        Optional list defining the ordering of group prefixes.

    Returns
    -------
    pandas.io.formats.style.Styler or pandas.DataFrame
        When ``normalize`` is ``True``, returns a styled percentage table;
        otherwise a styled count table with an extra ``\"Total\"`` column.
    """
    if isinstance(group_col, str):
        group_col = [group_col]

    category_counts = pd.crosstab(
        adata.obs[category_col],
        [adata.obs[col] for col in group_col],
        colnames=group_col,
    )

    if isinstance(category_counts.columns, pd.MultiIndex):
        category_counts.columns = [
            "_".join(map(str, col)) for col in category_counts.columns
        ]

    if selected_groups is not None:
        category_counts = category_counts[
            [
                col
                for col in category_counts.columns
                if any(col.startswith(prefix) for prefix in selected_groups)
            ]
        ]

    if sort_order is not None:
        category_counts = category_counts[
            sorted(
                category_counts.columns,
                key=lambda x: (
                    (
                        sort_order.index(x.split("_")[0])
                        if x.split("_")[0] in sort_order
                        else float("inf")
                    ),
                    x,
                ),
            )
        ]

    if normalize:
        category_counts = category_counts.div(category_counts.sum(axis=0), axis=1) * 100
        total_row = pd.DataFrame(category_counts.sum(axis=0)).T
        total_row.index = ["Total"]
        category_counts = pd.concat([category_counts, total_row])

        return (
            category_counts.round(1)
            .style.format("{:.1f}")
            .background_gradient(cmap="coolwarm", axis=1)
        )

    category_counts["Total"] = category_counts.sum(axis=1)
    styled_df = category_counts.style.background_gradient(
        cmap="coolwarm", axis=1, subset=category_counts.columns[:-1]
    )
    return styled_df
