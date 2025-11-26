"""
Helpers for extracting and reshaping data from AnnData.

These functions are generic utilities and do not perform any plotting.
"""

from __future__ import annotations

from typing import Callable, List, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc


def create_cell_mask(adata: sc.AnnData, conditions: dict) -> np.ndarray:
    """
    Create a boolean mask for cells based on multiple flexible conditions.

    Parameters
    ----------
    adata
        The AnnData object.
    conditions
        A mapping of ``obs`` column name to a condition. The condition can be:

        - a single value (``str``, ``int``, ``float``)
        - a list/tuple/array/Series of allowed values
        - a callable that takes a Series and returns a boolean Series

    Returns
    -------
    np.ndarray
        Boolean mask over observations.
    """

    def apply_callable(column_data: pd.Series, condition: Callable) -> pd.Series:
        return condition(column_data)

    def apply_list_like(
        column_data: pd.Series,
        condition: Union[List, Tuple, np.ndarray, pd.Series],
    ) -> pd.Series:
        return column_data.isin(condition)

    def apply_single_value(
        column_data: pd.Series,
        condition: Union[str, int, float],
    ) -> pd.Series:
        return column_data == condition

    condition_handlers = {
        Callable: apply_callable,
        (list, tuple, np.ndarray, pd.Series): apply_list_like,
        (str, int, float): apply_single_value,
    }

    def apply_condition(column_data: pd.Series, condition: object) -> pd.Series:
        for types, handler in condition_handlers.items():
            if isinstance(condition, types):
                return handler(column_data, condition)
        raise ValueError(f"Unsupported condition type: {type(condition)}")

    masks = []
    for column, condition in conditions.items():
        if column not in adata.obs.columns:
            raise ValueError(f"Column '{column}' not found in adata.obs")

        mask = apply_condition(adata.obs[column], condition)
        masks.append(mask)

    return np.all(masks, axis=0)


def subset_adata_by_mask(
    adata: sc.AnnData,
    cell_mask: np.ndarray,
    gene_list: list,
    use_raw: bool = True,
) -> sc.AnnData:
    """
    Subset an AnnData object using a cell mask and gene list.

    Parameters
    ----------
    adata
        The original AnnData object.
    cell_mask
        Boolean mask for selecting cells.
    gene_list
        Genes to keep.
    use_raw
        Whether to use ``adata.raw`` for gene subsetting.

    Returns
    -------
    AnnData
        A new subsetted AnnData object.
    """
    if use_raw:
        if adata.raw is None:
            raise ValueError(
                "use_raw is True but the AnnData object "
                "does not have a .raw attribute."
            )
        adata_subset = adata[cell_mask, :].raw.to_adata()
    else:
        adata_subset = adata[cell_mask, :]

    gene_mask = adata_subset.var_names.isin(gene_list)
    return adata_subset[:, gene_mask]


def extract_cell_cluster(adata: sc.AnnData, col_obs: str, keys) -> sc.AnnData | None:
    """
    Extract a subset of cells from an AnnData object based on values in ``.obs``.

    Parameters
    ----------
    adata
        AnnData object.
    col_obs
        Column in ``adata.obs`` used for filtering.
    keys
        Single key or list of keys to keep.

    Returns
    -------
    AnnData or None
        Subsetted AnnData or ``None`` if no cells match.
    """
    if col_obs not in adata.obs:
        raise ValueError(f"{col_obs} is not a valid column in adata.obs")

    if not isinstance(keys, list):
        keys = [keys]

    missing_keys = [key for key in keys if key not in adata.obs[col_obs].unique()]
    if missing_keys:
        print(f"The following keys are not in {col_obs}: {missing_keys}")

    mask = adata.obs[col_obs].isin(keys).to_numpy()
    if not mask.any():
        print("No cells match the specified keys. Check your keys and try again.")
        return None

    return adata[mask].copy()


def get_mat_df(
    adata: sc.AnnData,
    keys: list[str] | None = None,
    use_raw: bool = True,
) -> pd.DataFrame:
    """
    Extract expression data as a DataFrame.

    Parameters
    ----------
    adata
        AnnData object.
    keys
        Variable names to extract. If ``None``, all variables are used.
    use_raw
        If ``True``, use ``adata.raw``; otherwise use ``adata.X``.

    Returns
    -------
    pandas.DataFrame
        Expression matrix with observations as rows.
    """
    if use_raw and adata.raw is None:
        raise ValueError(
            "The provided AnnData object does not contain a .raw attribute."
        )

    if keys is not None and not isinstance(keys, list):
        raise ValueError("`keys` should be a list of strings or None.")

    if keys is None:
        keys = adata.raw.var_names if use_raw else adata.var_names

    data = sc.get.obs_df(adata, keys=keys, use_raw=use_raw)
    return data


def get_full_df(
    adata: sc.AnnData,
    keys: list[str] | None = None,
    scale: bool = True,
) -> pd.DataFrame:
    """
    Extract expression data from ``adata.raw`` and optionally z-score scale it.

    Parameters
    ----------
    adata
        AnnData object.
    keys
        Variable names to extract. If ``None``, all variables are used.
    scale
        If ``True``, apply z-score column-wise.

    Returns
    -------
    pandas.DataFrame
        Expression matrix (optionally scaled).
    """
    if adata.raw is None:
        raise ValueError(
            "The provided AnnData object does not contain a .raw attribute."
        )

    if keys is None:
        keys = list(adata.raw.var_names)

    data = sc.get.obs_df(adata, keys=keys, use_raw=True)

    if scale:
        data = scale_df(data)

    return data


def make_annotation_df(adata_dict: dict) -> pd.DataFrame:
    """
    Combine metadata from multiple AnnData objects into a single DataFrame.

    Currently expects keys ``\"age\"`` and ``\"rapa\"`` and normalizes the
    column name to ``\"group\"``.
    """
    age_df = adata_dict["age"].obs[["age"]].rename(columns={"age": "group"})
    rapa_df = (
        adata_dict["rapa"].obs[["treatment"]].rename(columns={"treatment": "group"})
    )
    combined_df = pd.concat([age_df, rapa_df], axis=0)
    return combined_df


def make_mat_dfs(adata_dict: dict, common_genes: list[str]) -> dict[str, pd.DataFrame]:
    """
    Extract scaled expression matrices for a list of common genes.

    Parameters
    ----------
    adata_dict
        Mapping of dataset name to AnnData.
    common_genes
        Genes to extract.

    Returns
    -------
    dict
        Mapping from dataset name to a transposed expression DataFrame.
    """
    dfs: dict[str, pd.DataFrame] = {}
    for name, adata in adata_dict.items():
        dfs[name] = get_full_df(adata, keys=common_genes, scale=True).T
    return dfs


def scale_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply z-score normalization to each column of a DataFrame.

    Parameters
    ----------
    df
        Input DataFrame.

    Returns
    -------
    pandas.DataFrame
        DataFrame with standardized columns (mean=0, std=1).
    """
    return (df - df.mean()) / df.std(ddof=0)
