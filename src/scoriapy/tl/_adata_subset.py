import scanpy as sc
from anndata import AnnData
from typing import Union, List


def subset_and_clean_adata(
    adata: AnnData, 
    obs_col: str, 
    keys: Union[str, List[str]]
) -> Union[AnnData, None]:
    """
    Extract a subset of cells from an AnnData object and clean categorical metadata.

    This function filters the input AnnData by one or more keys in a specified 
    ``.obs`` column. It automatically removes unused categories from all 
    categorical columns in the resulting subset to prevent downstream errors 
    during aggregation (e.g., PyDESeq2 or pseudobulking) and plotting.

    Parameters
    ----------
    adata
        The input AnnData object to be subsetted.
    obs_col
        The column name in ``adata.obs`` used for filtering (e.g., 'cluster_name').
    keys
        A single string or a list of strings representing the categories to keep.

    Returns
    -------
    Optional[AnnData]
        A copied subset of the original AnnData containing only the specified cells.
        Returns ``None`` if no cells match the provided keys.

    Raises
    ------
    ValueError
        If ``obs_col`` is not found in ``adata.obs``.

    Examples
    --------
    >>> # Extract a single cell type
    >>> myo_adata = extract_cell_cluster_pro(adata, "cluster_name", "Myo")
    >>>
    >>> # Combine multiple cluster IDs into one subset
    >>> epi_adata = extract_cell_cluster_pro(adata, "cluster-label", ["1", "5", "12"])
    """
    # 1. Validation
    if obs_col not in adata.obs:
        raise ValueError(f"'{obs_col}' is not a valid column in adata.obs")

    if not isinstance(keys, list):
        keys = [keys]

    # 2. Check for existence of keys in the data
    missing_keys = [key for key in keys if key not in adata.obs[obs_col].unique()]
    if missing_keys:
        print(f"Warning: The following keys were not found in '{obs_col}': {missing_keys}")

    # 3. Create boolean mask
    mask = adata.obs[obs_col].isin(keys)
    if not mask.any():
        print(f"No cells match the specified keys: {keys}. Check your input and try again.")
        return None

    # 4. Subset and Metadata Cleanup
    # We use .copy() to ensure the subset is an independent object
    subset = adata[mask].copy()
    
    # 2026 Best Practice: Remove 'ghost' categories that cause aggregation errors
    for col in subset.obs.select_dtypes(["category"]).columns:
        subset.obs[col] = subset.obs[col].cat.remove_unused_categories()
        
    print(f"Successfully created subset for {keys} ({subset.n_obs} cells).")
    return subset