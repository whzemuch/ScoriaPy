import logging
import scanpy as sc


def make_lean_adata_with_obs(
    adata,
    obs_keep,
    uns_keep=None,
    keep_X=True,
):
    """
    Create a reduced (lean) AnnData object by keeping only selected metadata,
    optionally retaining the expression matrix, and pruning bulky or
    non-serializable slots.

    This function is designed to shrink AnnData objects for:
    - exporting to R (zellkonverter / SingleCellExperiment)
    - saving compact .h5ad files
    - sharing results
    - keeping only DEG/GSEA summaries

    Behavior
    --------
    Always keeps:
        - selected ``obs`` columns
        - full ``var`` (gene metadata)

    Optionally keeps:
        - ``X`` (expression matrix) if ``keep_X=True``

    Always drops:
        - ``layers``
        - ``obsm`` (PCA/UMAP/embeddings)
        - ``varm``

    For ``uns``:
        - only keys listed in ``uns_keep`` are copied
        - if an entry is a Scanpy DE result containing ``{"df": ...}``,
          only the DataFrame is preserved and large ``params`` dictionaries
          are removed to ensure compatibility with zellkonverter

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.

    obs_keep : list[str]
        Column names in ``adata.obs`` to retain.

        Example:
        ["sample", "age", "cell_type", "condition"]

    uns_keep : list[str] or None, default=None
        Keys from ``adata.uns`` to copy into the new object.

        Special behavior:
        - DEG results produced by Scanpy wrappers that store
          ``uns[key]["df"]`` will keep only the DataFrame.

    keep_X : bool, default=True
        Whether to retain the expression matrix.

        - True  → keep counts/expression (analysis-ready)
        - False → drop X (metadata-only, very small object)

    Returns
    -------
    AnnData
        A pruned AnnData object containing only essential components.

    Notes
    -----
    Why drop ``params`` from ``uns``?

    Scanpy's ``rank_genes_groups`` stores Python-specific objects inside
    ``uns[key]["params"]`` that are not serializable by zellkonverter and
    may cause conversion failures. This function keeps only safe DataFrames.

    Dropping ``X`` typically reduces file size by 90–99% for large
    single-cell datasets.

    Examples
    --------
    Example 1 — metadata-only export for R / DESeq2 / clusterProfiler:

    >>> adata_meta = make_lean_adata_with_obs(
    ...     adata,
    ...     obs_keep=["sample", "age", "cell_type"],
    ...     uns_keep=["OvsY", "MvsY"],
    ...     keep_X=False
    ... )
    >>> adata_meta.write("metadata_only.h5ad")

    Example 2 — keep counts for continued Scanpy analysis:

    >>> adata_small = make_lean_adata_with_obs(
    ...     adata,
    ...     obs_keep=["sample", "age"],
    ...     keep_X=True
    ... )

    Example 3 — after running DE with Scanpy:

    >>> comparisons = {"OvsY": ("Old", "Young")}
    >>> run_and_store_pairwise_rg(adata, groupby="age", comparisons=comparisons)

    >>> adata_export = make_lean_adata_with_obs(
    ...     adata,
    ...     obs_keep=["age", "cell_type"],
    ...     uns_keep=["OvsY"],
    ...     keep_X=False
    ... )

    >>> adata_export.uns["OvsY"]["df"].head()

    Typical workflow
    ----------------
    DE → prune → save:

    >>> run_and_store_pairwise_rg(adata, "age", comparisons)
    >>> adata_out = make_lean_adata_with_obs(
    ...     adata,
    ...     obs_keep=["sample", "age"],
    ...     uns_keep=list(comparisons.keys()),
    ...     keep_X=False
    ... )
    >>> adata_out.write("results_only.h5ad")
    """
    import anndata as ad
    import copy

    X = adata.X if keep_X else None

    adata_new = ad.AnnData(
        X=X,
        obs=adata.obs[obs_keep].copy(),
        var=adata.var.copy(),
    )

    # drop heavy slots
    adata_new.layers = {}
    adata_new.obsm = {}
    adata_new.varm = {}

    # copy selected uns entries safely
    if uns_keep is not None:
        adata_new.uns = {}

        for k in uns_keep:
            if k not in adata.uns:
                continue

            v = adata.uns[k]

            if isinstance(v, dict) and "df" in v:
                adata_new.uns[k] = {"df": v["df"].copy()}
            else:
                adata_new.uns[k] = copy.deepcopy(v)

    return adata_new