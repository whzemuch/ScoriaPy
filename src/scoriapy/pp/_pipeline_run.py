import scanpy as sc
from ..utils.logging import setup_logger

def run_scanpy_basic_pipeline(
        adata,
        n_top_genes=2000,
        n_pcs=50,
        n_neighbors=15,
        resolution=0.5,
        hvg=True,
        scale=True,
        copy=True,
        logger=None,
    ):
    """
    Run standard Scanpy preprocessing with logging at each step:
    HVG → scale → PCA → neighbors → UMAP → Leiden.

    This pipeline automates the transition from raw/normalized counts to a 
    clustered embedding. It handles high-variable gene selection, data scaling 
    to unit variance, dimensionality reduction via PCA, and manifold learning 
    (UMAP/Leiden).

    Args:
        adata (AnnData): Annotated data matrix (typically normalized and log-transformed).
        n_top_genes (int): Number of highly variable genes to keep. Default is 2000.
        n_pcs (int): Number of principal components to compute and use for neighbors.
        n_neighbors (int): Size of the local neighborhood for the KNN graph.
        resolution (float): Resolution parameter for Leiden clustering (higher = more clusters).
        hvg (bool): If True, subsets the object to highly variable genes.
        scale (bool): If True, scales data to unit variance and zero mean.
        copy (bool): If True, returns a new AnnData object; otherwise modifies in-place.
        logger (logging.Logger, optional): Logger object for progress tracking. 
            If None, calls setup_logger().

    Returns:
        AnnData: The processed AnnData object with PCA, UMAP, and Leiden results.

    Example:
        >>> import scanpy as sc
        >>> # Load raw data and perform basic normalization
        >>> adata = sc.datasets.pbmc3k()
        >>> sc.pp.filter_cells(adata, min_genes=200)
        >>> sc.pp.filter_genes(adata, min_cells=3)
        >>> sc.pp.normalize_total(adata, target_sum=1e4)
        >>> sc.pp.log1p(adata)
        >>> 
        >>> # Run the automated pipeline
        >>> logger = setup_logger("preprocess.log")
        >>> adata_processed = run_scanpy_basic_pipeline(
        ...     adata, 
        ...     logger=logger
        ...     n_top_genes=2000, 
        ...     resolution=0.8
        ... )
        >>> 
        >>> # Visualize the results
        >>> sc.pl.umap(adata_processed, color='leiden')
    """

    if logger is None:
        logger = setup_logger()

    if copy:
        adata = adata.copy()

    # 1. HVG
    if hvg:
        logger.info("Computing highly variable genes...")
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
        adata = adata[:, adata.var["highly_variable"]]
        logger.info(f"Selected top {n_top_genes} HVGs.")
    else:
        logger.info("Skipping HVG selection.")

    # 2. Scaling
    if scale:
        logger.info("Scaling data...")
        sc.pp.scale(adata, max_value=10)
        logger.info("Scaling complete.")
    else:
        logger.info("Skipping scaling.")

    # 3. PCA
    logger.info(f"Running PCA (n_comps={n_pcs})...")
    sc.pp.pca(adata, n_comps=n_pcs)
    logger.info("PCA complete.")

    # 4. KNN graph
    logger.info(f"Computing neighbors (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    logger.info("Neighbors graph complete.")

    # 5. UMAP
    logger.info("Computing UMAP embedding...")
    sc.tl.umap(adata)
    logger.info("UMAP complete.")

    # 6. Leiden
    logger.info(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)
    adata.obs["leiden"] = adata.obs["leiden"].astype("category")
    logger.info("Leiden clustering complete.")

    logger.info("✨ Scanpy pipeline finished successfully!")

    return adata