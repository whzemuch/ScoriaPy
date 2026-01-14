from joblib import Parallel, delayed
import scanpy as sc

def run_resolution_sweep_parallel(adata, resolutions=[0.2, 0.4, 0.6, 0.8, 1.0], n_jobs=-1):
    """
    Perform a parallelized Leiden clustering sweep across multiple resolutions.

    This function utilizes multi-core processing to accelerate the computation of 
    multiple clustering resolutions. It uses the 2026 standard 'igraph' flavor 
    to ensure performance and compatibility. Each resolution is stored as a 
    new column in `adata.obs`.

    Args:
        adata (AnnData): Annotated data matrix. Must have a neighbors graph 
            pre-computed in `.obsp['connectivities']`.
        resolutions (list of float): Resolutions to evaluate. 
            Default is [0.2, 0.4, 0.6, 0.8, 1.0].
        n_jobs (int): Number of CPU cores to use. -1 uses all available processors.

    Returns:
        AnnData: The input object with added columns in `.obs` (e.g., 'leiden_0.2').

    Example:
        >>> import scanpy as sc
        >>> adata = sc.datasets.pbmc3k_processed()
        >>> # Ensure neighbors are computed
        >>> sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
        >>> # Run parallel sweep
        >>> resolutions = [0.1, 0.5, 1.0]
        >>> adata = run_resolution_sweep_parallel(adata, resolutions=resolutions)
        >>> # Visualize one of the results
        >>> sc.pl.umap(adata, color='leiden_0.5')
    """
    if 'connectivities' not in adata.obsp:
        raise ValueError("Adjacency matrix not found. Run sc.pp.neighbors(adata) first.")

    def _compute_leiden(res):
        # We explicitly pass the 2026 defaults to suppress warnings and ensure speed:
        # flavor="igraph", n_iterations=2, and directed=False
        return sc.tl.leiden(
            adata, 
            resolution=res, 
            copy=True, 
            flavor="igraph", 
            n_iterations=2, 
            directed=False
        ).obs['leiden']

    print(f"Running sweep across {len(resolutions)} resolutions using {n_jobs} cores...")
    results = Parallel(n_jobs=n_jobs)(delayed(_compute_leiden)(r) for r in resolutions)

    for r, labels in zip(resolutions, results):
        adata.obs[f"leiden_{r}"] = labels
    
    return adata