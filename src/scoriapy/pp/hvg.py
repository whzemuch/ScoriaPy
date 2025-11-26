import numpy as np
import scanpy as sc
from scoriapy.utils import setup_logger


def run_hvg(adata, n_top: int = 3000, flavor: str = "seurat", log=None):
    """
    Compute and subset highly variable genes (HVGs).

    This is a thin wrapper around :func:`scanpy.pp.highly_variable_genes` that
    logs progress and returns a view restricted to the selected genes. If the
    data has not been log-transformed yet (i.e. ``\"log1p\"`` is not present in
    ``adata.uns``), a ``log1p`` transform is applied before selecting HVGs.

    Parameters
    ----------
    adata
        Input AnnData object.
    n_top
        Number of highly variable genes to keep.
    flavor
        HVG selection flavor passed to Scanpy (e.g. ``\"seurat\"``).
    log
        Optional :class:`logging.Logger` instance. If ``None``, a default
        scoriapy logger is created.

    Returns
    -------
    AnnData
        Copy of ``adata`` restricted to the selected HVGs.
    """
    logger = log or setup_logger()

    # Ensure data is in a floating dtype to avoid issues inside Scanpy
    if not np.issubdtype(adata.X.dtype, np.floating):
        adata.X = adata.X.astype(float)

    # Apply log1p once if it has not been done already
    if "log1p" not in adata.uns:
        logger.info("Applying log1p transform before HVG selection.")
        sc.pp.log1p(adata)

    logger.info(f"Computing HVGs (flavor={flavor}, top={n_top})")
    sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=n_top)
    return adata[:, adata.var["highly_variable"]].copy()
