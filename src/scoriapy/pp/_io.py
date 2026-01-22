"""
Input/output helpers for reading scRNA-seq data into AnnData.

These functions are light wrappers around scanpy's 10x readers, with
some convenience columns added to ``.obs``.

Most users will access these via :mod:`scoriapy.pp`.
"""

from __future__ import annotations

import scanpy as sc
import anndata as ad
import pandas as pd
import logging
from pathlib import Path
from typing import Dict, Callable, Union, Optional
from tqdm import tqdm


def read_and_concat(
    sample_dict: Dict[str, Union[str, Path]],
    process_func: Callable[[Path, str], ad.AnnData],
    logger: Optional[logging.Logger] = None,
    **concat_kwargs
) -> ad.AnnData:
    """
    Load multiple samples using a custom processing function and concatenate them.

    This orchestrator applies a user-defined loading function to each sample path 
    in a dictionary and merges the results into a single AnnData object using 
    the modern ``anndata.concat`` API.

    Parameters
    ----------
    sample_dict
        Mapping of sample names to their respective file paths.
    process_func
        A callable that takes (Path, str) and returns an AnnData object. 
        Highly recommended to use ``functools.partial`` to pass extra 
        metadata or project settings.
    logger
        A :class:`logging.Logger` object for progress tracking. 
        If None, a basic console logger is initialized.
    **concat_kwargs
        Keyword arguments passed to :func:`anndata.concat`. 
        Commonly used: ``join='inner'``, ``uns_merge='unique'``.

    Returns
    -------
    ad.AnnData
        Concatenated AnnData object with ``.obs['sample_name']`` populated.

    Example
    -------
    >>> import pandas as pd
    >>> import scanpy as sc
    >>> from functools import partial
    >>> # 1. Define metadata and sample paths
    >>> meta = pd.DataFrame({"age": ["6m", "24m"]}, index=["S1", "S2"])
    >>> paths = {"S1": "data/s1.h5", "S2": "data/s2.h5"}
    >>> 
    >>> # 2. Define a loader that requires the metadata table
    >>> def my_loader(path, name, info_df):
    ...     adata = sc.read_10x_h5(path)
    ...     adata.var_names_make_unique()
    ...     if name in info_df.index:
    ...         for col in info_df.columns:
    ...             adata.obs[col] = info_df.loc[name, col]
    ...     return adata
    >>> 
    >>> # 3. Create a partial function to 'bake in' the metadata
    >>> loader_with_meta = partial(my_loader, info_df=meta)
    >>> 
    >>> # 4. Run orchestrator
    >>> adata = read_and_concat(paths, loader_with_meta, join="inner")
    """
    if logger is None:
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
        logger = logging.getLogger(__name__)

    adatas = {}
    
    # 1. Sequential Loading with Progress Bar
    for name, path in tqdm(sample_dict.items(), desc="Processing Samples"):
        logger.info(f"Loading {name} from {path}")
        
        try:
            # Execute the (potentially partial) function
            adata = process_func(Path(path), name)
            
            # Validation check
            if not isinstance(adata, ad.AnnData):
                raise TypeError(f"process_func must return AnnData, got {type(adata)}")
                
            logger.info(f"Success: {name} ({adata.n_obs} cells)")
            adatas[name] = adata
            
        except Exception as e:
            logger.error(f"Failed to process sample '{name}': {str(e)}")
            raise

    # 2. Modern Concatenation
    logger.info("Concatenating samples...")
    combined = ad.concat(
        adatas, 
        label="sample_name", 
        index_unique="_", 
        **concat_kwargs
    )
    
    logger.info(f"Final object: {combined.n_obs} cells x {combined.n_vars} genes.")
    
    return combined







def read_10x_h5_concat(
    samples: Iterable[Tuple[str, Path]],
    sample_info: pd.DataFrame,
) -> sc.AnnData:
    """
    Read multiple 10x HDF5 matrices and concatenate them into one AnnData.

    Parameters
    ----------
    samples
        Iterable of ``(sample_name, path_to_filtered_feature_bc_matrix.h5)``.
    sample_info
        DataFrame indexed by ``sample_name`` that contains at least a
        ``\"treatment\"`` column (additional columns are kept in the index
        and can be joined later if needed).

    Returns
    -------
    AnnData
        Concatenated object with ``sample_name`` and ``treatment`` in ``.obs``.
    """
    samples_iter = iter(samples)
    first_sample_name, first_sample_path = next(samples_iter)

    concatenated = sc.read_10x_h5(first_sample_path)
    concatenated.obs["sample_name"] = first_sample_name
    concatenated.obs["treatment"] = sample_info.loc[first_sample_name, "treatment"]
    concatenated.var_names_make_unique()

    for sample_name, sample_path in samples_iter:
        adata = sc.read_10x_h5(sample_path)
        adata.obs["sample_name"] = sample_name
        adata.obs["treatment"] = sample_info.loc[sample_name, "treatment"]
        adata.var_names_make_unique()

        concatenated = concatenated.concatenate(adata)

    return concatenated


def read_10x_mtx_concat(
    main_folder_path: Path | str,
    batch_key: str = "batch",
) -> sc.AnnData:
    """
    Read 10x MTX folders inside ``main_folder_path`` and concatenate them.

    Each direct subdirectory is assumed to contain a matrix in 10x MTX format.

    Parameters
    ----------
    main_folder_path
        Path to the main folder containing subfolders with 10x MTX data.
    batch_key
        Name of the column in ``.obs`` indicating the batch/subfolder.

    Returns
    -------
    AnnData
        Concatenated object with batches annotated by ``batch_key``.
    """
    main_folder = Path(main_folder_path)
    subfolders = [
        subfolder for subfolder in main_folder.iterdir() if subfolder.is_dir()
    ]

    if not subfolders:
        raise ValueError(f"No subfolders found in {main_folder!s}")

    adatas = []
    for subfolder in subfolders:
        adata = sc.read_10x_mtx(subfolder, var_names="gene_symbols", cache=True)
        adatas.append(adata)

    concatenated = adatas[0].concatenate(*adatas[1:], batch_key=batch_key)
    return concatenated
