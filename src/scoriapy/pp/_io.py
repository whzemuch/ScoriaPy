"""
Input/output helpers for reading scRNA-seq data into AnnData.

These functions are light wrappers around scanpy's 10x readers, with
some convenience columns added to ``.obs``.

Most users will access these via :mod:`scoriapy.pp`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd
import scanpy as sc


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
