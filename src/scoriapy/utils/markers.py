import numpy as np
from typing import List, Optional, Literal


def subset_by_markers(
    adata,
    markers: List[str],
    strategy: Literal["any", "all", "count", "threshold"] = "any",
    threshold: float = 0.0,
    require_n_markers: int = 1,
    min_count: int = 1,
    use_raw: bool = False,
    add_obs_key: Optional[str] = None,
):
    """
    Subset cells based on expression of a set of marker genes.

    Parameters
    ----------
    adata
        AnnData object containing the expression matrix.
    markers
        List of marker gene names. Genes not found in ``adata.var_names`` are
        silently ignored.
    strategy
        Strategy used to define positive cells:

        - ``\"any\"`` – keep cells expressing *any* marker above ``threshold``.
        - ``\"all\"`` – keep cells expressing *all* markers above ``threshold``.
        - ``\"count\"`` – keep cells with at least ``min_count`` markers above
          ``threshold``.
        - ``\"threshold\"`` – keep cells whose **mean** marker expression
          exceeds ``threshold``.
    threshold
        Expression cutoff to consider a marker as expressed.
    require_n_markers
        Minimum number of markers required to be present in
        ``adata.var_names``. Currently this is advisory only and does not alter
        behaviour, but can be used by callers for validation.
    min_count
        Minimum number of expressed markers when ``strategy='count'``.
    use_raw
        If ``True`` and ``adata.raw`` is present, use ``adata.raw.X`` instead
        of ``adata.X``.
    add_obs_key
        Optional name of an ``obs`` column to store the boolean mask.

    Returns
    -------
    AnnData
        Copy of ``adata`` restricted to cells that satisfy the marker rule.
    """

    X = adata.raw.X if (use_raw and adata.raw is not None) else adata.X

    gene_idx = [adata.var_names.get_loc(g) for g in markers if g in adata.var_names]
    expr = X[:, gene_idx]

    if strategy == "any":
        mask = (expr > threshold).any(axis=1)
    elif strategy == "all":
        mask = (expr > threshold).all(axis=1)
    elif strategy == "count":
        mask = (expr > threshold).sum(axis=1) >= min_count
    elif strategy == "threshold":
        mask = expr.mean(axis=1) > threshold
    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    if add_obs_key:
        adata.obs[add_obs_key] = mask

    return adata[mask].copy()
