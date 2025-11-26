import numpy as np
import pandas as pd
import pandas.testing as pdt

from scoriapy.utils import (
    create_cell_mask,
    subset_adata_by_mask,
    get_mat_df,
    get_full_df,
    make_annotation_df,
    make_mat_dfs,
    scale_df,
)


def test_create_cell_mask_simple(adata_base):
    conditions = {"age": "Y"}
    mask = create_cell_mask(adata_base, conditions)
    assert mask.dtype == bool
    assert mask.shape[0] == adata_base.n_obs
    assert set(adata_base.obs["age"][mask]) == {"Y"}


def test_create_cell_mask_list_and_callable(adata_base):
    conditions = {
        "age": ["Y", "M"],
        "group": lambda s: s == "A",
    }
    mask = create_cell_mask(adata_base, conditions)
    assert mask.any()
    # All selected rows should satisfy both conditions
    sub = adata_base.obs[mask]
    assert set(sub["age"]).issubset({"Y", "M"})
    assert set(sub["group"]) == {"A"}


def test_subset_adata_by_mask_and_use_raw(adata_base):
    mask = np.array([True] * adata_base.n_obs)
    genes = ["gene0", "gene1"]
    subset = subset_adata_by_mask(adata_base, mask, genes, use_raw=True)
    assert subset.n_vars == len(genes)
    # Ensure original is unchanged
    assert adata_base.n_vars == 10


def test_get_mat_df_and_get_full_df(adata_base):
    keys = ["gene0", "gene1", "gene2"]
    mat = get_mat_df(adata_base, keys=keys, use_raw=True)
    assert list(mat.columns) == keys

    full = get_full_df(adata_base, keys=keys, scale=True)
    assert list(full.columns) == keys
    # Rough z-score property
    np.testing.assert_allclose(full.mean(axis=0).values, 0.0, atol=1e-7)


def test_make_annotation_df_and_make_mat_dfs(adata_base):
    # Build a minimal dict with the expected keys
    adata_age = adata_base.copy()
    adata_age.obs["age"] = adata_age.obs["age"].astype(str)

    adata_rapa = adata_base.copy()
    adata_rapa.obs["treatment"] = ["Ctrl"] * adata_rapa.n_obs

    adata_dict = {"age": adata_age, "rapa": adata_rapa}

    ann = make_annotation_df(adata_dict)
    assert "group" in ann.columns
    assert ann.shape[0] == adata_age.n_obs + adata_rapa.n_obs

    common_genes = ["gene0", "gene1"]
    mat_dfs = make_mat_dfs(adata_dict, common_genes)
    assert set(mat_dfs.keys()) == {"age", "rapa"}
    for df in mat_dfs.values():
        assert list(df.index) == common_genes


def test_scale_df():
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [4.0, 5.0, 6.0]})
    scaled = scale_df(df)
    np.testing.assert_allclose(scaled.mean(axis=0).values, 0.0, atol=1e-7)
    np.testing.assert_allclose(scaled.std(axis=0, ddof=0).values, 1.0, atol=1e-7)
