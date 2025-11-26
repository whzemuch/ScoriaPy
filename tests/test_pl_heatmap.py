import numpy as np
import pandas as pd
import pandas.testing as pdt

from scoriapy import pl


def test_subset_adata_filters_cells_and_genes(adata_for_heatmap):
    genes_map = {g: "A" for g in adata_for_heatmap.var_names[:3]}
    subset = pl.subset_adata(
        adata_for_heatmap, group_col="age", keys=["Y"], genes_map=genes_map
    )

    assert set(subset.var_names) == set(genes_map.keys()) & set(
        adata_for_heatmap.raw.var_names
    )
    assert set(subset.obs["age"]) == {"Y"}
    assert "gene_cat" in subset.var


def test_sort_cells_by_metadata_returns_new_object(adata_for_heatmap):
    ad = adata_for_heatmap.copy()
    obs_before = ad.obs.copy()

    sorted_ad = pl.sort_cells_by_metadata(
        ad,
        primary_sort_col="age",
        primary_categories=["Y", "M", "O"],
        second_sort_col="phase",
        second_categories=["G1", "S", "G2M"],
    )

    # Original unchanged
    pdt.assert_frame_equal(ad.obs, obs_before)
    # Sorted object has same index but in sorted order
    assert set(sorted_ad.obs.index) == set(ad.obs.index)


def test_get_heatmap_dfs_shapes(adata_for_heatmap):
    res = pl.get_heatmap_dfs(
        adata_for_heatmap,
        obs_cols=["age", "phase"],
        var_cols=["gene_cat"],
        scale=False,
    )
    assert set(res.keys()) == {"df_mat", "col_df", "row_df"}
    df_mat = res["df_mat"]
    col_df = res["col_df"]
    row_df = res["row_df"]

    assert df_mat.shape[1] == adata_for_heatmap.n_obs
    assert col_df.shape[0] == adata_for_heatmap.n_obs
    assert row_df.shape[0] == adata_for_heatmap.n_vars


def test_plot_distribution_smoke(adata_for_heatmap):
    res = pl.get_heatmap_dfs(
        adata_for_heatmap,
        obs_cols=["age"],
        var_cols=["gene_cat"],
        scale=True,
    )
    pl.plot_distribution(res["df_mat"], n=10, x_intercept=-1.0)


def test_plot_cluster_heatmap_smoke(adata_for_heatmap):
    res = pl.get_heatmap_dfs(
        adata_for_heatmap,
        obs_cols=["age", "phase"],
        var_cols=["gene_cat"],
        scale=False,
    )
    pl.plot_cluster_heatmap(res["df_mat"], res["col_df"], res["row_df"], figsize=(4, 4))
