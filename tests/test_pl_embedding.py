import matplotlib.figure as mpl_fig
import pandas.testing as pdt

from scoriapy import pl


def test_plot_umap_by_categories_no_side_effects(adata_with_umap):
    ad = adata_with_umap.copy()
    obs_before = ad.obs.copy()

    pl.plot_umap_by_categories(ad, obs_columns=["age"], basis="X_umap", color="leiden")

    # Shape and obs unchanged
    assert ad.shape == adata_with_umap.shape
    pdt.assert_frame_equal(ad.obs, obs_before)


def test_calculate_median_position_keys_match_groups(adata_with_umap):
    positions = pl.calculate_median_position(
        adata_with_umap, basis="X_umap", overlay_attribute="age"
    )
    # Same labels as unique ages
    assert set(positions.keys()) == set(adata_with_umap.obs["age"].unique())
    # Each position is length-2
    assert all(len(v) == 2 for v in positions.values())


def test_plot_embeddings_with_labels_smoke(adata_with_umap):
    pl.plot_embeddings_with_labels(
        adata_with_umap,
        basis="X_umap",
        color_attribute="leiden",
        overlay_attribute="leiden",
        group_col=None,
    )


def test_plot_umap_by_group_smoke(adata_with_umap):
    pl.plot_umap_by_group(adata_with_umap, group_column="age")


def test_plot_pca_for_samples_smoke(adata_with_umap):
    # Ensure sample is categorical for fixture
    pl.plot_pca_for_samples(adata_with_umap, sample_list=["s1", "s2"])


def test_plot_umap_grid_by_sample_returns_figure(adata_with_umap):
    fig = pl.plot_umap_grid_by_sample(adata_with_umap, col_per_row=2)
    assert isinstance(fig, mpl_fig.Figure)


def test_plot_umap_param_grid_smoke(adata_with_umap):
    # Use very small grid to keep it light
    pl.plot_umap_param_grid(adata_with_umap, min_dists=(0.1,), spreads=(0.5,))
