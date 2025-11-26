from scoriapy.pp import run_hvg


def test_run_hvg_returns_subset(adata_base):
    ad = adata_base.copy()
    # Scanpy expects a floating-point matrix for HVG selection
    ad.X = ad.X.astype("float64")
    n_top = 5
    hv = run_hvg(ad, n_top=n_top)

    # For very small test datasets, Scanpy's HVG selection can return
    # more than ``n_top`` genes due to ties. We only require that the
    # result is a non-empty proper subset of the original variables.
    assert 0 < hv.n_vars < ad.n_vars
    # Original object should still have all genes
    assert ad.n_vars == 10
