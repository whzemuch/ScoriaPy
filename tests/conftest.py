import os
import sys

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

# Ensure src/ is on sys.path so `import scoriapy` works without setting PYTHONPATH.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC_DIR = os.path.join(ROOT_DIR, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)


@pytest.fixture
def adata_base():
    """Synthetic AnnData with basic obs/var structure for tests."""
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, size=(50, 10))  # 50 cells, 10 genes
    var_names = [f"gene{i}" for i in range(10)]

    obs = pd.DataFrame(
        {
            "group": pd.Categorical(rng.choice(["A", "B"], size=50)),
            "age": pd.Categorical(rng.choice(["Y", "M", "O"], size=50)),
            "phase": pd.Categorical(rng.choice(["G1", "S", "G2M"], size=50)),
            "sample": pd.Categorical(rng.choice(["s1", "s2"], size=50)),
            "sample_name": pd.Categorical(rng.choice(["s1", "s2"], size=50)),
        },
        index=[f"cell{i}" for i in range(50)],
    )

    adata = sc.AnnData(X=X, obs=obs, var=pd.DataFrame(index=var_names))
    adata.raw = adata.copy()
    return adata


@pytest.fixture
def adata_with_umap(adata_base):
    """AnnData with neighbors, UMAP embedding, and synthetic cluster labels."""
    ad = adata_base.copy()
    sc.pp.normalize_total(ad)
    sc.pp.log1p(ad)
    sc.pp.pca(ad)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)

    # Create a synthetic categorical 'leiden' label without relying on igraph
    rng = np.random.default_rng(1)
    clusters = rng.integers(0, 3, size=ad.n_obs).astype(str)
    ad.obs["leiden"] = pd.Categorical(clusters)

    categories = ad.obs["leiden"].cat.categories
    # Simple deterministic color palette
    palette = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
    ]
    ad.uns["leiden_colors"] = palette[: len(categories)]
    return ad


@pytest.fixture
def adata_for_heatmap(adata_base):
    """AnnData with gene categories for heatmap utilities."""
    ad = adata_base.copy()
    ad.raw = ad.copy()
    ad.var["gene_cat"] = np.random.default_rng(1).choice(["A", "B"], size=ad.n_vars)
    return ad
