import logging

import numpy as np
import pandas as pd
import pandas.testing as pdt

from scoriapy.utils import subset_by_markers, setup_logger


def test_subset_by_markers_any_strategy(adata_base):
    # Make one gene clearly expressed in a subset
    ad = adata_base.copy()
    ad.X[:, 0] = 0
    ad.X[:10, 0] = 10  # first 10 cells express gene0 strongly

    obs_before = ad.obs.copy()
    subset = subset_by_markers(ad, markers=["gene0"], strategy="any", threshold=1.0)

    assert subset.n_obs == 10
    # Original obs remains unchanged
    pdt.assert_frame_equal(ad.obs, obs_before)


def test_subset_by_markers_add_obs_key(adata_base):
    ad = adata_base.copy()
    ad.X[:, 0] = 0
    ad.X[:5, 0] = 5

    subset = subset_by_markers(
        ad,
        markers=["gene0"],
        strategy="any",
        threshold=1.0,
        add_obs_key="marker_mask",
    )

    assert "marker_mask" in subset.obs  # mask present on returned object
    assert subset.n_obs == 5


def test_setup_logger_idempotent(tmp_path, monkeypatch):
    log_file = tmp_path / "test.log"
    logger1 = setup_logger(log_file=str(log_file), level=logging.DEBUG)
    n_handlers_first = len(logger1.handlers)

    # Second call should not add new handlers
    logger2 = setup_logger(log_file=str(log_file), level=logging.INFO)
    assert logger1 is logger2
    assert len(logger2.handlers) == n_handlers_first
