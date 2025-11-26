from itertools import combinations

import numpy as np
import pandas as pd

from scoriapy import tl


def test_get_pairwise_deg_smoke(adata_base):
    groups = list(adata_base.obs["group"].unique())
    if len(groups) < 2:
        return
    pair = (groups[0], groups[1])
    res = tl.get_pairwise_deg(adata_base, [pair], groupby_col="group")
    assert len(res) == 1
    df = list(res.values())[0]
    assert {"names", "logfoldchanges"}.issubset(df.columns)


def test_run_rank_genes_groups_and_get_de_genes(adata_base):
    ad = adata_base.copy()
    tl.run_rank_genes_groups(ad, groupby_attribute="group")
    genes = tl.get_de_genes(ad, group_col="group")
    assert isinstance(genes, dict)
    assert set(genes.keys()) == set(ad.obs["group"].unique())


def test_get_deg_df_smoke(adata_base):
    ad = adata_base.copy()
    test_df, filter_df = tl.get_deg_df(
        ad,
        cluster_key="group",
        groupby_attribute="group",
        rank_key="group_rank",
        query="pvals_adj < 1.0",
    )
    assert isinstance(test_df, pd.DataFrame)
    assert isinstance(filter_df, pd.DataFrame)


def test_detect_gene_flips():
    one_down = ["g1", "g2"]
    one_up = ["g3"]
    two_down = ["g3"]
    two_up = ["g1"]
    df = tl.detect_gene_flips(one_down, one_up, two_down, two_up, return_df=True)
    assert set(df["Gene"]) == {"g1", "g3"}


def test_aggregate_single_cell_data_and_groups(adata_base):
    agg = tl.aggregate_single_cell_data(adata_base, group_cols=["group"])
    assert isinstance(agg, pd.DataFrame)
    assert len(agg.columns) == len(adata_base.obs["group"].unique())

    group_dicts = tl.create_group_dicts(adata_base, group_cols=["group"])
    assert all("group" in d for d in group_dicts)


def test_summarize_category_distribution_smoke(adata_base):
    styled = tl.summarize_category_distribution(
        adata_base,
        category_col="phase",
        group_col="age",
        normalize=True,
        sort_order=["Y", "M", "O"],
    )
    # styled is a Styler; we just ensure it was created
    assert hasattr(styled, "data")


def test_compare_cluster_phase_smoke(adata_base):
    cross = pd.crosstab(
        index=[adata_base.obs["group"], adata_base.obs["age"]],
        columns=adata_base.obs["phase"],
        normalize="index",
    )
    cluster_id = adata_base.obs["group"].unique()[0]
    res = tl.compare_cluster_phase(
        cross_table=cross, cluster_id=cluster_id, level=0, group_level="age"
    )
    assert set(res.keys()).issubset({"G1", "S", "G2M"})


def test_compare_gene_expression_smoke(adata_base):
    # Smoke test: ensure it runs; p-value may or may not be significant
    tl.compare_gene_expression(
        adata_base,
        clustering_key="group",
        cluster_ids=list(adata_base.obs["group"].unique()),
        gene="gene0",
        groupby_col="age",
    )


def test_compare_gene_expression_two_smoke(adata_base):
    res = tl.compare_gene_expression_two(
        adata_base,
        clustering_key="group",
        cluster_ids=list(adata_base.obs["group"].unique()),
        gene="gene0",
        groupby_col="group",
    )
    assert {"gene", "statistic", "p_value"} <= set(res.keys())
