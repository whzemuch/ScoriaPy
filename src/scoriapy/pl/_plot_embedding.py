import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
# Now support basis using X_umap or umap

def plot_umap_by_categories(
    adata,
    obs_columns,
    basis="X_umap",
    color="leiden",
    ncol=3,
    figsize=(15, 5),
    sort_order=None,
    cluster_to_cell_mapping=None,
    legend_location="on_data",
    point_size=None,
    return_fig=False,
):
    """
    Plots embeddings (e.g., UMAP) for each unique combination of categories in specified columns of adata.obs,
    ensuring consistent color mapping, optional sorting, and a single global legend.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    obs_columns : list[str]
        Columns in adata.obs to group by (e.g., ["age", "sample"]).
    basis : str
        Embedding basis in adata.obsm (e.g., "umap" or "X_umap").
    color : str
        Observation column used for coloring (e.g., "leiden" or "cluster_name").
    ncol : int
        Number of columns in the subplot layout.
    figsize : tuple
        Base figure size (width, height).
    sort_order : list[str], optional
        Custom order for sorting subsets (e.g., ["Y", "M", "O"]).
    cluster_to_cell_mapping : dict, optional
        Mapping from cluster names to cell type labels.
    legend_location : str
        Controls legend placement:
            - "on_data" → Cluster labels overlaid on embedding.
            - "right" → Individual legends on each subplot.
            - "none" → Only a single combined legend below all subplots.
    point_size : float, optional
        Point size for scatter plot. Defaults to Scanpy’s size formula.
    return_fig : bool
        If True, returns (fig, axes) instead of displaying the figure.
    """

    # --- normalize basis ---
    if basis.startswith("X_"):
        basis_key = basis          # for adata.obsm access
        basis_plot = basis[2:]     # for sc.pl.embedding
    else:
        basis_key = f"X_{basis}"
        basis_plot = basis

    # Ensure the color column is string-based and categorical
    adata.obs[color] = adata.obs[color].astype(str)

    # Validate color mapping
    if f"{color}_colors" not in adata.uns:
        raise ValueError(f"Color mapping not found in adata.uns['{color}_colors']. Check if colors were set.")

    original_categories = adata.obs[color].astype("category").cat.categories
    original_colors = adata.uns[f"{color}_colors"]

    if len(original_colors) < len(original_categories):
        raise ValueError(f"Colors for '{color}' are missing or incomplete in adata.uns['{color}_colors'].")

    # Default point size scaling
    if point_size is None:
        point_size = 120000 / adata.n_obs

    # --- Prepare grouping combinations ---
    unique_combinations = adata.obs[obs_columns].drop_duplicates().copy()
    unique_combinations["combined"] = unique_combinations.apply(lambda row: "_".join(map(str, row)), axis=1)

    # Sort panels if sort_order provided
    if sort_order is not None:
        unique_combinations = unique_combinations.sort_values(
            by="combined",
            key=lambda x: [
                sort_order.index(i.split("_")[0]) if i.split("_")[0] in sort_order else float('inf') for i in x
            ],
        )

    sorted_combinations = unique_combinations.drop(columns=["combined"])
    num_plots = len(sorted_combinations)
    nrow = (num_plots + ncol - 1) // ncol

    # --- Create subplot grid ---
    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(figsize[0], nrow * figsize[1] // ncol),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    axes = axes.flat

    for ax in axes[num_plots:]:
        ax.set_visible(False)

    # Track all clusters that appear in any subset
    all_present_clusters = set()

    # --- Iterate over subsets ---
    for idx, row in enumerate(sorted_combinations.itertuples(index=False)):
        subset_query = " & ".join([f"{col} == '{getattr(row, col)}'" for col in obs_columns])
        adata_sub = adata[adata.obs.eval(subset_query), :].copy()

        title = " | ".join([f"{col}: {getattr(row, col)}" for col in obs_columns])

        # Ensure consistent category order
        adata_sub.obs[color] = adata_sub.obs[color].astype("category")
        adata_sub.obs[color] = adata_sub.obs[color].cat.set_categories(original_categories)

        # Find clusters present in this subset
        present_clusters = [cat for cat in original_categories if cat in adata_sub.obs[color].values]
        all_present_clusters.update(present_clusters)

        # Assign filtered colors (in correct order)
        filtered_colors = [
            original_colors[i] for i, cat in enumerate(original_categories) if cat in present_clusters
        ]
        adata_sub.uns[f"{color}_colors"] = filtered_colors

        # Plot each panel
        sc.pl.embedding(
            adata_sub,
            basis=basis_plot,
            color=color,
            title=title,
            legend_loc="none",
            show=False,
            ax=axes[idx],
            size=point_size,
        )

        # Add cluster names directly on the UMAP
        if legend_location == "on_data":
            umap_coords = adata_sub.obsm[basis_key]
            clusters = adata_sub.obs[color].astype(str)
            cluster_medians = {
                c: np.median(umap_coords[clusters == c], axis=0) for c in present_clusters
            }
            for cluster, (x, y) in cluster_medians.items():
                axes[idx].text(
                    x, y, str(cluster),
                    fontsize=10, weight="bold", ha="center", va="center", color="black",
                )

    # --- Global legend ---
    handles = [
        plt.Line2D(
            [0], [0],
            marker="o", color="w",
            markerfacecolor=original_colors[i],
            markersize=10,
            label=(
                f"{cat} - {cluster_to_cell_mapping.get(str(cat), 'Unknown')}"
                if cluster_to_cell_mapping
                else f"{cat}"
            ),
        )
        for i, cat in enumerate(original_categories)
        if cat in all_present_clusters
    ]

    if legend_location != "right":
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            fontsize=8,
            ncol=min(ncol, len(handles)),
        )

    plt.tight_layout()
    if return_fig:
        return fig, axes
    plt.show()


def preserve_leiden_categories(adata_sub, adata_full, leiden_column="leiden"):
    """
    Ensures that the leiden column in a subset AnnData object retains the full category order
    and color mapping from the original AnnData object.

    Parameters:
    - adata_sub: Subset AnnData object.
    - adata_full: Full AnnData object (before subsetting).
    - leiden_column: The name of the leiden clustering column in `adata.obs`.

    Returns:
    - Modified `adata_sub` with preserved leiden categories and color mapping.
    """
    # Ensure leiden is categorical in both full and subset data
    if leiden_column not in adata_full.obs:
        raise ValueError(f"Column '{leiden_column}' not found in full AnnData object.")
    
    if leiden_column not in adata_sub.obs:
        raise ValueError(f"Column '{leiden_column}' not found in subset AnnData object.")

    # Convert leiden column to categorical and preserve full category order
    if not pd.api.types.is_categorical_dtype(adata_full.obs[leiden_column]):
        adata_full.obs[leiden_column] = adata_full.obs[leiden_column].astype("category")

    adata_sub.obs[leiden_column] = adata_sub.obs[leiden_column].astype("category")
    
    # Set categories to match the full dataset
    adata_sub.obs[leiden_column] = adata_sub.obs[leiden_column].cat.set_categories(
        adata_full.obs[leiden_column].cat.categories
    )

    # Preserve the original leiden color mapping
    if f"{leiden_column}_colors" in adata_full.uns:
        adata_sub.uns[f"{leiden_column}_colors"] = adata_full.uns[f"{leiden_column}_colors"]

    return adata_sub

 