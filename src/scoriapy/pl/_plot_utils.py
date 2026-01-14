import matplotlib.pyplot as plt
import scanpy as sc

def add_ncells_annotation(
    loc=(0.98, 0.02),
    fontsize=9,
    color="black",
    bbox=True
):
    """
    Decorator that adds the total cell count (n = XXX) to a Scanpy plotting function.

    This decorator automatically extracts the observation count from the AnnData 
    object and places a text annotation on the resulting axis. It forces 
    `show=False` on the wrapped function to allow for the annotation to be 
    added before the figure is displayed.

    Args:
        loc (tuple): XY coordinates in axes-relative units (0 to 1). 
            Default is (0.98, 0.02) (bottom right).
        fontsize (int): Size of the annotation text.
        color (str): Color of the annotation text.
        bbox (bool): If True, adds a semi-transparent white background box 
            behind the text for readability.

    Returns:
        function: The wrapped plotting function.

    Example:
        >>> # Standard usage as a decorator
        >>> @add_ncells_annotation(loc=(0.05, 0.05), fontsize=10)
        >>> def plot_custom_umap(adata, **kwargs):
        ...     return sc.pl.umap(adata, **kwargs)
        ...
        >>> # Using it directly on an existing Scanpy function
        >>> annotated_umap = add_ncells_annotation()(sc.pl.umap)
        >>> annotated_umap(adata, color='cell_type')
    """
    def decorator(plot_func):
        def wrapper(adata, *args, **kwargs):
            # Extract ax if provided
            ax = kwargs.get("ax", None)

            # Ensure show=False to prevent immediate rendering
            kwargs = dict(kwargs)
            kwargs["show"] = False

            # Call original plotting function
            result = plot_func(adata, *args, **kwargs)

            # Get the current axis if none was provided
            if ax is None:
                ax = plt.gca()

            # Prepare annotation styling
            txt_kwargs = dict(
                ha="right",
                va="bottom",
                transform=ax.transAxes,
                fontsize=fontsize,
                color=color,
            )

            if bbox:
                txt_kwargs["bbox"] = dict(
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.6,
                    pad=1
                )

            # Add the "n = ..." text
            ax.text(
                loc[0],
                loc[1],
                f"n = {adata.n_obs:,}",
                **txt_kwargs
            )

            return result
        return wrapper
    return decorator
