from pathlib import Path
import matplotlib.pyplot as plt
import datetime
import pandas as pd # Optional, for clean table display
class FigureRegistry:
    """
    A lightweight container for managing, annotating, displaying, and batch-saving matplotlib figures in Jupyter workflows.

    Features:
    - Tracks figures by name (and optional tags) with creation timestamps.
    - Batch-saving of all figures (e.g., as PDFs) to a chosen output folder.
    - Convenient `.show()` gallery for Jupyter (with names/tags above each figure).
    - Automatic rasterization of collections to control output file size.
    - Easy repr for auditing what's stored.

    Parameters
    ----------
    output_dir : str, optional
        Directory to save figures (default "figures").
    dpi : int, optional
        DPI for figure saving (default 300).

    Methods
    -------
    add(name, fig=None, rasterize=True, tags=None):
        Add a figure (current or supplied) to the registry, with optional tags and rasterization.
    save_all(overwrite=True):
        Save all registered figures to disk (default: PDF).
    show():
        Display all figures in a Jupyter notebook cell
        (with labels and tags).
    list_figures():
        List all tracked figure names.
    __repr__():
        Pretty print a summary table of tracked figures.

    Example
    -------
    >>> from pathlib import Path
    >>> import matplotlib.pyplot as plt
    >>> reg = FigureRegistry(output_dir="myfigs", dpi=200)
    >>> # Add a figure
    >>> fig, ax = plt.subplots()
    >>> ax.plot([1,2,3], [1,4,9])
    >>> reg.add("parabola", fig, tags=["demo", "math"])
    >>> # Show all figures in Jupyter
    >>> reg.show()
    >>> # Save all to PDF in "myfigs" directory
    >>> reg.save_all()
    >>> print(reg)
    Name                      | Tags                 | Created
    ----------------------------------------------------------------------
    parabola                  | demo, math           | 2026-01-07 16:24:33
    """

    
    def __init__(self, output_dir="figures", dpi=300):
        self.output_path = Path(output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.dpi = dpi
        self.registry = {}

    def add(self, name, fig=None, rasterize=True, tags=None):
        """Add figure and metadata to container."""
        fig = fig or plt.gcf()
        if rasterize:
            for ax in fig.axes:
                for c in ax.collections:
                    c.set_rasterized(True)
        
        self.registry[name] = {
            "fig": fig,
            "tags": tags or [],
            "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        return fig

    def __repr__(self):
        """Text summary of what is currently in the registry."""
        if not self.registry:
            return "FigureRegistry is empty."
        summary = [f"{'Name':<25} | {'Tags':<20} | {'Created'}"]
        summary.append("-" * 70)
        for name, meta in self.registry.items():
            tags_str = ", ".join(meta['tags'])
            summary.append(f"{name:<25} | {tags_str:<20} | {meta['timestamp']}")
        return "\n".join(summary)

    def show(self):
        """Visual gallery for Jupyter: displays every figure with its name and tags."""
        from IPython.display import display, HTML
        for name, meta in self.registry.items():
            display(HTML(f"<h3>{name} <small style='color:gray;'>({', '.join(meta['tags'])})</small></h3>"))
            display(meta["fig"])

    def save_all(self, overwrite=True):
        """Saves all figures as PDFs."""
        for name, meta in self.registry.items():
            fig = meta["fig"]
            file_path = self.output_path / f"{name}.pdf"
            fig.savefig(file_path, dpi=self.dpi, bbox_inches="tight")
            print(f"Saved: {file_path}")