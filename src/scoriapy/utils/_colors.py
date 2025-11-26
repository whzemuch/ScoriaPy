"""
Utility functions for working with color maps for AnnData plots.
"""

from __future__ import annotations

from typing import Dict, Iterable

from IPython.display import HTML, display


def map_colors(
    reference_dict: Dict[str, str],
    target_dict: Dict[str, str],
    mapping_dict: Dict[str, Iterable[str]],
) -> Dict[str, str]:
    """
    Map colors from a reference dictionary to a target dictionary.

    Parameters
    ----------
    reference_dict
        Mapping from reference category to color (e.g. cluster ID → hex code).
    target_dict
        Mapping whose keys define the target categories to color.
        Existing values are ignored; only the keys are used.
    mapping_dict
        Mapping from reference keys to iterable(s) of target keys. Each
        target key inherits the color of its reference key.

    Returns
    -------
    dict
        New mapping from target keys to colors. Keys not specified in
        ``mapping_dict`` are assigned gray (``\"#808080\"``).
    """
    new_dict = {key: "#808080" for key in target_dict.keys()}

    reversed_mapping = {
        target_key: ref_key
        for ref_key, target_keys in mapping_dict.items()
        for target_key in target_keys
    }

    for target_key, ref_key in reversed_mapping.items():
        if ref_key in reference_dict:
            new_dict[target_key] = reference_dict[ref_key]

    return new_dict


def display_colors(color_list):
    """
    Display colors as inline blocks in a notebook.

    Parameters
    ----------
    color_list
        Iterable of color specifications (e.g. hex strings) to display.
    """
    html_str = "".join(
        "<div style='width:100px; height:50px; line-height:50px; "
        f"color:white; font-weight:bold; background-color:{color}; "
        "display:inline-block; text-align:center;'>"
        f"{color}</div>"
        for color in color_list
    )
    display(HTML(html_str))


def display_colors_with_keys(color_map: Dict[str, str]):
    """
    Display a mapping of keys to colors as inline blocks in a notebook.

    Parameters
    ----------
    color_map
        Mapping from label to color code to visualize.
    """
    html_str = "".join(
        "<div style='display:inline-block; text-align:center; margin:10px'>"
        f"<div style='margin-bottom:5px'>{key}</div>"
        "<div style='width:100px; height:50px; line-height:50px; color:white; "
        f"background-color:{color}; display:inline-block; text-align:center;'>"
        f"{color}</div>"
        "</div>"
        for key, color in color_map.items()
    )
    display(HTML(html_str))


def confirm_and_set_colors(adata, color_attribute: str, color_map: Dict[str, str]):
    """
    Interactively confirm and apply a color map for a categorical obs column.

    Parameters
    ----------
    adata
        AnnData object.
    color_attribute
        Column in ``adata.obs`` used for coloring.
    color_map
        Mapping from category to color code.
    """
    current_color_key = f"{color_attribute}_colors"

    def _update_colors():
        display_colors_with_keys(color_map)
        if input("Confirm color change? (yes/no): ").lower().startswith("y"):
            adata.uns[current_color_key] = [
                color_map[cat] for cat in adata.obs[color_attribute].cat.categories
            ]
            print("Color map updated.")
        else:
            print("Color update canceled.")

    if current_color_key in adata.uns:
        print("Current default colors:")
        display_colors_with_keys(
            dict(
                zip(
                    adata.obs[color_attribute].cat.categories,
                    adata.uns[current_color_key],
                )
            )
        )
        print("The new colors:")
        _update_colors()
    else:
        print("Color key not found. Setting new colors.")
        _update_colors()
