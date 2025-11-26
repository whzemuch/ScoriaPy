"""
Sphinx configuration for the ScoriaPy documentation.

This is a minimal setup tailored for Read the Docs. It uses autodoc and
autosummary to build API pages from the scoriapy package.
"""

from __future__ import annotations

import os
import sys
from datetime import datetime


# -- Path setup --------------------------------------------------------------

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC_DIR = os.path.join(ROOT_DIR, "src")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)


# -- Project information -----------------------------------------------------

project = "ScoriaPy"
author = "ScoriaPy contributors"
copyright = f"{datetime.now():%Y}, {author}"

# The full version, including alpha/beta/rc tags
try:
    import scoriapy

    release = getattr(scoriapy, "__version__", "0.0.0")
except Exception:
    release = "0.0.0"


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "imported-members": True,
    "show-inheritance": False,
}
napoleon_google_docstring = False
napoleon_numpy_docstring = True

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
