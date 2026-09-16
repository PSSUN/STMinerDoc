# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


import os
import sys
from pathlib import Path
from unittest.mock import Mock

MOCK_MODULES = [
    "anndata",
    "bioservices",
    "matplotlib",
    "networkx",
    "numba",
    "numpy",
    "pandas",
    "Pillow",
    "plotly",
    "POT",
    "scanpy",
    "sklearn.cluster",
    "scikit_learn",
    "scipy",
    "seaborn",
    "setuptools",
    "tifffile",
    "tqdm",
    "umap_learn",
    "scikit-misc",
]

# for mod_name in MOCK_MODULES:
#     sys.modules[mod_name] = Mock()


# Prefer an explicitly selected checkout, then the sibling development checkout.
# Do not shadow the installed package with the historical bundled source copy.
docs_root = Path(__file__).resolve().parents[1]
source_override = os.environ.get("STMINER_SOURCE")
source_root = Path(source_override).resolve() if source_override else docs_root.parent / "STMiner"
if source_override and not (source_root / "STMiner" / "__init__.py").is_file():
    raise RuntimeError("STMINER_SOURCE must point to the STMiner repository root.")
if (source_root / "STMiner" / "__init__.py").is_file():
    sys.path.insert(0, str(source_root))


# -- Project information -----------------------------------------------------

project = "STMiner"
copyright = "2025, Peisen Sun"
author = "Peisen Sun"

# The full version, including alpha/beta/rc tags
import STMiner

release = STMiner.__version__
if not hasattr(STMiner.SPFinder, "enrich_patterns"):
    raise RuntimeError(
        "Documentation requires STMiner with enrich_patterns(); install the matching "
        "source checkout or set STMINER_SOURCE to its repository root."
    )


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx_markdown_tables",
    "sphinx_design",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon"
]
myst_enable_extensions = ["colon_fence"]

autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
master_doc = "index"
