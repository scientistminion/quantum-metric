# Configuration file for the Sphinx documentation builder.
# For the full list of options, see: https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime

# Make the package importable by Sphinx for autodoc.
sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information -----------------------------------------------------
project = "quantum-metric"
author = "Your Name"
copyright = f"{datetime.now().year}, {author}"

# Pull the version from the package
try:
    from quantum_metric._version import __version__
    release = __version__
    version = __version__
except ImportError:
    release = "0.1.0"
    version = "0.1.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",         # pull docstrings into docs
    "sphinx.ext.autosummary",     # auto-generate API summary tables
    "sphinx.ext.napoleon",        # support NumPy / Google style docstrings
    "sphinx.ext.viewcode",        # add links to source code
    "sphinx.ext.intersphinx",     # cross-link to numpy/pandas docs
    "myst_parser",                # allow Markdown alongside reStructuredText
    "sphinx_copybutton",          # copy button on code blocks
    "sphinx_click",               # document Click/Typer CLIs (optional)
]

# Autogenerate stub pages for API reference
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Allow both .md and .rst source files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_book_theme"
html_title = f"{project} v{version}"
html_static_path = ["_static"]

html_theme_options = {
    "repository_url": "https://github.com/yourusername/quantum-metric",
    "repository_branch": "main",
    "path_to_docs": "docs/source",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_issues_button": True,
    "use_download_button": True,
    "home_page_in_toc": True,
    "show_toc_level": 2,
    "show_navbar_depth": 2,
}

# Intersphinx: link to other projects' docs
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# MyST options (lets you write $math$ etc. in Markdown)
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
]
