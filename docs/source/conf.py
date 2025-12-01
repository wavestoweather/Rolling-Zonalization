# -- Project information -----------------------------------------------------

project = "rwguide"
copyright = "2025, Christopher Polster"
author = "Christopher Polster"


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon"
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

html_theme = "nature"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = False
html_show_copyright = False
html_show_sphinx = False

