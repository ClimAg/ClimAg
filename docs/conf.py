# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys

sys.path.append("..")

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "ClimAg"
copyright = "2022-2024, Nithiya Streethran"
author = "Nithiya Streethran"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc"]

# disable sorting of functions by alphabetical order
autodoc_member_order = "bysource"

# templates_path = ["_templates"]

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    ".ipynb_checkpoints",
    ".pytest_cache",
    "__pycache__",
    ".coverage",
    "environment.yml",
    "convert-notebooks.sh",
    ".gitignore",
    "README.md",
    "LICENCE.txt",
    "requirements.txt",
    ".readthedocs.yaml",
    "methods",
    "other",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "alabaster"
# html_static_path = ["_static"]
