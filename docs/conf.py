# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------
import datetime

project = "Ersilia"

now = datetime.datetime.now()
copyright = "{0}, Ersilia Open Source Initiative".format(now.year)
author = "Miquel Duran-Frigola"

# The short X.Y version
version = ""
release = ""

extensions = [
    "sphinx.ext.intersphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
]

autosummary_generate = True

templates_path = ["_templates"]

source_suffix = ".rst"
master_doc = "index"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_theme_options = {"collapse_navigation": True, "navigation_depth": 4}

html_static_path = []
htmlhelp_basename = "ersilia_doc"

pygments_style = "sphinx"

latex_documents = [
    (master_doc, "ersilia.tex", "Ersilia Documentation", author, "manual"),
]

man_pages = [(master_doc, "ersilia", "Ersilia Documentation", [author], 1)]

texinfo_documents = [
    (
        master_doc,
        "ersilia",
        "Ersilia Documentation",
        author,
        "ersilia",
        "One line description of the project.",
        "Miscellaneous",
    )
]

epub_title = project
epub_exclude_files = ["search.html"]
