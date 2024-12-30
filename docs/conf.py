import datetime

project = "Ersilia"
copyright = f"{datetime.datetime.now().year}, Ersilia Open Source Initiative"
author = "Miquel Duran-Frigola and Abel Legese"

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
