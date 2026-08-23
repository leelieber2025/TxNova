"""Sphinx configuration for the TxNova documentation (built on Read the Docs)."""

from __future__ import annotations

import os
from datetime import datetime

project = "TxNova"
author = "Zhao Li (李钊)"
copyright = f"{datetime.now():%Y}, {author}"
repository_url = "https://github.com/leelieber2025/TxNova"
default_branch = "main"
release = "0.1.8"
version = release

html_context = {
    "display_github": True,
    "github_user": "leelieber2025",
    "github_repo": "TxNova",
    "github_version": default_branch,
    "conf_py_path": "/docs/",
}

extensions = [
    "myst_nb",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "sphinx_design",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "substitution",
]
myst_heading_anchors = 3
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}
nb_execution_mode = "off"
nb_merge_streams = True

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "DEVELOPMENT.md",
    "DATASETS.md",
    "LOCUS_HYPOTHESIS.md",
    "draft",
    "draft_new",
    "requirements.txt",
    ".ipynb_checkpoints",
]
needs_sphinx = "5.0"
nitpicky = False

html_theme = "sphinx_book_theme"
html_title = project

html_theme_options = {
    "repository_url": repository_url,
    "repository_branch": default_branch,
    "path_to_docs": "docs",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_source_button": True,
    "use_issues_button": True,
    "use_download_button": True,
    "use_fullscreen_button": True,
    "home_page_in_toc": True,
    "show_navbar_depth": 1,
    "navigation_with_keys": True,
}

pygments_style = "tango"
pygments_dark_style = "monokai"

html_static_path = ["_static"]
html_extra_path = ["_extra"]
html_css_files = ["css/custom.css"]
html_show_sphinx = False

if os.environ.get("READTHEDOCS"):
    html_baseurl = "https://txnova.readthedocs.io/"
