"""Sphinx configuration for the consolidated HRIBO documentation."""

from __future__ import annotations

import os
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parent.parent

project = "HRIBO"
author = "HRIBO contributors"
copyright = "2020–2026, HRIBO contributors"
release = (REPOSITORY / "VERSION").read_text(encoding="utf-8").strip()
version = release.split("-", 1)[0].rsplit(".", 1)[0]

extensions = ["sphinxcontrib.bibtex"]
bibtex_bibfiles = ["references.bib"]

language = "en"
exclude_patterns = ["_build"]
nitpicky = True
root_doc = "index"

# Some publishers return HTTP 403 to automated clients after resolving valid
# DOI links.  Keep those stable identifiers in the bibliography, but exclude
# only the exact affected articles from Sphinx's network link check so a new or
# mistyped DOI cannot be hidden by a publisher-wide pattern.
linkcheck_ignore = [
    r"^https://doi\.org/10\.1093/bioinformatics/btaa959$",
    r"^https://doi\.org/10\.1093/bioinformatics/btq351$",
    r"^https://doi\.org/10\.1093/bioinformatics/bts480$",
    r"^https://doi\.org/10\.1093/bioinformatics/btu146$",
    r"^https://doi\.org/10\.1093/bioinformatics/btx047$",
    r"^https://doi\.org/10\.1093/nar/gkx758$",
    r"^https://doi\.org/10\.1093/nar/gkz061$",
    r"^https://doi\.org/10\.1002/cpmb\.108$",
]

html_theme = "sphinx_rtd_theme"
html_baseurl = os.environ.get(
    "READTHEDOCS_CANONICAL_URL", "https://hribo.readthedocs.io/"
)
html_title = f"HRIBO {release} documentation"

# Point each hosted page back to the branch or tag that produced it.  Pull
# request builds expose only a numeric Git identifier, so their checked-out
# commit is the stable GitHub target.  Local builds default to the 2.0.0 tag.
if os.environ.get("READTHEDOCS_VERSION_TYPE") == "external":
    github_version = os.environ.get("READTHEDOCS_GIT_COMMIT_HASH", "2.0.0")
else:
    github_version = os.environ.get("READTHEDOCS_GIT_IDENTIFIER", "2.0.0")

html_context = {
    "display_github": True,
    "github_user": "RickGelhausen",
    "github_repo": "HRIBO",
    "github_version": github_version,
    "conf_py_path": "/docs/",
}
