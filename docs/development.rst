Development and verification
============================

HRIBO is maintained as a repository-run Snakemake workflow.  Changes should
keep the workflow reproducible, make source dependencies visible to Snakemake,
and distinguish software integration tests from biological validation.

Development environment
-----------------------

Create the exact development environment from the explicit Linux pin:

.. code-block:: console

   $ micromamba create --name hribo-dev \
       --file environment-dev.linux-64.pin.txt
   $ micromamba activate hribo-dev

``environment-dev.yaml`` is the readable source manifest.  The explicit file
also includes the validators and runtime libraries needed to execute every
collected test without skips.  It is intentionally broader than the production
launcher environment; scientific commands still run in their rule-specific
environments or containers.

Do not use ``pip install .``: the repository has no supported Python package
interface.  Do not fix a test by installing an undeclared package into an
existing environment.  Add an approved dependency to the appropriate YAML,
regenerate its explicit pin, and review the solver diff.  Documentation is the
one separate Python toolchain.  Its direct dependencies live in
``docs/requirements.in`` and its complete Linux/Python 3.12 build closure is
pinned in ``docs/requirements.txt``.

Required checks
---------------

Run the complete Python suite from the repository root:

.. code-block:: console

   $ pytest -q

CI records JUnit output, verifies the expected collection count, and rejects
every skip.  If collection changes intentionally, inspect ``pytest
--collect-only -q`` and update the CI baseline only after accounting for every
added or removed test.

Run the maintained static checks as well:

.. code-block:: console

   $ ruff check .github/scripts tests workflow/scripts
   $ yamllint -c .yamllint.yaml \
       environment.yaml environment-dev.yaml .readthedocs.yaml CITATION.cff \
       config workflow/conda-lock workflow/envs workflow/schemas \
       .github/workflows
   $ python .github/scripts/check_config.py
   $ python .github/scripts/check_conda_pins.py
   $ git diff --check

Shell and R sources are syntax-checked independently of workflow execution:

.. code-block:: bash

   mapfile -t shell_files < <(
       find . workflow/scripts workflow/envs -maxdepth 2 \
           -type f -name '*.sh' -print
   )
   bash -n "${shell_files[@]}"
   shellcheck "${shell_files[@]}"

   mapfile -t r_files < <(
       find workflow/scripts -maxdepth 1 -type f -name '*.R' -print
   )
   for r_file in "${r_files[@]}"; do
       Rscript -e 'parse(file=commandArgs(trailingOnly=TRUE)[1])' "$r_file"
   done

Build the documentation with the separately pinned documentation environment.
Warnings, missing references, and malformed citations are errors:

.. code-block:: console

   $ python3.12 -m venv /tmp/hribo-docs
   $ /tmp/hribo-docs/bin/python -m pip install \
       --requirement docs/requirements.txt
   $ make -C docs strict \
       SPHINXBUILD=/tmp/hribo-docs/bin/sphinx-build

Before publishing documentation that changes external references, also check
its links (this step requires network access):

.. code-block:: console

   $ make -C docs linkcheck \
       SPHINXBUILD=/tmp/hribo-docs/bin/sphinx-build

When changing a documentation dependency, edit ``docs/requirements.in`` and
regenerate the full lock from a clean Python 3.12 environment with the recorded
compiler version.  Review all transitive changes before accepting them:

.. code-block:: console

   $ python3.12 -m venv /tmp/hribo-doc-lock
   $ /tmp/hribo-doc-lock/bin/python -m pip install pip-tools==7.6.1
   $ /tmp/hribo-doc-lock/bin/python -m piptools compile \
       --upgrade --strip-extras docs/requirements.in
   $ git diff -- docs/requirements.in docs/requirements.txt

Rule and script dependencies
----------------------------

Repository scripts executed by a rule are inputs to that rule, not invisible
implementation details.  Declare an entry-point script as a named input and
invoke it through its quoted ``{input.<name>:q}`` value.  Declare the complete
closure of repository-local Python imports as additional inputs.  For a
container rule, resolve repository code with ``workflow.source_path`` so it is
available inside the container.

``tests/test_rule_script_dependencies.py`` enforces these conventions.  They
matter for correctness: a helper edit must invalidate the jobs whose results it
can change, while dependency-only files must not be appended accidentally to a
tool's positional arguments.

Likewise, declare every ordinary deliverable in ``output``.  When an external
tool produces an indivisible group that Snakemake cannot publish safely as
separate files, use a checked transactional runner and a completion receipt.
Do not add a new untracked scientific side effect.

Tests and golden outputs
------------------------

Keep a regression test close to the boundary being changed:

* unit tests cover parsing, coordinate geometry, normalization, and error
  handling;
* DAG and execution tests cover stage selection, dependencies, quoting, and
  stable no-op reruns;
* golden spreadsheet snapshots expose workbook values as reviewable CSV;
* golden GFF3 files are checked both semantically and with a strict GFF3
  validator; and
* production-container smokes exercise the actual pinned deltaTE, REPARATION,
  and DeepRibo boundaries.

Regenerate golden files only for a deliberate output-contract or scientific
change:

.. code-block:: console

   $ python tests/regenerate_golden.py
   $ python tests/regenerate_golden_gff.py
   $ git diff -- tests/golden tests/golden_gff

Read the complete diff before accepting it.  A snapshot update is not a
substitute for an assertion that explains the intended new behavior.

Container smokes
----------------

The three smoke programs under ``.github/scripts/`` run bounded fixtures
through production image digests and validate their result contracts:

* ``deltate_container_smoke.py`` covers input adapters, the patched engine,
  result tables, report pages, and a no-op rerun;
* ``reparation_container_smoke.py`` covers annotation adaptation, transactional
  publication, GFF3 conversion, and a no-op rerun; and
* ``deepribo_container_smoke.py`` covers BAM-derived signals, parsing, cutoff
  publication, real pinned-model inference, GFF3 conversion, and a no-op
  rerun.

They require Linux/amd64, Apptainer, substantial temporary space, network or
pre-populated caches, and the pinned development environment.  Their command
line requires isolated work, Conda-prefix, and Apptainer-prefix directories;
use ``--help`` and the corresponding CI job as the authoritative invocation.
Do not point a smoke at a biological result directory.

These fixtures establish executable integration contracts, not biological
sensitivity or equivalence.  A release that changes scientific code must also
complete the representative real-data protocol described in
:doc:`real-data-validation`.

Changing dependencies
---------------------

Edit the appropriate readable manifest first: ``environment.yaml`` for the
launcher, ``environment-dev.yaml`` for maintenance tools, or one file under
``workflow/envs/`` for a rule.  Keep channels in the approved order and use an
explicit direct version.  Then regenerate all Linux pins with the repository's
approved lock tool:

.. code-block:: console

   $ python .github/scripts/update_conda_pins.py
   $ python .github/scripts/check_conda_pins.py

Pin regeneration resolves packages and therefore needs network access.  Review
every changed URL, version, and checksum, run the full suite, and run any
container smoke or real-data comparison affected by the dependency boundary.

.. _documentation-source-and-legacy-archive:

Documentation source and legacy archive
---------------------------------------

The canonical documentation is the ``docs/`` tree in this repository.  The
root ``.readthedocs.yaml`` builds that tree with the dependencies pinned in
``docs/requirements.txt``; documentation changes should therefore accompany
the workflow change they describe and pass the strict build above.

The former ``HRIBO_ReadTheDocs`` repository was audited at commit
``dcab179968a07cb83617429b972d35837e23d8f4`` before consolidation.  Its
complete 93-commit history is retained on the namespaced branch
``archive/hribo-readthedocs``.  Its four release tags are retained as
``archive/hribo-readthedocs/1.4.0``, ``archive/hribo-readthedocs/1.4.2``,
``archive/hribo-readthedocs/1.4.4``, and
``archive/hribo-readthedocs/1.5.1``.  This archive is provenance, not a second
documentation source.

The active documentation was reviewed and rewritten rather than copied
verbatim.  The legacy pages map as follows:

.. list-table:: Legacy documentation disposition
   :header-rows: 1
   :widths: 30 35 35

   * - Former page
     - Current page
     - Disposition
   * - ``overview``
     - :doc:`getting-started`, :doc:`samples`, :doc:`outputs`, and this page
     - Replaced by the current workflow, installation, and output contracts.
   * - ``workflow-configuration``
     - :doc:`configuration`, :doc:`samples`, and :doc:`stages`
     - Rewritten against the validated 2.0 schemas.
   * - ``analysis-results``
     - :doc:`outputs` and :doc:`table-reference`
     - Rebuilt from current declared outputs and tested workbook schemas.
   * - ``metagene-profiling``
     - :doc:`metagene-profiling`
     - Rewritten for the corrected coordinate and normalization contracts.
   * - ``example-workflow`` and ``extended-workflow``
     - :doc:`historical-example-data`, :doc:`tutorials/minimal`, and
       :doc:`tutorials/full`
     - Dataset provenance retained; obsolete commands and unvalidated results
       omitted.
   * - ``faq``
     - :doc:`troubleshooting`
     - Obsolete Singularity-specific advice removed.

The legacy NCBI-interface screenshots, workflow diagram, and metagene figures
are intentionally available only through the archive branch.  The first two no
longer reflect the current interfaces or workflow; the metagene figures
predate the corrected strand/anchor geometry.  Empty scratch files, stale
download recipes, and the old dependency workaround were likewise not moved
into the maintained tree.  The GPL-3.0 license in both repositories was
identical, so the root ``LICENSE`` remains authoritative.

The maintained tree is deliberately lean.  Every ``.rst`` page is reachable
from an ``index.rst`` toctree.  Its only non-page inputs are ``conf.py``,
``requirements.in``, ``requirements.txt``, ``references.bib``, and this local
``Makefile``; each has a current build or maintenance role.  Generated
``_build/`` and ``__pycache__/`` directories are ignored.  Do not copy files or
images from the archive into ``docs/`` unless a maintained page actually uses
them and their content still matches the current workflow.

Read the Docs cutover
---------------------

After the documentation commit and archive refs are pushed, update the
existing Read the Docs project rather than creating a second public project:

#. In the project's administration settings, change the repository URL from
   ``HRIBO_ReadTheDocs`` to this ``HRIBO`` repository.
#. Activate and build the intended branch.  Confirm in the build log that Read
   the Docs found the root ``.readthedocs.yaml`` and used ``docs/conf.py``.
#. Open the main pages, internal links, version selector, edit links, and
   project badge from the hosted site.
#. If the former repository is to be deleted, do so only after the active pages
   and the ``archive/hribo-readthedocs`` branch are visible on GitHub.  Preserve
   any repository-level issue or settings metadata that is still needed; Git
   history does not contain it.

Changing the service setting and deleting the former repository require
project-owner access and are therefore intentionally separate from the source
migration.
