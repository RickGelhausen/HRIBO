HRIBO |release|
================

HRIBO is a Snakemake workflow for reproducible processing and analysis of
bacterial ribosome-profiling data.  It combines read processing and quality
control with coverage tracks, metagene profiles, ORF prediction, matched
RNA/Ribo differential analysis, and consolidated result workbooks.

This is the canonical documentation and is maintained with the workflow
source.  It was consolidated and rewritten from the former
``HRIBO_ReadTheDocs`` repository, which was maintained from 2020 to 2024 by
Rick Gelhausen with contributions by Florian Eggenhofer.  Its audited head was
``dcab179968a07cb83617429b972d35837e23d8f4``; the complete history is retained
on the ``archive/hribo-readthedocs`` branch.  See
:ref:`documentation-source-and-legacy-archive` for the migration record.

.. important::

   |release| is under development.  Its automated suite and production
   container boundaries are validated, but the final comparison against a
   representative biological dataset remains a release gate.  See
   :doc:`real-data-validation` for the exact status and protocol.

Using HRIBO
-----------

.. toctree::
   :maxdepth: 2

   getting-started
   samples
   configuration
   stages
   outputs
   table-reference
   metagene-profiling
   tis-advisor

Guides
------

.. toctree::
   :maxdepth: 2

   tutorials/minimal
   tutorials/full
   historical-example-data
   migration-1.8-to-2.0
   real-data-validation
   troubleshooting

Project information
-------------------

.. toctree::
   :maxdepth: 1

   development
   references

Release history and support
---------------------------

See the repository `changelog
<https://github.com/RickGelhausen/HRIBO/blob/development/ChangeLog.md>`_ for
version history.  Report reproducible defects or documentation gaps in the
`HRIBO issue tracker <https://github.com/RickGelhausen/HRIBO/issues>`_ and
include the HRIBO commit, configuration, first failing rule, and relevant log.

License
-------

HRIBO is distributed under the GNU General Public License version 3.  When
publishing work that uses HRIBO, cite the HRIBO paper listed in
:doc:`references`; machine-readable citation metadata is provided in
``CITATION.cff``.
