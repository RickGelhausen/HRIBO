Upgrade from HRIBO 1.8 to 2.0
==============================

HRIBO 2.0 changes how the workflow is launched, configured, and organized.
Run it in a new analysis directory; do not point 2.0 at an existing 1.8 result
tree.  Keeping the old results unchanged makes comparisons possible and avoids
old intermediate files being mistaken for current work.

Before starting, save the 1.8 configuration and sample sheet and record the
exact genome, annotation, FASTQ files, and HRIBO version used.

Install and launch 2.0
----------------------

HRIBO is now run from a standard Snakemake workflow layout.  The workflow entry
point is ``workflow/Snakefile`` and the user templates are in ``config/``.
Create the launcher environment from the pinned file and use the supplied
script from the new analysis directory:

.. code-block:: console

   $ micromamba create --name hribo \
       --file /path/to/HRIBO/environment.linux-64.pin.txt
   $ micromamba activate hribo
   $ cd /path/to/new-analysis
   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --cores 20 \
       --dry-run all

The 2.0 execution path uses both Conda and Apptainer.  Old commands based on
``--use-conda``, ``--use-singularity``, or legacy cluster submission flags do
not describe the current launchers.  See :doc:`getting-started` for local and
SLURM commands.

Start with new templates
------------------------

Copy the 2.0 templates and transfer project values deliberately:

.. code-block:: console

   $ mkdir -p config
   $ cp /path/to/HRIBO/config/config.yaml config/config.yaml
   $ cp /path/to/HRIBO/config/samples.tsv config/samples.tsv

Do not copy the complete 1.8 YAML over the new file.  HRIBO 2.0 rejects unknown
settings inside its configuration sections.  The main user-visible changes
are:

* ``workflowSettings.workflow`` is replaced by
  ``workflowSettings.stages``.  It accepts ``full``, ``preprocessing``, or a
  list of names from :doc:`stages`.
* Differential analysis is selected with the ``differential_expression``
  stage rather than an on/off key.  The ``full`` preset includes it.
* ``differentialExpressionSettings.contrasts`` is a YAML list such as
  ``["Treated-Control"]``.  Positive effects mean higher signal in the
  condition on the left.
* ``predictionSettings.deepriboASiteOffset`` explicitly sets the distance from
  the 3' read end to the A-site.  It is not the P-site offset reported by the
  TIS advisor.
* Metagene and TIS-advisor read lengths, mapping methods, filters, and plotting
  choices are explicit in the configuration.

The sample-sheet headers remain:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2

Validation is stricter.  Methods are case-sensitive, condition names are
alphanumeric, replicate numbers cannot have leading zeros, and
``method``/``condition``/``replicate`` combinations must be unique.  Keep the
``fastqFile2`` header and leave its cell empty for single-end data.  See
:doc:`samples` for the complete format.

Result changes to account for
-----------------------------

Update downstream scripts to use the current paths and column names described
in :doc:`outputs` and :doc:`table-reference`.  In particular:

* final BAMs and indexes are exposed through ``maplink/``;
* the main combined results are
  ``auxiliary/overview.xlsx`` and ``auxiliary/overview.tsv``;
* browser-ready GFF files are in ``tracks/`` and ``auxiliary/``;
* differential workbooks use consistent effect and adjusted-p-value column
  names and include ``all`` plus filtered up/down sheets; and
* the old result-archive helper is no longer used.  Generate an optional native
  Snakemake report with ``--report report.html`` after the selected targets
  finish.

Do not load count tables below ``readcounts/`` into a genome browser merely
because their historical suffix is ``.gff`` or ``.gtf``.  Use
``tracks/updated_annotation.gff``, ``auxiliary/overview.gff``, or the motif
tracks instead.

Why numerical results may differ
--------------------------------

Several corrections in 2.0 can change values even when the biological inputs
are identical:

* RPKM, CPM, and normalized coverage use one mapped-read denominator summed
  across all contigs in a library.  This matters especially for genomes with
  plasmids or multiple chromosomes.
* Metagene windows are consistently oriented in transcript direction for both
  strands.  Start and stop axes therefore have clear but different zero
  positions, as described in :doc:`metagene-profiling`.
* Coverage follows aligned CIGAR blocks rather than treating deletions,
  reference skips, or soft clips as observed reference bases.
* Prediction parsing, coordinates, and combined GFF publication have been
  tightened, so candidate membership or boundaries can change.

These differences mean that a successful 2.0 run should not be expected to be
byte-identical to a 1.8 result.  Review MultiQC, compare feature counts and
coverage, inspect changed ORF coordinates in a genome browser, and re-evaluate
thresholds that were chosen using 1.8-normalized values.

Migration checklist
-------------------

Before switching an analysis or downstream pipeline to 2.0, confirm that:

* the old and new analyses use the same biological inputs and sample labels;
* every 1.8 option has been transferred to its current setting intentionally;
* the 2.0 dry-run passes before computation starts;
* expected stages and primary results are present;
* contrast directions and normalization choices are understood; and
* downstream parsers select columns by name rather than by a fixed position.
